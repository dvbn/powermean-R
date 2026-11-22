#define R_NO_REMAP
#include <R.h>
#include <R_ext/Lapack.h>
#include <R_ext/Utils.h>
#include <Rinternals.h>

#include "native_calls.h"

#include <algorithm>
#include <cstddef>
#include <numeric>

namespace {

#ifdef La_INT
using LapackInt = La_INT;
#else
using LapackInt = int;
#endif

}  // namespace

SEXP bbp_covariance_blocks_call(SEXP transformed, SEXP split, SEXP owners) {
    const SEXP shape = Rf_getAttrib(transformed, R_DimSymbol);
    if (TYPEOF(transformed) != REALSXP || TYPEOF(shape) != INTSXP ||
        Rf_xlength(shape) != 2 || INTEGER(shape)[0] != INTEGER(shape)[1] ||
        TYPEOF(split) != INTSXP || Rf_xlength(split) != 2 || TYPEOF(owners) != INTSXP) {
        return R_NilValue;
    }
    const LapackInt d = INTEGER(shape)[0];
    const LapackInt p = INTEGER(split)[0];
    const LapackInt cut = INTEGER(split)[1];
    if (p < 0 || p >= cut || cut >= d || Rf_xlength(owners) != d - p) return R_NilValue;
    const LapackInt nuisance = d - p;
    const int* country = INTEGER(owners);
    if (std::any_of(country, country + nuisance, [](int value) { return value < 1; })) return R_NilValue;

    auto* order = reinterpret_cast<LapackInt*>(R_alloc(nuisance, sizeof(LapackInt)));
    std::iota(order, order + nuisance, p);
    std::sort(order, order + nuisance, [country, p](LapackInt a, LapackInt b) {
        return country[a - p] == country[b - p] ? a < b : country[a - p] < country[b - p];
    });
    const auto end_country = [&](LapackInt start) {
        LapackInt end = start + 1;
        while (end < nuisance && country[order[end] - p] == country[order[start] - p]) ++end;
        return end;
    };

    const LapackInt left = cut - p;
    const LapackInt right = d - cut;
    double entries = static_cast<double>(p) * (d - (p - 1.0) / 2) +
        static_cast<double>(left) * right;
    for (LapackInt start = 0, end; start < nuisance; start = end) {
        end = end_country(start);
        if (order[start] < cut && order[end - 1] >= cut) return R_NilValue;
        const double count = end - start;
        entries += count * (count + 1) / 2;
    }
    // A dense symmetric crossproduct is cheaper when little work can be skipped.
    if (entries > 0.7 * d * (d + 1.0) / 2) return R_NilValue;

    // A centered row and its group curvature only combine slopes and
    // coordinates of one country from each role. The group score meat
    // has already entered the transformed bread at this point.
    SEXP result = PROTECT(Rf_allocMatrix(REALSXP, d, d));
    double* covariance = REAL(result);
    std::fill_n(covariance, Rf_xlength(result), 0.0);
    const double* matrix = REAL(transformed);
    const auto column = [d](LapackInt index) { return static_cast<std::size_t>(d) * index; };
    const double one = 1.0;
    const double zero = 0.0;
    const LapackInt stride = 1;

    R_CheckUserInterrupt();
    F77_CALL(dgemm)("T", "N", &left, &right, &d, &one, matrix + column(p), &d,
                    matrix + column(cut), &d, &zero, covariance + column(cut) + p, &d FCONE FCONE);
    for (LapackInt j = cut; j < d; ++j) {
        for (LapackInt i = p; i < cut; ++i) covariance[j + column(i)] = covariance[i + column(j)];
    }
    if (p) {
        F77_CALL(dgemm)("T", "N", &p, &d, &d, &one, matrix, &d,
                        matrix, &d, &zero, covariance, &d FCONE FCONE);
        for (LapackInt j = 0; j < d; ++j) {
            for (LapackInt i = 0; i < std::min(p, j); ++i) {
                covariance[j + column(i)] = covariance[i + column(j)];
            }
        }
    }
    for (LapackInt start = 0, end; start < nuisance; start = end) {
        end = end_country(start);
        for (LapackInt a = start; a < end; ++a) {
            const LapackInt i = order[a];
            for (LapackInt b = a; b < end; ++b) {
                const LapackInt j = order[b];
                const double dot = F77_CALL(ddot)(&d, matrix + column(i), &stride,
                                                  matrix + column(j), &stride);
                covariance[i + column(j)] = covariance[j + column(i)] = dot;
            }
        }
    }
    UNPROTECT(1);
    return result;
}
