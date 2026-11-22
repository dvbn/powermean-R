#define R_NO_REMAP
#include <R.h>
#include <R_ext/Lapack.h>
#include <Rinternals.h>

#include "native_calls.h"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace {

constexpr int max_columns = 2048;

#ifdef La_INT
using LapackInt = La_INT;
#else
using LapackInt = int;
#endif

}  // namespace

SEXP cancelled_csc_lu_call(SEXP columns, SEXP rows, SEXP entries, SEXP row_sum_scaling) {
    if (TYPEOF(columns) != INTSXP || TYPEOF(rows) != INTSXP || TYPEOF(entries) != REALSXP ||
        TYPEOF(row_sum_scaling) != LGLSXP || Rf_xlength(row_sum_scaling) != 1 ||
        LOGICAL(row_sum_scaling)[0] == NA_LOGICAL ||
        Rf_xlength(columns) < 2 || Rf_xlength(columns) > max_columns + 1 ||
        Rf_xlength(rows) != Rf_xlength(entries)) {
        return R_NilValue;
    }
    const bool use_sums = LOGICAL(row_sum_scaling)[0];
    const auto width = static_cast<int>(Rf_xlength(columns) - 1);
    const auto size = static_cast<std::size_t>(width);
    const LapackInt order = width;
    const int* offsets = INTEGER(columns);
    const int* indices = INTEGER(rows);
    const double* values = REAL(entries);
    if (offsets[0] != 0 || offsets[width] != Rf_xlength(entries)) {
        return R_NilValue;
    }
    auto* scale = reinterpret_cast<double*>(R_alloc(size, sizeof(double)));
    std::fill_n(scale, size, 0.0);
    auto* column_sum = use_sums ? reinterpret_cast<double*>(R_alloc(size, sizeof(double))) : nullptr;
    for (int column = 0; column < width; ++column) {
        if (offsets[column + 1] < offsets[column] || offsets[column + 1] > offsets[width]) {
            return R_NilValue;
        }
        int previous = -1;
        double total = 0;
        for (int at = offsets[column]; at < offsets[column + 1]; ++at) {
            const int row = indices[at];
            if (row <= previous || row >= width || !std::isfinite(values[at])) return R_NilValue;
            previous = row;
            const double amount = std::abs(values[at]);
            if (use_sums) {
                scale[row] += amount;
                total += amount;
            } else {
                scale[row] = std::max(scale[row], amount);
                scale[column] = std::max(scale[column], amount);
            }
        }
        if (use_sums) column_sum[column] = total;
    }
    auto* inverse_scale = use_sums ? reinterpret_cast<double*>(R_alloc(size, sizeof(double))) : nullptr;
    for (int column = 0; column < width; ++column) {
        if (use_sums) {
            scale[column] = std::sqrt(std::sqrt(scale[column]) * std::sqrt(column_sum[column]));
            if (!std::isfinite(scale[column]) || scale[column] <= 0) return R_NilValue;
            inverse_scale[column] = 1 / scale[column];
        } else {
            if (scale[column] == 0) return R_NilValue;
            scale[column] = std::sqrt(scale[column]);
        }
    }

    SEXP factor = PROTECT(Rf_allocMatrix(REALSXP, width, width));
    double* lu = REAL(factor);
    std::fill_n(lu, size * size, 0.0);
    // Use the same diagonal scaling on both sides. Save the one-norm
    // before LAPACK overwrites the matrix with its pivoted LU factors.
    double norm = 0;
    for (int column = 0; column < width; ++column) {
        double* output = lu + column * size;
        double sum = 0;
        for (int at = offsets[column]; at < offsets[column + 1]; ++at) {
            const int row = indices[at];
            output[row] = use_sums ?
                (inverse_scale[row] * values[at]) * inverse_scale[column] :
                values[at] / scale[row] / scale[column];
            sum += std::abs(output[row]);
        }
        norm = std::max(norm, sum);
    }
    if (!std::isfinite(norm)) {
        UNPROTECT(1);
        return R_NilValue;
    }
    auto* pivots = reinterpret_cast<LapackInt*>(R_alloc(size, sizeof(LapackInt)));
    LapackInt info = 0;
    F77_CALL(dgetrf)(&order, &order, lu, &order, pivots, &info);
    if (info != 0) {
        UNPROTECT(1);
        return R_NilValue;
    }
    SEXP pivot = PROTECT(Rf_allocVector(INTSXP, width));
    std::copy_n(pivots, size, INTEGER(pivot));
    SEXP units = PROTECT(Rf_allocVector(REALSXP, width));
    std::copy_n(scale, size, REAL(units));
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 4));
    SET_VECTOR_ELT(result, 0, factor);
    SET_VECTOR_ELT(result, 1, pivot);
    SET_VECTOR_ELT(result, 2, units);
    SET_VECTOR_ELT(result, 3, Rf_ScalarReal(norm));
    UNPROTECT(4);
    return result;
}

SEXP checked_dense_lu_solve(SEXP factor, SEXP pivot, SEXP right, SEXP transpose, int maximum_width) {
    if (TYPEOF(factor) != REALSXP || !Rf_isMatrix(factor) ||
        TYPEOF(right) != REALSXP || !Rf_isMatrix(right) || TYPEOF(pivot) != INTSXP ||
        TYPEOF(transpose) != LGLSXP || Rf_xlength(transpose) != 1 ||
        LOGICAL(transpose)[0] == NA_LOGICAL) {
        Rf_error("Invalid dense Jacobian factor or right-hand side.");
    }
    const int width = Rf_nrows(factor);
    if (width <= 0 || width > maximum_width || Rf_ncols(factor) != width ||
        Rf_xlength(pivot) != width || Rf_nrows(right) != width) {
        Rf_error("The dense Jacobian factor and right-hand side do not align.");
    }
    const LapackInt order = width;
    const LapackInt count = Rf_ncols(right);
    auto* pivots = reinterpret_cast<LapackInt*>(R_alloc(width, sizeof(LapackInt)));
    for (int column = 0; column < width; ++column) {
        const int row = INTEGER(pivot)[column];
        if (row <= column || row > width) Rf_error("Invalid dense Jacobian pivot.");
        pivots[column] = row;
    }
    SEXP result = PROTECT(Rf_duplicate(right));
    const char* operation = LOGICAL(transpose)[0] ? "T" : "N";
    LapackInt info = 0;
    F77_CALL(dgetrs)(operation, &order, &count, REAL(factor), &order, pivots, REAL(result),
                   &order, &info FCONE);
    if (info != 0) Rf_error("Dense Jacobian solve failed.");
    UNPROTECT(1);
    return result;
}

SEXP cancelled_dense_lu_solve_call(SEXP factor, SEXP pivot, SEXP right, SEXP transpose) {
    return checked_dense_lu_solve(factor, pivot, right, transpose, max_columns);
}
