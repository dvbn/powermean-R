#define R_NO_REMAP
#include <R.h>
#include <R_ext/Utils.h>
#include <Rinternals.h>

#include "native_calls.h"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace {

struct RoleShape {
    int groups;
    int dates;
    int countries;
};

RoleShape role_shape(SEXP matrix, SEXP country, int countries) {
    const SEXP dimensions = Rf_getAttrib(matrix, R_DimSymbol);
    if (TYPEOF(matrix) != REALSXP || TYPEOF(dimensions) != INTSXP || Rf_xlength(dimensions) != 2 ||
        TYPEOF(country) != INTSXP || countries < 1) Rf_error("Invalid WZ role inputs.");
    const int groups = INTEGER(dimensions)[0];
    const int dates = INTEGER(dimensions)[1];
    if (groups < 1 || dates < 2 || Rf_xlength(country) != groups) Rf_error("Misaligned WZ role inputs.");
    const int* codes = INTEGER(country);
    for (int group = 0; group < groups; ++group) {
        if (codes[group] < 1 || codes[group] > countries) Rf_error("Invalid WZ role country.");
    }
    return {groups, dates, countries};
}

void require_matrix_length(SEXP values, R_xlen_t expected) {
    if (TYPEOF(values) != REALSXP || Rf_xlength(values) != expected) Rf_error("Misaligned WZ role values.");
}

}  // namespace

SEXP wz_role_moments_call(SEXP mean, SEXP probability, SEXP score, SEXP country, SEXP countries) {
    if (TYPEOF(countries) != INTSXP || Rf_xlength(countries) != 1) Rf_error("Invalid WZ country count.");
    const auto shape = role_shape(mean, country, INTEGER(countries)[0]);
    require_matrix_length(probability, Rf_xlength(mean));
    require_matrix_length(score, Rf_xlength(mean));
    SEXP expected = PROTECT(Rf_alloc3DArray(REALSXP, shape.dates, shape.dates, shape.countries));
    SEXP meat = PROTECT(Rf_alloc3DArray(REALSXP, shape.dates, shape.dates, shape.countries));
    double* a = REAL(expected);
    double* b = REAL(meat);
    std::fill_n(a, Rf_xlength(expected), 0.0);
    std::fill_n(b, Rf_xlength(meat), 0.0);
    const double* mu = REAL(mean);
    const double* pi = REAL(probability);
    const double* residual = REAL(score);
    const int* codes = INTEGER(country);
    const std::size_t stride = static_cast<std::size_t>(shape.dates) * shape.dates;
    for (int group = 0; group < shape.groups; ++group) {
        if ((group & 8191) == 0) R_CheckUserInterrupt();
        const std::size_t block = stride * (codes[group] - 1);
        for (int col = 0; col < shape.dates; ++col) {
            const std::size_t j = group + static_cast<std::size_t>(col) * shape.groups;
            const std::size_t target = block + static_cast<std::size_t>(col) * shape.dates;
            for (int row = 0; row < shape.dates; ++row) {
                const std::size_t i = group + static_cast<std::size_t>(row) * shape.groups;
                a[target + row] += mu[i] * (static_cast<double>(col == row) - pi[j]);
                b[target + row] += residual[i] * residual[j];
            }
        }
    }
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(result, 0, expected);
    SET_VECTOR_ELT(result, 1, meat);
    UNPROTECT(3);
    return result;
}

SEXP wz_role_scores_call(SEXP probability, SEXP score, SEXP country, SEXP inverse, SEXP variance, SEXP observed) {
    const SEXP dimensions = Rf_getAttrib(inverse, R_DimSymbol);
    if (TYPEOF(dimensions) != INTSXP || Rf_xlength(dimensions) != 3) Rf_error("Invalid WZ role inverse.");
    const auto shape = role_shape(probability, country, INTEGER(dimensions)[2]);
    const std::size_t stride = static_cast<std::size_t>(shape.dates) * shape.dates;
    require_matrix_length(inverse, stride * shape.countries);
    require_matrix_length(variance, Rf_xlength(inverse));
    require_matrix_length(observed, shape.groups);
    require_matrix_length(score, Rf_xlength(probability));
    const double* lifted = REAL(inverse);
    const double* cov = REAL(variance);
    const double* pi = REAL(probability);
    const double* residual = REAL(score);
    const double* mass = REAL(observed);
    const int* codes = INTEGER(country);
    auto* first = reinterpret_cast<double*>(R_alloc(shape.dates, sizeof(double)));
    auto* second = reinterpret_cast<double*>(R_alloc(shape.dates, sizeof(double)));
    SEXP result = PROTECT(Rf_allocVector(REALSXP, Rf_xlength(probability)));
    double* output = REAL(result);
    for (int group = 0; group < shape.groups; ++group) {
        if ((group & 8191) == 0) R_CheckUserInterrupt();
        const std::size_t block = stride * (codes[group] - 1);
        double first_center = 0;
        double second_center = 0;
        for (int row = 0; row < shape.dates; ++row) {
            double a = 0;
            double v = 0;
            for (int col = 0; col < shape.dates; ++col) {
                const std::size_t j = group + static_cast<std::size_t>(col) * shape.groups;
                const std::size_t at = block + row + static_cast<std::size_t>(col) * shape.dates;
                a += lifted[at] * residual[j];
                v += 2 * cov[at] * pi[j];
            }
            v -= cov[block + row + static_cast<std::size_t>(row) * shape.dates];
            first[row] = a;
            second[row] = v;
            const double weight = pi[group + static_cast<std::size_t>(row) * shape.groups];
            first_center += weight * a;
            second_center += weight * v;
        }
        double* current = output + static_cast<std::size_t>(group) * shape.dates;
        for (int row = 0; row < shape.dates; ++row) {
            current[row] = mass[group] * pi[group + static_cast<std::size_t>(row) * shape.groups] *
                (-(first[row] - first_center) + 0.5 * (second[row] - second_center));
        }
    }
    UNPROTECT(1);
    return result;
}
