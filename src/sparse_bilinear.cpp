#define R_NO_REMAP
#include <R.h>
#include <R_ext/Utils.h>
#include <Rinternals.h>

#include "native_calls.h"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace {

struct SparseRows {
    const int* offsets;
    const int* columns;
    const double* values;
    int n;
};

SparseRows sparse_rows(SEXP storage, int width) {
    if (TYPEOF(storage) != VECSXP || Rf_xlength(storage) != 3) Rf_error("Invalid sparse row storage.");
    const SEXP offsets = VECTOR_ELT(storage, 0);
    const SEXP indices = VECTOR_ELT(storage, 1);
    const SEXP values = VECTOR_ELT(storage, 2);
    if (TYPEOF(offsets) != INTSXP || TYPEOF(indices) != INTSXP || TYPEOF(values) != REALSXP ||
        Rf_xlength(offsets) < 2 || Rf_xlength(offsets) > INT_MAX ||
        Rf_xlength(indices) != Rf_xlength(values) || Rf_xlength(values) >= INT_MAX) {
        Rf_error("Sparse bilinear forms need compressed numeric rows.");
    }
    const int n = static_cast<int>(Rf_xlength(offsets) - 1);
    const int entries = static_cast<int>(Rf_xlength(indices));
    SparseRows result{INTEGER(offsets), INTEGER(indices), REAL(values), n};
    if (result.offsets[0] != 0 || result.offsets[n] != entries) {
        Rf_error("Invalid offsets for sparse bilinear forms.");
    }
    for (int row = 0; row < n; ++row) {
        if (result.offsets[row] < 0 || result.offsets[row + 1] < result.offsets[row] ||
            result.offsets[row + 1] > entries) {
            Rf_error("Invalid sparse bilinear row.");
        }
        int previous = -1;
        for (int p = result.offsets[row]; p < result.offsets[row + 1]; ++p) {
            if (result.columns[p] <= previous || result.columns[p] >= width || !std::isfinite(result.values[p])) {
                Rf_error("Invalid sparse bilinear entry.");
            }
            previous = result.columns[p];
        }
    }
    return result;
}

int centered_row(const SparseRows& design, int row, const SparseRows& mean, int group,
                 int* columns, double* values) {
    int lhs = design.offsets[row];
    int rhs = mean.offsets[group];
    const int end_lhs = design.offsets[row + 1];
    const int end_rhs = mean.offsets[group + 1];
    int count = 0;
    while (lhs < end_lhs || rhs < end_rhs) {
        const int left = lhs < end_lhs ? design.columns[lhs] : INT_MAX;
        const int right = rhs < end_rhs ? mean.columns[rhs] : INT_MAX;
        const int column = std::min(left, right);
        const double a = left == column ? design.values[lhs++] : 0;
        const double b = right == column ? mean.values[rhs++] : 0;
        const double value = a - b;
        if (value == 0) continue;
        columns[count] = column;
        values[count] = value;
        ++count;
    }
    return count;
}

double metric_row(const double* matrix, int dimension, int column, const SparseRows& other, int row) {
    double total = 0;
    for (int j = other.offsets[row]; j < other.offsets[row + 1]; ++j) {
        total += matrix[column + static_cast<std::size_t>(dimension) * other.columns[j]] * other.values[j];
    }
    return total;
}

}  // namespace

SEXP bbp_mixed_products_call(SEXP z, SEXP mean, SEXP scores, SEXP curvature,
                            SEXP bread, SEXP covariance, SEXP group, SEXP sample) {
    const SEXP shape = Rf_getAttrib(bread, R_DimSymbol);
    const SEXP shape_cov = Rf_getAttrib(covariance, R_DimSymbol);
    if (TYPEOF(bread) != REALSXP || TYPEOF(shape) != INTSXP || Rf_xlength(shape) != 2 ||
        INTEGER(shape)[0] < 1 || INTEGER(shape)[0] != INTEGER(shape)[1] ||
        TYPEOF(covariance) != REALSXP || TYPEOF(shape_cov) != INTSXP || Rf_xlength(shape_cov) != 2 ||
        INTEGER(shape_cov)[0] != INTEGER(shape)[0] || INTEGER(shape_cov)[1] != INTEGER(shape)[0] ||
        TYPEOF(group) != INTSXP || TYPEOF(sample) != INTSXP) {
        Rf_error("Mixed BBP products need aligned sparse rows and square numerical metrics.");
    }
    const int d = INTEGER(shape)[0];
    const auto design = sparse_rows(z, d);
    const auto group_mean = sparse_rows(mean, d);
    const auto score_rows = sparse_rows(scores, d);
    const bool coupled = curvature != R_NilValue;
    const auto curvature_rows = coupled ? sparse_rows(curvature, d) : SparseRows{};
    if (Rf_xlength(group) != design.n || Rf_xlength(sample) != design.n ||
        (coupled && curvature_rows.n != group_mean.n)) {
        Rf_error("Mixed BBP product indices do not align.");
    }
    const int* groups = INTEGER(group);
    const int* samples = INTEGER(sample);
    const double* inverse = REAL(bread);
    const double* variance = REAL(covariance);
    auto* columns = reinterpret_cast<int*>(R_alloc(d, sizeof(int)));
    auto* stamp = reinterpret_cast<int*>(R_alloc(d, sizeof(int)));
    auto* values = reinterpret_cast<double*>(R_alloc(d, sizeof(double)));
    auto* impact = reinterpret_cast<double*>(R_alloc(d, sizeof(double)));
    auto* shared_action = reinterpret_cast<double*>(R_alloc(d, sizeof(double)));
    std::fill_n(stamp, d, 0);
    SEXP displacement = PROTECT(Rf_allocVector(REALSXP, design.n));
    SEXP own = PROTECT(Rf_allocVector(REALSXP, design.n));
    SEXP shared = PROTECT(Rf_allocVector(REALSXP, design.n));
    double* a = REAL(displacement);
    double* b = REAL(own);
    double* c = REAL(shared);
    int previous_group = -1;
    int previous_sample = -1;
    int epoch = 0;
    for (int row = 0; row < design.n; ++row) {
        if ((row & 16383) == 0) R_CheckUserInterrupt();
        if (groups[row] < 1 || groups[row] > group_mean.n || samples[row] < 1 || samples[row] > score_rows.n) {
            Rf_error("Invalid mixed BBP product row index.");
        }
        const int g = groups[row] - 1;
        const int s = samples[row] - 1;
        if (g != previous_group || s != previous_sample) {
            previous_group = g;
            previous_sample = s;
            ++epoch;
        }
        const int length = centered_row(design, row, group_mean, g, columns, values);
        double step = 0;
        double row_variance = 0;
        double between = 0;
        for (int i = 0; i < length; ++i) {
            const int col = columns[i];
            if (stamp[col] != epoch) {
                impact[col] = metric_row(inverse, d, col, score_rows, s);
                if (coupled) shared_action[col] = metric_row(variance, d, col, curvature_rows, g);
                stamp[col] = epoch;
            }
            double diagonal = 0;
            for (int j = 0; j < length; ++j) {
                diagonal += variance[col + static_cast<std::size_t>(d) * columns[j]] * values[j];
            }
            step += values[i] * impact[col];
            row_variance += values[i] * diagonal;
            if (coupled) between += values[i] * shared_action[col];
        }
        if (!std::isfinite(step) || !std::isfinite(row_variance) || !std::isfinite(between)) {
            UNPROTECT(3);
            return R_NilValue;
        }
        a[row] = step;
        b[row] = row_variance;
        c[row] = between;
    }
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 3));
    SET_VECTOR_ELT(result, 0, displacement);
    SET_VECTOR_ELT(result, 1, own);
    SET_VECTOR_ELT(result, 2, shared);
    UNPROTECT(4);
    return result;
}
