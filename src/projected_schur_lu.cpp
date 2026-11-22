#define R_NO_REMAP
#include <R.h>
#include <R_ext/Lapack.h>
#include <Rinternals.h>

#include "native_calls.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize("fp-contract=off")
#endif

namespace {

constexpr std::size_t workspace_bytes = 128 * 1024 * 1024;
constexpr int max_dense_width = 4096;  // floor(sqrt(workspace_bytes / sizeof(double)))

#ifdef La_INT
using LapackInt = La_INT;
#else
using LapackInt = int;
#endif

struct CscView {
    const int* column;
    const int* row;
    const double* value;
    int entries;
};

bool csc_view(SEXP supplied, int rows, int columns, CscView& view, bool triangle = false) {
    if (TYPEOF(supplied) != VECSXP || XLENGTH(supplied) != 3) return false;
    const SEXP column = VECTOR_ELT(supplied, 0);
    const SEXP row = VECTOR_ELT(supplied, 1);
    const SEXP value = VECTOR_ELT(supplied, 2);
    if (TYPEOF(column) != INTSXP || TYPEOF(row) != INTSXP || TYPEOF(value) != REALSXP ||
        XLENGTH(column) != static_cast<R_xlen_t>(columns) + 1 ||
        XLENGTH(row) != XLENGTH(value) || XLENGTH(row) > std::numeric_limits<int>::max()) {
        return false;
    }
    view = CscView{INTEGER(column), INTEGER(row), REAL(value), static_cast<int>(XLENGTH(row))};
    if (view.column[0] != 0 || view.column[columns] != view.entries) return false;
    bool above = false;
    bool below = false;
    for (int c = 0; c < columns; ++c) {
        if (view.column[c + 1] < view.column[c] || view.column[c + 1] > view.entries) return false;
        int previous = -1;
        for (int k = view.column[c]; k < view.column[c + 1]; ++k) {
            const int i = view.row[k];
            if (i <= previous || i >= rows || !std::isfinite(view.value[k])) return false;
            previous = i;
            if (triangle) {
                above = above || i < c;
                below = below || i > c;
                if (above && below) return false;
            }
        }
    }
    return true;
}

bool budget_fits(int rows, int width, int entries) {
    // Width bounds the matrix term. The other terms are checked before
    // multiplying, including on hosts using 32-bit allocation sizes.
    std::size_t left = workspace_bytes - sizeof(double) *
        static_cast<std::size_t>(width) * static_cast<std::size_t>(width);
    const auto charge = [&left](std::size_t count, std::size_t size) {
        if (count > left / size) return false;
        left -= count * size;
        return true;
    };
    return charge(static_cast<std::size_t>(rows) + 1U, 2U * sizeof(std::size_t)) &&
        charge(static_cast<std::size_t>(entries), sizeof(int) + sizeof(double)) &&
        charge(static_cast<std::size_t>(width), sizeof(LapackInt) + sizeof(int));
}

}  // namespace

SEXP projected_schur_lu_call(SEXP coupling, SEXP diagonal, SEXP remainder) {
    if (TYPEOF(diagonal) != REALSXP || XLENGTH(diagonal) < 1 ||
        XLENGTH(diagonal) > std::numeric_limits<int>::max() ||
        TYPEOF(remainder) != VECSXP || XLENGTH(remainder) != 3) {
        return R_NilValue;
    }
    const SEXP pointer = VECTOR_ELT(remainder, 0);
    if (TYPEOF(pointer) != INTSXP || XLENGTH(pointer) < 2 ||
        XLENGTH(pointer) > max_dense_width + 1) return R_NilValue;
    const int width = static_cast<int>(XLENGTH(pointer) - 1);
    const int m = static_cast<int>(XLENGTH(diagonal));
    const double* d = REAL(diagonal);
    CscView b{}, g{};
    if (!csc_view(coupling, m, width, b) || !csc_view(remainder, width, width, g, true) ||
        !budget_fits(m, width, b.entries)) return R_NilValue;
    for (int i = 0; i < m; ++i) {
        if (!std::isfinite(d[i]) || d[i] <= 0) return R_NilValue;
    }

    // Convert the coupling once to ascending remainder coordinates within
    // each eliminated row. No full sparse Schur product is materialized.
    const auto rows = static_cast<std::size_t>(m);
    auto* offset = reinterpret_cast<std::size_t*>(R_alloc(rows + 1U, sizeof(std::size_t)));
    auto* next = reinterpret_cast<std::size_t*>(R_alloc(rows + 1U, sizeof(std::size_t)));
    auto* index = reinterpret_cast<int*>(R_alloc(b.entries, sizeof(int)));
    auto* values = reinterpret_cast<double*>(R_alloc(b.entries, sizeof(double)));
    std::fill_n(offset, rows + 1U, std::size_t{0});
    for (int k = 0; k < b.entries; ++k) ++offset[static_cast<std::size_t>(b.row[k]) + 1U];
    for (std::size_t i = 0; i < rows; ++i) offset[i + 1U] += offset[i];
    std::copy_n(offset, rows, next);
    for (int c = 0; c < width; ++c) {
        for (int k = b.column[c]; k < b.column[c + 1]; ++k) {
            const auto slot = next[b.row[k]]++;
            index[slot] = c;
            values[slot] = b.value[k];
        }
    }

    SEXP factor = PROTECT(Rf_allocMatrix(REALSXP, width, width));
    double* lu = REAL(factor);
    const auto dimension = static_cast<std::size_t>(width);
    const auto cells = dimension * dimension;
    std::fill_n(lu, cells, 0.0);
    for (std::size_t i = 0; i < rows; ++i) {
        // CSC Schur product visits its contributing eliminated coordinates in
        // this order for each output. Divide the left term before multiplying.
        for (auto q = offset[i]; q != offset[i + 1U]; ++q) {
            const double left = values[q] / d[i];
            if (!std::isfinite(left)) {
                UNPROTECT(1);
                return R_NilValue;
            }
            for (auto k = offset[i]; k != offset[i + 1U]; ++k) {
                const auto target = static_cast<std::size_t>(index[k]) * dimension + index[q];
                lu[target] += left * values[k];
            }
        }
    }
    for (std::size_t i = 0; i < cells; ++i) lu[i] = -lu[i];
    for (int c = 0; c < width; ++c) {
        for (int k = g.column[c]; k < g.column[c + 1]; ++k) {
            const int i = g.row[k];
            const auto position = static_cast<std::size_t>(c) * dimension + i;
            lu[position] = g.value[k] - (-lu[position]);
            if (i != c) {
                const auto mirror = static_cast<std::size_t>(i) * dimension + c;
                lu[mirror] = g.value[k] - (-lu[mirror]);
            }
        }
    }
    for (std::size_t i = 0; i < cells; ++i) {
        if (!std::isfinite(lu[i])) {
            UNPROTECT(1);
            return R_NilValue;
        }
        if (lu[i] == 0) lu[i] = 0;  // canonical omitted zero in the sparse counterpart
    }

    const LapackInt order = width;
    auto* pivots = reinterpret_cast<LapackInt*>(R_alloc(dimension, sizeof(LapackInt)));
    LapackInt info = 0;
    F77_CALL(dgetrf)(&order, &order, lu, &order, pivots, &info);
    if (info != 0 || !std::all_of(lu, lu + cells, [](double value) { return std::isfinite(value); })) {
        UNPROTECT(1);
        return R_NilValue;
    }
    SEXP pivot = PROTECT(Rf_allocVector(INTSXP, width));
    std::copy_n(pivots, dimension, INTEGER(pivot));
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(result, 0, factor);
    SET_VECTOR_ELT(result, 1, pivot);
    UNPROTECT(3);
    return result;
}

SEXP projected_schur_lu_solve_call(SEXP factor, SEXP pivot, SEXP right, SEXP transpose) {
    return checked_dense_lu_solve(factor, pivot, right, transpose, max_dense_width);
}
