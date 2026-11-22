#define R_NO_REMAP
#include <R.h>
#include <R_ext/Lapack.h>
#include <R_ext/Utils.h>
#include <Rinternals.h>

#include "native_calls.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>

namespace {

// A compact, nearly dense Jacobian can use a dense LAPACK factor. The
// temporary row support stays small and the public coordinates stay fixed.
constexpr int max_columns = 2048;
constexpr double max_products = 40e6;

#ifdef La_INT
using LapackInt = La_INT;
#else
using LapackInt = int;
#endif

struct GammaRows {
    int width;
    const int* offsets;
    const int* columns;
    const double* values;
    const double* weight;
};

double subset_weight(std::uint32_t mask, const double* weight) {
    double sum = 0;
    for (int row = 0; mask != 0; ++row, mask >>= 1) {
        if ((mask & 1U) != 0) sum += weight[row];
    }
    return sum;
}

bool add_binary_group(std::size_t width, int size, const int* columns,
                      const std::uint32_t* mask, const double* average,
                      const double* weight, int count, double* mass, double* hessian) {
    double group_weight = 0;
    for (int row = 0; row < count; ++row) group_weight += weight[row];
    if (group_weight > std::numeric_limits<double>::max() / 4) return false;
    for (int j = 0; j < size; ++j) mass[j] = subset_weight(mask[columns[j]], weight);
    // Expand each centered binary crossproduct using its intersection and
    // first moments. Use two nonnegative terms on the diagonal.
    for (int j = 0; j < size; ++j) {
        const int column = columns[j];
        const double center = average[column];
        if (center == 1) continue;
        double* target = hessian + static_cast<std::size_t>(column) * width;
        const double complement = 1 - center;
        target[column] += mass[j] * complement * complement + (group_weight - mass[j]) * center * center;
        for (int k = 0; k < j; ++k) {
            const int other = columns[k];
            const double other_center = average[other];
            if (other_center == 1) continue;
            const double joint = subset_weight(mask[column] & mask[other], weight);
            target[other] += joint - center * mass[k] - other_center * mass[j] +
                center * other_center * group_weight;
        }
    }
    return true;
}

bool add_groups(const GammaRows& rows, const int* counts, int groups,
                bool use_binary, double* hessian) {
    const auto width = static_cast<std::size_t>(rows.width);
    auto* columns = reinterpret_cast<int*>(R_alloc(width, sizeof(int)));
    auto* position = reinterpret_cast<int*>(R_alloc(width, sizeof(int)));
    auto* total = reinterpret_cast<double*>(R_alloc(width, sizeof(double)));
    auto* difference = reinterpret_cast<double*>(R_alloc(width, sizeof(double)));
    auto* mask = reinterpret_cast<std::uint32_t*>(R_alloc(width, sizeof(std::uint32_t)));
    std::fill_n(position, width, -1);
    double products = 0;
    int begin = 0;
    for (int group = 0; group < groups; ++group) {
        if ((group & 1023) == 0) R_CheckUserInterrupt();
        const int count = counts[group];
        const int end = begin + count;
        if (count < 2) {
            begin = end;
            continue;
        }
        int size = 0;
        bool binary = use_binary && count <= 16;
        for (int i = begin; i < end; ++i) {
            for (int a = rows.offsets[i]; a < rows.offsets[i + 1]; ++a) {
                const int column = rows.columns[a];
                if (position[column] < 0) {
                    position[column] = size;
                    columns[size++] = column;
                    total[column] = 0;
                    mask[column] = 0;
                }
                total[column] += rows.values[a];
                if (binary) {
                    const auto bit = std::uint32_t{1} << (i - begin);
                    binary = rows.values[a] == 1 && (mask[column] & bit) == 0;
                    mask[column] |= bit;
                }
            }
        }
        products += static_cast<double>(count) * size * (size + 1) / 2;
        if (products > max_products) return false;
        std::sort(columns, columns + size);
        for (int j = 0; j < size; ++j) {
            position[columns[j]] = j;
            total[columns[j]] /= count;
        }
        if (!binary || !add_binary_group(width, size, columns, mask, total, rows.weight + begin,
                                         count, difference, hessian)) {
            for (int i = begin; i < end; ++i) {
                if (rows.weight[i] == 0) continue;
                for (int j = 0; j < size; ++j) difference[j] = -total[columns[j]];
                for (int a = rows.offsets[i]; a < rows.offsets[i + 1]; ++a) {
                    difference[position[rows.columns[a]]] += rows.values[a];
                }
                for (int j = 0; j < size; ++j) {
                    const double value = rows.weight[i] * difference[j];
                    auto* target = hessian + static_cast<std::size_t>(columns[j]) * width;
                    for (int k = 0; k <= j; ++k) target[columns[k]] += value * difference[k];
                }
            }
        }
        for (int j = 0; j < size; ++j) position[columns[j]] = -1;
        begin = end;
    }
    for (std::size_t j = 0; j < width; ++j) {
        for (std::size_t i = 0; i <= j; ++i) {
            const double entry = hessian[i + j * width];
            if (!std::isfinite(entry)) return false;
            hessian[j + i * width] = entry;
        }
    }
    return true;
}

}  // namespace

SEXP cancel_gamma_hessian_call(SEXP offsets, SEXP columns, SEXP values, SEXP counts,
                               SEXP weight, SEXP dimension) {
    if (TYPEOF(offsets) != INTSXP || TYPEOF(columns) != INTSXP || TYPEOF(values) != REALSXP ||
        TYPEOF(counts) != INTSXP || TYPEOF(weight) != REALSXP || TYPEOF(dimension) != INTSXP ||
        (Rf_xlength(dimension) != 1 && Rf_xlength(dimension) != 2) ||
        Rf_xlength(weight) < 1 || Rf_xlength(weight) >= INT_MAX ||
        Rf_xlength(offsets) != Rf_xlength(weight) + 1 || Rf_xlength(counts) < 1 ||
        Rf_xlength(counts) > Rf_xlength(weight) || Rf_xlength(columns) != Rf_xlength(values) ||
        Rf_xlength(values) >= INT_MAX) {
        Rf_error("Invalid Gamma Jacobian arrays.");
    }
    const int width = INTEGER(dimension)[0];
    if (width < 1 || width > max_columns) return R_NilValue;
    const int mode = Rf_xlength(dimension) == 2 ? INTEGER(dimension)[1] : 1;
    if (mode != 0 && mode != 1) Rf_error("Invalid Gamma Jacobian mode.");
    const int n = static_cast<int>(Rf_xlength(weight));
    const int groups = static_cast<int>(Rf_xlength(counts));
    const int entries = static_cast<int>(Rf_xlength(values));
    const int* sizes = INTEGER(counts);
    const GammaRows rows{width, INTEGER(offsets), INTEGER(columns), REAL(values), REAL(weight)};
    int passed = 0;
    for (int group = 0; group < groups; ++group) {
        if (sizes[group] <= 0 || sizes[group] > n - passed) Rf_error("Invalid Gamma Jacobian group size.");
        passed += sizes[group];
    }
    if (passed != n || rows.offsets[0] != 0 || rows.offsets[n] != entries) {
        Rf_error("Gamma Jacobian arrays do not align.");
    }
    for (int row = 0; row < n; ++row) {
        if (rows.offsets[row] > rows.offsets[row + 1] || rows.offsets[row + 1] > entries ||
            !std::isfinite(rows.weight[row]) || rows.weight[row] < 0) return R_NilValue;
        for (int a = rows.offsets[row]; a < rows.offsets[row + 1]; ++a) {
            if (rows.columns[a] < 0 || rows.columns[a] >= width || !std::isfinite(rows.values[a])) {
                return R_NilValue;
            }
        }
    }
    SEXP result = PROTECT(Rf_allocMatrix(REALSXP, width, width));
    auto* output = REAL(result);
    std::fill_n(output, static_cast<std::size_t>(width) * width, 0.0);
    const bool accepted = add_groups(rows, sizes, groups, mode == 1, output);
    UNPROTECT(1);
    return accepted ? result : R_NilValue;
}

SEXP cancel_gamma_dense_factor_call(SEXP hessian) {
    const SEXP dimensions = Rf_getAttrib(hessian, R_DimSymbol);
    if (TYPEOF(hessian) != REALSXP || TYPEOF(dimensions) != INTSXP ||
        Rf_xlength(dimensions) != 2 || Rf_getAttrib(hessian, R_DimNamesSymbol) != R_NilValue) {
        return R_NilValue;
    }
    const int dimension = INTEGER(dimensions)[0];
    if (dimension < 1 || dimension > max_columns || INTEGER(dimensions)[1] != dimension) {
        return R_NilValue;
    }
    const auto n = static_cast<std::size_t>(dimension);
    auto* sums = reinterpret_cast<long double*>(R_alloc(n, sizeof(long double)));
    std::fill_n(sums, n, 0.0L);
    const double* values = REAL(hessian);
    // Match rowSums(abs(H)) in each row's column order while reading H
    // contiguously. Preserve H for the subsequent residual and covariance.
    for (std::size_t col = 0; col < n; ++col) {
        const double* source = values + col * n;
        for (std::size_t row = 0; row < n; ++row) sums[row] += std::abs(source[row]);
    }
    const double smallest = std::sqrt(std::numeric_limits<double>::min());
    const double largest = std::sqrt(std::numeric_limits<double>::max() / 4);
    SEXP units = PROTECT(Rf_allocVector(REALSXP, dimension));
    double* scale = REAL(units);
    for (std::size_t row = 0; row < n; ++row) {
        scale[row] = std::sqrt(static_cast<double>(sums[row]));
        if (!(scale[row] >= smallest && scale[row] <= largest)) {
            UNPROTECT(1);
            return R_NilValue;
        }
    }
    SEXP factor = PROTECT(Rf_allocMatrix(REALSXP, dimension, dimension));
    double* triangle = REAL(factor);
    // LAPACK only reads the upper triangle. Form the scaled entries in
    // its output buffer instead of materializing the scaling outer product.
    for (std::size_t col = 0; col < n; ++col) {
        const double* source = values + col * n;
        double* output = triangle + col * n;
        const double unit = scale[col];
        for (std::size_t row = 0; row <= col; ++row) output[row] = source[row] / (scale[row] * unit);
        std::fill(output + col + 1, output + n, 0.0);
    }
    const auto order = static_cast<LapackInt>(dimension);
    LapackInt info = 0;
    F77_CALL(dpotrf)("U", &order, triangle, &order, &info FCONE);
    for (std::size_t col = 0; col < n && !info; ++col) {
        if (!std::isfinite(triangle[col * n + col])) info = -1;
    }
    if (info) {
        UNPROTECT(2);
        return R_NilValue;
    }
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(result, 0, factor);
    SET_VECTOR_ELT(result, 1, units);
    UNPROTECT(3);
    return result;
}
