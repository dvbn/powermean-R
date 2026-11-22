// Exact binary-rational fallbacks for coefficient targets and restriction tests.
#define R_NO_REMAP

#include <R.h>
#include <Rinternals.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <vector>

#include "binary_round.h"
#include "native_calls.h"

namespace {

using boost::multiprecision::cpp_rational;
using powermean::dyadic;

double slope_scale(double value, double unit, double order, bool report) {
    if (std::isnan(value)) return value;
    if (value == 0 || std::isinf(value)) return order < 0 ? -value : value;
    int value_exp = 0, unit_exp = 0, order_exp = 0;
    const double value_mant = std::frexp(value, &value_exp);
    const double unit_mant = std::frexp(unit, &unit_exp);
    const double order_mant = std::frexp(order, &order_exp);
    if (unit_mant == 0.5) {
        const double mantissa = report ? 2 * (value_mant / order_mant) :
                                         0.5 * (value_mant * order_mant);
        const int shift = report ? value_exp - unit_exp - order_exp :
                                   value_exp + unit_exp + order_exp;
        const double normal = std::ldexp(mantissa, shift);
        if (!report) return std::isfinite(normal) && normal != 0 ? normal : R_NaN;
        const double magnitude = std::abs(normal);
        if (magnitude > std::numeric_limits<double>::min() &&
            magnitude < std::numeric_limits<double>::max()) return normal;
    }

    // Round once when the final value lies at a binary64 range boundary.
    const cpp_rational scale = dyadic(unit) * dyadic(order);
    cpp_rational exact = dyadic(value);
    if (report) {
        exact /= scale;
    } else {
        exact *= scale;
    }
    try {
        const double rounded = powermean::nearest_binary64(exact);
        if (rounded != 0) return rounded;
    } catch (const std::overflow_error&) {
        // The caller reports an unavailable slope coordinate.
    }
    return R_NaN;
}

bool finite_real(SEXP values) {
    if (TYPEOF(values) != REALSXP) return false;
    auto* pointer = REAL(values);
    return std::all_of(pointer, pointer + XLENGTH(values), [](double v) {
        return std::isfinite(v);
    });
}

struct Pivot {
    int column;
    std::vector<cpp_rational> row;
};

bool rank_is(const double* matrix, int rows, int columns, int expected) {
    if (expected < 0 || expected > std::min(rows, columns)) return false;
    std::vector<Pivot> basis;
    std::vector<cpp_rational> row(columns);
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < columns; ++j) {
            row[j] = dyadic(matrix[i + static_cast<R_xlen_t>(j) * rows]);
        }
        for (const auto& pivot : basis) {
            const cpp_rational multiplier = row[pivot.column];
            if (multiplier == 0) continue;
            for (int j = pivot.column; j < columns; ++j) {
                row[j] -= multiplier * pivot.row[j];
            }
        }
        const auto found = std::find_if(row.begin(), row.end(),
                                       [](const cpp_rational& x) { return x != 0; });
        if (found == row.end()) continue;
        if (basis.size() == static_cast<std::size_t>(expected)) return false;
        const auto column = static_cast<int>(found - row.begin());
        const cpp_rational pivot = *found;
        for (int j = column; j < columns; ++j) row[j] /= pivot;
        basis.push_back({column, row});
        if (basis.size() == static_cast<std::size_t>(columns)) return true;
    }
    return basis.size() == static_cast<std::size_t>(expected);
}

double bilinear(const double* left, const double* matrix, const double* right,
                int rows, int columns) {
    cpp_rational total(0);
    std::vector<cpp_rational> a(rows);
    for (int i = 0; i < rows; ++i) a[i] = dyadic(left[i]);
    for (int j = 0; j < columns; ++j) {
        if (right[j] == 0) continue;
        cpp_rational column(0);
        for (int i = 0; i < rows; ++i) {
            const auto entry = matrix[i + static_cast<R_xlen_t>(j) * rows];
            if (entry != 0 && left[i] != 0) column += a[i] * dyadic(entry);
        }
        total += column * dyadic(right[j]);
    }
    return total.convert_to<double>();
}

double quotient_derivative(const double* values, const double* derivative) {
    const auto x = dyadic(values[0]);
    const auto y = dyadic(values[1]);
    const auto dx = dyadic(derivative[0]);
    const auto dy = dyadic(derivative[1]);
    cpp_rational exact = (dx * y - x * dy) / (y * y);
    return exact.convert_to<double>();
}

double power_derivative(const double* values, const double* derivative) {
    const auto powered = dyadic(values[2]);
    cpp_rational exact(0);
    if (derivative[0] != 0) {
        exact += powered * dyadic(values[1]) * dyadic(derivative[0]) / dyadic(values[0]);
    }
    if (derivative[1] != 0) {
        exact += powered * dyadic(values[3]) * dyadic(derivative[1]);
    }
    return exact.convert_to<double>();
}

bool inexact_tiny_product(double left, double right, double rounded) {
    constexpr double normal_min = std::numeric_limits<double>::min();
    const double magnitude = std::abs(rounded);
    if (magnitude > normal_min || left == 0 || right == 0) return false;
    if (magnitude == 0) return true;
    const cpp_rational exact = dyadic(std::abs(left)) * dyadic(std::abs(right));
    return exact < dyadic(normal_min) && exact != dyadic(magnitude);
}

}  // namespace

SEXP pm_result_exact_rank_is(SEXP matrix, SEXP expected) {
    const auto shape = Rf_getAttrib(matrix, R_DimSymbol);
    if (TYPEOF(shape) != INTSXP || XLENGTH(shape) != 2 || !finite_real(matrix) ||
        TYPEOF(expected) != INTSXP || XLENGTH(expected) != 1) {
        Rf_error("Invalid exact restriction rank arguments.");
    }
    bool answer = false;
    try {
        answer = rank_is(REAL(matrix), INTEGER(shape)[0], INTEGER(shape)[1], INTEGER(expected)[0]);
    } catch (...) {
        Rf_error("Exact restriction rank could not be formed.");
    }
    return Rf_ScalarLogical(answer);
}

SEXP pm_result_exact_bilinear(SEXP left, SEXP matrix, SEXP right) {
    if (!finite_real(left) || !finite_real(right)) {
        Rf_error("Invalid exact coefficient contraction operands.");
    }
    double answer = 0;
    if (matrix == R_NilValue) {
        const auto values = XLENGTH(left);
        if (XLENGTH(right) != 2 ||
            (values == 2 && REAL(left)[1] == 0) ||
            (values == 4 && REAL(left)[0] == 0) ||
            (values != 2 && values != 4)) {
            Rf_error("Invalid exact derivative operands.");
        }
        try {
            answer = values == 2 ? quotient_derivative(REAL(left), REAL(right)) :
                                   power_derivative(REAL(left), REAL(right));
        } catch (...) {
            Rf_error("Exact coefficient derivative could not be formed.");
        }
    } else {
        const auto shape = Rf_getAttrib(matrix, R_DimSymbol);
        if (TYPEOF(shape) != INTSXP || XLENGTH(shape) != 2 || !finite_real(matrix) ||
            XLENGTH(left) != INTEGER(shape)[0] || XLENGTH(right) != INTEGER(shape)[1]) {
            Rf_error("Invalid exact coefficient contraction shape.");
        }
        try {
            answer = bilinear(REAL(left), REAL(matrix), REAL(right), INTEGER(shape)[0], INTEGER(shape)[1]);
        } catch (...) {
            Rf_error("Exact coefficient contraction could not be formed.");
        }
    }
    return Rf_ScalarReal(answer);
}

SEXP pm_exact_dot_sign(SEXP left, SEXP right) {
    if (!finite_real(left) || !finite_real(right) || XLENGTH(left) != XLENGTH(right)) {
        Rf_error("Invalid exact dot product operands.");
    }
    int answer = 0;
    try {
        cpp_rational total(0);
        const double* a = REAL(left);
        const double* b = REAL(right);
        for (R_xlen_t i = 0; i < XLENGTH(left); ++i) {
            if (a[i] != 0 && b[i] != 0) total += dyadic(a[i]) * dyadic(b[i]);
        }
        answer = (total > 0) - (total < 0);
    } catch (...) {
        Rf_error("Exact dot product could not be formed.");
    }
    return Rf_ScalarInteger(answer);
}

SEXP pm_prediction_exact_shift(SEXP index, SEXP changed, SEXP original, SEXP beta, SEXP multiplier) {
    const auto shape = Rf_getAttrib(changed, R_DimSymbol);
    const auto old_shape = Rf_getAttrib(original, R_DimSymbol);
    const auto n = XLENGTH(index);
    const auto p = XLENGTH(beta);
    if (!finite_real(index) || !finite_real(beta) || !finite_real(multiplier) ||
        XLENGTH(multiplier) != 1 || !finite_real(changed) || !finite_real(original) ||
        TYPEOF(shape) != INTSXP || XLENGTH(shape) != 2 ||
        TYPEOF(old_shape) != INTSXP || XLENGTH(old_shape) != 2 ||
        INTEGER(shape)[0] != n || INTEGER(shape)[1] != p ||
        INTEGER(old_shape)[0] != n || INTEGER(old_shape)[1] != p) {
        Rf_error("Invalid exact prediction operands.");
    }
    SEXP result = PROTECT(Rf_allocVector(REALSXP, n));
    try {
        const auto scale = dyadic(REAL(multiplier)[0]);
        std::vector<cpp_rational> weights;
        weights.reserve(p);
        for (R_xlen_t j = 0; j < p; ++j) weights.push_back(scale * dyadic(REAL(beta)[j]));
        for (R_xlen_t i = 0; i < n; ++i) {
            cpp_rational found = dyadic(REAL(index)[i]);
            for (R_xlen_t j = 0; j < p; ++j) {
                const auto k = i + j * n;
                const double old = REAL(original)[k];
                const double next = REAL(changed)[k];
                if (old != next && weights[j] != 0) found += (dyadic(next) - dyadic(old)) * weights[j];
            }
            try {
                REAL(result)[i] = powermean::nearest_binary64(found);
            } catch (const std::overflow_error&) {
                REAL(result)[i] = R_NaN;
            }
        }
    } catch (...) {
        Rf_error("Exact prediction index could not be formed.");
    }
    UNPROTECT(1);
    return result;
}

SEXP pm_index_within_forcing(SEXP operands) {
    if (!finite_real(operands) || XLENGTH(operands) != 3 ||
        std::any_of(REAL(operands), REAL(operands) + 3, [](double x) { return x < 0; })) {
        Rf_error("Index forcing requires three finite nonnegative operands.");
    }
    const double* terms = REAL(operands);
    if (terms[0] == 0) return Rf_ScalarLogical(TRUE);
    if (terms[1] == 0 || terms[2] == 0) return Rf_ScalarLogical(FALSE);
    bool answer = false;
    try {
        answer = dyadic(terms[0]) <= dyadic(terms[1]) * dyadic(terms[2]);
    } catch (...) {
        Rf_error("The index forcing comparison could not be formed.");
    }
    return Rf_ScalarLogical(answer);
}

// Recovery uses unit row weights. Only a potentially inexact subnormal
// term needs an exact product check. Exact subnormal terms remain valid.
SEXP pm_index_product_range(SEXP left, SEXP right, SEXP rounded) {
    const auto n = XLENGTH(rounded);
    if (TYPEOF(left) != REALSXP || TYPEOF(right) != REALSXP || TYPEOF(rounded) != REALSXP ||
        XLENGTH(left) != n || XLENGTH(right) != n) {
        Rf_error("Index product checks require three aligned numeric vectors.");
    }
    const double* a = REAL(left);
    const double* b = REAL(right);
    const double* terms = REAL(rounded);
    int status = 0;
    try {
        for (R_xlen_t i = 0; i < n; ++i) {
            const double term = terms[i];
            if (!std::isfinite(term)) {
                status = 1;
                break;
            }
            if (status == 0 && inexact_tiny_product(a[i], b[i], term)) status = 2;
        }
    } catch (...) {
        Rf_error("The subnormal index product could not be checked.");
    }
    return Rf_ScalarInteger(status);
}

namespace {

// Product range takes precedence over range in a later row-order reduction.
SEXP linear_slope_equations(SEXP design, SEXP residual, bool retain) {
    const SEXP shape = Rf_getAttrib(design, R_DimSymbol);
    if (TYPEOF(design) != REALSXP || TYPEOF(shape) != INTSXP || XLENGTH(shape) != 2 ||
        TYPEOF(residual) != REALSXP || XLENGTH(residual) != INTEGER(shape)[0]) {
        Rf_error("Linear equation checks require an aligned design and residual.");
    }
    const auto rows = XLENGTH(residual);
    const int columns = INTEGER(shape)[1];
    const double* x = REAL(design);
    const double* e = REAL(residual);
    SEXP values = PROTECT(retain ? Rf_allocVector(REALSXP, columns) : R_NilValue);
    bool lost = false;
    bool prefix_overflow = false;
    int status = 0;
    try {
        for (int j = 0; j < columns && status == 0; ++j) {
            double sum = 0;
            for (R_xlen_t i = 0; i < rows; ++i, ++x) {
                // Use the rounded product as the reduction input.
                const volatile double product = *x * e[i];
                if (!std::isfinite(product)) {
                    status = 1;
                    break;
                }
                if (!lost) lost = inexact_tiny_product(*x, e[i], product);
                sum += product;
            }
            if (retain) REAL(values)[j] = sum;
            prefix_overflow = prefix_overflow || !std::isfinite(sum);
        }
    } catch (...) {
        Rf_error("The subnormal linear equation product could not be checked.");
    }
    if (status == 0) status = lost ? 2 : prefix_overflow ? 3 : 0;
    if (!retain || status != 0) values = Rf_ScalarInteger(status);
    UNPROTECT(1);
    return values;
}

}  // namespace

SEXP pm_linear_slope_range(SEXP design, SEXP residual) {
    return linear_slope_equations(design, residual, false);
}

SEXP pm_linear_slope_sums(SEXP design, SEXP residual) {
    return linear_slope_equations(design, residual, true);
}

// The R caller screens normally represented quotients. Reconstruct tiny
// means in exact binary arithmetic to detect inexact subnormal division.
SEXP pm_linear_mean_loss(SEXP sums, SEXP totals, SEXP means) {
    const auto levels = XLENGTH(sums);
    if (!finite_real(sums) || !finite_real(totals) || !finite_real(means) ||
        XLENGTH(totals) != levels || XLENGTH(means) != levels) {
        Rf_error("Group equation checks require aligned finite numeric vectors.");
    }
    const double* a = REAL(sums);
    const double* b = REAL(totals);
    const double* value = REAL(means);
    constexpr double normal_min = std::numeric_limits<double>::min();
    try {
        for (R_xlen_t i = 0; i < levels; ++i) {
            if (a[i] == 0 || std::abs(value[i]) > normal_min) continue;
            const cpp_rational exact = dyadic(std::abs(a[i])) / dyadic(b[i]);
            if (exact < dyadic(normal_min) && exact != dyadic(std::abs(value[i]))) {
                return Rf_ScalarLogical(true);
            }
        }
    } catch (...) {
        Rf_error("The subnormal level mean could not be checked.");
    }
    return Rf_ScalarLogical(false);
}

SEXP pm_result_covariance_scale(SEXP matrix, SEXP units, SEXP order) {
    const SEXP shape = Rf_getAttrib(matrix, R_DimSymbol);
    if (TYPEOF(matrix) != REALSXP || TYPEOF(shape) != INTSXP || XLENGTH(shape) != 2 ||
        INTEGER(shape)[0] != INTEGER(shape)[1] || TYPEOF(units) != REALSXP ||
        XLENGTH(units) != INTEGER(shape)[0] || TYPEOF(order) != REALSXP || XLENGTH(order) != 1 ||
        !std::isfinite(REAL(order)[0]) || REAL(order)[0] == 0) {
        Rf_error("Invalid reported covariance scale.");
    }
    const int size = INTEGER(shape)[0];
    auto* fractions = reinterpret_cast<double*>(R_alloc(size, sizeof(double)));
    auto* exponents = reinterpret_cast<int*>(R_alloc(size, sizeof(int)));
    for (int j = 0; j < size; ++j) {
        const double unit = REAL(units)[j];
        if (!std::isfinite(unit) || unit <= 0) Rf_error("Invalid reported slope unit.");
        fractions[j] = std::frexp(unit, &exponents[j]);
    }
    int order_exponent = 0;
    const double order_fraction = std::frexp(REAL(order)[0], &order_exponent);
    const double order_square = order_fraction * order_fraction;
    SEXP result = PROTECT(Rf_duplicate(matrix));
    for (int j = 0; j < size; ++j) {
        for (int i = 0; i < size; ++i) {
            const R_xlen_t index = i + static_cast<R_xlen_t>(j) * size;
            const double value = REAL(matrix)[index];
            // Undefined covariance entries are already indicated by NA.
            if (std::isnan(value) || value == 0) continue;
            int exponent = 0;
            const double fraction = std::frexp(value, &exponent);
            const double quotient = (fraction / (fractions[i] * fractions[j])) / order_square;
            const double reported = std::ldexp(quotient, exponent - exponents[i] -
                                                           exponents[j] - 2 * order_exponent);
            if (!std::isfinite(reported) || reported == 0) {
                UNPROTECT(1);
                return R_NilValue;
            }
            REAL(result)[index] = reported;
        }
    }
    UNPROTECT(1);
    return result;
}

SEXP pm_result_slope_scale(SEXP values, SEXP units, SEXP order, SEXP report) {
    if (TYPEOF(values) != REALSXP || !finite_real(units) || XLENGTH(values) != XLENGTH(units) ||
        !finite_real(order) || XLENGTH(order) != 1 || REAL(order)[0] == 0 ||
        TYPEOF(report) != LGLSXP || XLENGTH(report) != 1 || LOGICAL(report)[0] == NA_LOGICAL) {
        Rf_error("Invalid reported slope scale.");
    }
    if (std::any_of(REAL(units), REAL(units) + XLENGTH(units), [](double unit) { return unit <= 0; })) {
        Rf_error("Invalid reported slope unit.");
    }
    SEXP scaled = PROTECT(Rf_duplicate(values));
    try {
        for (R_xlen_t i = 0; i < XLENGTH(values); ++i) {
            REAL(scaled)[i] = slope_scale(REAL(values)[i], REAL(units)[i], REAL(order)[0], LOGICAL(report)[0]);
        }
    } catch (...) {
        Rf_error("Reported slope scale could not be formed.");
    }
    UNPROTECT(1);
    return scaled;
}
