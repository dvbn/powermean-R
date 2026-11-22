// Exact sign checks for the original rows of the cancelled Gamma quotient.
#define R_NO_REMAP

#include <boost/multiprecision/cpp_int.hpp>

#include <R.h>
#include <Rinternals.h>

#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <limits>
#include <map>
#include <new>
#include <utility>
#include <vector>

#include "native_calls.h"

namespace {

using boost::multiprecision::cpp_int;
using boost::multiprecision::cpp_rational;
using RationalRow = std::vector<cpp_rational>;

constexpr int digits = std::numeric_limits<double>::digits;
static_assert(std::numeric_limits<double>::is_iec559 && digits == 53,
              "Cancelled Gamma support requires IEEE binary64.");

struct Dyadic {
    std::int64_t mantissa;
    int exponent;
};

Dyadic dyadic(double value) {
    int exponent = 0;
    const double fraction = std::frexp(value, &exponent);
    return {static_cast<std::int64_t>(std::ldexp(fraction, digits)), exponent - digits};
}

cpp_rational rational(double value) {
    if (value == 0) return 0;
    const auto part = dyadic(value);
    cpp_rational exact(part.mantissa);
    if (part.exponent >= 0) exact *= cpp_int(1) << part.exponent;
    else exact /= cpp_int(1) << -part.exponent;
    return exact;
}

struct SupportRows {
    int n;
    int p;
    const int* start;
    const int* row;
    const double* value;
    const int* group;
    const int* positive;
    std::vector<int> count;

    int signs(const RationalRow& direction) const {
        // Clear coefficient denominators once, and choose a shared dyadic
        // exponent for the nonzero raw matrix entries actually used.
        cpp_int common = 1;
        int exponent = 0;
        for (int j = 0; j < p; ++j) {
            if (direction[j] == 0) continue;
            common = boost::multiprecision::lcm(common, denominator(direction[j]));
            for (int a = start[j]; a < start[j + 1]; ++a) {
                if (value[a] != 0) exponent = std::min(exponent, dyadic(value[a]).exponent);
            }
        }
        std::vector<cpp_int> action(n);
        for (int j = 0; j < p; ++j) {
            if (direction[j] == 0) continue;
            const cpp_int coefficient = numerator(direction[j]) * (common / denominator(direction[j]));
            for (int a = start[j]; a < start[j + 1]; ++a) {
                if (value[a] == 0) continue;
                const auto part = dyadic(value[a]);
                cpp_int contribution = coefficient * part.mantissa;
                contribution <<= part.exponent - exponent;
                action[row[a]] += contribution;
            }
        }
        std::vector<cpp_int> sum(count.size());
        for (int i = 0; i < n; ++i) sum[group[i]] += action[i];
        int signs = 0;
        for (int i = 0; i < n && signs != 3; ++i) {
            if (!positive[i]) continue;
            const cpp_int centered = count[group[i]] * action[i] - sum[group[i]];
            if (centered > 0) signs |= 1;
            if (centered < 0) signs |= 2;
        }
        return signs;
    }

    std::map<int, RationalRow> sums(const int* face, int length) const {
        std::map<int, RationalRow> result;
        for (int i = 0; i < length; ++i) result.try_emplace(group[face[i]], p);
        for (int j = 0; j < p; ++j) {
            for (int a = start[j]; a < start[j + 1]; ++a) {
                auto where = result.find(group[row[a]]);
                if (where != result.end()) where->second[j] += rational(value[a]);
            }
        }
        return result;
    }

    RationalRow centered_row(int index, const RationalRow& sum) const {
        RationalRow result(p);
        const int size = count[group[index]];
        for (int j = 0; j < p; ++j) {
            result[j] = -sum[j];
            const auto end = row + start[j + 1];
            const auto at = std::lower_bound(row + start[j], end, index);
            if (at != end && *at == index) result[j] += size * rational(value[at - row]);
        }
        return result;
    }
};

bool add_equality(std::vector<RationalRow>& basis, RationalRow& rhs,
                  RationalRow equation, cpp_rational right) {
    const auto p = basis.size();
    for (std::size_t j = 0; j < p; ++j) {
        const cpp_rational pivot = equation[j];
        if (pivot == 0) continue;
        if (basis[j].empty()) {
            for (std::size_t k = j; k < p; ++k) equation[k] /= pivot;
            basis[j] = std::move(equation);
            rhs[j] = right / pivot;
            return true;
        }
        for (std::size_t k = j; k < p; ++k) equation[k] -= pivot * basis[j][k];
        right -= pivot * rhs[j];
    }
    return right == 0;
}

bool repair_face(const SupportRows& rows, RationalRow& direction, const int* face, int length) {
    constexpr int max_columns = 64;
    constexpr double max_repair_work = 2000000;
    const int p = rows.p;
    // Bound the number and size of retained pivots and the elimination work.
    // Declining the repair leaves the proposal unverified.
    if (!length || p > max_columns ||
        static_cast<double>(p) * p * (p + static_cast<double>(length)) > max_repair_work) {
        return false;
    }
    const auto sums = rows.sums(face, length);
    const int anchor = face[length - 1];
    auto selected = rows.centered_row(anchor, sums.at(rows.group[anchor]));
    cpp_rational anchor_value = 0;
    for (int j = 0; j < p; ++j) anchor_value += selected[j] * direction[j];
    if (anchor_value <= 0) return false;

    std::vector<RationalRow> basis(p);
    RationalRow rhs(p);
    for (int i = 0; i < length - 1; ++i) {
        const int index = face[i];
        if (!add_equality(basis, rhs, rows.centered_row(index, sums.at(rows.group[index])), 0)) {
            return false;
        }
    }
    if (!add_equality(basis, rhs, std::move(selected), 1)) return false;
    for (auto& value : direction) value /= anchor_value;
    for (int j = p; j-- > 0;) {
        if (basis[j].empty()) continue;
        direction[j] = rhs[j];
        for (int k = j + 1; k < p; ++k) direction[j] -= basis[j][k] * direction[k];
    }
    return rows.signs(direction) == 1;
}

} // namespace

SEXP cancel_gamma_support_call(SEXP design, SEXP group, SEXP positive, SEXP units,
                               SEXP direction, SEXP face) {
    if (TYPEOF(design) != VECSXP || XLENGTH(design) != 3 ||
        TYPEOF(group) != INTSXP || TYPEOF(positive) != LGLSXP || TYPEOF(face) != INTSXP ||
        TYPEOF(units) != REALSXP || TYPEOF(direction) != REALSXP ||
        XLENGTH(positive) != XLENGTH(group) || XLENGTH(direction) != XLENGTH(units) ||
        XLENGTH(group) > INT_MAX || XLENGTH(face) > INT_MAX || XLENGTH(units) >= INT_MAX) {
        Rf_error("Invalid cancelled Gamma support inputs.");
    }
    const SEXP offsets = VECTOR_ELT(design, 0);
    const SEXP index = VECTOR_ELT(design, 1);
    const SEXP value = VECTOR_ELT(design, 2);
    const int n = static_cast<int>(XLENGTH(group));
    const int p = static_cast<int>(XLENGTH(units));
    if (TYPEOF(offsets) != INTSXP || XLENGTH(offsets) != p + 1 ||
        TYPEOF(index) != INTSXP || TYPEOF(value) != REALSXP || XLENGTH(index) != XLENGTH(value)) {
        Rf_error("Invalid cancelled Gamma support matrix.");
    }
    const auto* starts = INTEGER(offsets);
    const auto* observed = LOGICAL(positive);
    const auto* groups = INTEGER(group);
    const auto* constraint = INTEGER(face);
    const auto* slope_unit = REAL(units);
    const auto* proposed = REAL(direction);
    const auto* positions = INTEGER(index);
    const auto* entries = REAL(value);
    if (starts[0] != 0 || starts[p] != XLENGTH(index)) {
        Rf_error("Invalid cancelled Gamma support offsets.");
    }
    int n_groups = 0;
    for (int i = 0; i < n; ++i) {
        if (groups[i] < 0 || groups[i] >= n || observed[i] == NA_LOGICAL) {
            Rf_error("Invalid cancelled Gamma support row.");
        }
        n_groups = std::max(n_groups, groups[i] + 1);
    }
    for (int j = 0; j < p; ++j) {
        if (starts[j + 1] < starts[j] || starts[j + 1] > starts[p] ||
            !std::isfinite(slope_unit[j]) || slope_unit[j] <= 0 || !std::isfinite(proposed[j])) {
            Rf_error("Invalid cancelled Gamma support column.");
        }
        int prior = -1;
        for (int a = starts[j]; a < starts[j + 1]; ++a) {
            if (positions[a] <= prior || positions[a] >= n || !std::isfinite(entries[a])) {
                Rf_error("Invalid cancelled Gamma support coefficient.");
            }
            prior = positions[a];
        }
    }
    const int face_length = static_cast<int>(XLENGTH(face));
    for (int j = 0; j < face_length; ++j) {
        if (constraint[j] < 0 || constraint[j] >= n || !observed[constraint[j]]) {
            Rf_error("Invalid cancelled Gamma support face.");
        }
    }
    int signs = NA_INTEGER;
    try {
        SupportRows raw{n, p, starts, positions, entries, groups, observed,
                        std::vector<int>(n_groups)};
        for (int i = 0; i < n; ++i) ++raw.count[groups[i]];
        RationalRow candidate(p);
        for (int j = 0; j < p; ++j) {
            if (proposed[j] != 0) candidate[j] = rational(proposed[j]) / rational(slope_unit[j]);
        }
        signs = raw.signs(candidate);
        if (signs != 1 && face_length && repair_face(raw, candidate, constraint, face_length)) signs = 1;
    } catch (const std::bad_alloc&) {
        // Resource exhaustion proves no sign pattern.
    }
    return Rf_ScalarInteger(signs);
}
