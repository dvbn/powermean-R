#define R_NO_REMAP
#include <R.h>
#include <R_ext/Utils.h>
#include <Rinternals.h>

#include "native_calls.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>

namespace {

struct Rows {
    int n;
    int d;
    int groups;
    const int* offset;
    const int* column;
    const double* value;
    const int* group;
    const double* inverse_count;
    const int* ordered_counts;

    double product(int row, const double* beta) const {
        double value_sum = 0;
        for (int a = offset[row]; a < offset[row + 1]; ++a) {
            value_sum += value[a] * beta[column[a]];
        }
        return value_sum;
    }

    void add_row(int row, double weight, double* result) const {
        for (int a = offset[row]; a < offset[row + 1]; ++a) {
            result[column[a]] += value[a] * weight;
        }
    }

    void project(const double* beta, double* row_values, double* total) const {
        std::fill_n(total, groups, 0.0);
        for (int row = 0; row < n; ++row) {
            row_values[row] = product(row, beta);
            total[group[row]] += row_values[row];
        }
        for (int row = 0; row < n; ++row) {
            const int g = group[row];
            row_values[row] -= total[g] * inverse_count[g];
        }
    }

    void project_transpose(const double* row_values, double* result, double* total) const {
        std::fill_n(total, groups, 0.0);
        std::fill_n(result, d, 0.0);
        for (int row = 0; row < n; ++row) total[group[row]] += row_values[row];
        for (int row = 0; row < n; ++row) {
            const int g = group[row];
            add_row(row, row_values[row] - total[g] * inverse_count[g], result);
        }
    }
};

double inner(const double* left, const double* right, int size) {
    double total = 0;
    for (int j = 0; j < size; ++j) total += left[j] * right[j];
    return total;
}

double maximum(const double* data, int size) {
    double result = 0;
    for (int j = 0; j < size; ++j) {
        if (!std::isfinite(data[j])) return std::numeric_limits<double>::infinity();
        result = std::max(result, std::abs(data[j]));
    }
    return result;
}

void apply(const Rows& design, const double* weight, const double* source, double* result,
           double* row_values, double* total) {
    if (design.ordered_counts) {
        std::fill_n(result, design.d, 0.0);
        int first = 0;
        for (int g = 0; g < design.groups; ++g) {
            const int count = design.ordered_counts[g];
            const int end = first + count;
            if (count == 2) {
                const double delta = (design.product(first, source) - design.product(first + 1, source)) *
                    (0.25 * weight[first] + 0.25 * weight[first + 1]);
                design.add_row(first, delta, result);
                design.add_row(first + 1, -delta, result);
            } else if (count > 2) {
                double sum = 0;
                for (int row = first; row < end; ++row) {
                    row_values[row] = design.product(row, source);
                    sum += row_values[row];
                }
                const double mean = sum * design.inverse_count[g];
                double weighted_sum = 0;
                for (int row = first; row < end; ++row) {
                    row_values[row] = (row_values[row] - mean) * weight[row];
                    weighted_sum += row_values[row];
                }
                const double weighted_mean = weighted_sum * design.inverse_count[g];
                for (int row = first; row < end; ++row) {
                    design.add_row(row, row_values[row] - weighted_mean, result);
                }
            }
            first = end;
        }
        return;
    }
    design.project(source, row_values, total);
    for (int row = 0; row < design.n; ++row) row_values[row] *= weight[row];
    design.project_transpose(row_values, result, total);
}

bool direction(const Rows& design, const double* weight, const double* score,
               const double* diagonal, double* step, double* work,
               double* row_values, double* total, double relative) {
    const int d = design.d;
    double* residual = work;
    double* scaled = residual + d;
    double* search = scaled + d;
    double* image = search + d;
    std::copy_n(score, d, residual);
    std::fill_n(step, d, 0.0);
    const double largest = maximum(score, d);
    if (!std::isfinite(largest) || largest == 0) return largest == 0;
    const double goal = largest * relative;
    for (int j = 0; j < d; ++j) {
        if (!(diagonal[j] > 0) || !std::isfinite(diagonal[j])) return false;
        search[j] = scaled[j] = residual[j] / diagonal[j];
    }
    double norm = inner(residual, scaled, d);
    for (int iteration = 0; iteration < 400; ++iteration) {
        apply(design, weight, search, image, row_values, total);
        const double curvature = inner(search, image, d);
        if (!(curvature > 0) || !std::isfinite(curvature) || !(norm > 0)) return false;
        const double size = norm / curvature;
        for (int j = 0; j < d; ++j) {
            step[j] += size * search[j];
            residual[j] -= size * image[j];
        }
        if (maximum(residual, d) <= goal) {
            apply(design, weight, step, image, row_values, total);
            for (int j = 0; j < d; ++j) image[j] -= score[j];
            return maximum(image, d) <= std::max(largest * 1e-9, 2 * goal);
        }
        for (int j = 0; j < d; ++j) scaled[j] = residual[j] / diagonal[j];
        const double next = inner(residual, scaled, d);
        const double ratio = next / norm;
        for (int j = 0; j < d; ++j) search[j] = scaled[j] + ratio * search[j];
        norm = next;
    }
    return false;
}

bool evaluate(const Rows& design, const double* logs, const double* shift, const double* beta,
              double* weight, double* score, double* diagonal, double* row_values, double* total) {
    design.project(beta, row_values, total);
    std::fill_n(diagonal, design.d, 0.0);
    for (int row = 0; row < design.n; ++row) {
        const double index = shift[row] + row_values[row];
        if (!std::isfinite(index)) return false;
        const double t = std::exp(logs[row] - index);
        if (!std::isfinite(t) || (std::isfinite(logs[row]) && t == 0)) return false;
        weight[row] = t;
        // The raw weighted diagonal preconditions the centered operator. The
        // actual residual check uses the centered operator in every iteration.
        for (int a = design.offset[row]; a < design.offset[row + 1]; ++a) {
            const double value = design.value[a];
            diagonal[design.column[a]] += t * value * value;
        }
    }
    design.project_transpose(weight, score, total);
    return true;
}

bool solve(const Rows& design, const double* logs, const double* shift, double* beta,
           double* weight, double* changes, double* score, double* diagonal,
           double* step, double* work, double* total, int limit, double tolerance,
           bool adaptive, int& iteration) {
    double starting_score = 0;
    for (iteration = 0; iteration <= limit; ++iteration) {
        R_CheckUserInterrupt();
        if (!evaluate(design, logs, shift, beta, weight, score, diagonal, changes, total)) return false;
        double relative = 1e-11;
        if (adaptive) {
            const double size = maximum(score, design.d);
            if (!std::isfinite(size)) return false;
            if (iteration == 0) starting_score = size;
            // Save linear iterations away from the root. The outer decrease
            // uses the original objective and the final direction is tight.
            if (starting_score > 0 && size > starting_score * 1e-8) {
                relative = 0.01 * std::sqrt(std::min(1.0, size / starting_score));
            }
        }
        if (!direction(design, weight, score, diagonal, step, work, changes, total, relative)) return false;
        design.project(step, changes, total);
        double max_change = maximum(changes, design.n);
        if (max_change <= tolerance && relative > 1e-11) {
            if (!direction(design, weight, score, diagonal, step, work, changes, total, 1e-11)) return false;
            design.project(step, changes, total);
            max_change = maximum(changes, design.n);
        }
        if (max_change <= tolerance) return true;
        long double energy = 0;
        if (adaptive) {
            for (int column = 0; column < design.d; ++column) {
                energy += static_cast<long double>(score[column]) * step[column];
            }
        } else {
            for (int row = 0; row < design.n; ++row) {
                const double delta = changes[row];
                energy += static_cast<long double>(weight[row]) * delta * delta;
            }
        }
        if (iteration == limit || !std::isfinite(max_change) || !(energy > 0)) return false;
        double fraction = 1;
        bool accepted = false;
        for (int backtrack = 0; backtrack < 40; ++backtrack, fraction *= 0.5) {
            long double loss = 0;
            for (int row = 0; row < design.n; ++row) {
                if (weight[row] != 0) {
                    loss += static_cast<long double>(weight[row]) * std::expm1(-fraction * changes[row]);
                }
            }
            if (std::isfinite(loss) && loss <= -1e-4L * fraction * energy) {
                for (int j = 0; j < design.d; ++j) beta[j] += fraction * step[j];
                accepted = true;
                break;
            }
        }
        if (!accepted) return false;
    }
    return false;
}

double* workspace(std::size_t size) {
    return reinterpret_cast<double*>(R_alloc(size, sizeof(double)));
}

// Reproduce the frozen criterion's binary64 reduction. In particular the
// late Armijo comparison subtracts two losses separated by a few ulps.
double trial_mass_sum(const double* mass, std::size_t rows) {
    if (rows == 0) return 0.0;
    if (rows < 8) {
        double total = -0.0;
        for (std::size_t i = 0; i < rows; ++i) total += mass[i];
        return total;
    }
    if (rows > 128) {
        const std::size_t split = (rows / 16) * 8;
        return trial_mass_sum(mass, split) + trial_mass_sum(mass + split, rows - split);
    }
    double lanes[8];
    std::copy_n(mass, 8, lanes);
    const std::size_t full = rows - rows % 8;
    for (std::size_t first = 8; first < full; first += 8) {
        for (std::size_t lane = 0; lane < 8; ++lane) lanes[lane] += mass[first + lane];
    }
    double total = ((lanes[0] + lanes[1]) + (lanes[2] + lanes[3])) +
        ((lanes[4] + lanes[5]) + (lanes[6] + lanes[7]));
    for (std::size_t i = full; i < rows; ++i) total += mass[i];
    return total;
}

}  // namespace

SEXP cancel_gamma_sum_call(SEXP mass) {
    if (TYPEOF(mass) != REALSXP) Rf_error("Gamma trial mass must be numeric.");
    return Rf_ScalarReal(trial_mass_sum(REAL(mass), static_cast<std::size_t>(Rf_xlength(mass))));
}

SEXP cancel_gamma_call(SEXP offsets, SEXP columns, SEXP values, SEXP logs,
                       SEXP layout, SEXP start, SEXP controls) {
    if (TYPEOF(layout) != VECSXP || Rf_xlength(layout) != 3) Rf_error("Invalid Gamma group layout.");
    const SEXP offset = VECTOR_ELT(layout, 0);
    const SEXP group = VECTOR_ELT(layout, 1);
    const SEXP counts = VECTOR_ELT(layout, 2);
    if (TYPEOF(offsets) != INTSXP || TYPEOF(columns) != INTSXP || TYPEOF(values) != REALSXP ||
        TYPEOF(logs) != REALSXP || TYPEOF(offset) != REALSXP || TYPEOF(start) != REALSXP ||
        TYPEOF(controls) != REALSXP || (Rf_xlength(controls) != 2 && Rf_xlength(controls) != 3) ||
        Rf_xlength(logs) < 1 || Rf_xlength(logs) >= INT_MAX ||
        Rf_xlength(start) < 1 || Rf_xlength(start) >= INT_MAX ||
        Rf_xlength(offset) != Rf_xlength(logs) || Rf_xlength(offsets) != Rf_xlength(logs) + 1 ||
        TYPEOF(group) != INTSXP || Rf_xlength(group) != Rf_xlength(logs) ||
        TYPEOF(counts) != INTSXP || Rf_xlength(counts) < 1 || Rf_xlength(counts) > Rf_xlength(logs) ||
        Rf_xlength(columns) != Rf_xlength(values) || Rf_xlength(values) >= INT_MAX) {
        Rf_error("Invalid sparse Gamma solver arrays.");
    }
    const int n = static_cast<int>(Rf_xlength(logs));
    const int d = static_cast<int>(Rf_xlength(start));
    const int groups = static_cast<int>(Rf_xlength(counts));
    double* inverse_count = workspace(groups);
    double* total = workspace(groups);
    std::fill_n(total, groups, 0.0);
    Rows design{n, d, groups, INTEGER(offsets), INTEGER(columns), REAL(values),
                INTEGER(group), inverse_count, nullptr};
    const int entries = static_cast<int>(Rf_xlength(values));
    if (design.offset[0] != 0 || design.offset[n] != entries) Rf_error("Invalid sparse row offsets.");
    bool ordered = true;
    int last_group = -1;
    for (int row = 0; row < n; ++row) {
        if (design.offset[row] < 0 || design.offset[row] > design.offset[row + 1] ||
            design.offset[row + 1] > entries || std::isnan(REAL(logs)[row]) ||
            design.group[row] < 0 || design.group[row] >= groups) {
            Rf_error("Invalid sparse Gamma row.");
        }
        ordered = ordered && design.group[row] >= last_group;
        last_group = design.group[row];
        ++total[last_group];
    }
    for (int g = 0; g < groups; ++g) {
        if (total[g] < 1 || total[g] != INTEGER(counts)[g]) Rf_error("Invalid Gamma group count.");
        inverse_count[g] = 1.0 / total[g];
    }
    if (ordered) design.ordered_counts = INTEGER(counts);
    for (int a = 0; a < entries; ++a) {
        if (design.column[a] < 0 || design.column[a] >= d || !std::isfinite(design.value[a])) {
            Rf_error("Invalid sparse Gamma column.");
        }
    }
    const double budget = REAL(controls)[0];
    const double tolerance = REAL(controls)[1];
    const double mode = Rf_xlength(controls) == 3 ? REAL(controls)[2] : 0;
    if (!std::isfinite(budget) || budget < 0 || budget >= INT_MAX ||
        budget != std::floor(budget) || !std::isfinite(tolerance) || tolerance <= 0 ||
        (mode != 0 && mode != 1)) {
        Rf_error("Invalid sparse Gamma solver controls.");
    }
    double* score = workspace(d);
    double* diagonal = workspace(d);
    double* step = workspace(d);
    double* work = workspace(static_cast<std::size_t>(4) * d);
    double* weight = workspace(n);
    double* changes = workspace(n);
    SEXP beta = PROTECT(Rf_duplicate(start));
    int iteration = 0;
    bool accepted = solve(design, REAL(logs), REAL(offset), REAL(beta), weight, changes,
                          score, diagonal, step, work, total,
                          static_cast<int>(budget), tolerance, mode == 1, iteration);
    if (!accepted) {
        UNPROTECT(1);
        return R_NilValue;
    }
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(result, 0, beta);
    SET_VECTOR_ELT(result, 1, Rf_ScalarInteger(iteration));
    UNPROTECT(2);
    return result;
}

SEXP cancelled_log_recover_call(SEXP group, SEXP counts, SEXP response, SEXP raw, SEXP centered) {
    if (TYPEOF(group) != INTSXP || TYPEOF(counts) != INTSXP || TYPEOF(response) != REALSXP ||
        TYPEOF(raw) != REALSXP || TYPEOF(centered) != REALSXP || Rf_xlength(group) >= INT_MAX ||
        Rf_xlength(counts) < 1 || Rf_xlength(counts) > Rf_xlength(group) ||
        Rf_xlength(group) != Rf_xlength(response) || Rf_xlength(group) != Rf_xlength(raw) ||
        Rf_xlength(group) != Rf_xlength(centered)) {
        Rf_error("Invalid cancelled recovery arrays.");
    }
    const int n = static_cast<int>(Rf_xlength(group));
    const int groups = static_cast<int>(Rf_xlength(counts));
    const double missing = -std::numeric_limits<double>::infinity();
    const int* g = INTEGER(group);
    const double* outcome = REAL(response);
    const double* eta[] = {REAL(raw), REAL(centered)};
    double* maxima[] = {workspace(groups), workspace(groups)};
    int* seen = reinterpret_cast<int*>(R_alloc(groups, sizeof(int)));
    std::fill_n(seen, groups, 0);
    for (auto* top : maxima) std::fill_n(top, groups, missing);
    for (int row = 0; row < n; ++row) {
        if (g[row] < 0 || g[row] >= groups) Rf_error("Invalid cancelled recovery group.");
        if (!std::isfinite(outcome[row]) || outcome[row] < 0) Rf_error("Invalid cancelled response.");
        ++seen[g[row]];
        if (outcome[row] == 0) continue;
        const double log_response = std::log(outcome[row]);
        for (int j = 0; j < 2; ++j) {
            if (!std::isfinite(eta[j][row])) Rf_error("Invalid cancelled recovery index.");
            maxima[j][g[row]] = std::max(maxima[j][g[row]], log_response - eta[j][row]);
        }
    }
    for (int i = 0; i < groups; ++i) {
        if (seen[i] == 0 || seen[i] != INTEGER(counts)[i]) Rf_error("Invalid cancelled recovery count.");
    }
    SEXP result = PROTECT(Rf_allocVector(VECSXP, 2));
    SEXP absolute = PROTECT(Rf_allocVector(REALSXP, groups));
    SEXP relative = PROTECT(Rf_allocVector(REALSXP, groups));
    double* sums[] = {REAL(absolute), REAL(relative)};
    for (auto* sum : sums) std::fill_n(sum, groups, 0.0);
    for (int row = 0; row < n; ++row) {
        if (outcome[row] == 0) continue;
        const double log_response = std::log(outcome[row]);
        for (int j = 0; j < 2; ++j) {
            const double top = maxima[j][g[row]];
            if (std::isfinite(top)) sums[j][g[row]] += std::exp(log_response - eta[j][row] - top);
        }
    }
    for (int j = 0; j < 2; ++j) {
        for (int i = 0; i < groups; ++i) {
            const double top = maxima[j][i];
            sums[j][i] = std::isfinite(top) ? std::log(sums[j][i] / seen[i]) + top : top;
        }
    }
    SET_VECTOR_ELT(result, 0, absolute);
    SET_VECTOR_ELT(result, 1, relative);
    UNPROTECT(3);
    return result;
}
