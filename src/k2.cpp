#include "k2.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace powermean {
namespace {

using Size = std::size_t;
using Status = ProjectionStatus;

template <class T>
constexpr Size element_limit() noexcept {
    return static_cast<Size>(std::numeric_limits<std::ptrdiff_t>::max()) / sizeof(T);
}

template <class T>
bool extent(Size rows, Size columns, Size stride, Size& result) noexcept {
    const Size limit = element_limit<T>();
    if (rows == 0 || columns == 0 || stride < rows || rows > limit ||
        columns - 1 > (limit - rows) / stride) {
        return false;
    }
    result = (columns - 1) * stride + rows;
    return true;
}

struct Layout {
    ProjectionSizes sizes;
    Size fractions;
    Size max_levels;
};

Status plan(ProjectionShape shape, ProjectionBuffer<const Size> offsets,
            Layout& layout) noexcept {
    const auto [rows, blocks, columns, code_stride, value_stride] = shape;
    if (blocks >= element_limit<Size>() ||
        columns > element_limit<ProjectionDiagnostic>() ||
        !extent<std::int32_t>(rows, blocks, code_stride, layout.sizes.codes) ||
        !extent<double>(rows, columns, value_stride, layout.sizes.values) ||
        !extent<double>(rows, blocks, rows, layout.fractions)) {
        return Status::invalid_dimensions;
    }
    if (!offsets.data || offsets.size < blocks + 1) return Status::invalid_buffer;
    const Size limit = element_limit<double>();
    if (offsets.data[0] != 0 || offsets.data[blocks] > limit) {
        return Status::invalid_dimensions;
    }
    layout.max_levels = 0;
    for (Size block = 0; block < blocks; ++block) {
        const Size first = offsets.data[block];
        const Size last = offsets.data[block + 1];
        if (last <= first || last - first > UINT64_C(2147483648)) {
            return Status::invalid_dimensions;
        }
        layout.max_levels = std::max(layout.max_levels, last - first);
    }
    Size remaining = limit - offsets.data[blocks];
    if (layout.max_levels > remaining) return Status::invalid_dimensions;
    remaining -= layout.max_levels;
    if (layout.fractions > remaining) return Status::invalid_dimensions;
    remaining -= layout.fractions;
    if (rows > remaining / 5) return Status::invalid_dimensions;
    layout.sizes.workspace = limit - remaining + 5 * rows;
    layout.sizes.diagnostics = columns;
    return Status::converged;
}

template <class T>
bool has(ProjectionBuffer<T> buffer, Size count) noexcept {
    return buffer.data && buffer.size >= count;
}

struct Groups {
    ProjectionShape shape;
    const Size* offsets;
    const std::int32_t* codes;
    double* weights;
    double* fractions;
    double* totals;
    double* means;
};

Status normalize_weights(ProjectionBuffer<const double> input, Size rows,
                         double* weights) noexcept {
    if (!input.data) {
        std::fill_n(weights, rows, 1.0);
        return Status::converged;
    }
    double largest = 0.0;
    for (Size row = 0; row < rows; ++row) {
        const double weight = input.data[row];
        if (!(weight > 0.0) || !std::isfinite(weight)) return Status::invalid_weight;
        largest = std::max(largest, weight);
    }
    for (Size row = 0; row < rows; ++row) {
        weights[row] = input.data[row] / largest;
        if (!std::isnormal(weights[row])) return Status::numerical_failure;
    }
    return Status::converged;
}

Status prepare_groups(const Groups& groups) noexcept {
    const Size rows = groups.shape.rows;
    for (Size block = 0; block < groups.shape.blocks; ++block) {
        const auto* code = groups.codes + block * groups.shape.code_stride;
        const Size levels = groups.offsets[block + 1] - groups.offsets[block];
        auto* totals = groups.totals + groups.offsets[block];
        auto* fractions = groups.fractions ? groups.fractions + block * rows : nullptr;
        std::fill_n(totals, levels, 0.0);
        for (Size row = 0; row < rows; ++row) {
            if (code[row] < 0 || static_cast<Size>(code[row]) >= levels) {
                return Status::invalid_code;
            }
            totals[code[row]] += groups.weights[row];
        }
        if (fractions) {
            for (Size row = 0; row < rows; ++row) {
                fractions[row] = groups.weights[row] / totals[code[row]];
            }
        }
    }
    return Status::converged;
}

struct GroupBlock {
    Size rows;
    Size levels;
    const std::int32_t* code;
    const double* fractions;
    const double* totals;
    double* means;
};

GroupBlock group_block(const Groups& groups, Size block) noexcept {
    const Size first = groups.offsets[block];
    return {groups.shape.rows, groups.offsets[block + 1] - first,
            groups.codes + block * groups.shape.code_stride,
            groups.fractions ? groups.fractions + block * groups.shape.rows : nullptr,
            groups.totals + first, groups.means};
}

void group_means(const GroupBlock& block, const double* values) noexcept {
    std::fill_n(block.means, block.levels, 0.0);
    if (block.fractions) {
        for (Size row = 0; row < block.rows; ++row) {
            block.means[block.code[row]] += block.fractions[row] * values[row];
        }
    } else {
        for (Size row = 0; row < block.rows; ++row) block.means[block.code[row]] += values[row];
        for (Size level = 0; level < block.levels; ++level) {
            if (block.totals[level] > 0) block.means[level] /= block.totals[level];
        }
    }
}

double certificate(const Groups& groups, const double* residual) noexcept {
    double maximum = 0.0;
    for (Size j = 0; j < groups.shape.blocks; ++j) {
        const auto block = group_block(groups, j);
        group_means(block, residual);
        for (Size level = 0; level < block.levels; ++level) {
            const double value = std::abs(block.means[level]);
            if (!std::isfinite(value)) return value;
            maximum = std::max(maximum, value);
        }
    }
    return maximum;
}

void sweep_block(const GroupBlock& block, double* work, double* action) noexcept {
    for (Size row = 0; row < block.rows; ++row) {
        const double mean = block.means[block.code[row]];
        work[row] -= mean;
        action[row] += mean;
    }
}

// Apply H = I - F*F in the weighted row inner product. Accumulate the
// removed group means directly, avoiding subtraction of two row vectors.
template <bool capture = false>
double apply(const Groups& groups, const double* input, double* action,
             double* work, double* effects = nullptr) noexcept {
    const Size rows = groups.shape.rows;
    std::copy_n(input, rows, work);
    std::fill_n(action, rows, 0.0);
    if constexpr (capture) std::fill_n(effects, groups.offsets[groups.shape.blocks], 0.0);
    double energy = 0.0;
    for (Size j = 0; j < groups.shape.blocks; ++j) {
        const auto block = group_block(groups, j);
        group_means(block, work);
        for (Size level = 0; level < block.levels; ++level) {
            const double mean = block.means[level];
            energy += (block.totals[level] * mean) * mean;
            if constexpr (capture) effects[groups.offsets[j] + level] += mean;
        }
        sweep_block(block, work, action);
    }
    for (Size j = groups.shape.blocks - 1; j > 0; --j) {
        const auto block = group_block(groups, j - 1);
        group_means(block, work);
        if constexpr (capture) {
            for (Size level = 0; level < block.levels; ++level) {
                effects[groups.offsets[j - 1] + level] += block.means[level];
            }
        }
        sweep_block(block, work, action);
    }
    return energy;
}

struct Cg {
    double* direction;
    double* action;
    double* defect;
    double* work;
    double norm;
};

struct EffectRecurrence {
    Size levels;
    double* increment;
    double* direction;
    double* action;
    double* defect;
};

double squared_norm(const Groups& groups, const double* values) noexcept {
    double result = 0.0;
    for (Size row = 0; row < groups.shape.rows; ++row) {
        result += (groups.weights[row] * values[row]) * values[row];
    }
    return result;
}

template <bool capture = false>
void restart(const Groups& groups, Cg& cg, const double* residual,
             EffectRecurrence* effects = nullptr) noexcept {
    apply<capture>(groups, residual, cg.defect, cg.work, capture ? effects->defect : nullptr);
    std::copy_n(cg.defect, groups.shape.rows, cg.direction);
    if constexpr (capture) std::copy_n(effects->defect, effects->levels, effects->direction);
    cg.norm = squared_norm(groups, cg.defect);
}

template <bool capture = false>
bool advance_direction(const Groups& groups, Cg& cg, const double* residual,
                       bool check_recurrence, EffectRecurrence* effects = nullptr) noexcept {
    const double next = squared_norm(groups, cg.defect);
    const bool recheck = check_recurrence || !(next > 0.0) || !std::isfinite(next);
    if (recheck) {
        apply<capture>(groups, residual, cg.action, cg.work, capture ? effects->action : nullptr);
        for (Size row = 0; row < groups.shape.rows; ++row) {
            cg.work[row] = cg.action[row] - cg.defect[row];
        }
        const double actual = squared_norm(groups, cg.action);
        const double error = squared_norm(groups, cg.work);
        if (!(next > 0.0) || !std::isfinite(error) || error > 0.01 * actual) {
            std::copy_n(cg.action, groups.shape.rows, cg.defect);
            std::copy_n(cg.action, groups.shape.rows, cg.direction);
            cg.norm = actual;
            if constexpr (capture) {
                std::copy_n(effects->action, effects->levels, effects->defect);
                std::copy_n(effects->action, effects->levels, effects->direction);
            }
            return true;
        }
    }
    const double beta = next / cg.norm;
    for (Size row = 0; row < groups.shape.rows; ++row) {
        cg.direction[row] = cg.defect[row] + beta * cg.direction[row];
    }
    if constexpr (capture) {
        for (Size level = 0; level < effects->levels; ++level) {
            effects->direction[level] = effects->defect[level] + beta * effects->direction[level];
        }
    }
    cg.norm = next;
    return recheck;
}

template <bool capture = false>
bool cg_step(const Groups& groups, Cg& cg, double* residual,
             EffectRecurrence* effects = nullptr) noexcept {
    const double energy = apply<capture>(groups, cg.direction, cg.action, cg.work,
                                         capture ? effects->action : nullptr);
    if (!(energy > 0.0) || !std::isfinite(energy)) return false;
    const double alpha = cg.norm / energy;
    if (!(alpha > 0.0) || !std::isfinite(alpha)) return false;
    for (Size row = 0; row < groups.shape.rows; ++row) {
        residual[row] -= alpha * cg.direction[row];
        cg.defect[row] -= alpha * cg.action[row];
        if (!std::isfinite(residual[row])) return false;
    }
    if constexpr (capture) {
        for (Size level = 0; level < effects->levels; ++level) {
            effects->increment[level] += alpha * effects->direction[level];
            effects->defect[level] -= alpha * effects->action[level];
        }
    }
    return true;
}

bool forward_step(const Groups& groups, double* residual, double* work) noexcept {
    const Size rows = groups.shape.rows;
    std::copy_n(residual, rows, work);
    for (Size j = 0; j < groups.shape.blocks; ++j) {
        const auto block = group_block(groups, j);
        group_means(block, residual);
        for (Size row = 0; row < rows; ++row) {
            residual[row] -= block.means[block.code[row]];
        }
    }
    return !std::equal(residual, residual + rows, work);
}

bool scale_rhs(const double* source, double* work, Size rows, int exponent) noexcept {
    for (Size row = 0; row < rows; ++row) {
        const double scaled = std::scalbn(source[row], -exponent);
        if (std::abs(scaled) < std::numeric_limits<double>::min() &&
            std::scalbn(scaled, exponent) != source[row]) {
            return false;
        }
        work[row] = scaled;
    }
    return true;
}

ProjectionDiagnostic solve(const Groups& groups, Cg& cg, double* residual,
                           double largest, const ProjectionControls& controls) noexcept {
    const Size rows = groups.shape.rows;
    int exponent = 0;
    const double scaled_max = std::frexp(largest, &exponent);
    const double limit = std::scalbn(controls.absolute_tolerance, -exponent) +
                         controls.relative_tolerance * scaled_max;
    if (!scale_rhs(residual, cg.work, rows, exponent)) {
        // Check an already valid column without rounding its smallest entries.
        const double mean = certificate(groups, residual);
        const double unscaled_limit =
            controls.absolute_tolerance + controls.relative_tolerance * largest;
        const auto status = std::isfinite(mean) && mean <= unscaled_limit
            ? Status::converged : Status::numerical_failure;
        return {status, 0, largest, mean};
    }

    double mean = certificate(groups, cg.work);
    ProjectionDiagnostic result{Status::iteration_limit, 0, largest, std::scalbn(mean, exponent)};
    if (mean <= limit) {
        result.status = Status::converged;
        return result;
    }
    const bool forward = controls.method == ProjectionMethod::forward_sweeps ||
        (groups.shape.blocks == 1 && controls.method == ProjectionMethod::accelerated_sweeps);
    const bool count_operators = !forward && controls.method == ProjectionMethod::accelerated_sweeps;
    const std::uint64_t cost = count_operators ? 2 : 1;
    if (controls.max_iterations < cost) return result;

    std::copy_n(cg.work, rows, residual);
    if (!forward) {
        if (count_operators) result.iterations = cost;
        restart(groups, cg, residual);
    }
    std::uint64_t updates = 0;
    while (mean > limit && controls.max_iterations - result.iterations >= cost) {
        if (count_operators || forward) result.iterations += cost;
        if (forward) {
            if (!forward_step(groups, residual, cg.work)) {
                result.status = Status::numerical_failure;
                break;
            }
        } else if (!cg_step(groups, cg, residual)) {
            result.status = Status::numerical_failure;
            break;
        }
        if (!count_operators && !forward) ++result.iterations;
        ++updates;
        mean = certificate(groups, residual);
        if (!std::isfinite(mean)) break;
        if (!forward && mean > limit && controls.max_iterations - result.iterations >= cost) {
            const bool periodic = updates % 32 == 0 &&
                (!count_operators || controls.max_iterations - result.iterations >= 2 * cost);
            const bool rechecked = advance_direction(groups, cg, residual, periodic);
            if (count_operators && rechecked) result.iterations += cost;
        }
    }

    // Certify the values representable in the caller's original units.
    for (Size row = 0; row < rows; ++row) {
        residual[row] = std::scalbn(residual[row], exponent);
        cg.work[row] = std::scalbn(residual[row], -exponent);
    }
    const double returned_mean = certificate(groups, cg.work);
    result.max_abs_group_mean = std::scalbn(returned_mean, exponent);
    if (!std::isfinite(mean) || !std::isfinite(result.max_abs_group_mean) ||
        (mean <= limit && returned_mean > limit)) {
        result.status = Status::numerical_failure;
    } else if (returned_mean <= limit && result.status != Status::numerical_failure) {
        result.status = Status::converged;
    }
    return result;
}

bool valid_controls(const ProjectionControls& controls) noexcept {
    return std::isfinite(controls.absolute_tolerance) &&
           std::isfinite(controls.relative_tolerance) &&
           controls.absolute_tolerance >= 0.0 && controls.relative_tolerance >= 0.0;
}

bool validate_rhs(ProjectionShape shape, const double* values) noexcept {
    for (Size column = 0; column < shape.columns; ++column) {
        const double* rhs = values + column * shape.value_stride;
        for (Size row = 0; row < shape.rows; ++row) {
            if (!std::isfinite(rhs[row])) return false;
        }
    }
    return true;
}

double max_abs(const double* values, Size rows) noexcept {
    double result = 0.0;
    for (Size row = 0; row < rows; ++row) result = std::max(result, std::abs(values[row]));
    return result;
}

bool scaled_effects(const EffectRecurrence& recurrence, const double* initial,
                    int exponent, double* candidate) noexcept {
    for (Size level = 0; level < recurrence.levels; ++level) {
        const double increment = recurrence.increment[level];
        const double physical = std::scalbn(increment, exponent);
        const double value = initial ? initial[level] + physical : physical;
        if (!std::isfinite(value) || (physical == 0.0 && increment != 0.0)) return false;
        candidate[level] = value;
    }
    return true;
}

bool hold_effect_candidate(const Groups& groups, const double* rhs, const double* candidate,
                           double* work, double* residual, double* effects,
                           double& maximum, bool scaled_warm_check = false) noexcept {
    const Size rows = groups.shape.rows;
    const Size stride = groups.shape.code_stride;
    for (Size row = 0; row < rows; ++row) {
        work[row] = candidate[groups.codes[row]];
    }
    for (Size block = 1; block < groups.shape.blocks; ++block) {
        const auto* code = groups.codes + block * stride;
        const auto* column = candidate + groups.offsets[block];
        for (Size row = 0; row < rows; ++row) work[row] += column[code[row]];
    }
    for (Size row = 0; row < rows; ++row) {
        work[row] = rhs[row] - work[row];
        if (!std::isfinite(work[row])) return false;
    }
    std::copy_n(work, rows, residual);
    std::copy_n(candidate, groups.offsets[groups.shape.blocks], effects);
    maximum = certificate(groups, work);
    if (scaled_warm_check && !std::isfinite(maximum)) {
        // A warm start may be computed in a smaller exact binary chart even
        // when group sums of its initial physical residual would overflow.
        // The caller's requested absorption checks the actual physical rows.
        int shift = 0;
        std::frexp(max_abs(residual, rows), &shift);
        if (scale_rhs(residual, work, rows, shift)) {
            maximum = std::scalbn(certificate(groups, work), shift);
        }
    }
    return true;
}

bool forward_effect_step(const Groups& groups, const double* residual, const double* effects,
                         double* work, double* candidate) noexcept {
    const Size rows = groups.shape.rows;
    std::copy_n(residual, rows, work);
    for (Size j = 0; j < groups.shape.blocks; ++j) {
        const auto block = group_block(groups, j);
        group_means(block, work);
        for (Size level = 0; level < block.levels; ++level) {
            const Size index = groups.offsets[j] + level;
            candidate[index] = effects[index] + block.means[level];
            if (!std::isfinite(candidate[index])) return false;
        }
        for (Size row = 0; row < rows; ++row) {
            work[row] -= block.means[block.code[row]];
            if (!std::isfinite(work[row])) return false;
        }
    }
    return true;
}

ProjectionEffectDiagnostic solve_effect_candidates(
    const Groups& groups, Cg& cg, const double* rhs, const double* initial,
    double* residual, double* effects, const ProjectionControls& controls,
    double* working_residual, EffectRecurrence& recurrence, double* candidate) noexcept {
    constexpr double unavailable = std::numeric_limits<double>::quiet_NaN();
    ProjectionEffectDiagnostic result{
        {Status::candidate_range, 0, unavailable, unavailable}, false
    };
    ProjectionDiagnostic& diagnostic = result.projection;
    const bool scaled_warm_check = controls.method == ProjectionMethod::native_iterations;
    std::fill_n(recurrence.increment, recurrence.levels, 0.0);
    if (!scaled_effects(recurrence, initial, 0, candidate) ||
        !hold_effect_candidate(groups, rhs, candidate, cg.work, residual, effects,
                               diagnostic.max_abs_group_mean, scaled_warm_check)) return result;
    result.formed = true;
    const Size rows = groups.shape.rows;
    const double largest = max_abs(residual, rows);
    diagnostic.rhs_max_abs = largest;
    if (!std::isfinite(diagnostic.max_abs_group_mean)) return result;

    int exponent = 0;
    const double scaled_max = std::frexp(largest, &exponent);
    const double limit = std::scalbn(controls.absolute_tolerance, -exponent) +
                         controls.relative_tolerance * scaled_max;
    const auto within = [&] {
        return std::scalbn(diagnostic.max_abs_group_mean, -exponent) <= limit;
    };
    if (within()) {
        diagnostic.status = Status::converged;
        return result;
    }
    diagnostic.status = Status::iteration_limit;
    const bool forward = controls.method == ProjectionMethod::forward_sweeps ||
        (groups.shape.blocks == 1 && controls.method == ProjectionMethod::accelerated_sweeps);
    if (forward) {
        while (diagnostic.iterations < controls.max_iterations) {
            ++diagnostic.iterations;
            // A level update can move even when each individual row subtraction
            // rounds away. Start the next sweep from the newly lifted candidate.
            if (!forward_effect_step(groups, residual, effects, cg.work, candidate) ||
                !hold_effect_candidate(groups, rhs, candidate, cg.work, residual, effects,
                                       diagnostic.max_abs_group_mean) ||
                !std::isfinite(diagnostic.max_abs_group_mean)) {
                diagnostic.status = Status::candidate_range;
                break;
            }
            if (within()) {
                diagnostic.status = Status::converged;
                break;
            }
        }
        return result;
    }
    const bool count_operators = controls.method == ProjectionMethod::accelerated_sweeps;
    const std::uint64_t cost = count_operators ? 2 : 1;
    if (controls.max_iterations < cost) return result;
    if (!scale_rhs(residual, working_residual, rows, exponent)) {
        diagnostic.status = Status::candidate_range;
        return result;
    }

    if (count_operators) diagnostic.iterations = cost;
    restart<true>(groups, cg, working_residual, &recurrence);
    std::uint64_t updates = 0;
    while (controls.max_iterations - diagnostic.iterations >= cost) {
        if (count_operators) diagnostic.iterations += cost;
        if (!cg_step<true>(groups, cg, working_residual, &recurrence)) {
            diagnostic.status = Status::numerical_failure;
            break;
        }
        if (!count_operators) ++diagnostic.iterations;
        ++updates;

        // The row recurrence supplies the search direction. The candidate's
        // coefficients, original lift and freshly checked residual own acceptance.
        if (!scaled_effects(recurrence, initial, exponent, candidate) ||
            !hold_effect_candidate(groups, rhs, candidate, cg.work, residual, effects,
                                   diagnostic.max_abs_group_mean, scaled_warm_check) ||
            !std::isfinite(diagnostic.max_abs_group_mean)) {
            diagnostic.status = Status::candidate_range;
            break;
        }
        if (within()) {
            diagnostic.status = Status::converged;
            break;
        }
        if (controls.max_iterations - diagnostic.iterations >= cost) {
            const bool periodic = updates % 32 == 0 &&
                (!count_operators || controls.max_iterations - diagnostic.iterations >= 2 * cost);
            const bool rechecked = advance_direction<true>(groups, cg, working_residual,
                                                          periodic, &recurrence);
            if (count_operators && rechecked) diagnostic.iterations += cost;
        }
    }
    return result;
}

Status solve_columns(const Groups& groups, Cg& cg, double* values,
                     const ProjectionControls& controls,
                     ProjectionDiagnostic* diagnostics) noexcept {
    Status status = Status::converged;
    for (Size column = 0; column < groups.shape.columns; ++column) {
        double* rhs = values + column * groups.shape.value_stride;
        const double largest = max_abs(rhs, groups.shape.rows);
        const auto result = solve(groups, cg, rhs, largest, controls);
        diagnostics[column] = result;
        if (status == Status::converged || result.status == Status::numerical_failure) {
            status = result.status;
        }
    }
    return status;
}

}  // namespace

ProjectionStatus projection_sizes(ProjectionShape shape,
                                  ProjectionBuffer<const Size> offsets,
                                  ProjectionSizes& sizes) noexcept {
    Layout layout{};
    const auto status = plan(shape, offsets, layout);
    if (status == Status::converged) sizes = layout.sizes;
    return status;
}

ProjectionStatus project_fixed_effects(
    ProjectionShape shape, ProjectionBuffer<const Size> offsets,
    ProjectionBuffer<const std::int32_t> codes, ProjectionBuffer<const double> weights,
    ProjectionBuffer<double> values, ProjectionControls controls,
    ProjectionBuffer<double> workspace, ProjectionBuffer<ProjectionDiagnostic> diagnostics) noexcept {
    Layout layout{};
    const auto status = plan(shape, offsets, layout);
    if (status != Status::converged) return status;
    if (!has(codes, layout.sizes.codes) || !has(values, layout.sizes.values) ||
        !has(diagnostics, layout.sizes.diagnostics) ||
        (weights.data ? weights.size < shape.rows : weights.size != 0)) {
        return Status::invalid_buffer;
    }
    if (!has(workspace, layout.sizes.workspace)) return Status::insufficient_workspace;
    if (!valid_controls(controls)) return Status::invalid_controls;
    if (!validate_rhs(shape, values.data)) return Status::invalid_rhs;

    const Size rows = shape.rows;
    double* fractions = workspace.data + 5 * rows;
    double* totals = fractions + layout.fractions;
    Groups groups{shape, offsets.data, codes.data, workspace.data, weights.data ? fractions : nullptr,
                  totals, totals + offsets.data[shape.blocks]};
    const auto weight_status = normalize_weights(weights, rows, groups.weights);
    if (weight_status != Status::converged) return weight_status;
    const auto group_status = prepare_groups(groups);
    if (group_status != Status::converged) return group_status;

    Cg cg{workspace.data + rows, workspace.data + 2 * rows,
          workspace.data + 3 * rows, workspace.data + 4 * rows, 0.0};
    return solve_columns(groups, cg, values.data, controls, diagnostics.data);
}

ProjectionStatus projection_effect_sizes(
    ProjectionShape shape, ProjectionBuffer<const Size> offsets,
    ProjectionEffectSizes& sizes) noexcept {
    Layout layout{};
    const auto status = plan(shape, offsets, layout);
    if (status != Status::converged) return status;
    const Size levels = offsets.data[shape.blocks];
    const Size limit = element_limit<double>();
    if (levels > limit / shape.columns ||
        shape.columns > element_limit<ProjectionEffectDiagnostic>() ||
        shape.rows > limit - layout.sizes.workspace ||
        levels > (limit - layout.sizes.workspace - shape.rows) / 5) {
        return Status::invalid_dimensions;
    }
    sizes = {layout.sizes, levels * shape.columns,
             layout.sizes.workspace + shape.rows + 5 * levels};
    return Status::converged;
}

ProjectionStatus project_fixed_effect_candidates(
    ProjectionShape shape, ProjectionBuffer<const Size> offsets,
    ProjectionBuffer<const std::int32_t> codes, ProjectionBuffer<const double> weights,
    ProjectionBuffer<const double> rhs, ProjectionBuffer<const double> initial,
    ProjectionBuffer<double> values, ProjectionBuffer<double> effects,
    ProjectionControls controls, ProjectionBuffer<double> workspace,
    ProjectionBuffer<ProjectionEffectDiagnostic> diagnostics) noexcept {
    constexpr double unavailable = std::numeric_limits<double>::quiet_NaN();
    const auto fail = [&](Status status) {
        if (has(diagnostics, shape.columns)) {
            std::fill_n(diagnostics.data, shape.columns,
                        ProjectionEffectDiagnostic{{status, 0, unavailable, unavailable}, false});
        }
        return status;
    };
    ProjectionEffectSizes sizes{};
    const auto status = projection_effect_sizes(shape, offsets, sizes);
    if (status != Status::converged) return fail(status);
    const auto& base = sizes.projection;
    if (!has(codes, base.codes) || !has(rhs, base.values) || !has(values, base.values) ||
        !has(effects, sizes.effects) || !has(diagnostics, base.diagnostics) ||
        (initial.data ? initial.size < sizes.effects : initial.size != 0) ||
        (weights.data ? weights.size < shape.rows : weights.size != 0)) {
        return fail(Status::invalid_buffer);
    }
    if (!has(workspace, sizes.workspace)) return fail(Status::insufficient_workspace);
    if (!valid_controls(controls)) return fail(Status::invalid_controls);
    if (!validate_rhs(shape, rhs.data)) return fail(Status::invalid_rhs);

    const Size rows = shape.rows;
    const Size levels = offsets.data[shape.blocks];
    double* fractions = workspace.data + 5 * rows;
    double* totals = fractions + rows * shape.blocks;
    Groups groups{shape, offsets.data, codes.data, workspace.data, weights.data ? fractions : nullptr,
                  totals, totals + levels};
    const auto weight_status = normalize_weights(weights, rows, groups.weights);
    if (weight_status != Status::converged) return fail(weight_status);
    const auto group_status = prepare_groups(groups);
    if (group_status != Status::converged) return fail(group_status);
    Cg cg{workspace.data + rows, workspace.data + 2 * rows,
          workspace.data + 3 * rows, workspace.data + 4 * rows, 0.0};

    double* working_residual = workspace.data + base.workspace;
    double* coordinates = working_residual + rows;
    EffectRecurrence recurrence{levels, coordinates, coordinates + levels,
                                coordinates + 2 * levels, coordinates + 3 * levels};
    double* candidate = coordinates + 4 * levels;
    Status overall = Status::converged;
    for (Size col = 0; col < shape.columns; ++col) {
        const Size row_start = col * shape.value_stride;
        const Size level_start = col * levels;
        const double* seed = initial.data ? initial.data + level_start : nullptr;
        const auto result = solve_effect_candidates(groups, cg, rhs.data + row_start, seed,
            values.data + row_start, effects.data + level_start, controls,
            working_residual, recurrence, candidate);
        diagnostics.data[col] = result;
        const Status ending = result.projection.status;
        if (static_cast<int>(ending) > static_cast<int>(overall)) overall = ending;
    }
    return overall;
}

}  // namespace powermean
