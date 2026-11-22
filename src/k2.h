#ifndef POWERMEAN_K2_H
#define POWERMEAN_K2_H

#include <cstddef>
#include <cstdint>

namespace powermean {

enum class ProjectionStatus {
    converged,
    invalid_dimensions,
    invalid_buffer,
    invalid_code,
    invalid_weight,
    invalid_rhs,
    invalid_controls,
    insufficient_workspace,
    iteration_limit,
    numerical_failure,
    candidate_range
};

template <class T>
struct ProjectionBuffer {
    T* data;
    std::size_t size;
};

struct ProjectionShape {
    std::size_t rows;
    std::size_t blocks;
    std::size_t columns;
    std::size_t code_stride;
    std::size_t value_stride;
};

struct ProjectionSizes {
    std::size_t codes;
    std::size_t values;
    std::size_t workspace;
    std::size_t diagnostics;
};

enum class ProjectionMethod {
    native_iterations,
    forward_sweeps,
    accelerated_sweeps
};

struct ProjectionControls {
    double absolute_tolerance;
    double relative_tolerance;
    std::uint64_t max_iterations;
    ProjectionMethod method = ProjectionMethod::native_iterations;
};

struct ProjectionDiagnostic {
    ProjectionStatus status;
    std::uint64_t iterations;
    double rhs_max_abs;
    double max_abs_group_mean;
};

struct ProjectionEffectSizes {
    ProjectionSizes projection;
    std::size_t effects;
    std::size_t workspace;
};

struct ProjectionEffectDiagnostic {
    ProjectionDiagnostic projection;
    bool formed;
};

// All sizes and strides are in elements. Block offsets are cumulative level
// counts starting at zero. Each code column uses its own local zero-based IDs.
ProjectionStatus projection_sizes(ProjectionShape shape,
                                  ProjectionBuffer<const std::size_t> offsets,
                                  ProjectionSizes& sizes) noexcept;

// Replaces active values in place. Use {nullptr, 0} for unit weights.
// All buffer ranges must be disjoint. See README.md for output on failure.
ProjectionStatus project_fixed_effects(
    ProjectionShape shape, ProjectionBuffer<const std::size_t> offsets,
    ProjectionBuffer<const std::int32_t> codes,
    ProjectionBuffer<const double> weights, ProjectionBuffer<double> values,
    ProjectionControls controls, ProjectionBuffer<double> workspace,
    ProjectionBuffer<ProjectionDiagnostic> diagnostics) noexcept;

ProjectionStatus projection_effect_sizes(
    ProjectionShape shape, ProjectionBuffer<const std::size_t> offsets,
    ProjectionEffectSizes& sizes) noexcept;

// The RHS is immutable. Initial and returned effect columns contain all levels
// consecutively in block order. Use {nullptr, 0} for a zero initial state or unit
// weights. A formed column holds its last finite effects and RHS minus their
// lift. A nonfinite group diagnostic means that candidate's check refused.
// All buffers are disjoint, including the original and returned row columns.
ProjectionStatus project_fixed_effect_candidates(
    ProjectionShape shape, ProjectionBuffer<const std::size_t> offsets,
    ProjectionBuffer<const std::int32_t> codes, ProjectionBuffer<const double> weights,
    ProjectionBuffer<const double> rhs, ProjectionBuffer<const double> initial,
    ProjectionBuffer<double> values, ProjectionBuffer<double> effects,
    ProjectionControls controls, ProjectionBuffer<double> workspace,
    ProjectionBuffer<ProjectionEffectDiagnostic> diagnostics) noexcept;

}  // namespace powermean

#endif
