// Adapted from Powermean 0.1.0, BSD-3-Clause. See inst/LICENSE.K4.
// Frozen source: 7b22dcf01c6fea74176823e24a4b18df559cd386
// src/powermean/_native/iols/bind.cpp
// SHA-256 7d12b4308bbff9f4ee6f3dc3f72f77c06874462c7e762d6bc7c8ad707060e4e5.

#include "k4.h"
#include "k4_core.h"

namespace powermean {
namespace {

using namespace k4;

void validate_direct_nofe_inputs_core(
    const Eigen::Ref<const Eigen::VectorXd>& endog,
    const Eigen::Ref<const RowMatrixXd>& exog_slopes,
    double tau) {
    if (exog_slopes.rows() != endog.size()) {
        throw std::invalid_argument(
            "endog and exog_slopes must have the same number of rows");
    }
    if (exog_slopes.cols() < 1) {
        throw std::invalid_argument(
            "exog_slopes must have at least one slope column; the intercept is "
            "profiled");
    }
    bool any_positive = false;
    for (Eigen::Index row = 0; row < endog.size(); ++row) {
        if (!std::isfinite(endog(row)) || endog(row) < 0.0) {
            throw std::invalid_argument("endog must be finite and nonnegative");
        }
        any_positive = any_positive || endog(row) > 0.0;
    }
    if (!any_positive) {
        throw std::invalid_argument("endog must contain at least one positive value");
    }
    for (Eigen::Index row = 0; row < exog_slopes.rows(); ++row) {
        for (Eigen::Index col = 0; col < exog_slopes.cols(); ++col) {
            if (!std::isfinite(exog_slopes(row, col))) {
                throw std::invalid_argument("exog_slopes must be finite");
            }
        }
    }
    for (Eigen::Index col = 0; col < exog_slopes.cols(); ++col) {
        if (exog_slopes.col(col).maxCoeff() == exog_slopes.col(col).minCoeff()) {
            throw std::invalid_argument(
                "exog_slopes column " + std::to_string(col) +
                " is constant; pass slopes only, the intercept is profiled");
        }
    }
    if (!std::isfinite(tau) || tau < 0.0 || tau > 3.0) {
        throw std::invalid_argument("tau must be finite and in [0, 3]");
    }
}

void validate_direct_nofe_control_core(const FullrecoverControlCore& control) {
    if (control.delta_grid_size < 1) {
        throw std::invalid_argument(
            "direct fixed-tau route requires delta_grid_size >= 1 so Phase 1 runs");
    }
    if (control.anderson_m < 0) {
        throw std::invalid_argument("anderson_m must be nonnegative");
    }
    if (!std::isfinite(control.eps) || control.eps <= 0.0) {
        throw std::invalid_argument("eps must be finite and positive");
    }
    if (control.stopcrit < 2 || control.max_phase1_iterations < 2) {
        throw std::invalid_argument(
            "stopcrit and max_phase1_iterations must be at least 2");
    }
    if (control.max_phase2_evaluations < 1) {
        throw std::invalid_argument("max_phase2_evaluations must be at least 1");
    }
    if (control.direct_restart_start != kDirectRestartStartPhase1 &&
        control.direct_restart_start != kDirectRestartStartBestScore) {
        throw std::invalid_argument(
            "direct_restart_start must be one of ('phase1_start', "
            "'best_score_iterate')");
    }
}


// The public adapter changes the PPML response unit only outside the
// default grid's scaling span, and only when every original value is restored.
int condition_response(Eigen::VectorXd& response, double tau,
                       const DenseIolsOptions& options) {
    if (tau != 1 || !options.pilot_deltas.empty() ||
        options.rho_log_min != -7 || options.rho_log_max != 14) return 0;
    const double largest = response.maxCoeff();
    const double average = largest * (response / largest).mean();
    const double span = (1 + std::exp(options.rho_log_max)) /
                        (1 + std::exp(options.rho_log_min));
    if (!std::isfinite(average) || average <= 0 || (average >= 1 / span && average <= span)) {
        return 0;
    }
    const int exponent = std::ilogb(average);
    Eigen::VectorXd adjusted(response.size());
    for (Eigen::Index i = 0; i < response.size(); ++i) {
        adjusted(i) = std::scalbn(response(i), -exponent);
        if (std::scalbn(adjusted(i), exponent) != response(i)) return 0;
    }
    response.swap(adjusted);
    return exponent;
}

} // namespace

DenseIolsResult dense_iols(const double* y, const double* slopes, int n, int p,
                          double tau, const DenseIolsOptions& options) {
    if (y == nullptr || slopes == nullptr || n < 3 || p < 1 || p >= n - 1 ||
        !std::isfinite(tau) || tau < 1 || tau > 2) {
        throw std::invalid_argument("Invalid no-fixed-effect IOLS dimensions or variance");
    }
    const Eigen::Map<const Eigen::VectorXd> outcome(y, n);
    const Eigen::Map<const Eigen::MatrixXd> columns(slopes, n, p);
    const RowMatrixXd exog = columns;

    FullrecoverControlCore control;
    control.eps = options.eps;
    control.anderson_m = options.anderson_memory;
    control.stopcrit = options.stopcrit;
    control.max_phase1_iterations = options.max_phase1_iterations;
    control.max_phase2_evaluations = options.max_phase2_evaluations;
    control.rho_log_min = options.rho_log_min;
    control.rho_log_max = options.rho_log_max;
    control.interrupt = options.interrupt;
    if (!options.pilot_deltas.empty()) {
        for (const double delta : options.pilot_deltas) {
            if (!std::isfinite(delta) || delta <= 0) {
                throw std::invalid_argument("IOLS pilot shifts must be positive and finite");
            }
        }
        control.pilot_deltas = options.pilot_deltas;
    }
    validate_direct_nofe_inputs_core(outcome, exog, tau);
    validate_direct_nofe_control_core(control);

    Eigen::VectorXd kernel_outcome = outcome;
    const int exponent = condition_response(kernel_outcome, tau, options);
    const std::vector<std::int64_t> groups(static_cast<std::size_t>(n), 0);
    DirectRhoGridHookCore grid;
    grid.policy = kDirectRhoGridCurvatureStep;
    const auto state = fullrecover_invariant_core(exog, groups, 1);
    const auto result = concentrated_fullrecover_profiled_core(
        kernel_outcome, exog, groups, 1, tau, control, state,
        std::nullopt, std::nullopt, std::nullopt, true, &grid
    );
    const auto& trace = result.iteration;
    const double missing = std::numeric_limits<double>::quiet_NaN();
    double intercept = missing;
    std::vector<double> index;
    if (all_finite_vec(result.beta)) {
        Eigen::VectorXd base = exog * result.beta;
        intercept = profile_intercept_adjustment_core(outcome, base, tau);
        base.array() += intercept;
        if (all_finite_vec(base)) index.assign(base.data(), base.data() + base.size());
    }
    const std::vector<double> beta(result.beta.data(), result.beta.data() + result.beta.size());
    return {trace.converged, trace.iterations, trace.phase, trace.restarts, trace.phase1_attempts,
            trace.phase1_accepted, trace.phase1_rejected, trace.phase1_evaluations,
            trace.phase2_evaluations, exponent, intercept, trace.post_convergence_residual,
            trace.post_convergence_score_residual, trace.post_convergence_fe_score_residual,
            trace.has_final_rho ? trace.final_rho : missing,
            trace.has_final_criterion ? trace.final_criterion : missing,
            trace.failure_reason, trace.score_validation_status, beta, index};
}

} // namespace powermean
