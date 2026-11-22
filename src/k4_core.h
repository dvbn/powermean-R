// Adapted from Powermean 0.1.0, BSD-3-Clause. See inst/LICENSE.K4.
// Frozen source: 7b22dcf01c6fea74176823e24a4b18df559cd386
// src/powermean/_native/iols/bind.cpp
// SHA-256 7d12b4308bbff9f4ee6f3dc3f72f77c06874462c7e762d6bc7c8ad707060e4e5.

#ifndef POWERMEAN_K4_CORE_H
#define POWERMEAN_K4_CORE_H

#if defined(__FAST_MATH__) || (defined(__FINITE_MATH_ONLY__) && __FINITE_MATH_ONLY__ > 0)
#error K4 requires finite-value checks and the frozen floating-point evaluation order.
#endif
#if defined(__clang__)
#pragma clang fp contract(off)
#elif defined(__GNUC__)
#pragma GCC optimize ("fp-contract=off")
#endif

#include <Eigen/Core>

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <limits>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace powermean::k4 {

using RowMatrixXd =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;

struct FullrecoverControlCore {
    void (*interrupt)() = nullptr;
    int anderson_m = 5;
    double eps = 1e-6;
    int stopcrit = 10000;
    int tau_path_steps = 20;
    bool has_min_dtau = false;
    double min_dtau = 0.0;
    int max_path_evaluations = 10000;
    int max_phase1_iterations = 10000;
    int max_phase2_evaluations = 50000;
    int fe_demean_maxiter = 10000;
    double fe_demean_tol = 1e-10;
    double delta_log_min = -5.0;
    double delta_log_max = 5.0;
    int delta_grid_size = 5;
    std::optional<std::vector<double>> pilot_deltas;
    double rho_log_min = -7.0;
    double rho_log_max = 14.0;
    int rho_grid_size = 150;
    // runtime_v3: Phase 2 restart start policy, read ONLY by the direct no-FE
    // fixed-tau route ("phase1_start" = runtime_v2 behaviour,
    // "best_score_iterate" = restart from / return the best preconditioned
    // score iterate). Legacy callers never read it.
    std::string direct_restart_start = "phase1_start";
    bool research_tiered_fe_tol = false;
    double research_inner_fe_demean_tol = 1e-8;
    bool speed_fullrecover_joint_anderson = false;
    bool speed_fullrecover_accept_score_warning_path_steps = false;
    bool speed_small_manyway_exact_demean = false;
    bool speed_compact_fixed_tau_details = false;
    int speed_small_manyway_exact_demean_max_rows = 1200;
    int speed_small_manyway_exact_demean_max_levels = 512;
};

// runtime_v3: one Phase 2 rho attempt of the direct no-FE fixed-tau route
// (twin of the NumPy phase2_attempt_trace entries).
struct Phase2AttemptTraceCore {
    int attempt = 0;
    double rho = 0.0;
    int iterations = 0;
    int evaluations = 0;
    std::string failure_reason;  // empty = converged
    bool has_final_criterion = false;
    double final_criterion = 0.0;
    double best_criterion = std::numeric_limits<double>::quiet_NaN();
    double kappa_last = std::numeric_limits<double>::quiet_NaN();
    std::string start_policy = "phase1_start";
    // runtime_v4 progress-monitor fields (twins of the NumPy trace keys).
    double final_progress_residual = std::numeric_limits<double>::quiet_NaN();
    double best_progress_residual = std::numeric_limits<double>::quiet_NaN();
    int anderson_m_start = 0;
    int anderson_m_end = 0;
    int anderson_accepted = 0;
};

struct FullrecoverIterationCore {
    bool converged = false;
    std::string phase = "not_started";
    int iterations = 0;
    int restarts = 0;
    bool has_final_rho = false;
    double final_rho = 0.0;
    std::vector<double> rho_path;
    bool has_final_criterion = false;
    double final_criterion = 0.0;
    int anderson_m_initial = 5;
    int anderson_m_final = 5;
    std::string failure_reason;
    int phase1_attempts = 0;
    int phase1_accepted = 0;
    int phase1_rejected = 0;
    int phase1_evaluations = 0;  // runtime_v2: profiled core only
    int phase2_evaluations = 0;
    int max_phase2_evaluations = 0;
    double last_finite_criterion = std::numeric_limits<double>::quiet_NaN();
    double best_finite_criterion = std::numeric_limits<double>::quiet_NaN();
    std::string rho_policy;
    double fullcancel_initial_scale = std::numeric_limits<double>::quiet_NaN();
    int rho_candidate_count = 0;
    bool has_rho_hint = false;
    double rho_hint = std::numeric_limits<double>::quiet_NaN();
    double rho_log_min = std::numeric_limits<double>::quiet_NaN();
    double rho_log_max = std::numeric_limits<double>::quiet_NaN();
    int rho_grid_size = 0;
    bool has_rho_min = false;
    double rho_min = std::numeric_limits<double>::quiet_NaN();
    bool has_rho_max = false;
    double rho_max = std::numeric_limits<double>::quiet_NaN();
    std::string algorithm;
    std::string index_policy;
    std::string fe_policy;
    std::string fe_state_policy;
    std::string anderson_policy;
    double qhat_last = std::numeric_limits<double>::quiet_NaN();
    int qhat_invalid_count = 0;
    double kappa_monitor_last = std::numeric_limits<double>::quiet_NaN();
    double kappa_monitor_max = std::numeric_limits<double>::quiet_NaN();
    int kappa_monitor_checks = 0;
    double post_convergence_residual = std::numeric_limits<double>::quiet_NaN();
    double post_convergence_kappa_est = std::numeric_limits<double>::quiet_NaN();
    double post_convergence_score_residual = std::numeric_limits<double>::quiet_NaN();
    double post_convergence_fe_score_residual = std::numeric_limits<double>::quiet_NaN();
    double score_tolerance = std::numeric_limits<double>::quiet_NaN();
    std::string score_validation_policy = "diagnostic";
    std::string score_validation_status = "not_checked";
    std::vector<std::string> score_warning_codes;
    std::vector<std::string> score_warning_messages;
    bool post_convergence_score_exceeds_tolerance = false;
    bool post_convergence_fe_score_exceeds_tolerance = false;
    int joint_hdfe_phi_calls = 0;
    int joint_hdfe_phi_demean_sweeps = 0;
    bool research_tiered_fe_tol = false;
    double research_inner_fe_demean_tol =
        std::numeric_limits<double>::quiet_NaN();
    int phase2_rho_restart_events = 0;
    int phase2_rho_attempts = 0;
    // runtime_v3 direct no-FE fixed-tau route only (profiled core with
    // direct_nofe_fixed_tau = true); defaults are untouched on legacy routes.
    std::string direct_restart_start;
    std::string nonconverged_iterate_policy;
    double kappa_coordinate_row_max = std::numeric_limits<double>::quiet_NaN();
    double best_score_norm = std::numeric_limits<double>::infinity();
    double best_score_rho = std::numeric_limits<double>::quiet_NaN();
    int phase2_warm_restarts = 0;
    int phase2_attempts = 0;
    std::vector<Phase2AttemptTraceCore> phase2_attempt_trace;
    int post_convergence_invalid_map_count = 0;  // runtime_v4, Phase 2
    int phase1_post_convergence_invalid_map_count = 0;  // runtime_v4, Phase 1
};

struct FullrecoverConcentratedCore {
    Eigen::VectorXd beta;
    Eigen::VectorXd linear;
    bool has_linear = false;
    FullrecoverIterationCore iteration;
};

struct FullrecoverInvariantCore {
    Eigen::MatrixXd x_demeaned;
    Eigen::MatrixXd projection;
    Eigen::MatrixXd manyway_q_basis;
};

inline constexpr const char* kDirectRestartStartPhase1 = "phase1_start";
inline constexpr const char* kDirectRestartStartBestScore = "best_score_iterate";
// runtime_v4: renamed because the windowed rate now reads the plain-map
// residual history instead of accepted moves.
inline constexpr const char* kDirectKappaMonitorPolicy = "windowed_plain_residual_ratio_10";
inline constexpr const char* kDirectProgressMonitorPolicy =
    "plain_map_residual_maxabs_star_at_current_iterate";
inline constexpr const char* kDirectPostCheckPolicy =
    "invalid_verification_map_clears_convergence";
inline constexpr const char* kDirectNofeRuntimeVersion = "direct_nofe_runtime_v4_2026_10_03";

inline constexpr int kDirectKappaWindow = 10;
inline constexpr std::size_t kProfiledBetaPathRetained = 4;
inline constexpr std::size_t kProfiledKappaMinEntriesLegacy = 3;
inline constexpr std::size_t kProfiledKappaMinEntriesNumpyParity = 4;

inline constexpr const char* kDirectRhoGridHookVersion = "rho_grid_hook_v1_2026_10_05";
inline constexpr const char* kDirectRhoGridLegacy = "legacy";
inline constexpr const char* kDirectRhoGridStep = "step_grid";
inline constexpr const char* kDirectRhoGridCurvatureStep = "curvature_step_grid";

struct DirectRhoGridHookCore {
    std::string policy = kDirectRhoGridLegacy;
    // Outputs, filled once at the Phase-2 handoff.
    bool reached = false;
    int grid_count = 0;
    std::vector<double> values;
    double alpha0 = std::numeric_limits<double>::quiet_NaN();
    double alpha_min = std::numeric_limits<double>::quiet_NaN();
    bool clipped = false;
    bool has_curvature = false;
    double curvature_L0 = std::numeric_limits<double>::quiet_NaN();
    double curvature_lambda_min = std::numeric_limits<double>::quiet_NaN();
    std::int64_t curvature_wall_ns = 0;
};



double stable_l2_norm(
    const Eigen::Ref<const Eigen::VectorXd>& values);

Eigen::MatrixXd min_norm_lstsq(
    const Eigen::Ref<const Eigen::MatrixXd>& a,
    const Eigen::Ref<const Eigen::MatrixXd>& b);

Eigen::MatrixXd min_norm_lstsq_anderson(
    const Eigen::Ref<const Eigen::MatrixXd>& a,
    const Eigen::Ref<const Eigen::MatrixXd>& b);

bool all_finite_vec(const Eigen::Ref<const Eigen::VectorXd>& values);

double exp_guard_scalar(double value);

bool is_singleton_intercept_group(
    const std::vector<std::int64_t>& groups,
    std::int64_t n_groups);

double profile_intercept_adjustment_core(
    const Eigen::Ref<const Eigen::VectorXd>& endog,
    const Eigen::Ref<const Eigen::VectorXd>& linear,
    double tau);

Eigen::MatrixXd demean_profiled_core(
    const Eigen::Ref<const Eigen::MatrixXd>& values,
    const std::vector<std::int64_t>& groups,
    std::int64_t n_groups);

Eigen::VectorXd initial_linear_profiled_core(
    const Eigen::Ref<const Eigen::VectorXd>& endog,
    const Eigen::Ref<const RowMatrixXd>& exog,
    const Eigen::Ref<const Eigen::VectorXd>& beta,
    const std::vector<std::int64_t>& groups,
    std::int64_t n_groups,
    double tau);

FullrecoverInvariantCore fullrecover_invariant_core(
    const Eigen::Ref<const RowMatrixXd>& exog,
    const std::vector<std::int64_t>& groups,
    std::int64_t n_groups);

bool phi_step_profiled_core(
    const Eigen::Ref<const Eigen::VectorXd>& endog,
    const Eigen::Ref<const RowMatrixXd>& exog,
    const Eigen::Ref<const Eigen::VectorXd>& beta,
    const std::vector<std::int64_t>& groups,
    std::int64_t n_groups,
    double tau,
    const Eigen::Ref<const Eigen::MatrixXd>& projection,
    int phase,
    double param,
    Eigen::VectorXd& beta_new);

double profiled_slope_score_residual_core(
    const Eigen::Ref<const Eigen::VectorXd>& endog,
    const Eigen::Ref<const RowMatrixXd>& exog,
    const Eigen::Ref<const Eigen::VectorXd>& beta,
    const std::vector<std::int64_t>& groups,
    std::int64_t n_groups,
    const Eigen::Ref<const Eigen::MatrixXd>& x_demeaned,
    double tau);

double profiled_fe_score_residual_core(
    const Eigen::Ref<const Eigen::VectorXd>& endog,
    const Eigen::Ref<const RowMatrixXd>& exog,
    const Eigen::Ref<const Eigen::VectorXd>& beta,
    const std::vector<std::int64_t>& groups,
    std::int64_t n_groups,
    double tau);

std::vector<double> rho_path_core(
    const FullrecoverControlCore& control,
    std::optional<double> rho_hint);

std::vector<double> delta_grid_core(const FullrecoverControlCore& control);

std::vector<double> direct_step_rho_grid_core(double alpha0, double alpha_min, bool& clipped);

bool phase2_initialization_usable_core(
    const Eigen::Ref<const Eigen::VectorXd>& endog,
    const Eigen::Ref<const RowMatrixXd>& exog,
    const Eigen::Ref<const Eigen::VectorXd>& beta,
    const std::vector<std::int64_t>& groups,
    std::int64_t n_groups,
    double tau,
    const Eigen::Ref<const Eigen::MatrixXd>& projection,
    const FullrecoverControlCore& control,
    bool policy_grid);

double max_abs_change_core(
    const Eigen::Ref<const Eigen::VectorXd>& left,
    const Eigen::Ref<const Eigen::VectorXd>& right);

bool record_criterion_progress_core(
    std::vector<double>& history,
    double criterion,
    int phase,
    double tau,
    double eps,
    int max_phase2_evaluations,
    double& last_finite_criterion,
    double& best_finite_criterion,
    std::string& failure_reason,
    const std::vector<double>* stagnation_history = nullptr);

int adapt_anderson_m_core(int aa_m, double kappa_val);

double windowed_contraction_rate_core(
    const std::vector<double>& history,
    int window);

void set_direct_progress_trace_fields(
    Phase2AttemptTraceCore& record,
    const std::vector<double>& progress_history,
    int aa_m_start,
    int aa_m_end,
    int aa_accepted);

FullrecoverIterationCore not_started_iteration_core(
    const FullrecoverControlCore& control);

void push_profiled_beta_path_core(
    std::vector<Eigen::VectorXd>& beta_path,
    const Eigen::VectorXd& beta,
    std::size_t& entries_seen,
    double& kappa_row_max,
    std::size_t min_entries = kProfiledKappaMinEntriesLegacy);

void validate_direct_rho_grid_policy_core(const std::string& policy);

std::vector<double> direct_policy_rho_path_core(
    const Eigen::Ref<const Eigen::VectorXd>& endog,
    const Eigen::Ref<const RowMatrixXd>& exog,
    const Eigen::Ref<const Eigen::VectorXd>& beta,
    const std::vector<std::int64_t>& groups,
    std::int64_t n_groups,
    double tau,
    const FullrecoverControlCore& control,
    const FullrecoverInvariantCore& state,
    DirectRhoGridHookCore& hook);

FullrecoverConcentratedCore concentrated_fullrecover_profiled_core(
    const Eigen::Ref<const Eigen::VectorXd>& endog,
    const Eigen::Ref<const RowMatrixXd>& exog,
    const std::vector<std::int64_t>& groups,
    std::int64_t n_groups,
    double tau,
    const FullrecoverControlCore& control,
    const FullrecoverInvariantCore& state,
    std::optional<Eigen::VectorXd> beta_init,
    std::optional<Eigen::VectorXd> linear_init,
    std::optional<double> rho_hint,
    bool direct_nofe_fixed_tau = false,
    DirectRhoGridHookCore* rho_hook = nullptr);

} // namespace powermean::k4
#endif
