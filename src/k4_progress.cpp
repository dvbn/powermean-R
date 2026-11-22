// Adapted from Powermean 0.1.0, BSD-3-Clause. See inst/LICENSE.K4.
// Frozen source: 7b22dcf01c6fea74176823e24a4b18df559cd386
// src/powermean/_native/iols/bind.cpp
// SHA-256 7d12b4308bbff9f4ee6f3dc3f72f77c06874462c7e762d6bc7c8ad707060e4e5.

#include "k4_core.h"

namespace powermean::k4 {

bool phase2_initialization_usable_core(
    const Eigen::Ref<const Eigen::VectorXd>& endog,
    const Eigen::Ref<const RowMatrixXd>& exog,
    const Eigen::Ref<const Eigen::VectorXd>& beta,
    const std::vector<std::int64_t>& groups,
    std::int64_t n_groups,
    double tau,
    const Eigen::Ref<const Eigen::MatrixXd>& projection,
    const FullrecoverControlCore& control,
    bool policy_grid) {
    if (!all_finite_vec(beta)) {
        return false;
    }
    std::vector<double> rho_path;
    if (policy_grid && control.rho_log_max <= control.rho_log_min) {
        // Both step policies clip to the upper rho bound in this case.
        const double alpha_min = 1.0 / (1.0 + std::exp(control.rho_log_max));
        bool clipped = false;
        rho_path = direct_step_rho_grid_core(alpha_min, alpha_min, clipped);
    } else {
        rho_path = rho_path_core(control, std::nullopt);
    }
    for (const double rho : rho_path) {
        Eigen::VectorXd beta_test(beta.size());
        if (phi_step_profiled_core(
                endog,
                exog,
                beta,
                groups,
                n_groups,
                tau,
                projection,
                1,
                rho,
                beta_test)) {
            return true;
        }
    }
    return false;
}

double max_abs_change_core(
    const Eigen::Ref<const Eigen::VectorXd>& left,
    const Eigen::Ref<const Eigen::VectorXd>& right) {
    if (left.size() == 0 && right.size() == 0) {
        return 0.0;
    }
    return (left - right).cwiseAbs().maxCoeff();
}

bool criterion_stagnated_core(
    const std::vector<double>& history,
    double eps,
    int window = 50,
    int min_iterations = 100,
    double min_relative_drop = 0.01) {
    if (history.size() < static_cast<std::size_t>(std::max(window * 2, min_iterations))) {
        return false;
    }
    const double current = history.back();
    const double previous = history[history.size() - static_cast<std::size_t>(window)];
    if (!(std::isfinite(current) && std::isfinite(previous))) {
        return false;
    }
    if (current <= eps * 10.0) {
        return false;
    }
    const double relative_drop = (previous - current) / std::max(std::abs(previous), eps);
    return relative_drop < min_relative_drop;
}

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
    const std::vector<double>* stagnation_history) {
    // runtime_v4: the direct route passes its plain-map residual history as
    // stagnation_history; every other caller keeps the criterion history.
    if (!std::isfinite(criterion)) {
        return false;
    }
    history.push_back(criterion);
    last_finite_criterion = criterion;
    best_finite_criterion =
        std::isfinite(best_finite_criterion)
            ? std::min(best_finite_criterion, criterion)
            : criterion;
    if (phase == 1 && tau > 2.0 && history.size() >= 20) {
        const double criterion_ref = history[history.size() - 20];
        if (std::isfinite(criterion_ref) && criterion_ref > eps &&
            criterion > std::max(criterion_ref * 100.0, 1e6)) {
            failure_reason = "criterion_explosion";
            return true;
        }
    }
    const int no_progress_min_iterations =
        std::min(std::max(max_phase2_evaluations / 2500, 20), 200);
    const int no_progress_window =
        std::min(std::max(no_progress_min_iterations / 2, 10), 50);
    if (phase == 1 &&
        criterion_stagnated_core(
            stagnation_history == nullptr ? history : *stagnation_history,
            eps,
            no_progress_window,
            no_progress_min_iterations,
            1e-4)) {
        failure_reason = "criterion_stagnation";
        return true;
    }
    return false;
}

int adapt_anderson_m_core(int aa_m, double kappa_val) {
    if (!std::isfinite(kappa_val)) {
        return aa_m;
    }
    if (kappa_val >= 0.95) {
        return std::max(aa_m, 20);
    }
    if (kappa_val >= 0.85) {
        return std::max(aa_m, 12);
    }
    if (kappa_val >= 0.70) {
        return std::max(aa_m, 8);
    }
    return aa_m;
}

// runtime_v3 direct-route kappa statistic (twin of the NumPy
// _windowed_contraction_rate): geometric-mean per-iteration step contraction
// over the last `window` recorded finite criteria,
// (c_last / c_{last-window})^(1/window). NaN when not computable.

double windowed_contraction_rate_core(
    const std::vector<double>& history,
    int window) {
    if (window < 1 ||
        history.size() < static_cast<std::size_t>(window) + 1) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    const double current = history.back();
    const double previous =
        history[history.size() - 1 - static_cast<std::size_t>(window)];
    if (!(std::isfinite(current) && std::isfinite(previous))) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    if (previous <= 0.0 || current < 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return std::pow(current / previous, 1.0 / static_cast<double>(window));
}

// runtime_v4 per-attempt progress-monitor fields (twin of the NumPy
// _direct_progress_trace_fields); the v3 criterion fields keep their meaning.
void set_direct_progress_trace_fields(
    Phase2AttemptTraceCore& record,
    const std::vector<double>& progress_history,
    int aa_m_start,
    int aa_m_end,
    int aa_accepted) {
    if (!progress_history.empty()) {
        record.final_progress_residual = progress_history.back();
        record.best_progress_residual =
            *std::min_element(progress_history.begin(), progress_history.end());
    }
    record.anderson_m_start = aa_m_start;
    record.anderson_m_end = aa_m_end;
    record.anderson_accepted = aa_accepted;
}

FullrecoverIterationCore not_started_iteration_core(
    const FullrecoverControlCore& control) {
    FullrecoverIterationCore iteration;
    iteration.anderson_m_initial = control.anderson_m;
    iteration.anderson_m_final = control.anderson_m;
    iteration.failure_reason = "not_started";
    return iteration;
}

void push_profiled_beta_path_core(
    std::vector<Eigen::VectorXd>& beta_path,
    const Eigen::VectorXd& beta,
    std::size_t& entries_seen,
    double& kappa_row_max,
    std::size_t min_entries) {
    beta_path.push_back(beta);
    entries_seen += 1;
    if (entries_seen >= min_entries) {
        const std::size_t row = beta_path.size() - 1;
        const Eigen::VectorXd d_now =
            (beta_path[row] - beta_path[row - 1]).cwiseAbs();
        const Eigen::VectorXd d_prev =
            (beta_path[row - 1] - beta_path[row - 2]).cwiseAbs();
        double mean_ratio = 0.0;
        int count = 0;
        for (Eigen::Index idx = 0; idx < d_now.size(); ++idx) {
            const double ratio = d_now(idx) / d_prev(idx);
            if (!std::isnan(ratio)) {
                mean_ratio += ratio;
                count += 1;
            }
        }
        if (count > 0) {
            kappa_row_max = std::max(kappa_row_max, mean_ratio / count);
        }
    }
    if (beta_path.size() > kProfiledBetaPathRetained) {
        beta_path.erase(beta_path.begin());
    }
}

// rho_grid_hook_v1 (research diagnostic, twin of the NumPy hook): opt-in
// Phase-2 rho list for the direct no-FE route. A null hook or the "legacy"
// policy keeps rho_path_core untouched; no other control flow changes.

} // namespace powermean::k4
