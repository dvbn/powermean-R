// Adapted from Powermean 0.1.0, BSD-3-Clause. See inst/LICENSE.K4.
// Frozen source: 7b22dcf01c6fea74176823e24a4b18df559cd386
// src/powermean/_native/iols/bind.cpp
// SHA-256 7d12b4308bbff9f4ee6f3dc3f72f77c06874462c7e762d6bc7c8ad707060e4e5.

#include "k4_core.h"
#include <Eigen/Cholesky>

namespace powermean::k4 {

struct EigenVectorPairRingCore {
    std::vector<Eigen::VectorXd> x;
    std::vector<Eigen::VectorXd> f;
    std::size_t capacity = 0;
    std::size_t start = 0;
    std::size_t count = 0;

    // Callers consume references immediately; do not hold x_at/f_at/f_back
    // references across another push because a full ring overwrites slots.
    void reset(std::size_t new_capacity) {
        capacity = new_capacity;
        start = 0;
        count = 0;
        x.clear();
        f.clear();
        x.reserve(capacity);
        f.reserve(capacity);
    }

    void clear() {
        start = 0;
        count = 0;
        x.clear();
        f.clear();
    }

    // runtime_v2 (direct fixed-tau route only): enlarge the ring in place when
    // the Anderson memory adapts upward, keeping the logical order of the
    // retained pairs. The NumPy reference trims its history lists to the
    // adapted aa_m + 1, so after adaptation it uses up to aa_m lags; a fixed
    // ring silently kept the initial memory. No-op unless new_capacity is larger.
    void grow(std::size_t new_capacity) {
        if (new_capacity <= capacity) {
            return;
        }
        std::vector<Eigen::VectorXd> x_ordered;
        std::vector<Eigen::VectorXd> f_ordered;
        x_ordered.reserve(new_capacity);
        f_ordered.reserve(new_capacity);
        for (std::size_t logical = 0; logical < count; ++logical) {
            x_ordered.emplace_back(x_at(logical));
            f_ordered.emplace_back(f_at(logical));
        }
        x = std::move(x_ordered);
        f = std::move(f_ordered);
        capacity = new_capacity;
        start = 0;
    }

    void push(const Eigen::VectorXd& x_value, const Eigen::VectorXd& f_value) {
        if (capacity == 0) {
            return;
        }
        if (count < capacity) {
            x.emplace_back(x_value);
            f.emplace_back(f_value);
            ++count;
            return;
        }
        x[start] = x_value;
        f[start] = f_value;
        start = (start + 1) % capacity;
    }

    std::size_t size() const {
        return count;
    }

    std::size_t physical_index(std::size_t logical_index) const {
        assert(capacity > 0);
        assert(logical_index < count);
        return (start + logical_index) % capacity;
    }

    const Eigen::VectorXd& x_at(std::size_t logical_index) const {
        return x[physical_index(logical_index)];
    }

    const Eigen::VectorXd& f_at(std::size_t logical_index) const {
        return f[physical_index(logical_index)];
    }

    const Eigen::VectorXd& f_back() const {
        assert(count > 0);
        return f_at(count - 1);
    }
};


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
    bool direct_nofe_fixed_tau,
    DirectRhoGridHookCore* rho_hook) {
    // rho_grid_hook_v1: opt-in policy is direct-route only.
    if (rho_hook != nullptr) {
        validate_direct_rho_grid_policy_core(rho_hook->policy);
        if (rho_hook->policy != kDirectRhoGridLegacy && !direct_nofe_fixed_tau) {
            throw std::invalid_argument(
                "rho_grid_policy other than 'legacy' requires direct_nofe_fixed_tau.");
        }
    }
    const bool policy_grid = direct_nofe_fixed_tau && rho_hook != nullptr &&
                             rho_hook->policy != kDirectRhoGridLegacy;
    if (direct_nofe_fixed_tau &&
        (beta_init.has_value() || !is_singleton_intercept_group(groups, n_groups))) {
        throw std::invalid_argument(
            "Direct fixed-tau Phase 1 + Phase 2 is authorized only for the "
            "noFE singleton-intercept route without beta_init");
    }
    if (!beta_init.has_value() && std::abs(tau - 2.0) > 1e-10 &&
        !direct_nofe_fixed_tau) {
        throw std::invalid_argument(
            "Non-anchor iOLS path steps require beta_init from tau=2");
    }
    const std::vector<double> delta_grid =
        beta_init.has_value() ? std::vector<double>{} : delta_grid_core(control);
    Eigen::VectorXd beta;
    Eigen::VectorXd linear;
    int n_outer = 0;
    if (beta_init.has_value()) {
        beta = *beta_init;
        linear = linear_init.has_value()
                     ? *linear_init
                     : initial_linear_profiled_core(
                           endog, exog, beta, groups, n_groups, tau);
        n_outer = 1;
    } else {
        Eigen::VectorXd log_y(endog.size());
        for (Eigen::Index row = 0; row < endog.size(); ++row) {
            log_y(row) = std::log(endog(row) + 1.0);
        }
        const Eigen::MatrixXd log_y_mat = log_y;
        const Eigen::VectorXd y_init =
            demean_profiled_core(log_y_mat, groups, n_groups).col(0);
        beta = state.projection * y_init;
        linear = initial_linear_profiled_core(endog, exog, beta, groups, n_groups, tau);
        n_outer = static_cast<int>(delta_grid.size()) + 1;
    }

    const Eigen::VectorXd beta_initial = beta;
    const Eigen::VectorXd linear_initial = linear;
    Eigen::VectorXd phase1_accepted_beta = beta;
    Eigen::VectorXd phase1_accepted_linear = linear;
    int phase1_attempts = 0;
    int phase1_accepted = 0;
    int phase1_rejected = 0;
    // runtime_v2: total Phase 1 map evaluations (plain steps plus Anderson
    // trial maps over every delta attempt); twin of the NumPy counter.
    int phase1_evaluations = 0;
    // runtime_v2 direct route: NumPy-parity kappa window and Anderson history
    // that grows with the adapted memory. Legacy callers keep their semantics.
    const std::size_t kappa_min_entries = direct_nofe_fixed_tau
                                              ? kProfiledKappaMinEntriesNumpyParity
                                              : kProfiledKappaMinEntriesLegacy;
    FullrecoverIterationCore final_diag = not_started_iteration_core(control);
    // runtime_v3 direct-route state (declared at function scope for C++
    // scoping; read and written only inside direct_nofe_fixed_tau blocks).
    bool v3_restart_from_best = false;
    double v3_best_score_norm = std::numeric_limits<double>::infinity();
    double v3_best_score_rho = std::numeric_limits<double>::quiet_NaN();
    bool v3_has_best_score_beta = false;
    Eigen::VectorXd v3_best_score_beta;
    int v3_phase2_warm_restarts = 0;
    std::vector<Phase2AttemptTraceCore> v3_phase2_attempt_trace;
    double v3_kappa_coordinate_row_max = std::numeric_limits<double>::quiet_NaN();
    double v3_kappa_monitor_last = std::numeric_limits<double>::quiet_NaN();
    double v3_kappa_monitor_max = std::numeric_limits<double>::quiet_NaN();
    int v3_kappa_monitor_checks = 0;
    int v3_attempt_evaluations_start = 0;
    std::string v3_attempt_start_policy = kDirectRestartStartPhase1;
    // runtime_v4: final verification maps that came back invalid, by
    // internal phase (0 = Phase 1 delta grid, 1 = Phase 2 rho grid).
    int v4_post_convergence_invalid_map_counts[2] = {0, 0};
    if (direct_nofe_fixed_tau) {
        v3_restart_from_best =
            control.direct_restart_start == kDirectRestartStartBestScore;
    }

    for (int outer_index = 0; outer_index < n_outer; ++outer_index) {
        if (control.interrupt) control.interrupt();
        if (!all_finite_vec(beta)) {
            beta = beta_initial;
            linear = linear_initial;
        }
        Eigen::VectorXd beta_outer_initial = beta;
        Eigen::VectorXd linear_outer_initial = linear;
        int phase = 0;
        double param = 0.0;
        double eps_current = control.eps;
        int stopcrit_current = control.stopcrit;
        std::vector<double> rho_path;
        if (outer_index < static_cast<int>(delta_grid.size())) {
            phase = 0;
            param = delta_grid[static_cast<std::size_t>(outer_index)];
            eps_current = std::max(control.eps * 10.0, 1e-5);
            stopcrit_current = std::min(control.stopcrit, control.max_phase1_iterations);
        } else {
            phase = 1;
            if (!beta_init.has_value()) {
                beta = phase1_accepted_beta;
                linear = phase1_accepted_linear;
            }
            if (policy_grid) {
                // rho_grid_hook_v1: step grid; curvature L0 at the Phase-2 pilot.
                rho_path = direct_policy_rho_path_core(
                    endog, exog, beta, groups, n_groups, tau, control, state, *rho_hook);
            } else {
                rho_path = rho_path_core(control, rho_hint);
            }
            if (direct_nofe_fixed_tau && rho_hook != nullptr) {
                rho_hook->reached = true;
                rho_hook->grid_count = static_cast<int>(rho_path.size());
            }
            if (rho_path.empty()) {
                FullrecoverConcentratedCore out;
                out.beta = beta;
                out.linear = linear;
                out.has_linear = all_finite_vec(linear);
                out.iteration = not_started_iteration_core(control);
                out.iteration.phase = "P2";
                out.iteration.failure_reason = "empty_finite_rho_path";
                out.iteration.phase1_attempts = phase1_attempts;
                out.iteration.phase1_accepted = phase1_accepted;
                out.iteration.phase1_rejected = phase1_rejected;
                out.iteration.phase1_evaluations = phase1_evaluations;
                return out;
            }
            param = rho_path.front();
            beta_outer_initial = beta;
            linear_outer_initial = linear;
        }

        int restart_index = 0;
        bool verified = false;
        int phase2_evaluations = 0;
        bool phase2_budget_exhausted = false;
        int iteration_count = 0;
        bool has_final_criterion = false;
        double final_criterion = 0.0;
        int aa_m_final = control.anderson_m;
        std::string failure_reason;
        double last_finite_criterion = std::numeric_limits<double>::quiet_NaN();
        double best_finite_criterion = std::numeric_limits<double>::quiet_NaN();
        double post_convergence_residual =
            std::numeric_limits<double>::quiet_NaN();
        double post_convergence_kappa_est =
            std::numeric_limits<double>::quiet_NaN();
        double post_convergence_score_residual =
            std::numeric_limits<double>::quiet_NaN();
        double post_convergence_fe_score_residual =
            std::numeric_limits<double>::quiet_NaN();
        double score_tol = std::numeric_limits<double>::quiet_NaN();
        std::string score_validation_status = "not_checked";
        std::vector<std::string> score_warning_codes;
        std::vector<std::string> score_warning_messages;
        bool score_exceeds_tolerance = false;
        bool fe_score_exceeds_tolerance = false;

        while (!verified) {
            int k = 1;
            double criterion = 100.0;
            std::vector<double> criterion_history;
            int kappa_bad_count = 0;
            std::vector<Eigen::VectorXd> beta_path;
            EigenVectorPairRingCore aa_hist;
            beta_path.reserve(kProfiledBetaPathRetained + 1);
            beta_path.push_back(beta);
            std::size_t beta_path_entries = 1;
            double kappa_row_max = -std::numeric_limits<double>::infinity();
            int aa_m = control.anderson_m;
            if (direct_nofe_fixed_tau && phase == 1) {
                v3_attempt_evaluations_start = phase2_evaluations;
            }
            criterion_history.reserve(
                static_cast<std::size_t>(std::max(stopcrit_current, 0)));
            if (aa_m > 0) {
                aa_hist.reset(static_cast<std::size_t>(aa_m + 1));
            }
            // runtime_v4 direct-route Phase 2 progress monitor (twin of the
            // NumPy progress_history): max-abs plain-map residual at the
            // current iterate, recorded at the same iterations as the
            // criterion. It drives stagnation, kappa and the Anderson memory
            // adaptation; the eps stop still uses `criterion`.
            std::vector<double> progress_history;
            double plain_residual_maxabs = std::numeric_limits<double>::quiet_NaN();
            int attempt_aa_m_start = aa_m;
            int attempt_aa_accepted = 0;

            while (criterion > eps_current && k < stopcrit_current &&
                   !std::isnan(criterion)) {
                if (control.interrupt && k % 64 == 0) control.interrupt();
                const Eigen::VectorXd beta0 = beta;
                if (phase == 1) {
                    if (phase2_evaluations >= control.max_phase2_evaluations) {
                        criterion = std::numeric_limits<double>::quiet_NaN();
                        failure_reason = "phase2_evaluation_budget_exhausted";
                        phase2_budget_exhausted = true;
                        break;
                    }
                    phase2_evaluations += 1;
                } else {
                    phase1_evaluations += 1;
                }
                Eigen::VectorXd beta_phi(beta.size());
                const bool valid = phi_step_profiled_core(
                    endog,
                    exog,
                    beta,
                    groups,
                    n_groups,
                    tau,
                    state.projection,
                    phase,
                    param,
                    beta_phi);
                if (!valid) {
                    criterion = std::numeric_limits<double>::quiet_NaN();
                    failure_reason = "invalid_phi_step";
                    break;
                }
                if (direct_nofe_fixed_tau && phase == 1) {
                    // runtime_v3: rho-free preconditioned score norm of the
                    // current iterate, (1+rho)*||phi(beta)-beta||_2 =
                    // ||P score_obs(beta)||_2 (twin of the NumPy tracker).
                    const Eigen::VectorXd v3_plain_step = beta_phi - beta;
                    const double v3_step_score_norm =
                        stable_l2_norm(v3_plain_step) *
                        (1.0 + param);
                    if (std::isfinite(v3_step_score_norm) &&
                        v3_step_score_norm < v3_best_score_norm) {
                        v3_best_score_norm = v3_step_score_norm;
                        v3_best_score_rho = param;
                        v3_best_score_beta = beta;
                        v3_has_best_score_beta = true;
                    }
                    // runtime_v4: same operands as the plain-step criterion,
                    // so it equals `criterion` whenever no Anderson step is
                    // accepted; no extra map evaluation.
                    plain_residual_maxabs = max_abs_change_core(beta_phi, beta);
                }

                bool accepted_aa_step = false;
                if (aa_m > 0 && k > aa_m + 1) {
                    const Eigen::VectorXd f_k = beta_phi - beta;
                    aa_hist.push(beta, f_k);
                    const int m_cur = static_cast<int>(aa_hist.size()) - 1;
                    if (m_cur >= 1) {
                        Eigen::MatrixXd f_mat(beta.size(), m_cur);
                        for (int j = 1; j <= m_cur; ++j) {
                            f_mat.col(j - 1) =
                                aa_hist.f_back() -
                                aa_hist.f_at(aa_hist.size() - static_cast<std::size_t>(j + 1));
                        }
                        Eigen::VectorXd gamma =
                            min_norm_lstsq_anderson(
                                f_mat,
                                aa_hist.f_back())
                                .col(0);
                        Eigen::VectorXd beta_aa = beta_phi;
                        for (int j = 0; j < m_cur; ++j) {
                            const std::size_t idx =
                                aa_hist.size() - static_cast<std::size_t>(j + 2);
                            const Eigen::VectorXd phi_prev =
                                aa_hist.x_at(idx) + aa_hist.f_at(idx);
                            beta_aa += gamma(j) * (phi_prev - beta_phi);
                        }
                        if (all_finite_vec(beta_aa)) {
                            if (phase == 1) {
                                if (phase2_evaluations >=
                                    control.max_phase2_evaluations) {
                                    criterion =
                                        std::numeric_limits<double>::quiet_NaN();
                                    failure_reason =
                                        "phase2_evaluation_budget_exhausted";
                                    phase2_budget_exhausted = true;
                                    break;
                                }
                                phase2_evaluations += 1;
                            } else {
                                phase1_evaluations += 1;
                            }
                            Eigen::VectorXd beta_aa_phi(beta.size());
                            const bool aa_valid = phi_step_profiled_core(
                                endog,
                                exog,
                                beta_aa,
                                groups,
                                n_groups,
                                tau,
                                state.projection,
                                phase,
                                param,
                                beta_aa_phi);
                            if (aa_valid) {
                                const Eigen::VectorXd plain_residual = beta_phi - beta;
                                const Eigen::VectorXd aa_residual =
                                    beta_aa_phi - beta_aa;
                                const double res_plain =
                                    stable_l2_norm(plain_residual);
                                const double res_aa =
                                    stable_l2_norm(aa_residual);
                                if (res_aa < res_plain * 1.1) {
                                    beta = beta_aa;
                                    push_profiled_beta_path_core(
                                        beta_path,
                                        beta,
                                        beta_path_entries,
                                        kappa_row_max,
                                        kappa_min_entries);
                                    if (k >= 2) {
                                        criterion = max_abs_change_core(beta, beta0);
                                    }
                                    accepted_aa_step = true;
                                    attempt_aa_accepted += 1;
                                }
                            }
                        }
                    }
                } else if (aa_m > 0) {
                    aa_hist.push(beta, beta_phi - beta);
                }

                if (!accepted_aa_step) {
                    beta = beta_phi;
                    push_profiled_beta_path_core(
                        beta_path,
                        beta,
                        beta_path_entries,
                        kappa_row_max,
                        kappa_min_entries);
                }
                if (k >= 2) {
                    if (!accepted_aa_step) {
                        criterion = max_abs_change_core(beta, beta0);
                    }
                    const bool direct_progress = direct_nofe_fixed_tau && phase == 1;
                    if (direct_progress && std::isfinite(plain_residual_maxabs)) {
                        progress_history.push_back(plain_residual_maxabs);
                    }
                    if (record_criterion_progress_core(
                            criterion_history,
                            criterion,
                            phase,
                            tau,
                            eps_current,
                            control.max_phase2_evaluations,
                            last_finite_criterion,
                            best_finite_criterion,
                            failure_reason,
                            direct_progress ? &progress_history : nullptr)) {
                        criterion = std::numeric_limits<double>::quiet_NaN();
                    }
                }

                if (phase == 1 && std::abs(tau - 2.0) > 1e-10 && k >= 30 &&
                    k % 10 == 0 && beta_path_entries >= kappa_min_entries) {
                    if (direct_nofe_fixed_tau) {
                        // runtime_v3: the lifetime coordinate-mean ratio
                        // maximum is kept as a diagnostic only (+inf rows
                        // included: they are the pollution it exists to show).
                        if (!std::isnan(kappa_row_max) &&
                            kappa_row_max > -std::numeric_limits<double>::infinity()) {
                            v3_kappa_coordinate_row_max =
                                std::isnan(v3_kappa_coordinate_row_max)
                                    ? kappa_row_max
                                    : std::max(v3_kappa_coordinate_row_max,
                                               kappa_row_max);
                        }
                    }
                    // runtime_v3 direct route: windowed geometric-mean
                    // contraction (twin of _windowed_contraction_rate);
                    // legacy routes keep the lifetime coordinate statistic.
                    // runtime_v4: taken over the plain-map residual history,
                    // not accepted moves.
                    const double kappa_val = direct_nofe_fixed_tau
                                                 ? windowed_contraction_rate_core(
                                                       progress_history,
                                                       kDirectKappaWindow)
                                                 : kappa_row_max;
                    if (direct_nofe_fixed_tau) {
                        v3_kappa_monitor_last = kappa_val;
                        v3_kappa_monitor_checks += 1;
                        if (std::isfinite(kappa_val)) {
                            v3_kappa_monitor_max =
                                std::isnan(v3_kappa_monitor_max)
                                    ? kappa_val
                                    : std::max(v3_kappa_monitor_max, kappa_val);
                        }
                    }
                    const int aa_m_new = adapt_anderson_m_core(aa_m, kappa_val);
                    if (aa_m_new > aa_m) {
                        aa_m = aa_m_new;
                        kappa_bad_count = 0;
                        if (direct_nofe_fixed_tau && aa_m > 0) {
                            // NumPy parity: the history now holds up to aa_m lags.
                            aa_hist.grow(static_cast<std::size_t>(aa_m + 1));
                        }
                    } else if (std::isfinite(kappa_val) && kappa_val >= 1.0) {
                        kappa_bad_count += 1;
                    } else {
                        kappa_bad_count = 0;
                    }
                    if (kappa_bad_count >= 3) {
                        criterion = std::numeric_limits<double>::quiet_NaN();
                        failure_reason = "noncontractive";
                    }
                }

                k += 1;
                iteration_count = k;

                if (std::isnan(criterion)) {
                    restart_index += 1;
                    if (phase == 1 && !phase2_budget_exhausted &&
                        restart_index < static_cast<int>(rho_path.size())) {
                        if (direct_nofe_fixed_tau) {
                            Phase2AttemptTraceCore v3_record;
                            // restart_index was already advanced: the attempt
                            // that just ended is the previous one.
                            v3_record.attempt = restart_index - 1;
                            v3_record.rho = param;
                            v3_record.iterations = k;
                            v3_record.evaluations =
                                phase2_evaluations - v3_attempt_evaluations_start;
                            v3_record.failure_reason = failure_reason;
                            v3_record.has_final_criterion =
                                std::isfinite(last_finite_criterion);
                            v3_record.final_criterion = last_finite_criterion;
                            v3_record.best_criterion =
                                criterion_history.empty()
                                    ? std::numeric_limits<double>::quiet_NaN()
                                    : *std::min_element(criterion_history.begin(),
                                                        criterion_history.end());
                            v3_record.kappa_last = v3_kappa_monitor_last;
                            v3_record.start_policy = v3_attempt_start_policy;
                            set_direct_progress_trace_fields(
                                v3_record,
                                progress_history,
                                attempt_aa_m_start,
                                aa_m,
                                attempt_aa_accepted);
                            v3_phase2_attempt_trace.push_back(v3_record);
                        }
                        param = rho_path[static_cast<std::size_t>(restart_index)];
                        criterion = 100.0;
                        beta = beta_outer_initial;
                        if (direct_nofe_fixed_tau) {
                            // runtime_v3 restart start policy.
                            v3_attempt_start_policy = kDirectRestartStartPhase1;
                            if (v3_restart_from_best && v3_has_best_score_beta) {
                                beta = v3_best_score_beta;
                                v3_phase2_warm_restarts += 1;
                                v3_attempt_start_policy = kDirectRestartStartBestScore;
                            }
                            v3_attempt_evaluations_start = phase2_evaluations;
                            v3_kappa_monitor_last = std::numeric_limits<double>::quiet_NaN();
                        }
                        k = 1;
                        criterion_history.clear();
                        kappa_bad_count = 0;
                        beta_path.clear();
                        beta_path.push_back(beta);
                        beta_path_entries = 1;
                        kappa_row_max = -std::numeric_limits<double>::infinity();
                        aa_hist.clear();
                        // runtime_v4: aa_m itself is still carried over (v3).
                        progress_history.clear();
                        plain_residual_maxabs = std::numeric_limits<double>::quiet_NaN();
                        attempt_aa_m_start = aa_m;
                        attempt_aa_accepted = 0;
                    } else if (phase == 0) {
                        break;
                    }
                }
            }

            bool converged = !std::isnan(criterion) && criterion <= eps_current;
            has_final_criterion = !std::isnan(criterion);
            final_criterion = criterion;
            if (converged) {
                Eigen::VectorXd beta_check(beta.size());
                bool valid = phi_step_profiled_core(
                    endog,
                    exog,
                    beta,
                    groups,
                    n_groups,
                    tau,
                    state.projection,
                    phase,
                    param,
                    beta_check);
                if (valid) {
                    const double residual = max_abs_change_core(beta_check, beta);
                    post_convergence_residual = residual;
                    if (residual > eps_current * 10.0) {
                        converged = false;
                        failure_reason = "post_convergence_residual";
                    }
                } else if (direct_nofe_fixed_tau) {
                    // runtime_v4: an invalid final verification map clears
                    // native convergence (v3 kept it). No residual exists, so
                    // an earlier attempt's post-check values are not reported.
                    converged = false;
                    failure_reason = "post_convergence_invalid_map";
                    v4_post_convergence_invalid_map_counts[phase] += 1;
                    post_convergence_residual = std::numeric_limits<double>::quiet_NaN();
                    post_convergence_kappa_est = std::numeric_limits<double>::quiet_NaN();
                }
                if (converged && phase == 1) {
                    const double score_residual = profiled_slope_score_residual_core(
                        endog,
                        exog,
                        beta,
                        groups,
                        n_groups,
                        state.x_demeaned,
                        tau);
                    const double fe_score_residual = profiled_fe_score_residual_core(
                        endog,
                        exog,
                        beta,
                        groups,
                        n_groups,
                        tau);
                    score_tol = std::max(eps_current * 10.0, 1e-8);
                    post_convergence_score_residual = score_residual;
                    score_exceeds_tolerance =
                        !std::isfinite(score_residual) || score_residual > score_tol;
                    if (score_exceeds_tolerance) {
                        score_warning_codes.push_back(
                            "post_convergence_score_residual");
                        score_warning_messages.push_back(
                            "Slope score residual exceeds tolerance at this "
                            "tau-path state.");
                        score_validation_status = "diagnostic_warning";
                    } else if (score_validation_status == "not_checked") {
                        score_validation_status = "passed";
                    }
                    post_convergence_fe_score_residual = fe_score_residual;
                    fe_score_exceeds_tolerance =
                        !std::isfinite(fe_score_residual) ||
                        fe_score_residual > score_tol;
                    if (fe_score_exceeds_tolerance) {
                        score_warning_codes.push_back(
                            "post_convergence_fe_score_residual");
                        score_warning_messages.push_back(
                            "Recovered-FE score residual exceeds tolerance at "
                            "this tau-path state.");
                        score_validation_status = "diagnostic_warning";
                    } else if (score_validation_status == "not_checked") {
                        score_validation_status = "passed";
                    }
                }
                if (converged && beta_path.size() >= 4 &&
                    !(phase == 1 && std::abs(tau - 2.0) <= 1e-10)) {
                    const double d0 =
                        max_abs_change_core(beta_path[beta_path.size() - 1],
                                            beta_path[beta_path.size() - 2]);
                    const double d1 =
                        max_abs_change_core(beta_path[beta_path.size() - 2],
                                            beta_path[beta_path.size() - 3]);
                    const double d2 =
                        max_abs_change_core(beta_path[beta_path.size() - 3],
                                            beta_path[beta_path.size() - 4]);
                    if (d1 > 1e-15 && d2 > 1e-15) {
                        const double kappa_est = d0 / d1;
                        post_convergence_kappa_est = kappa_est;
                        if (kappa_est >= 1.0) {
                            converged = false;
                            failure_reason = "post_convergence_noncontractive";
                        }
                    }
                }
            }

            if (direct_nofe_fixed_tau && phase == 1) {
                Phase2AttemptTraceCore v3_record;
                // restart_index can equal rho_path.size() when the last
                // attempt was killed in-loop with no rho left.
                v3_record.attempt =
                    std::min(restart_index, static_cast<int>(rho_path.size()) - 1);
                v3_record.rho = param;
                // k (not iteration_count) so an attempt that ended on its
                // first map does not inherit the previous count.
                v3_record.iterations = k;
                v3_record.evaluations =
                    phase2_evaluations - v3_attempt_evaluations_start;
                // a finite criterion above eps here means the per-rho cap
                // (stopcrit) ended the attempt; failure_reason is per phase
                // and would be stale.
                v3_record.failure_reason =
                    converged ? std::string()
                              : ((!std::isnan(criterion) && criterion > eps_current)
                                     ? std::string("stopcrit_reached")
                                     : failure_reason);
                v3_record.has_final_criterion = has_final_criterion;
                v3_record.final_criterion = final_criterion;
                v3_record.best_criterion =
                    criterion_history.empty()
                        ? std::numeric_limits<double>::quiet_NaN()
                        : *std::min_element(criterion_history.begin(),
                                            criterion_history.end());
                v3_record.kappa_last = v3_kappa_monitor_last;
                v3_record.start_policy = v3_attempt_start_policy;
                set_direct_progress_trace_fields(
                    v3_record,
                    progress_history,
                    attempt_aa_m_start,
                    aa_m,
                    attempt_aa_accepted);
                v3_phase2_attempt_trace.push_back(v3_record);
            }
            if (!converged && phase == 1 && !phase2_budget_exhausted &&
                restart_index < static_cast<int>(rho_path.size()) - 1) {
                restart_index = std::min(
                    restart_index + std::max(restart_index, 5),
                    static_cast<int>(rho_path.size()) - 1);
                param = rho_path[static_cast<std::size_t>(restart_index)];
                beta = beta_outer_initial;
                if (direct_nofe_fixed_tau) {
                    // runtime_v3 restart start policy.
                    v3_attempt_start_policy = kDirectRestartStartPhase1;
                    if (v3_restart_from_best && v3_has_best_score_beta) {
                        beta = v3_best_score_beta;
                        v3_phase2_warm_restarts += 1;
                        v3_attempt_start_policy = kDirectRestartStartBestScore;
                    }
                    v3_kappa_monitor_last = std::numeric_limits<double>::quiet_NaN();
                }
            } else {
                verified = true;
            }
            aa_m_final = aa_m;

            if (phase == 0) {
                phase1_attempts += 1;
                const bool candidate_usable =
                    converged &&
                    phase2_initialization_usable_core(
                        endog,
                        exog,
                        beta,
                        groups,
                        n_groups,
                        tau,
                        state.projection,
                        control,
                        policy_grid);
                if (candidate_usable) {
                    phase1_accepted += 1;
                    phase1_accepted_beta = beta;
                    phase1_accepted_linear = initial_linear_profiled_core(
                        endog, exog, beta, groups, n_groups, tau);
                } else {
                    phase1_rejected += 1;
                    beta = phase1_accepted_beta;
                    if (converged && failure_reason.empty()) {
                        failure_reason = "phase1_candidate_unusable";
                    }
                }
                final_diag.converged = candidate_usable;
            } else {
                if (!converged && !direct_nofe_fixed_tau) {
                    // Legacy tau-path semantics: hand the anchor/previous
                    // state back to the caller. The runtime_v2 direct route
                    // keeps the final iterate of the last rho attempt so that
                    // the external exact-score audit can still accept it.
                    beta = beta_outer_initial;
                }
                if (direct_nofe_fixed_tau) {
                    final_diag.nonconverged_iterate_policy =
                        "last_attempt_final_iterate";
                    if (v3_restart_from_best) {
                        final_diag.nonconverged_iterate_policy =
                            kDirectRestartStartBestScore;
                        if (!converged && v3_has_best_score_beta) {
                            // runtime_v3: hand back the smallest
                            // preconditioned-score iterate of the whole
                            // Phase 2; the external audit decides.
                            beta = v3_best_score_beta;
                        }
                    }
                    final_diag.direct_restart_start = control.direct_restart_start;
                    final_diag.kappa_coordinate_row_max = v3_kappa_coordinate_row_max;
                    final_diag.kappa_monitor_last = v3_kappa_monitor_last;
                    final_diag.kappa_monitor_max = v3_kappa_monitor_max;
                    final_diag.kappa_monitor_checks = v3_kappa_monitor_checks;
                    final_diag.best_score_norm = v3_best_score_norm;
                    final_diag.best_score_rho = v3_best_score_rho;
                    final_diag.phase2_warm_restarts = v3_phase2_warm_restarts;
                    final_diag.phase2_attempts =
                        static_cast<int>(v3_phase2_attempt_trace.size());
                    final_diag.phase2_attempt_trace = v3_phase2_attempt_trace;
                    final_diag.post_convergence_invalid_map_count =
                        v4_post_convergence_invalid_map_counts[1];
                    final_diag.phase1_post_convergence_invalid_map_count =
                        v4_post_convergence_invalid_map_counts[0];
                }
                final_diag.converged = converged;
            }

            final_diag.phase = phase == 0 ? "P1" : "P2";
            final_diag.iterations = iteration_count;
            final_diag.restarts = restart_index;
            final_diag.has_final_rho = phase != 0;
            final_diag.final_rho = param;
            final_diag.rho_path.clear();
            if (phase != 0) {
                const int upto =
                    std::min(restart_index + 1, static_cast<int>(rho_path.size()));
                for (int idx = 0; idx < upto; ++idx) {
                    final_diag.rho_path.push_back(rho_path[static_cast<std::size_t>(idx)]);
                }
            }
            final_diag.has_final_criterion = has_final_criterion;
            final_diag.final_criterion = final_criterion;
            final_diag.anderson_m_initial = control.anderson_m;
            final_diag.anderson_m_final = aa_m_final;
            final_diag.failure_reason =
                final_diag.converged
                    ? ""
                    : (!failure_reason.empty() ? failure_reason : "not_converged");
            final_diag.phase1_attempts = phase1_attempts;
            final_diag.phase1_accepted = phase1_accepted;
            final_diag.phase1_rejected = phase1_rejected;
            final_diag.phase1_evaluations = phase1_evaluations;
            final_diag.phase2_evaluations = phase2_evaluations;
            final_diag.max_phase2_evaluations = control.max_phase2_evaluations;
            final_diag.last_finite_criterion = last_finite_criterion;
            final_diag.best_finite_criterion = best_finite_criterion;
            final_diag.fe_state_policy = "analytical_profile";
            final_diag.anderson_policy = "beta_anderson";
            final_diag.post_convergence_residual = post_convergence_residual;
            final_diag.post_convergence_kappa_est = post_convergence_kappa_est;
            final_diag.post_convergence_score_residual =
                post_convergence_score_residual;
            final_diag.post_convergence_fe_score_residual =
                post_convergence_fe_score_residual;
            final_diag.score_tolerance = score_tol;
            final_diag.score_validation_policy = "diagnostic";
            final_diag.score_validation_status = score_validation_status;
            final_diag.score_warning_codes = score_warning_codes;
            final_diag.score_warning_messages = score_warning_messages;
            final_diag.post_convergence_score_exceeds_tolerance =
                score_exceeds_tolerance;
            final_diag.post_convergence_fe_score_exceeds_tolerance =
                fe_score_exceeds_tolerance;
        }
        if (phase == 1) {
            break;
        }
    }

    FullrecoverConcentratedCore out;
    out.beta = beta;
    out.linear = initial_linear_profiled_core(endog, exog, beta, groups, n_groups, tau);
    out.has_linear = all_finite_vec(out.linear);
    out.iteration = final_diag;
    return out;
}


// Direct fixed-tau noFE route, twin of the NumPy
// _solve_nofe_direct_fixed_tau. Messages and check order match the NumPy
// validators so both backends reject the same inputs the same way.

} // namespace powermean::k4
