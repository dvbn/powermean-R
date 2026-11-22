// Adapted from Powermean 0.1.0, BSD-3-Clause. See inst/LICENSE.K4.
// Frozen source: 7b22dcf01c6fea74176823e24a4b18df559cd386
// src/powermean/_native/iols/bind.cpp
// SHA-256 7d12b4308bbff9f4ee6f3dc3f72f77c06874462c7e762d6bc7c8ad707060e4e5.

#include "k4_core.h"

namespace powermean::k4 {

constexpr double kFloat64Max = std::numeric_limits<double>::max();
const double kFloat64LogMax = std::log(kFloat64Max);


struct DemeanWorkspaceCore {
    Eigen::MatrixXd scratch;
    Eigen::MatrixXd previous;
    Eigen::MatrixXd sums;
    Eigen::VectorXd denom;
};

void demean_oneway_matrix_into_core(
    const Eigen::Ref<const Eigen::MatrixXd>& values,
    const std::vector<std::int64_t>& groups,
    std::int64_t n_groups,
    const Eigen::VectorXd* weights,
    Eigen::MatrixXd& out,
    DemeanWorkspaceCore& workspace) {
    workspace.sums.resize(static_cast<Eigen::Index>(n_groups), values.cols());
    workspace.sums.setZero();
    workspace.denom.resize(static_cast<Eigen::Index>(n_groups));
    workspace.denom.setZero();
    for (Eigen::Index row = 0; row < values.rows(); ++row) {
        const auto group = static_cast<Eigen::Index>(
            groups[static_cast<std::size_t>(row)]);
        const double weight = weights != nullptr ? (*weights)(row) : 1.0;
        workspace.sums.row(group) += weight * values.row(row);
        workspace.denom(group) += weight;
    }

    for (Eigen::Index group = 0; group < workspace.sums.rows(); ++group) {
        if (weights != nullptr) {
            workspace.sums.row(group) /=
                std::max(workspace.denom(group), 1.0e-300);
        } else {
            const double divisor = std::max(workspace.denom(group), 1.0);
            workspace.sums.row(group) /= divisor;
        }
    }

    out.resize(values.rows(), values.cols());
    for (Eigen::Index row = 0; row < values.rows(); ++row) {
        out.row(row) = values.row(row) -
                       workspace.sums.row(groups[static_cast<std::size_t>(row)]);
    }
}

Eigen::MatrixXd demean_oneway_matrix(
    const Eigen::Ref<const Eigen::MatrixXd>& values,
    const std::vector<std::int64_t>& groups,
    std::int64_t n_groups,
    const std::optional<Eigen::VectorXd>& weights) {
    Eigen::MatrixXd out;
    DemeanWorkspaceCore workspace;
    const Eigen::VectorXd* weights_ptr = weights.has_value() ? &(*weights) : nullptr;
    demean_oneway_matrix_into_core(
        values,
        groups,
        n_groups,
        weights_ptr,
        out,
        workspace);
    return out;
}


bool all_finite_vec(const Eigen::Ref<const Eigen::VectorXd>& values) {
    for (Eigen::Index idx = 0; idx < values.size(); ++idx) {
        if (!std::isfinite(values(idx))) {
            return false;
        }
    }
    return true;
}

double exp_guard_scalar(double value) {
    if (std::isinf(value) && value > 0.0) {
        return kFloat64Max;
    }
    if (value > kFloat64LogMax) {
        return kFloat64Max;
    }
    return std::exp(value);
}

double logsumexp_all_core(const Eigen::Ref<const Eigen::VectorXd>& values) {
    if (values.size() == 0) {
        return -std::numeric_limits<double>::infinity();
    }
    bool has_posinf = false;
    bool has_finite = false;
    double max_value = -std::numeric_limits<double>::infinity();
    for (Eigen::Index idx = 0; idx < values.size(); ++idx) {
        const double value = values(idx);
        if (std::isinf(value) && value > 0.0) {
            has_posinf = true;
        }
        if (std::isfinite(value)) {
            has_finite = true;
            max_value = std::max(max_value, value);
        }
    }
    if (has_posinf) {
        return std::numeric_limits<double>::infinity();
    }
    if (!has_finite) {
        return -std::numeric_limits<double>::infinity();
    }
    double sum = 0.0;
    for (Eigen::Index idx = 0; idx < values.size(); ++idx) {
        const double value = values(idx);
        if (std::isfinite(value)) {
            sum += std::exp(value - max_value);
        }
    }
    return max_value + std::log(sum);
}

bool is_singleton_intercept_group(
    const std::vector<std::int64_t>& groups,
    std::int64_t n_groups) {
    if (n_groups != 1) {
        return false;
    }
    for (const std::int64_t group : groups) {
        if (group != 0) {
            return false;
        }
    }
    return true;
}

double profile_intercept_adjustment_core(
    const Eigen::Ref<const Eigen::VectorXd>& endog,
    const Eigen::Ref<const Eigen::VectorXd>& linear,
    double tau) {
    Eigen::VectorXd numerator =
        Eigen::VectorXd::Constant(linear.size(), -std::numeric_limits<double>::infinity());
    Eigen::VectorXd denominator(linear.size());
    for (Eigen::Index row = 0; row < linear.size(); ++row) {
        if (endog(row) > 0.0) {
            numerator(row) = std::log(endog(row)) + (1.0 - tau) * linear(row);
        }
        denominator(row) = (2.0 - tau) * linear(row);
    }
    const double log_num = logsumexp_all_core(numerator);
    const double log_den = logsumexp_all_core(denominator);
    if (std::isfinite(log_num) && std::isfinite(log_den)) {
        return log_num - log_den;
    }
    return 0.0;
}

Eigen::VectorXd group_logsumexp_core(
    const Eigen::Ref<const Eigen::VectorXd>& values,
    const std::vector<std::int64_t>& groups,
    std::int64_t n_groups) {
    Eigen::VectorXd max_by_group =
        Eigen::VectorXd::Constant(n_groups, -std::numeric_limits<double>::infinity());
    std::vector<bool> has_nan(static_cast<std::size_t>(n_groups), false);
    for (Eigen::Index row = 0; row < values.size(); ++row) {
        const auto group = static_cast<Eigen::Index>(groups[static_cast<std::size_t>(row)]);
        const double value = values(row);
        if (std::isnan(value)) {
            has_nan[static_cast<std::size_t>(group)] = true;
            continue;
        }
        if (value > max_by_group(group)) {
            max_by_group(group) = value;
        }
    }
    Eigen::VectorXd out =
        Eigen::VectorXd::Constant(n_groups, -std::numeric_limits<double>::infinity());
    Eigen::VectorXd sums = Eigen::VectorXd::Zero(n_groups);
    for (Eigen::Index row = 0; row < values.size(); ++row) {
        const auto group = static_cast<Eigen::Index>(groups[static_cast<std::size_t>(row)]);
        const double max_value = max_by_group(group);
        const double value = values(row);
        if (has_nan[static_cast<std::size_t>(group)] || !std::isfinite(max_value)) {
            continue;
        }
        if (std::isfinite(value)) {
            sums(group) += std::exp(value - max_value);
        }
    }
    for (Eigen::Index group = 0; group < n_groups; ++group) {
        if (has_nan[static_cast<std::size_t>(group)]) {
            continue;
        }
        if (std::isinf(max_by_group(group)) && max_by_group(group) > 0.0) {
            out(group) = std::numeric_limits<double>::infinity();
        } else if (std::isfinite(max_by_group(group)) && sums(group) > 0.0) {
            out(group) = max_by_group(group) + std::log(sums(group));
        }
    }
    return out;
}

Eigen::VectorXd profile_one_effect_core(
    const Eigen::Ref<const Eigen::VectorXd>& endog,
    const Eigen::Ref<const Eigen::VectorXd>& partial_linear,
    const std::vector<std::int64_t>& groups,
    std::int64_t n_groups,
    double tau) {
    Eigen::VectorXd numerator =
        Eigen::VectorXd::Constant(partial_linear.size(), -std::numeric_limits<double>::infinity());
    Eigen::VectorXd denominator(partial_linear.size());
    for (Eigen::Index row = 0; row < partial_linear.size(); ++row) {
        if (endog(row) > 0.0) {
            numerator(row) =
                std::log(endog(row)) + (1.0 - tau) * partial_linear(row);
        }
        denominator(row) = (2.0 - tau) * partial_linear(row);
    }
    const Eigen::VectorXd log_num = group_logsumexp_core(numerator, groups, n_groups);
    const Eigen::VectorXd log_den = group_logsumexp_core(denominator, groups, n_groups);
    Eigen::VectorXd effect = Eigen::VectorXd::Zero(n_groups);
    for (Eigen::Index group = 0; group < effect.size(); ++group) {
        if (std::isfinite(log_num(group)) && std::isfinite(log_den(group))) {
            effect(group) = log_num(group) - log_den(group);
        }
    }
    return effect;
}

Eigen::VectorXd profile_linear_profiled_core(
    const Eigen::Ref<const Eigen::VectorXd>& endog,
    const Eigen::Ref<const RowMatrixXd>& exog,
    const Eigen::Ref<const Eigen::VectorXd>& beta,
    const std::vector<std::int64_t>& groups,
    std::int64_t n_groups,
    double tau) {
    Eigen::VectorXd base = exog * beta;
    if (is_singleton_intercept_group(groups, n_groups)) {
        const double adjustment = profile_intercept_adjustment_core(endog, base, tau);
        return Eigen::VectorXd(base.array() + adjustment);
    }
    const Eigen::VectorXd effect =
        profile_one_effect_core(endog, base, groups, n_groups, tau);
    Eigen::VectorXd linear(base.size());
    for (Eigen::Index row = 0; row < base.size(); ++row) {
        linear(row) = base(row) + effect(groups[static_cast<std::size_t>(row)]);
    }
    return linear;
}

Eigen::MatrixXd demean_singleton_intercept_core(
    const Eigen::Ref<const Eigen::MatrixXd>& values) {
    if (values.size() == 0) {
        return Eigen::MatrixXd(values);
    }
    Eigen::RowVectorXd mean = values.colwise().mean();
    return values.rowwise() - mean;
}

Eigen::MatrixXd demean_profiled_core(
    const Eigen::Ref<const Eigen::MatrixXd>& values,
    const std::vector<std::int64_t>& groups,
    std::int64_t n_groups) {
    if (is_singleton_intercept_group(groups, n_groups)) {
        return demean_singleton_intercept_core(values);
    }
    return demean_oneway_matrix(values, groups, n_groups, std::nullopt);
}

FullrecoverInvariantCore fullrecover_invariant_core(
    const Eigen::Ref<const RowMatrixXd>& exog,
    const std::vector<std::int64_t>& groups,
    std::int64_t n_groups) {
    FullrecoverInvariantCore state;
    state.x_demeaned = demean_profiled_core(exog, groups, n_groups);
    if (state.x_demeaned.cols() == 0) {
        state.projection = Eigen::MatrixXd(0, state.x_demeaned.rows());
        return state;
    }
    const Eigen::MatrixXd& x = state.x_demeaned;
    const Eigen::MatrixXd gram = x.transpose() * x;
    bool range_loss = !gram.allFinite();
    for (Eigen::Index j = 0; !range_loss && j < x.cols(); ++j) {
        range_loss = gram(j, j) < std::numeric_limits<double>::min() &&
                     (x.col(j).array() != 0.0).any();
    }
    if (range_loss) {
        const double largest = x.cwiseAbs().maxCoeff();
        if (largest > 0.0 && std::isfinite(largest) && x.allFinite()) {
            int shift;
            std::frexp(largest, &shift);
            const auto scale = [&shift](double v) { return std::scalbn(v, -shift); };
            const auto reversible = [&x](int k) {
                return (x.unaryExpr([k](double v) {
                    return std::scalbn(std::scalbn(v, -k), k);
                }).array() == x.array()).all();
            };
            // Keep the smallest entries exact if a full common shift loses bits.
            if (shift > 0 && !reversible(shift)) {
                int low = 0, high = shift;
                while (high - low > 1) {
                    const int mid = (low + high) / 2;
                    if (reversible(mid)) low = mid;
                    else high = mid;
                }
                shift = low;
            }
            const Eigen::MatrixXd z = x.unaryExpr(scale);
            const Eigen::MatrixXd scaled_gram = z.transpose() * z;
            if (shift && scaled_gram.allFinite()) {
                const Eigen::MatrixXd scaled = min_norm_lstsq(
                    scaled_gram, z.transpose());
                state.projection = scaled.unaryExpr(scale);
                if (state.projection.allFinite()) {
                    return state;
                }
            }
        }
    }
    state.projection = min_norm_lstsq(gram, x.transpose());
    return state;
}

Eigen::VectorXd initial_linear_profiled_core(
    const Eigen::Ref<const Eigen::VectorXd>& endog,
    const Eigen::Ref<const RowMatrixXd>& exog,
    const Eigen::Ref<const Eigen::VectorXd>& beta,
    const std::vector<std::int64_t>& groups,
    std::int64_t n_groups,
    double tau) {
    return profile_linear_profiled_core(endog, exog, beta, groups, n_groups, tau);
}

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
    Eigen::VectorXd& beta_new) {
    const Eigen::VectorXd linear_current =
        profile_linear_profiled_core(endog, exog, beta, groups, n_groups, tau);
    if (!all_finite_vec(linear_current)) {
        beta_new = beta;
        return false;
    }
    const double phase1_log_denom = phase == 1 ? std::log1p(param) : 0.0;
    const double phase1_inv_denom =
        phase == 1 ? 1.0 / (1.0 + param) : std::numeric_limits<double>::quiet_NaN();
    const bool tau_is_one = std::abs(tau - 1.0) <= 1e-12;
    const bool tau_is_two = std::abs(tau - 2.0) <= 1e-12;
    Eigen::VectorXd pseudo(linear_current.size());
    for (Eigen::Index row = 0; row < linear_current.size(); ++row) {
        const double linear = linear_current(row);
        if (phase == 0) {
            const double mu = exp_guard_scalar(linear);
            const double u_ratio = endog(row) / mu;
            const double weight = exp_guard_scalar((2.0 - tau) * linear);
            const double delta = param;
            const double delta_tilde = delta / std::max(weight, 1e-300);
            pseudo(row) =
                linear + std::log1p((u_ratio - 1.0) / (1.0 + delta_tilde));
        } else {
            double positive_component = 0.0;
            if (endog(row) > 0.0) {
                if (tau_is_one) {
                    positive_component = endog(row) * phase1_inv_denom;
                } else {
                    const double positive_log =
                        std::log(endog(row)) + (1.0 - tau) * linear -
                        phase1_log_denom;
                    positive_component = exp_guard_scalar(positive_log);
                }
            }
            const double negative_component =
                tau_is_two ? phase1_inv_denom
                           : exp_guard_scalar((2.0 - tau) * linear -
                                              phase1_log_denom);
            pseudo(row) = linear + positive_component - negative_component;
        }
    }
    if (!all_finite_vec(pseudo)) {
        beta_new = beta;
        return false;
    }
    const Eigen::MatrixXd pseudo_mat = pseudo;
    const Eigen::VectorXd demeaned =
        demean_profiled_core(pseudo_mat, groups, n_groups).col(0);
    beta_new = projection * demeaned;
    if (!all_finite_vec(beta_new)) {
        beta_new = beta;
        return false;
    }
    return true;
}

double profiled_score_row_core(double y, double linear, double tau) {
    const double log_weight = (1.0 - tau) * linear;
    const double log_base = (2.0 - tau) * linear;
    const double base = std::exp(log_base);
    if (y == 0.0) {
        return -base;
    }
    const double weight = std::exp(log_weight);
    const double score = y * weight - base;
    if (weight >= std::numeric_limits<double>::min() && std::isfinite(score)) {
        return score;
    }
    // Factor the difference by its larger term if the direct product lost range.
    const double log_y = std::log(y);
    const double d = log_y - linear;
    if (d == 0.0) {
        return 0.0;
    }
    const double larger = d > 0.0 ? log_y + log_weight : log_base;
    return std::copysign(std::exp(larger + std::log(-std::expm1(-std::abs(d)))), d);
}

double profiled_slope_score_residual_core(
    const Eigen::Ref<const Eigen::VectorXd>& endog,
    const Eigen::Ref<const RowMatrixXd>& exog,
    const Eigen::Ref<const Eigen::VectorXd>& beta,
    const std::vector<std::int64_t>& groups,
    std::int64_t n_groups,
    const Eigen::Ref<const Eigen::MatrixXd>& x_demeaned,
    double tau) {
    const Eigen::VectorXd linear =
        profile_linear_profiled_core(endog, exog, beta, groups, n_groups, tau);
    Eigen::VectorXd score_obs(linear.size());
    for (Eigen::Index row = 0; row < linear.size(); ++row) {
        score_obs(row) = profiled_score_row_core(endog(row), linear(row), tau);
    }
    if (!all_finite_vec(score_obs)) {
        return std::numeric_limits<double>::infinity();
    }
    const Eigen::VectorXd score = x_demeaned.transpose() * score_obs;
    if (!all_finite_vec(score)) {
        return std::numeric_limits<double>::infinity();
    }
    if (score.size() == 0) {
        return 0.0;
    }
    return score.cwiseAbs().maxCoeff() / std::max<double>(endog.size(), 1.0);
}

double profiled_fe_score_residual_core(
    const Eigen::Ref<const Eigen::VectorXd>& endog,
    const Eigen::Ref<const RowMatrixXd>& exog,
    const Eigen::Ref<const Eigen::VectorXd>& beta,
    const std::vector<std::int64_t>& groups,
    std::int64_t n_groups,
    double tau) {
    const Eigen::VectorXd linear =
        profile_linear_profiled_core(endog, exog, beta, groups, n_groups, tau);
    Eigen::VectorXd score_obs(linear.size());
    for (Eigen::Index row = 0; row < linear.size(); ++row) {
        score_obs(row) = profiled_score_row_core(endog(row), linear(row), tau);
    }
    if (!all_finite_vec(score_obs)) {
        return std::numeric_limits<double>::infinity();
    }
    if (is_singleton_intercept_group(groups, n_groups)) {
        return std::abs(score_obs.sum()) / std::max<double>(endog.size(), 1.0);
    }
    Eigen::VectorXd sums = Eigen::VectorXd::Zero(n_groups);
    for (Eigen::Index row = 0; row < score_obs.size(); ++row) {
        sums(groups[static_cast<std::size_t>(row)]) += score_obs(row);
    }
    return sums.cwiseAbs().maxCoeff() / std::max<double>(endog.size(), 1.0);
}

} // namespace powermean::k4
