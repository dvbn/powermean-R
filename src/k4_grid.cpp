// Adapted from Powermean 0.1.0, BSD-3-Clause. See inst/LICENSE.K4.
// Frozen source: 7b22dcf01c6fea74176823e24a4b18df559cd386
// src/powermean/_native/iols/bind.cpp
// SHA-256 7d12b4308bbff9f4ee6f3dc3f72f77c06874462c7e762d6bc7c8ad707060e4e5.

#include "binary_round.h"
#include "k4_core.h"
#include <Eigen/Cholesky>
#include <Eigen/Eigenvalues>

namespace powermean::k4 {

std::vector<double> exp_linspace(double log_min, double log_max, int size) {
    std::vector<double> out;
    const int n = std::max(size, 0);
    if (n == 0) {
        return out;
    }
    if (n == 1) {
        out.push_back(std::exp(log_min));
        return out;
    }
    out.reserve(static_cast<std::size_t>(n));
    for (int idx = 0; idx < n; ++idx) {
        const double frac = static_cast<double>(idx) / static_cast<double>(n - 1);
        out.push_back(std::exp(log_min + frac * (log_max - log_min)));
    }
    return out;
}

std::vector<double> unique_positive_finite(std::vector<double> values) {
    std::sort(values.begin(), values.end());
    std::vector<double> out;
    for (const double value : values) {
        if (!std::isfinite(value) || value <= 0.0) {
            continue;
        }
        if (out.empty() || value != out.back()) {
            out.push_back(value);
        }
    }
    return out;
}

std::vector<double> focus_rho_path_core(
    const std::vector<double>& path,
    std::optional<double> rho_hint) {
    if (path.empty() || !rho_hint.has_value() || !std::isfinite(*rho_hint) ||
        *rho_hint <= 0.0) {
        return path;
    }
    const double low = path.front();
    const double high = path.back();
    const double hint = std::min(std::max(*rho_hint, low), high);
    const double mults[] = {0.5, 0.75, 1.0, 1.25, 1.5, 2.0, 3.0, 5.0};
    std::vector<double> ordered;
    auto push_unique = [&ordered](double value) {
        for (const double seen : ordered) {
            if (seen == value) {
                return;
            }
        }
        ordered.push_back(value);
    };
    for (const double mult : mults) {
        const double value = std::min(std::max(hint * mult, low), high);
        if (std::isfinite(value) && value > 0.0) {
            push_unique(value);
        }
    }
    for (const double value : path) {
        push_unique(value);
    }
    return ordered;
}

std::vector<double> rho_path_core(
    const FullrecoverControlCore& control,
    std::optional<double> rho_hint) {
    if (!(std::isfinite(control.rho_log_min) &&
          std::isfinite(control.rho_log_max)) ||
        control.rho_log_max <= control.rho_log_min) {
        throw std::invalid_argument("rho_log_max must be greater than rho_log_min");
    }
    std::vector<double> path = unique_positive_finite(exp_linspace(
        control.rho_log_min,
        control.rho_log_max,
        std::max(control.rho_grid_size, 2)));
    return focus_rho_path_core(path, rho_hint);
}

std::vector<double> delta_grid_core(const FullrecoverControlCore& control) {
    // Internal package adapter; absent input preserves the paper path below.
    if (control.pilot_deltas.has_value()) {
        return *control.pilot_deltas;
    }
    const int grid_size = std::max(control.delta_grid_size, 0);
    if (grid_size == 0) {
        return {};
    }
    if (!(std::isfinite(control.delta_log_min) &&
          std::isfinite(control.delta_log_max)) ||
        control.delta_log_max <= control.delta_log_min) {
        throw std::invalid_argument("delta_log_max must be greater than delta_log_min");
    }
    return exp_linspace(control.delta_log_min, control.delta_log_max, grid_size);
}

void validate_direct_rho_grid_policy_core(const std::string& policy) {
    if (policy != kDirectRhoGridLegacy && policy != kDirectRhoGridStep &&
        policy != kDirectRhoGridCurvatureStep) {
        throw std::invalid_argument(
            "rho_grid_policy must be one of ('legacy', 'step_grid', "
            "'curvature_step_grid')");
    }
}

// Halving alpha grid in plain double arithmetic (libm exp only in callers).
std::vector<double> direct_step_rho_grid_core(
    double alpha0,
    double alpha_min,
    bool& clipped) {
    if (!(std::isfinite(alpha0) && alpha0 > 0.0)) {
        throw std::runtime_error("direct_rho_grid_invalid: alpha0");
    }
    if (!(std::isfinite(alpha_min) && alpha_min > 0.0)) {
        throw std::runtime_error("direct_rho_grid_invalid: alpha_min");
    }
    std::vector<double> alphas;
    clipped = alpha0 <= alpha_min;
    if (!clipped) {
        double a = alpha0;
        while (a > alpha_min) {
            alphas.push_back(a);
            a = a * 0.5;
        }
    }
    alphas.push_back(alpha_min);
    std::vector<double> rho;
    rho.reserve(alphas.size());
    for (const double a : alphas) {
        const double value = 1.0 / a - 1.0;
        if (!(std::isfinite(value) && value > 0.0)) {
            throw std::runtime_error("direct_rho_grid_invalid: rho");
        }
        if (!rho.empty() && !(value > rho.back())) {
            throw std::runtime_error(
                "direct_rho_grid_invalid: not strictly increasing");
        }
        rho.push_back(value);
    }
    return rho;
}

// Anchor the weighted mean at a heaviest row to retain small displacements.
Eigen::MatrixXd anchored_curvature_core(
    Eigen::MatrixXd centered, const Eigen::VectorXd& h,
    const Eigen::VectorXd& mean_weights, Eigen::Index peak) {
    const Eigen::RowVectorXd anchor = centered.row(peak);
    centered.rowwise() -= anchor;
    const Eigen::RowVectorXd shift = mean_weights.transpose() * centered;
    centered.rowwise() -= shift;
    centered.array().colwise() *= h.array().sqrt();
    return centered.transpose() * centered;
}


// Translation of the frozen binary64 rank check in data/design.py.
bool exact_dyadic_full_rank(const Eigen::MatrixXd& z) {
    const Eigen::Index p = z.cols();
    if (z.rows() < p || !z.allFinite()) return false;
    constexpr std::int64_t prime = 65537;
    std::vector<std::vector<std::int64_t>> minor(p, std::vector<std::int64_t>(p));

    const auto mod = [](std::int64_t x) {
        const std::int64_t remainder = x % prime;
        return remainder < 0 ? remainder + prime : remainder;
    };
    for (Eigen::Index i = 0; i < p; ++i) {
        for (Eigen::Index j = 0; j < p; ++j) {
            int exponent;
            const double fraction = std::frexp(z(i, j), &exponent);
            const auto significand = static_cast<std::int64_t>(std::ldexp(fraction, 53));
            const int power = ((exponent - 53) % 32 + 32) % 32;
            minor[i][j] = mod(mod(significand) * ((std::int64_t{1} << power) % prime));
        }
    }
    bool full_modular_rank = true;
    for (Eigen::Index column = 0; column < p; ++column) {
        Eigen::Index pivot = column;
        while (pivot < p && minor[pivot][column] == 0) ++pivot;
        if (pivot == p) {
            full_modular_rank = false;
            break;
        }
        if (pivot != column) minor[pivot].swap(minor[column]);
        const auto pivot_value = minor[column][column];
        for (Eigen::Index i = column + 1; i < p; ++i) {
            const auto multiplier = minor[i][column];
            for (Eigen::Index j = column + 1; j < p; ++j) {
                minor[i][j] = mod(pivot_value * minor[i][j] -
                                  multiplier * minor[column][j]);
            }
        }
    }
    if (full_modular_rank) return true;

    using boost::multiprecision::cpp_rational;
    std::vector<std::pair<Eigen::Index, std::vector<cpp_rational>>> basis;
    for (Eigen::Index i = 0; i < z.rows(); ++i) {
        std::vector<cpp_rational> row;
        row.reserve(p);
        for (Eigen::Index j = 0; j < p; ++j) row.emplace_back(dyadic(z(i, j)));
        for (const auto& previous : basis) {
            const cpp_rational multiplier = row[previous.first];
            if (multiplier == 0) continue;
            for (Eigen::Index j = 0; j < p; ++j) row[j] -= multiplier * previous.second[j];
        }
        for (Eigen::Index j = 0; j < p; ++j) {
            if (row[j] == 0) continue;
            const cpp_rational leading = row[j];
            for (auto& entry : row) entry /= leading;
            basis.emplace_back(j, std::move(row));
            if (basis.size() == static_cast<std::size_t>(p)) return true;
            break;
        }
    }
    return false;
}


// Form the same generalized spectrum in orthonormal slope coordinates when
// an unscaled intermediate overflows. Only the weighted mean uses relative
// weights; small independent directions retain their original signed curvature.
Eigen::MatrixXd equilibrated_curvature_core(
    const Eigen::MatrixXd& z, const Eigen::VectorXd& h, Eigen::Index anchor) {
    Eigen::MatrixXd u = z;
    for (Eigen::Index j = 0; j < u.cols(); ++j) {
        const double largest = u.col(j).cwiseAbs().maxCoeff();
        if (largest != 0.0) {
            u.col(j) /= std::scalbn(1.0, std::ilogb(largest));
        }
    }
    const Eigen::LLT<Eigen::MatrixXd> llt(u.transpose() * u);
    if (llt.info() != Eigen::Success) {
        throw std::runtime_error("direct_rho_curvature_failed: gram_cholesky_failed");
    }
    // A positive rounded pivot is not proof of independence. Use the same
    // exact dyadic rank certificate as the NumPy fallback on the original rows.
    if (z.cols() > 1 && !exact_dyadic_full_rank(z)) {
        throw std::runtime_error("direct_rho_curvature_failed: gram_cholesky_failed");
    }
    Eigen::MatrixXd q = llt.matrixL().solve(u.transpose()).transpose();
    const Eigen::VectorXd relative = h / h.cwiseAbs().maxCoeff();
    const double total = relative.sum();
    if (!(total > 0.0)) {
        throw std::runtime_error("direct_rho_curvature_failed: nonpositive_curvature_sum");
    }
    if (anchor >= 0) {
        return anchored_curvature_core(std::move(q), h, relative / total, anchor);
    }
    const Eigen::RowVectorXd mean = (relative.transpose() * q) / total;
    q.rowwise() -= mean;
    q.array().colwise() *= h.cwiseAbs().cwiseSqrt().array();
    return q.transpose() * (q.array().colwise() * h.array().sign()).matrix();
}

Eigen::MatrixXd direct_curvature_matrix_core(
    const Eigen::MatrixXd& x_demeaned, const Eigen::VectorXd& h, double total) {
    const double n_double = static_cast<double>(x_demeaned.rows());
    Eigen::Index peak;
    const double largest = h.maxCoeff(&peak);
    const double smallest = h.minCoeff();
    // On each centered contrast the raw moment / profiled variance is at
    // most max(h)/min(h). Switch when its n*eps error bound can cost half
    // the significand. Keep the original signed Schur form for signed weights.
    const double relative_floor = n_double * std::sqrt(std::numeric_limits<double>::epsilon());
    const Eigen::Index anchor =
        smallest >= 0.0 && smallest < largest && smallest / largest < relative_floor
            ? peak
            : -1;
    if (!std::isfinite(total)) {
        return equilibrated_curvature_core(x_demeaned, h, anchor);
    }
    const Eigen::VectorXd g = x_demeaned.transpose() * h;
    Eigen::MatrixXd hess_raw;
    if (anchor >= 0 && (g.array() != 0.0).any()) {
        if (!g.allFinite() ||
            g.cwiseAbs().maxCoeff() > std::sqrt(std::numeric_limits<double>::max())) {
            return equilibrated_curvature_core(x_demeaned, h, anchor);
        }
        hess_raw = anchored_curvature_core(x_demeaned, h, h / total, anchor) / n_double;
    } else {
        const Eigen::MatrixXd weighted =
            (x_demeaned.array().colwise() * h.array()).matrix();
        hess_raw =
            (x_demeaned.transpose() * weighted - (g * g.transpose()) / total) / n_double;
    }
    const Eigen::MatrixXd hess = 0.5 * (hess_raw + hess_raw.transpose());
    const Eigen::MatrixXd gram_raw =
        (x_demeaned.transpose() * x_demeaned) / n_double;
    const Eigen::MatrixXd gram = 0.5 * (gram_raw + gram_raw.transpose());
    // Detect finite subnormal loss as well as overflow, including g*g
    // before division by the signed total in the Schur correction.
    constexpr double tiny = std::numeric_limits<double>::min();
    if (!hess.allFinite() || !gram.allFinite() ||
        (gram.diagonal().array() < tiny).any() ||
        (hess.diagonal().cwiseAbs().array() < tiny).any() ||
        ((g.array() != 0.0) && (g.cwiseAbs().array() < std::sqrt(tiny))).any()) {
        return equilibrated_curvature_core(x_demeaned, h, anchor);
    }
    const Eigen::LLT<Eigen::MatrixXd> llt(gram);
    if (llt.info() != Eigen::Success) {
        throw std::runtime_error(
            "direct_rho_curvature_failed: gram_cholesky_failed");
    }
    const Eigen::MatrixXd lower = llt.matrixL();
    const Eigen::MatrixXd left =
        lower.triangularView<Eigen::Lower>().solve(hess);
    const Eigen::MatrixXd cmat_raw =
        lower.triangularView<Eigen::Lower>().solve(left.transpose());
    const Eigen::MatrixXd cmat = 0.5 * (cmat_raw + cmat_raw.transpose());
    return cmat.allFinite() ? cmat : equilibrated_curvature_core(x_demeaned, h, anchor);
}

// Largest generalized eigenvalue of the observed intercept-profiled negative
// Hessian H against Q = Z'Z/n at beta (Z = centered slope design). Diagnostic
// only: the projection and the update direction are not touched.
void direct_profiled_curvature_core(
    const Eigen::Ref<const Eigen::VectorXd>& endog,
    const Eigen::Ref<const RowMatrixXd>& exog,
    const Eigen::Ref<const Eigen::VectorXd>& beta,
    const std::vector<std::int64_t>& groups,
    std::int64_t n_groups,
    double tau,
    const Eigen::MatrixXd& x_demeaned,
    DirectRhoGridHookCore& hook) {
    const auto start = std::chrono::steady_clock::now();
    const Eigen::VectorXd linear =
        initial_linear_profiled_core(endog, exog, beta, groups, n_groups, tau);
    if (!all_finite_vec(linear)) {
        throw std::runtime_error("direct_rho_curvature_failed: nonfinite_index");
    }
    const Eigen::Index n = x_demeaned.rows();
    const double c1 = tau - 1.0;
    const double c2 = 2.0 - tau;
    Eigen::VectorXd h = Eigen::VectorXd::Zero(n);
    for (Eigen::Index row = 0; row < n; ++row) {
        double value = 0.0;
        if (c1 != 0.0 && endog(row) > 0.0) {
            value += c1 * exp_guard_scalar(
                              std::log(endog(row)) + (1.0 - tau) * linear(row));
        }
        if (c2 != 0.0) {
            value += c2 * exp_guard_scalar((2.0 - tau) * linear(row));
        }
        h(row) = value;
    }
    if (!all_finite_vec(h)) {
        throw std::runtime_error(
            "direct_rho_curvature_failed: nonfinite_row_curvature");
    }
    const double total = h.sum();
    if (std::isfinite(total) && total <= 0.0) {
        throw std::runtime_error(
            "direct_rho_curvature_failed: nonpositive_curvature_sum");
    }
    const Eigen::MatrixXd cmat = direct_curvature_matrix_core(x_demeaned, h, total);
    if (!cmat.allFinite()) {
        throw std::runtime_error("direct_rho_curvature_failed: nonfinite_hessian");
    }
    const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> spectrum(
        cmat, Eigen::EigenvaluesOnly);
    if (spectrum.info() != Eigen::Success) {
        throw std::runtime_error("direct_rho_curvature_failed: eigensolver_failed");
    }
    if (!spectrum.eigenvalues().allFinite()) {
        throw std::runtime_error(
            "direct_rho_curvature_failed: nonfinite_eigenvalue");
    }
    hook.has_curvature = true;
    hook.curvature_L0 = spectrum.eigenvalues().maxCoeff();
    hook.curvature_lambda_min = spectrum.eigenvalues().minCoeff();
    hook.curvature_wall_ns = static_cast<std::int64_t>(
        std::chrono::duration_cast<std::chrono::nanoseconds>(
            std::chrono::steady_clock::now() - start)
            .count());
}

// Non-legacy rho list at the Phase-2 handoff (curvature computed once here).
std::vector<double> direct_policy_rho_path_core(
    const Eigen::Ref<const Eigen::VectorXd>& endog,
    const Eigen::Ref<const RowMatrixXd>& exog,
    const Eigen::Ref<const Eigen::VectorXd>& beta,
    const std::vector<std::int64_t>& groups,
    std::int64_t n_groups,
    double tau,
    const FullrecoverControlCore& control,
    const FullrecoverInvariantCore& state,
    DirectRhoGridHookCore& hook) {
    const double alpha_max = 1.0 / (1.0 + std::exp(control.rho_log_min));
    const double alpha_min = 1.0 / (1.0 + std::exp(control.rho_log_max));
    double alpha0 = alpha_max;
    if (hook.policy == kDirectRhoGridCurvatureStep) {
        direct_profiled_curvature_core(
            endog, exog, beta, groups, n_groups, tau, state.x_demeaned, hook);
        alpha0 = 1.0 / std::max(1.0 + std::exp(control.rho_log_min), hook.curvature_L0);
    } else if (hook.policy != kDirectRhoGridStep) {
        validate_direct_rho_grid_policy_core(hook.policy);
        throw std::invalid_argument("rho_grid_policy legacy has no policy path");
    }
    bool clipped = false;
    std::vector<double> rho = direct_step_rho_grid_core(alpha0, alpha_min, clipped);
    hook.values = rho;
    hook.alpha0 = alpha0;
    hook.alpha_min = alpha_min;
    hook.clipped = clipped;
    return rho;
}

} // namespace powermean::k4
