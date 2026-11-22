#define R_NO_REMAP
#include <R.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include <Rinternals.h>

#include "groups.h"
#include "native_calls.h"

namespace {

SEXP group_counts_call(SEXP codes, SEXP n_levels) {
    if (TYPEOF(codes) != INTSXP || TYPEOF(n_levels) != INTSXP ||
        Rf_xlength(n_levels) != 1 || INTEGER(n_levels)[0] < 0) {
        Rf_error("Group codes and the number of levels must be integers.");
    }
    const auto levels = static_cast<std::size_t>(INTEGER(n_levels)[0]);
    const auto rows = static_cast<std::size_t>(Rf_xlength(codes));
    auto* counts = reinterpret_cast<std::uint64_t*>(
        R_alloc(levels, sizeof(std::uint64_t))
    );
    if (powermean::group_counts(INTEGER(codes), rows, levels, counts) !=
        powermean::GroupStatus::ok) {
        Rf_error("Group codes must be present and between zero and the number of levels minus one.");
    }
    SEXP output = PROTECT(Rf_allocVector(REALSXP, levels));
    for (std::size_t level = 0; level < levels; ++level) {
        REAL(output)[level] = static_cast<double>(counts[level]);
    }
    UNPROTECT(1);
    return output;
}

const R_CallMethodDef call_methods[] = {
    {"pm_group_counts", reinterpret_cast<DL_FUNC>(&group_counts_call), 2},
    {"pm_project_fe", reinterpret_cast<DL_FUNC>(&project_call), 5},
    {"pm_project_fe_effects", reinterpret_cast<DL_FUNC>(&project_effect_call), 6},
    {"pm_support_ir_penalty", reinterpret_cast<DL_FUNC>(&support_ir_penalty_call), 2},
    {"pm_exact_basis", reinterpret_cast<DL_FUNC>(&basis_call), 2},
    {"pm_basis_action", reinterpret_cast<DL_FUNC>(&basis_action_call), 4},
    {"pm_gravity_proof", reinterpret_cast<DL_FUNC>(&gravity_proof_call), 6},
    {"pm_group_score", reinterpret_cast<DL_FUNC>(&group_score_call), 4},
    {"pm_group_weighted_sums", reinterpret_cast<DL_FUNC>(&group_weighted_sums_call), 5},
    {"pm_group_center", reinterpret_cast<DL_FUNC>(&group_center_call), 3},
    {"pm_sample_fingerprint", reinterpret_cast<DL_FUNC>(&pm_sample_fingerprint), 11},
    {"pm_jochmans_fingerprint", reinterpret_cast<DL_FUNC>(&pm_jochmans_fingerprint), 4},
    {"pm_dense_iols", reinterpret_cast<DL_FUNC>(&pm_dense_iols), 4},
    {"pm_wz_pair_meats", reinterpret_cast<DL_FUNC>(&wz_pair_meats_call), 8},
    {"pm_wz_pair_meats_controlled", reinterpret_cast<DL_FUNC>(&wz_pair_meats_controlled_call), 9},
    {"pm_wz_nuisance_matrix", reinterpret_cast<DL_FUNC>(&wz_nuisance_matrix_call), 5},
    {"pm_wz_role_moments", reinterpret_cast<DL_FUNC>(&wz_role_moments_call), 5},
    {"pm_wz_role_scores", reinterpret_cast<DL_FUNC>(&wz_role_scores_call), 6},
    {"pm_result_exact_rank_is", reinterpret_cast<DL_FUNC>(&pm_result_exact_rank_is), 2},
    {"pm_index_within_forcing", reinterpret_cast<DL_FUNC>(&pm_index_within_forcing), 1},
    {"pm_index_product_range", reinterpret_cast<DL_FUNC>(&pm_index_product_range), 3},
    {"pm_linear_slope_range", reinterpret_cast<DL_FUNC>(&pm_linear_slope_range), 2},
    {"pm_linear_slope_sums", reinterpret_cast<DL_FUNC>(&pm_linear_slope_sums), 2},
    {"pm_linear_mean_loss", reinterpret_cast<DL_FUNC>(&pm_linear_mean_loss), 3},
    {"pm_result_exact_bilinear", reinterpret_cast<DL_FUNC>(&pm_result_exact_bilinear), 3},
    {"pm_prediction_exact_shift", reinterpret_cast<DL_FUNC>(&pm_prediction_exact_shift), 5},
    {"pm_exact_dot_sign", reinterpret_cast<DL_FUNC>(&pm_exact_dot_sign), 2},
    {"pm_result_covariance_scale", reinterpret_cast<DL_FUNC>(&pm_result_covariance_scale), 3},
    {"pm_design_scale_status", reinterpret_cast<DL_FUNC>(&pm_design_scale_status), 1},
    {"pm_result_slope_scale", reinterpret_cast<DL_FUNC>(&pm_result_slope_scale), 4},
    {"pm_cancel_gamma", reinterpret_cast<DL_FUNC>(&cancel_gamma_call), 7},
    {"pm_cancel_gamma_sum", reinterpret_cast<DL_FUNC>(&cancel_gamma_sum_call), 1},
    {"pm_cancel_gamma_hessian", reinterpret_cast<DL_FUNC>(&cancel_gamma_hessian_call), 6},
    {"pm_cancel_gamma_dense_factor", reinterpret_cast<DL_FUNC>(&cancel_gamma_dense_factor_call), 1},
    {"pm_cancelled_csc_lu", reinterpret_cast<DL_FUNC>(&cancelled_csc_lu_call), 4},
    {"pm_cancelled_dense_lu_solve", reinterpret_cast<DL_FUNC>(&cancelled_dense_lu_solve_call), 4},
    {"pm_projected_schur_lu", reinterpret_cast<DL_FUNC>(&projected_schur_lu_call), 3},
    {"pm_projected_schur_lu_solve", reinterpret_cast<DL_FUNC>(&projected_schur_lu_solve_call), 4},
    {"pm_cancel_gamma_support", reinterpret_cast<DL_FUNC>(&cancel_gamma_support_call), 6},
    {"pm_cancelled_log_recover", reinterpret_cast<DL_FUNC>(&cancelled_log_recover_call), 5},
    {"pm_bbp_mixed_products", reinterpret_cast<DL_FUNC>(&bbp_mixed_products_call), 8},
    {"pm_bbp_covariance_blocks", reinterpret_cast<DL_FUNC>(&bbp_covariance_blocks_call), 3},
    {nullptr, nullptr, 0}
};

}  // namespace

extern "C" attribute_visible void R_init_powermean(DllInfo* dll) {
    R_registerRoutines(dll, nullptr, call_methods, nullptr, nullptr);
    R_useDynamicSymbols(dll, FALSE);
    R_forceSymbols(dll, TRUE);
}
