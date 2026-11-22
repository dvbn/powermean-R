#ifndef POWERMEAN_NATIVE_CALLS_H
#define POWERMEAN_NATIVE_CALLS_H

#include <Rinternals.h>

SEXP project_call(SEXP codes, SEXP levels, SEXP weight, SEXP rhs, SEXP control);
SEXP project_effect_call(SEXP codes, SEXP levels, SEXP weight, SEXP rhs, SEXP control, SEXP initial);
SEXP support_ir_penalty_call(SEXP zeros, SEXP threshold);
SEXP gravity_proof_call(SEXP codes, SEXP roles, SEXP maps, SEXP levels, SEXP y, SEXP x);
SEXP pm_dense_iols(SEXP outcome, SEXP slopes, SEXP tau, SEXP controls);
SEXP basis_call(SEXP codes, SEXP levels);
SEXP basis_action_call(SEXP codes, SEXP levels, SEXP effects, SEXP workspace);
SEXP group_center_call(SEXP codes, SEXP levels, SEXP values);
SEXP bbp_mixed_products_call(SEXP z, SEXP mean, SEXP scores, SEXP curvature,
                            SEXP bread, SEXP covariance, SEXP group, SEXP sample);
SEXP bbp_covariance_blocks_call(SEXP transformed, SEXP split, SEXP owners);
SEXP cancel_gamma_call(SEXP offsets, SEXP columns, SEXP values, SEXP logs,
                       SEXP offset, SEXP start, SEXP controls);
SEXP cancel_gamma_sum_call(SEXP mass);
SEXP cancel_gamma_hessian_call(SEXP offsets, SEXP columns, SEXP values, SEXP counts,
                               SEXP weight, SEXP dimension);
SEXP cancel_gamma_dense_factor_call(SEXP hessian);
SEXP cancelled_csc_lu_call(SEXP offsets, SEXP rows, SEXP values, SEXP row_sum_scaling);
SEXP cancelled_dense_lu_solve_call(SEXP factor, SEXP pivot, SEXP right, SEXP transpose);
SEXP projected_schur_lu_call(SEXP coupling, SEXP diagonal, SEXP remainder);
SEXP projected_schur_lu_solve_call(SEXP factor, SEXP pivot, SEXP right, SEXP transpose);
SEXP checked_dense_lu_solve(SEXP factor, SEXP pivot, SEXP right, SEXP transpose, int maximum_width);
SEXP cancel_gamma_support_call(SEXP design, SEXP group, SEXP positive, SEXP units,
                               SEXP direction, SEXP face);
SEXP cancelled_log_recover_call(SEXP group, SEXP counts, SEXP response, SEXP raw, SEXP centered);
SEXP wz_pair_meats_call(SEXP mean, SEXP probability, SEXP score, SEXP slopes,
                        SEXP mapped, SEXP inverse, SEXP origins, SEXP destinations);
SEXP wz_pair_meats_controlled_call(SEXP mean, SEXP probability, SEXP score,
                                   SEXP slope, SEXP mapped, SEXP nuisance,
                                   SEXP origin, SEXP destination, SEXP threshold);
SEXP wz_nuisance_matrix_call(SEXP mean, SEXP probability, SEXP left, SEXP right, SEXP dimensions);
SEXP wz_role_moments_call(SEXP mean, SEXP probability, SEXP score, SEXP country, SEXP countries);
SEXP wz_role_scores_call(SEXP probability, SEXP score, SEXP country, SEXP inverse, SEXP variance, SEXP observed);
SEXP group_score_call(SEXP codes, SEXP levels, SEXP score, SEXP scale);
SEXP group_weighted_sums_call(SEXP codes, SEXP levels, SEXP block, SEXP design, SEXP score);
SEXP pm_sample_fingerprint(SEXP y, SEXP x, SEXP positions, SEXP rows, SEXP columns,
                           SEXP names, SEXP weights, SEXP weight_kind, SEXP offset,
                           SEXP offset_kind, SEXP exhausted);
SEXP pm_jochmans_fingerprint(SEXP prefix, SEXP codes, SEXP y, SEXP x);

SEXP pm_result_exact_rank_is(SEXP matrix, SEXP expected);
SEXP pm_result_exact_bilinear(SEXP left, SEXP matrix, SEXP right);
SEXP pm_prediction_exact_shift(SEXP index, SEXP changed, SEXP original, SEXP beta, SEXP multiplier);
SEXP pm_exact_dot_sign(SEXP left, SEXP right);
SEXP pm_index_within_forcing(SEXP operands);
SEXP pm_index_product_range(SEXP left, SEXP right, SEXP rounded);
SEXP pm_linear_slope_range(SEXP design, SEXP residual);
SEXP pm_linear_slope_sums(SEXP design, SEXP residual);
SEXP pm_linear_mean_loss(SEXP sums, SEXP totals, SEXP means);
SEXP pm_result_covariance_scale(SEXP matrix, SEXP units, SEXP order);
SEXP pm_design_scale_status(SEXP design);
SEXP pm_result_slope_scale(SEXP values, SEXP units, SEXP order, SEXP report);

#endif
