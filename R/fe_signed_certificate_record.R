# Retain selected-point diagnostics without keeping the fitted design or the
# solver's factors. The incidence coordinates are in the native exact basis.
pm_fe_signed_certificate_record <- function(prepared, root, options, controls) {
    certificate <- root$certificate
    geometry <- root$certificate_geometry
    if (is.null(geometry)) return(certificate)

    p <- ncol(prepared$x)
    variance <- list(type = "PowerVariance", tau = options$tau)
    coefficients <- pm_certificate_coefficients(certificate$coefficients, controls$coef_tol)
    basis <- list(rank = prepared$basis$rank, rows = geometry$n_rows,
                  blocks = prepared$fe, levels = prepared$levels,
                  position = prepared$basis$position)
    point <- list(kind = "root", n_slopes = p, slopes = colnames(prepared$x), variance = variance,
                  theta = utils::tail(root$beta, p), coordinates = root$beta,
                  slope_units = prepared$slope_units, omega = options$omega,
                  inner = pm_fe_certificate_inner_controls(options,
                      !identical(options$absorbed_options$absorbed, "index")))
    rows <- if (!isFALSE(options$keep_influence)) {
        pm_fe_joint_certificate_rows(prepared, root$term$score, root$term$curvature)
    }
    measures <- list(rcond = certificate$rcond, conditioned = certificate$conditioned,
        positive_definite = certificate$positive_definite, scaled_score = certificate$scaled_score,
        index_correction = certificate$index_correction, score_metric = certificate$score_metric,
        correction = root$correction)
    effects <- list(maximum = certificate$group_score, unit = geometry$score_unit)
    work <- root$certificate_factor_work
    if (!is.null(work)) {
        measures[c("ordinary", "curvature", "inverse_norm_estimate")] <-
            list(work$ordinary, work$curvature, work$inverse_norm)
        step <- pm_fe_joint_work_step
        effect_step <- if (is.finite(effects$maximum)) step("effect equations", "completed") else
            step("effect equations", "refused", "The fixed-effect group equations are not representable.")
        origin <- root$certificate_origin_work %||% work
        assembly_step <- if (is.null(origin$assembly_refusal)) {
            step("assembly", "completed", counters = list("H stored entries" = origin$hessian_entries))
        } else step("assembly", "refused", origin$assembly_refusal)
        factor_steps <- function(value) c(list(value$equilibration), pm_fe_joint_factor_steps(value))
        certificate$work <- c(list(
            step("basis", "completed", "exact indicator coordinates", list(rows = geometry$n_rows)),
            step("rows", "completed"), effect_step, assembly_step), factor_steps(origin))
        for (retry in root$certificate_retry_work) {
            certificate$work <- c(certificate$work,
                if (is.environment(retry)) factor_steps(retry) else list(retry))
        }
        if (!is.null(root$certificate_origin_work)) {
            certificate$work <- c(certificate$work, factor_steps(work))
        }
    }
    inapplicable <- function(name, reason) {
        list(name = name, value = NULL, bounds = numeric(), passed = FALSE,
             reason = paste0(name, ": not applicable: an exact full-system root; ", reason),
             applicable = FALSE)
    }
    checks <- list(inapplicable("nuisance representation",
        paste0("the nuisance coordinates of E_J are unknowns of the sparse joint system, ",
               "not formed by an absorption, profile or recovery, so no achieved inner tolerance exists")))
    if (p) {
        checks[[2L]] <- inapplicable("auxiliary absorption",
            paste0("the full sparse joint Jacobian is factorized, no C-weighted slope absorption runs, ",
                   "and the observed curvature of a nonconvex equation is not an absorption weight"))
    }
    inner <- list(applicable = FALSE, passed = FALSE, checks = checks, reason = "",
                  not_applicable = vapply(checks, `[[`, "", "reason"), steps = list())
    not_applicable <- c("inner", "curvature", if (geometry$n_rows == geometry$dimension) "coefficient")
    certificate[c("controls", "point", "basis", "geometry", "rows", "effects", "measures", "coefficients", "inner",
                  "error_messages", "errors", "not_applicable", "require_positive_definite", "empty", "degenerate")] <-
        list(controls[c("tol", "tol_secondary", "coef_tol", "rcond_min")], point, basis, geometry, rows, effects,
             measures, coefficients, inner, list(), list(), not_applicable, FALSE, p == 0L,
             isTRUE(coefficients$degenerate))
    certificate$assessment <- pm_point_assessment(certificate)
    certificate
}
