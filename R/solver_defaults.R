solver_defaults <- function(solver = "newton") {
    if (!is.character(solver) || length(solver) != 1L || !is.null(dim(solver)) ||
        is.object(solver) || is.na(solver) || !solver %in% c("newton", "iols")) {
        pm_invalid_value_record("solver", solver, "'iols', 'newton'")
    }
    if (solver == "iols") {
        return(structure(pm_iols_defaults(), class = "powermean_iols_options", backend_explicit = FALSE))
    }
    structure(list(
        tol = 1e-8, tol_secondary = 1e-6, coef_tol = 1e-6, rcond_min = 1e-13,
        newton_provider = "bundled", gradient_tol = 1e-10, max_iter = 100000L,
        max_evaluations = 1000000L, armijo = 2^-11, backtrack_factor = 0.5,
        max_backtracks = 21L, retry = 1L, curvature = "observed"
    ), class = "powermean_solver_options")
}

as.list.powermean_solver_options <- function(x, ...) {
    attributes(x) <- list(names = names(x))
    x
}

pm_iols_defaults <- function() {
    list(
        eps = 1e-12, anderson_m = 5L, max_phase1_iterations = 100000L,
        stopcrit = 1000000L, max_phase2_evaluations = 5000000L,
        # Preserve the wheel's binary64 values independently of libm.
        delta_grid = c(0x1.b993fe00d5376p-8, 0x1.50385c094f425p-4, 1,
                       0x1.85d6fd931e0bbp+3, 0x1.28d389970338fp+7),
        rho_log_min = -7, rho_log_max = 14, iols_backend = "cpp"
    )
}

iols_options <- function(...) {
    values <- list(...)
    if (length(values) == 1L && is.null(names(values))) values <- values[[1L]]
    if (!pm_iols_is_mapping(values)) {
        pm_iols_control_invalid("solver_options", values, "mappings of solver option names to values")
    }
    structure(pm_iols_controls(values), class = "powermean_iols_options",
              backend_explicit = "iols_backend" %in% names(values))
}

`$.powermean_iols_options` <- function(x, name) {
    if (identical(name, "backend_explicit")) return(attr(x, "backend_explicit", exact = TRUE))
    if (identical(name, "effective_backend_request")) {
        explicit <- pm_iols_record_presence(x)
        pm_solver_backend_choice(x[["iols_backend"]], "cpp", "iols_backend")
        return(pm_iols_backend_request(explicit))
    }
    x[[name, exact = FALSE]]
}

as.list.powermean_iols_options <- function(x, ...) {
    attributes(x) <- list(names = names(x))
    x
}
