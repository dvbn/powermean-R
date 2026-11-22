power_variance <- function(tau = 2) {
    if (!is.numeric(tau) || is.complex(tau) || length(tau) != 1L ||
        !is.null(dim(tau)) || !is.finite(tau)) {
        stop(errorCondition("tau must be a finite real scalar", class = "powermean_value_error"))
    }
    structure(list(tau = if (tau == 0) 0 else as.double(tau)), class = "powermean_power_variance")
}

as.list.powermean_power_variance <- function(x, ...) {
    list(type = "PowerVariance", tau = x$tau)
}

as.list.powermean_nb2 <- function(x, ...) {
    list(type = "NB2Variance", kappa = x$kappa)
}

pm_variance_type <- function(variance) {
    if (!is.list(variance)) return(NULL)
    if (inherits(variance, "powermean_power_variance") ||
        (identical(variance[["type"]], "PowerVariance") &&
         identical(sort(names(variance)), c("tau", "type")))) return("PowerVariance")
    if (inherits(variance, "powermean_nb2") ||
        (identical(variance[["type"]], "NB2Variance") &&
         identical(sort(names(variance)), c("kappa", "type")))) return("NB2Variance")
    NULL
}

pm_variance_parameters <- function(variance) {
    type <- pm_variance_type(variance)
    if (identical(type, "PowerVariance")) {
        return(as.list.powermean_power_variance(power_variance(variance[["tau"]])))
    }
    if (identical(type, "NB2Variance")) {
        return(as.list.powermean_nb2(NB2(variance[["kappa"]])))
    }
    pm_prepare_invalid("variance", variance, "a power_variance() or NB2() record")
}

pm_fit_variance_admission <- function(variance, type, tau, tau_explicit, fixed_effects) {
    if (is.null(variance)) return(invisible(NULL))
    if (tau_explicit && !is.null(tau)) {
        pm_covariance_conflict(paste0("tau=", pm_condition_repr(tau)),
                               paste0("variance=", pm_named_fit_repr(variance)))
    }
    if (is.null(type)) {
        valid <- if (fixed_effects) "'NB2Variance', 'PowerVariance'" else
            "a working variance: NB2Variance (pm.NB2) or PowerVariance"
        label <- pm_named_fit_repr(variance)
        pm_stop("specification.invalid_value", list(argument = "variance", value = variance, valid = valid),
                sprintf("variance=%s is not accepted; the valid values are %s.", label, valid))
    }
    invisible(NULL)
}

pm_variance_array <- function(log_mean) {
    if ((!is.numeric(log_mean) && !is.logical(log_mean)) || is.complex(log_mean)) {
        stop(errorCondition("log_mean must be real-valued", class = "powermean_value_error"))
    }
    if (any(!is.finite(log_mean))) {
        stop(errorCondition("log_mean must contain only finite values", class = "powermean_value_error"))
    }
    storage.mode(log_mean) <- "double"
    log_mean
}

log_value <- function(variance, log_mean) {
    spec <- pm_variance_parameters(variance)
    q <- pm_variance_array(log_mean)
    value <- if (!is.null(spec$tau)) spec$tau * q else if (spec$kappa == 0) q else
        q + pm_softplus(log(spec$kappa) + q)
    if (any(!is.finite(value))) {
        stop(errorCondition("log variance is not representable in float64",
                            class = "powermean_floating_point_error"))
    }
    value
}

elasticity <- function(variance, log_mean) {
    spec <- pm_variance_parameters(variance)
    q <- pm_variance_array(log_mean)
    if (!is.null(spec$tau)) return(q * 0 + spec$tau)
    if (spec$kappa == 0) return(q * 0 + 1)
    1 + exp(-pm_softplus(-log(spec$kappa) - q))
}

elasticity_derivative <- function(variance, log_mean) {
    spec <- pm_variance_parameters(variance)
    q <- pm_variance_array(log_mean)
    if (!is.null(spec$tau) || spec$kappa == 0) return(abs(q * 0))
    tail <- exp(-abs(log(spec$kappa) + q))
    tail / (1 + tail)^2
}
