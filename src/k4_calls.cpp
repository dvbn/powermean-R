#include "k4.h"

#define R_NO_REMAP
#include <Rinternals.h>
#include <R_ext/Utils.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>

namespace {

constexpr const char* fields[] = {
    "slopes", "index", "intercept", "converged", "phase", "iterations", "restarts",
    "phase1_attempts", "phase1_accepted", "phase1_rejected", "phase1_evaluations",
    "phase2_evaluations", "response_unit_exponent", "residual", "score_residual",
    "intercept_score_residual", "final_rho", "final_criterion", "failure_reason",
    "score_validation_status"
};

void check_for_interrupt(void*) {
    R_CheckUserInterrupt();
}

void interrupt() {
    if (!R_ToplevelExec(check_for_interrupt, nullptr)) {
        throw std::runtime_error("IOLS iteration interrupted.");
    }
}

double scalar(SEXP x, const char* description) {
    if ((TYPEOF(x) != REALSXP && TYPEOF(x) != INTSXP) || XLENGTH(x) != 1) {
        throw std::invalid_argument(description);
    }
    const double value = Rf_asReal(x);
    if (!std::isfinite(value)) throw std::invalid_argument(description);
    return value;
}

int count(SEXP x, int minimum, const char* description) {
    const double value = scalar(x, description);
    if (value < minimum || value >= std::numeric_limits<int>::max() ||
        value != std::trunc(value)) {
        throw std::invalid_argument(description);
    }
    return static_cast<int>(value);
}

powermean::DenseIolsOptions parse_options(SEXP controls) {
    if (TYPEOF(controls) != VECSXP || XLENGTH(controls) != 8) {
        throw std::invalid_argument("Invalid IOLS control record.");
    }
    powermean::DenseIolsOptions options;
    options.eps = scalar(VECTOR_ELT(controls, 0), "Invalid IOLS eps.");
    if (options.eps <= 0) throw std::invalid_argument("Invalid IOLS eps.");
    options.anderson_memory = count(VECTOR_ELT(controls, 1), 0, "Invalid Anderson memory.");
    options.max_phase1_iterations = count(VECTOR_ELT(controls, 2), 2, "Invalid Phase 1 budget.");
    options.stopcrit = count(VECTOR_ELT(controls, 3), 2, "Invalid IOLS stop budget.");
    options.max_phase2_evaluations = count(VECTOR_ELT(controls, 4), 1, "Invalid Phase 2 budget.");
    SEXP shifts = VECTOR_ELT(controls, 5);
    if (shifts != R_NilValue) {
        if (TYPEOF(shifts) != REALSXP || XLENGTH(shifts) < 1 ||
            XLENGTH(shifts) >= std::numeric_limits<int>::max()) {
            throw std::invalid_argument("Invalid IOLS pilot shifts.");
        }
        const double* values = REAL(shifts);
        options.pilot_deltas.assign(values, values + XLENGTH(shifts));
    }
    options.rho_log_min = scalar(VECTOR_ELT(controls, 6), "Invalid IOLS lower rho bound.");
    options.rho_log_max = scalar(VECTOR_ELT(controls, 7), "Invalid IOLS upper rho bound.");
    options.interrupt = interrupt;
    return options;
}

} // namespace

SEXP pm_dense_iols(SEXP outcome, SEXP slopes, SEXP tau, SEXP controls) {
    SEXP dimensions = Rf_getAttrib(slopes, R_DimSymbol);
    if (TYPEOF(outcome) != REALSXP || TYPEOF(slopes) != REALSXP ||
        TYPEOF(dimensions) != INTSXP || XLENGTH(dimensions) != 2) {
        Rf_error("IOLS requires a double outcome and a double slope matrix.");
    }
    const int rows = INTEGER(dimensions)[0];
    const int cols = INTEGER(dimensions)[1];
    if (rows < 3 || cols < 1 || cols >= rows - 1 || XLENGTH(outcome) != rows ||
        XLENGTH(slopes) != static_cast<R_xlen_t>(rows) * cols) {
        Rf_error("The IOLS outcome and slope matrix are misaligned.");
    }
    constexpr int field_count = sizeof(fields) / sizeof(fields[0]);
    SEXP result = PROTECT(Rf_allocVector(VECSXP, field_count));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, field_count));
    SEXP beta = PROTECT(Rf_allocVector(REALSXP, cols));
    SEXP index = PROTECT(Rf_allocVector(REALSXP, rows));
    for (int j = 0; j < field_count; ++j) SET_STRING_ELT(names, j, Rf_mkChar(fields[j]));
    Rf_setAttrib(result, R_NamesSymbol, names);
    SET_VECTOR_ELT(result, 0, beta);

    bool has_index = false;
    bool converged = false;
    int counts[8];
    double values[6];
    char phase[32];
    char failure[192];
    char score[96];
    char error[384] = {};
    try {
        const auto fit = powermean::dense_iols(
            REAL(outcome), REAL(slopes), rows, cols, scalar(tau, "Invalid IOLS tau."),
            parse_options(controls)
        );
        if (fit.slopes.size() != static_cast<std::size_t>(cols)) {
            throw std::runtime_error("IOLS returned an incomplete candidate.");
        }
        std::copy(fit.slopes.begin(), fit.slopes.end(), REAL(beta));
        has_index = fit.index.size() == static_cast<std::size_t>(rows);
        if (has_index) std::copy(fit.index.begin(), fit.index.end(), REAL(index));
        converged = fit.converged;
        counts[0] = fit.iterations;
        counts[1] = fit.restarts;
        counts[2] = fit.phase1_attempts;
        counts[3] = fit.phase1_accepted;
        counts[4] = fit.phase1_rejected;
        counts[5] = fit.phase1_evaluations;
        counts[6] = fit.phase2_evaluations;
        counts[7] = fit.response_unit_exponent;
        values[0] = fit.intercept;
        values[1] = fit.residual;
        values[2] = fit.score_residual;
        values[3] = fit.intercept_score_residual;
        values[4] = fit.final_rho;
        values[5] = fit.final_criterion;
        std::snprintf(phase, sizeof(phase), "%s", fit.phase.c_str());
        std::snprintf(failure, sizeof(failure), "%s", fit.failure_reason.c_str());
        std::snprintf(score, sizeof(score), "%s", fit.score_validation_status.c_str());
    } catch (const std::exception& problem) {
        std::snprintf(error, sizeof(error), "%s", problem.what());
    } catch (...) {
        std::snprintf(error, sizeof(error), "%s", "Unexpected IOLS kernel failure.");
    }
    if (error[0]) Rf_error("%s", error);

    if (has_index) SET_VECTOR_ELT(result, 1, index);
    SET_VECTOR_ELT(result, 2, Rf_ScalarReal(values[0]));
    SET_VECTOR_ELT(result, 3, Rf_ScalarLogical(converged));
    SET_VECTOR_ELT(result, 4, Rf_mkString(phase));
    for (int j = 0; j < 8; ++j) SET_VECTOR_ELT(result, j + 5, Rf_ScalarInteger(counts[j]));
    for (int j = 1; j < 6; ++j) SET_VECTOR_ELT(result, j + 12, Rf_ScalarReal(values[j]));
    SET_VECTOR_ELT(result, 18, Rf_mkString(failure));
    SET_VECTOR_ELT(result, 19, Rf_mkString(score));
    UNPROTECT(4);
    return result;
}
