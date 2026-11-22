#include "k3.h"
#include "native_calls.h"

#include <algorithm>
#include <cmath>
#include <limits>

SEXP basis_action_call(SEXP codes, SEXP levels, SEXP effects, SEXP workspace) {
    if (TYPEOF(codes) != INTSXP || !Rf_isMatrix(codes) ||
        TYPEOF(levels) != INTSXP || TYPEOF(effects) != REALSXP ||
        Rf_getAttrib(effects, R_DimSymbol) != R_NilValue) {
        Rf_error("The basis action requires an integer code matrix, integer levels and double effects.");
    }
    const SEXP dimensions = Rf_getAttrib(codes, R_DimSymbol);
    const int rows = INTEGER(dimensions)[0];
    const int blocks = INTEGER(dimensions)[1];
    if (XLENGTH(levels) != blocks) {
        Rf_error("The basis action requires one level count per code column.");
    }
    if (XLENGTH(workspace) != 1 ||
        (TYPEOF(workspace) != REALSXP && TYPEOF(workspace) != INTSXP)) {
        Rf_error("The basis workspace must be a positive integer byte count or Inf.");
    }
    const double bytes = Rf_asReal(workspace);
    constexpr auto bits = std::numeric_limits<std::size_t>::digits;
    if (!(bytes > 0) || (std::isfinite(bytes) &&
        (bytes != std::trunc(bytes) || bytes >= std::ldexp(1.0, bits)))) {
        Rf_error("The basis workspace must be a representable positive integer byte count or Inf.");
    }
    const auto limit = std::isfinite(bytes) ? static_cast<std::size_t>(bytes)
                                          : std::numeric_limits<std::size_t>::max();
    const auto result = powermean::exact_fe_action(INTEGER(codes), rows, INTEGER(levels), blocks,
                                                    REAL(effects), XLENGTH(effects), limit);
    const SEXP output = PROTECT(Rf_allocVector(VECSXP, 3));
    SET_VECTOR_ELT(output, 0, Rf_allocVector(REALSXP, result.coefficients.size()));
    SET_VECTOR_ELT(output, 1, Rf_allocVector(INTSXP, result.position.size()));
    SET_VECTOR_ELT(output, 2, Rf_ScalarInteger(static_cast<int>(result.status)));
    std::copy(result.coefficients.begin(), result.coefficients.end(), REAL(VECTOR_ELT(output, 0)));
    std::copy(result.position.begin(), result.position.end(), INTEGER(VECTOR_ELT(output, 1)));
    const SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
    SET_STRING_ELT(names, 0, Rf_mkChar("coefficients"));
    SET_STRING_ELT(names, 1, Rf_mkChar("position"));
    SET_STRING_ELT(names, 2, Rf_mkChar("status"));
    Rf_setAttrib(output, R_NamesSymbol, names);
    UNPROTECT(2);
    return output;
}
