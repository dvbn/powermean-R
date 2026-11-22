#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include "native_calls.h"

#include <boost/multiprecision/cpp_int.hpp>
#include <cmath>
#include <cstdint>
#include <limits>

SEXP support_ir_penalty_call(SEXP zeros, SEXP threshold) {
    namespace mp = boost::multiprecision;
    // INT_MAX * 2^(2 * (53 + 1073)) needs at most 2283 bits.
    using Integer = mp::number<mp::cpp_int_backend<2304, 2304, mp::unsigned_magnitude,
                                                  mp::unchecked, void>, mp::et_off>;
    using Double = std::numeric_limits<double>;
    if (TYPEOF(zeros) != INTSXP || Rf_xlength(zeros) != 1 ||
        INTEGER(zeros)[0] <= 0 || TYPEOF(threshold) != REALSXP ||
        Rf_xlength(threshold) != 1 || !(REAL(threshold)[0] > 0 && REAL(threshold)[0] < 1)) {
        Rf_error("Invalid iterative rectifier penalty inputs.");
    }
    double result = NA_REAL;
    {
        int exponent = 0;
        const double fraction = std::frexp(REAL(threshold)[0], &exponent);
        const Integer mantissa = static_cast<std::uint64_t>(std::ldexp(fraction, Double::digits));
        const Integer denominator = mantissa * mantissa;
        const Integer numerator = Integer(INTEGER(zeros)[0]) << (2 * (Double::digits - exponent));
        const Integer penalty = (numerator + denominator - 1) / denominator;
        const auto bits = mp::msb(penalty) + 1;
        const auto trailing = mp::lsb(penalty);
        if (bits <= Double::max_exponent && bits - trailing <= Double::digits) {
            result = penalty.convert_to<double>();
        }
    }
    return Rf_ScalarReal(result);
}
