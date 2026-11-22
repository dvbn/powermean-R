#ifndef POWERMEAN_BINARY_ROUND_H
#define POWERMEAN_BINARY_ROUND_H

#include <boost/multiprecision/cpp_int.hpp>

#include <cmath>
#include <cstdint>
#include <stdexcept>

namespace powermean {

inline double nearest_binary64(const boost::multiprecision::cpp_rational& ratio) {
    using boost::multiprecision::cpp_int;
    cpp_int n = numerator(ratio);
    const cpp_int d = denominator(ratio);
    if (n == 0) return 0;
    const bool negative = n < 0;
    if (negative) n = -n;
    // Floor of log2(|n|/d) by comparing to an exact power of two.
    std::int64_t e = static_cast<std::int64_t>(boost::multiprecision::msb(n)) -
                     static_cast<std::int64_t>(boost::multiprecision::msb(d));
    if (e > 1024) throw std::overflow_error("binary64 ratio overflow");
    if (e < -1075) return negative ? -0.0 : 0.0;
    const bool below = e >= 0 ? n < (d << e) : (n << -e) < d;
    if (below) --e;
    if (e > 1023) throw std::overflow_error("binary64 ratio overflow");
    if (e < -1075) return negative ? -0.0 : 0.0;
    const int shift = e >= -1022 ? 52 - static_cast<int>(e) : 1074;
    if (shift >= 0) n <<= shift;
    cpp_int divisor = d;
    if (shift < 0) divisor <<= -shift;
    cpp_int rounded = n / divisor;
    const cpp_int twice_remainder = (n % divisor) << 1;
    if (twice_remainder > divisor ||
        (twice_remainder == divisor && boost::multiprecision::bit_test(rounded, 0))) {
        ++rounded;
    }
    const double magnitude = std::ldexp(rounded.convert_to<double>(), -shift);
    if (!std::isfinite(magnitude)) throw std::overflow_error("binary64 ratio overflow");
    return negative ? -magnitude : magnitude;
}

}  // namespace powermean

#endif
