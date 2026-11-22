#ifndef POWERMEAN_K3_FRACTION_H
#define POWERMEAN_K3_FRACTION_H

// Checked rational arithmetic derived from the frozen Powermean source.
// Copyright (c) 2026, David Benatia. BSD 3-Clause. See K3.md and LICENSE.K3.
#include <boost/multiprecision/cpp_int.hpp>

#include <cstdint>
#include <limits>
#include <numeric>

namespace powermean::k3 {

struct ArithmeticRange {};

constexpr std::int64_t I64_MAX = std::numeric_limits<std::int64_t>::max();
constexpr std::int64_t I64_MIN = std::numeric_limits<std::int64_t>::min();

inline std::int64_t checked_add(std::int64_t a, std::int64_t b) {
    if ((b > 0 && a > I64_MAX - b) || (b < 0 && a < I64_MIN - b)) throw ArithmeticRange();
    return a + b;
}

inline std::int64_t checked_neg(std::int64_t a) {
    if (a == I64_MIN) throw ArithmeticRange();
    return -a;
}

inline std::int64_t checked_mul(std::int64_t a, std::int64_t b) {
    if (a == 0 || b == 0) return 0;
    if (a > 0) {
        if (b > 0) {
            if (a > I64_MAX / b) throw ArithmeticRange();
        } else if (b < I64_MIN / a) {
            throw ArithmeticRange();
        }
    } else if (b > 0) {
        if (a < I64_MIN / b) throw ArithmeticRange();
    } else if (a == I64_MIN || b == I64_MIN || -a > I64_MAX / -b) {
        throw ArithmeticRange();
    }
    return a * b;
}

inline std::int64_t positive_gcd(std::int64_t a, std::int64_t b) {
    if (a == I64_MIN || b == I64_MIN) throw ArithmeticRange();
    return std::gcd(a < 0 ? -a : a, b < 0 ? -b : b);
}

// An exact rational in lowest terms with a positive denominator (Fraction's form).
struct Small {
    std::int64_t n;
    std::int64_t d = 1;
};

#if defined(__SIZEOF_INT128__)
using Wide = __int128;
using UWide = unsigned __int128;

inline UWide wide_gcd(UWide a, UWide b) {
    while (b != 0) {
        const UWide r = a % b;
        a = b;
        b = r;
    }
    return a;
}

// n / d in lowest terms; declines unless both fit 64 bits
inline Small normalized_wide(Wide n, Wide d) {
    if (d < 0) {
        n = -n;
        d = -d;
    }
    if (n == 0) return {0, 1};
    const UWide g = wide_gcd(n < 0 ? static_cast<UWide>(-n) : static_cast<UWide>(n),
                             static_cast<UWide>(d));
    n /= static_cast<Wide>(g);
    d /= static_cast<Wide>(g);
    if (n > static_cast<Wide>(I64_MAX) || n < static_cast<Wide>(I64_MIN)
            || d > static_cast<Wide>(I64_MAX)) {
        throw ArithmeticRange();
    }
    return {static_cast<std::int64_t>(n), static_cast<std::int64_t>(d)};
}

inline Small q_add(Small a, Small b) {
    if (a.d == 1 && b.d == 1) return {checked_add(a.n, b.n), 1};
    return normalized_wide(static_cast<Wide>(a.n) * b.d + static_cast<Wide>(b.n) * a.d,
                           static_cast<Wide>(a.d) * b.d);
}

inline Small q_mul(Small a, Small b) {
    if (a.d == 1 && b.d == 1) return {checked_mul(a.n, b.n), 1};
    return normalized_wide(static_cast<Wide>(a.n) * b.n, static_cast<Wide>(a.d) * b.d);
}

// value / lead exactly (_over): an integer when it divides, else the fraction.
inline Small q_div(Small a, Small b) {
    return normalized_wide(static_cast<Wide>(a.n) * b.d, static_cast<Wide>(a.d) * b.n);
}
#else
inline Small normalized(std::int64_t n, std::int64_t d) {
    if (d < 0) {
        n = checked_neg(n);
        d = checked_neg(d);
    }
    if (n == 0) return {0, 1};
    const std::int64_t g = positive_gcd(n, d);
    return {n / g, d / g};
}

inline Small q_add(Small a, Small b) {
    if (a.d == 1 && b.d == 1) return {checked_add(a.n, b.n), 1};
    return normalized(checked_add(checked_mul(a.n, b.d), checked_mul(b.n, a.d)),
                      checked_mul(a.d, b.d));
}

inline Small q_mul(Small a, Small b) {
    if (a.d == 1 && b.d == 1) return {checked_mul(a.n, b.n), 1};
    return normalized(checked_mul(a.n, b.n), checked_mul(a.d, b.d));
}

inline Small q_div(Small a, Small b) {
    return normalized(checked_mul(a.n, b.d), checked_mul(a.d, b.n));
}
#endif

inline Small q_neg(Small a) { return {checked_neg(a.n), a.d}; }

inline bool is_zero(Small value) { return value.n == 0; }

using Big = boost::multiprecision::cpp_rational;
inline Big q_add(const Big& a, const Big& b) { return a + b; }
inline Big q_mul(const Big& a, const Big& b) { return a * b; }
inline Big q_div(const Big& a, const Big& b) { return a / b; }
inline Big q_neg(const Big& a) { return -a; }
inline bool is_zero(const Big& value) { return value == 0; }

}  // namespace powermean::k3

#endif
