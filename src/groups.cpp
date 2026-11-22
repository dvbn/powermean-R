#include "groups.h"

#include <algorithm>

namespace powermean {

GroupStatus group_counts(const int* codes, std::size_t n_rows,
                         std::size_t n_levels, std::uint64_t* counts) noexcept {
    std::fill_n(counts, n_levels, std::uint64_t{0});
    for (std::size_t row = 0; row < n_rows; ++row) {
        const int code = codes[row];
        if (code < 0 || static_cast<std::size_t>(code) >= n_levels) {
            return GroupStatus::invalid_code;
        }
        ++counts[code];
    }
    return GroupStatus::ok;
}

GroupStatus group_weighted_sums(const int* codes, std::size_t n_rows, std::size_t n_levels,
                                const double* design, const double* score, std::size_t columns,
                                double* totals) noexcept {
    for (std::size_t row = 0; row < n_rows; ++row) {
        const int level = codes[row];
        if (level < 0 || static_cast<std::size_t>(level) >= n_levels) return GroupStatus::invalid_code;
    }
    for (std::size_t column = 0; column < columns; ++column) {
        double* target = totals + column * n_levels;
        const double* values = design + column * n_rows;
        std::fill_n(target, n_levels, 0.0);
        for (std::size_t row = 0; row < n_rows; ++row) target[codes[row]] += values[row] * score[row];
    }
    return GroupStatus::ok;
}

}  // namespace powermean
