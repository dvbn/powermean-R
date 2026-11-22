#ifndef POWERMEAN_GROUPS_H
#define POWERMEAN_GROUPS_H

#include <cstddef>
#include <cstdint>

namespace powermean {

enum class GroupStatus { ok, invalid_code };

// Codes are zero-based. The caller supplies n_levels writable counts.
// Counts are usable only on success.
GroupStatus group_counts(const int* codes, std::size_t n_rows,
                         std::size_t n_levels, std::uint64_t* counts) noexcept;

// Accumulate weighted matrix columns into a column-major levels-by-columns buffer.
GroupStatus group_weighted_sums(const int* codes, std::size_t n_rows, std::size_t n_levels,
                                const double* design, const double* score, std::size_t columns,
                                double* totals) noexcept;

}  // namespace powermean

#endif
