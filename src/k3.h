#ifndef POWERMEAN_K3_H
#define POWERMEAN_K3_H

#include <cstddef>
#include <vector>

namespace powermean {

enum class BasisStatus {
    ok,
    invalid_shape,
    invalid_code,
    missing_level,
    capacity,
    allocation,
    not_representable,
    resource_limit
};

struct BasisResult {
    BasisStatus status = BasisStatus::invalid_shape;
    int n_components = 0;
    bool multiprecision = false;
    // Flat by formula block, then by level code. E omits the reference
    // levels. J is the exact independent subset of E, in E-column order.
    std::vector<int> component_by_level;
    std::vector<int> column_by_level;    // E column, -1 for references
    std::vector<int> position_by_level;  // J position, -1 for omitted levels
    std::vector<int> selected;           // E columns in J; size is exact rank
};

// Codes are zero-based and column-major: codes[block * n_rows + row].
// level_counts has n_blocks entries. Every declared level must be observed.
// The caller keeps both buffers alive for the call. Results own their storage
// and are populated only on success. No dense dummy matrix is formed.
BasisResult exact_fe_basis(const int* codes, std::size_t n_rows,
                           const int* level_counts, std::size_t n_blocks);

struct BasisActionResult {
    BasisStatus status = BasisStatus::invalid_shape;
    std::vector<int> position;
    std::vector<double> coefficients;
};

// Apply the canonical quotient map to the full level effects, in block order.
// workspace_bytes constrains auxiliary work, excluding caller-owned inputs.
BasisActionResult exact_fe_action(const int* codes, std::size_t n_rows,
                                 const int* level_counts, std::size_t n_blocks,
                                 const double* effects, std::size_t n_effects,
                                 std::size_t workspace_bytes);

}  // namespace powermean

#endif
