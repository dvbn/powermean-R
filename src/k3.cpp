// This file maps retained Powermean FE codes to the E and J coordinates.
// The mapping follows the frozen certificate_basis.py. The exact reducer
// and source provenance are in k3_echelon.h and k3_fraction.h.
#include "k3.h"
#include "k3_echelon.h"
#include "k3_layout.h"

#include <limits>
#include <new>
#include <numeric>
#include <stdexcept>

namespace powermean {
namespace {

using k3::Grouping;

class Components {
public:
    explicit Components(int size) : parent(size), height(size, 0) {
        std::iota(parent.begin(), parent.end(), 0);
    }

    int root(int vertex) {
        while (parent[vertex] != vertex) {
            parent[vertex] = parent[parent[vertex]];
            vertex = parent[vertex];
        }
        return vertex;
    }

    void join(int a, int b) {
        a = root(a);
        b = root(b);
        if (a == b) return;
        if (height[a] < height[b]) {
            parent[a] = b;
        } else {
            parent[b] = a;
            if (height[a] == height[b]) ++height[a];
        }
    }

private:
    std::vector<int> parent;
    std::vector<unsigned char> height;
};

BasisResult failure(BasisStatus status) {
    BasisResult result;
    result.status = status;
    return result;
}

BasisStatus compute_offsets(Grouping& groups) {
    int total = 0;
    for (std::size_t k = 0; k < groups.blocks; ++k) {
        const int count = groups.levels[k];
        if (count <= 0) return BasisStatus::invalid_shape;
        if (static_cast<std::size_t>(count) > groups.rows) return BasisStatus::missing_level;
        if (count > std::numeric_limits<int>::max() - total) return BasisStatus::capacity;
        groups.offset[k] = total;
        total += count;
    }
    groups.offset.back() = total;
    return BasisStatus::ok;
}

BasisStatus compute_components(const Grouping& groups, BasisResult& result) {
    const int total = groups.offset.back();
    Components graph(total);
    std::vector<unsigned char> seen(total, 0);
    for (std::size_t i = 0; i < groups.rows; ++i) {
        const int head = groups.code(0, i);
        if (head < 0 || head >= groups.levels[0]) return BasisStatus::invalid_code;
        seen[head] = 1;
        for (std::size_t k = 1; k < groups.blocks; ++k) {
            const int code = groups.code(k, i);
            if (code < 0 || code >= groups.levels[k]) return BasisStatus::invalid_code;
            const int vertex = groups.offset[k] + code;
            seen[vertex] = 1;
            graph.join(head, vertex);
        }
    }
    for (unsigned char value : seen) {
        if (value == 0) return BasisStatus::missing_level;
    }
    std::vector<int> component(total, -1);
    for (std::size_t i = 0; i < groups.rows; ++i) {
        int& found = component[graph.root(groups.code(0, i))];
        if (found < 0) found = result.n_components++;
    }
    result.component_by_level.resize(total);
    for (int v = 0; v < total; ++v) {
        result.component_by_level[v] = component[graph.root(v)];
    }
    return BasisStatus::ok;
}

int map_columns(const Grouping& groups, BasisResult& result) {
    result.column_by_level.resize(groups.offset.back(), -1);
    std::vector<std::size_t> last_block(result.n_components, 0);
    int next = 0;
    for (std::size_t k = 0; k < groups.blocks; ++k) {
        for (int v = groups.offset[k]; v < groups.offset[k + 1]; ++v) {
            auto& last = last_block[result.component_by_level[v]];
            // The first level per component of each non-head block is its reference.
            if (k == 0 || last == k) result.column_by_level[v] = next++;
            last = k;
        }
    }
    return next;
}

}  // namespace

namespace k3 {

BasisStatus grouping_shape(const int* codes, std::size_t rows,
                           const int* levels, std::size_t blocks) {
    if (!codes || !levels || rows == 0 || blocks == 0) return BasisStatus::invalid_shape;
    const auto max_cells = std::numeric_limits<std::size_t>::max() / sizeof(int);
    if (rows > max_cells / blocks || blocks > std::numeric_limits<int>::max()) {
        return BasisStatus::capacity;
    }
    return BasisStatus::ok;
}

BasisStatus map_grouping(Grouping& groups, BasisResult& result, int& width) {
    const BasisStatus offsets = compute_offsets(groups);
    if (offsets != BasisStatus::ok) return offsets;
    const BasisStatus components = compute_components(groups, result);
    if (components != BasisStatus::ok) return components;
    width = map_columns(groups, result);
    return BasisStatus::ok;
}

void map_selected_positions(int width, BasisResult& result) {
    std::vector<int> positions(width, -1);
    for (std::size_t j = 0; j < result.selected.size(); ++j) {
        positions[result.selected[j]] = static_cast<int>(j);
    }
    result.position_by_level.resize(result.column_by_level.size(), -1);
    for (std::size_t v = 0; v < result.column_by_level.size(); ++v) {
        const int column = result.column_by_level[v];
        if (column >= 0) result.position_by_level[v] = positions[column];
    }
}

}  // namespace k3

namespace {

BasisResult build_basis(Grouping& groups) {
    BasisResult result;
    int width;
    const BasisStatus mapped = k3::map_grouping(groups, result, width);
    if (mapped != BasisStatus::ok) return failure(mapped);
    const int head_count = groups.levels[0];
    const auto column = [&](std::size_t block, std::size_t row) {
        return result.column_by_level[groups.vertex(block, row)];
    };
    if (groups.blocks <= 2 || width == head_count ||
        k3::full_tail_rank_mod2(groups.rows, groups.blocks, head_count, width, column)) {
        // Connectivity suffices for at most two blocks. For more blocks,
        // a full modular tail proof likewise keeps every canonical column.
        result.selected.resize(width);
        std::iota(result.selected.begin(), result.selected.end(), 0);
    } else {
        result.selected = k3::ordered_pivots(
            groups.rows, groups.blocks, head_count, width, column, result.multiprecision
        );
    }
    k3::map_selected_positions(width, result);
    result.status = BasisStatus::ok;
    return result;
}

}  // namespace

BasisResult exact_fe_basis(const int* codes, std::size_t n_rows,
                           const int* level_counts, std::size_t n_blocks) {
    const BasisStatus shape = k3::grouping_shape(codes, n_rows, level_counts, n_blocks);
    if (shape != BasisStatus::ok) return failure(shape);
    try {
        Grouping groups{codes, n_rows, level_counts, n_blocks,
                        std::vector<int>(n_blocks + 1)};
        return build_basis(groups);
    } catch (const std::bad_alloc&) {
        return failure(BasisStatus::allocation);
    } catch (const std::length_error&) {
        return failure(BasisStatus::capacity);
    }
}

}  // namespace powermean
