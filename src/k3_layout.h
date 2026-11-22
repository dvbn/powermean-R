#ifndef POWERMEAN_K3_LAYOUT_H
#define POWERMEAN_K3_LAYOUT_H

#include "k3.h"

namespace powermean::k3 {

struct Grouping {
    const int* codes;
    std::size_t rows;
    const int* levels;
    std::size_t blocks;
    std::vector<int> offset;

    int code(std::size_t block, std::size_t row) const {
        return codes[block * rows + row];
    }

    int vertex(std::size_t block, std::size_t row) const {
        return offset[block] + code(block, row);
    }
};

BasisStatus grouping_shape(const int* codes, std::size_t rows,
                           const int* levels, std::size_t blocks);
BasisStatus map_grouping(Grouping& groups, BasisResult& result, int& width);
void map_selected_positions(int width, BasisResult& result);

}  // namespace powermean::k3

#endif
