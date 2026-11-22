#ifndef POWERMEAN_K3_ECHELON_H
#define POWERMEAN_K3_ECHELON_H

// Derived from Powermean's _native/basis/echelon.cpp, 7b22dcf (see K3.md).
// Work on differences from each head's first row. Selecting the smallest
// remaining column gives the same ordered pivot set as full row reduction.
#include "k3_fraction.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <numeric>
#include <unordered_set>
#include <utility>
#include <vector>

namespace powermean::k3 {

inline unsigned lowest_set_bit(std::uint64_t value) {
#if defined(__GNUC__) || defined(__clang__)
    return static_cast<unsigned>(__builtin_ctzll(value));
#else
    unsigned bit = 0;
    while ((value & 1) == 0) {
        ++bit;
        value >>= 1;
    }
    return bit;
#endif
}

template <class Column>
bool full_tail_rank_mod2(std::size_t rows, std::size_t blocks, int heads, int width,
                         const Column& column,
                         std::size_t workspace = std::numeric_limits<std::size_t>::max()) {
    constexpr int max_tail = 512;
    constexpr int word_bits = 64;
    const int tail = width - heads;
    if (blocks < 3 || tail < 1 || tail > max_tail || heads > 8 * tail ||
        rows < static_cast<std::size_t>(width)) return false;

    const int words = (tail + word_bits - 1) / word_bits;
    const auto needed = static_cast<std::size_t>(heads) * sizeof(std::size_t) +
        (tail * words + max_tail / word_bits) * sizeof(std::uint64_t);
    if (needed > workspace) return false;
    const auto unseen = std::numeric_limits<std::size_t>::max();
    const std::size_t prefix = std::min(rows, static_cast<std::size_t>(heads + 8 * tail));
    std::vector<std::size_t> anchors(heads, unseen);
    std::vector<std::uint64_t> pivots(tail * words, 0);
    std::array<std::uint64_t, max_tail / word_bits> difference{};
    int retained = 0;
    for (std::size_t i = 0; i < prefix; ++i) {
        auto& anchor = anchors[column(0, i)];
        if (anchor == unseen) {
            anchor = i;
            continue;
        }
        std::fill_n(difference.begin(), words, 0);
        for (std::size_t block = 1; block < blocks; ++block) {
            const int next = column(block, i);
            const int previous = column(block, anchor);
            if (next == previous) continue;
            for (int mapped : {next, previous}) {
                if (mapped < 0) continue;
                const int bit = mapped - heads;
                difference[bit / word_bits] ^= std::uint64_t{1} << (bit % word_bits);
            }
        }
        bool inserted = false;
        for (int word = 0; word < words && !inserted; ++word) {
            while (difference[word] != 0) {
                const int pivot_index = word * word_bits + lowest_set_bit(difference[word]);
                auto* pivot = pivots.data() + pivot_index * words;
                if (pivot[word] == 0) {
                    std::copy_n(difference.begin() + word, words - word, pivot + word);
                    // The head indicators are independent. An odd full tail
                    // minor then proves that all canonical E columns survive.
                    if (++retained == tail) return true;
                    inserted = true;
                    break;
                }
                for (int j = word; j < words; ++j) difference[j] ^= pivot[j];
            }
        }
    }
    return false;
}

template <class Q>
using Row = std::vector<std::pair<int, Q>>;

struct UnrestrictedWorkspace {
    template <class T>
    void grow(std::vector<T>&, std::size_t) const {}

    void add_holder(std::unordered_set<int>& holders, int holder) const {
        holders.insert(holder);
    }

    void remove_holder(std::unordered_set<int>& holders, int holder) const {
        holders.erase(holder);
    }

    void clear_holders(std::unordered_set<int>& holders) const { holders.clear(); }
};

// O(1) column lookup for the one sparse row under modification.
template <class Q>
class Scratch {
public:
    explicit Scratch(int width) : position(width, -1) {}

    void load(const Row<Q>& row) {
        for (std::size_t i = 0; i < row.size(); ++i) {
            position[row[i].first] = static_cast<int>(i);
        }
    }

    void unload(const Row<Q>& row) {
        for (const auto& entry : row) position[entry.first] = -1;
    }

    bool has(int column) const { return position[column] >= 0; }

    template <class Workspace>
    void add(Row<Q>& row, int column, const Q& value, Workspace& workspace) {
        const int i = position[column];
        if (i < 0) {
            if (!is_zero(value)) {
                workspace.grow(row, row.size() + 1);
                position[column] = static_cast<int>(row.size());
                row.emplace_back(column, value);
            }
            return;
        }
        auto total = q_add(row[i].second, value);
        if (is_zero(total)) {
            erase(row, i);
        } else {
            row[i].second = std::move(total);
        }
    }

    Q pop(Row<Q>& row, int column) {
        const int i = position[column];
        Q value = std::move(row[i].second);
        erase(row, i);
        return value;
    }

private:
    void erase(Row<Q>& row, int i) {
        position[row[i].first] = -1;
        const std::size_t last = row.size() - 1;
        if (static_cast<std::size_t>(i) != last) {
            row[i] = std::move(row[last]);
            position[row[i].first] = i;
        }
        row.pop_back();
    }

    std::vector<int> position;
};

template <class Q, class Workspace = UnrestrictedWorkspace>
class Echelon : private Workspace {
public:
    explicit Echelon(int width, Workspace workspace = {})
        : Workspace(workspace), pivot(width, 0), tails(width), holders(width),
          scratch(width) {}

    template <class Column>
    void add_difference(std::size_t i, std::size_t anchor, std::size_t blocks,
                        int head_count, const Column& column) {
        // Mapped columns are disjoint by formula block, so each incidence
        // row contributes at most one of each column.
        for (std::size_t k = 1; k < blocks; ++k) {
            const int next = column(k, i);
            const int previous = column(k, anchor);
            if (next == previous) continue;
            if (next >= 0) scratch.add(difference, next - head_count, Q{1}, workspace());
            if (previous >= 0) scratch.add(difference, previous - head_count, Q{-1}, workspace());
        }
        reduce();
        if (difference.empty()) return;
        const int p = smallest_column();
        const Q lead = scratch.pop(difference, p);
        for (auto& entry : difference) entry.second = q_div(entry.second, lead);
        scratch.unload(difference);
        insert_pivot(p);
        difference.clear();
    }

    const std::vector<Row<Q>>& pivot_tails() const noexcept { return tails; }

    std::vector<int> selected(int head_count) const {
        std::vector<int> answer(head_count);
        std::iota(answer.begin(), answer.end(), 0);
        for (std::size_t p = 0; p < pivot.size(); ++p) {
            if (pivot[p]) answer.push_back(head_count + static_cast<int>(p));
        }
        return answer;
    }

private:
    Workspace& workspace() { return *this; }

    void reduce() {
        present.clear();
        this->grow(present, difference.size());
        for (const auto& entry : difference) {
            if (pivot[entry.first]) present.push_back(entry.first);
        }
        for (int p : present) {
            const Q minus = q_neg(scratch.pop(difference, p));
            for (const auto& entry : tails[p]) {
                scratch.add(difference, entry.first, q_mul(minus, entry.second), workspace());
            }
        }
    }

    int smallest_column() const {
        int p = difference.front().first;
        for (const auto& entry : difference) p = std::min(p, entry.first);
        return p;
    }

    void replace_tail(int p, int holder) {
        auto& old = tails[holder];
        scratch.load(old);
        const Q minus = q_neg(scratch.pop(old, p));
        for (const auto& entry : difference) {
            const int column = entry.first;
            const bool had = scratch.has(column);
            scratch.add(old, column, q_mul(minus, entry.second), workspace());
            const bool has = scratch.has(column);
            if (has && !had) this->add_holder(holders[column], holder);
            if (had && !has) this->remove_holder(holders[column], holder);
        }
        scratch.unload(old);
    }

    void insert_pivot(int p) {
        // The replacement tail contains no pivot columns, including p.
        for (int holder : holders[p]) replace_tail(p, holder);
        this->clear_holders(holders[p]);
        for (const auto& entry : difference) this->add_holder(holders[entry.first], p);
        this->grow(tails[p], difference.size());
        tails[p] = difference;
        pivot[p] = 1;
    }

    std::vector<unsigned char> pivot;
    std::vector<Row<Q>> tails;
    std::vector<std::unordered_set<int>> holders;
    Scratch<Q> scratch;
    Row<Q> difference;
    std::vector<int> present;
};

template <class State, class Column>
void reduce_incidence(State& state, std::vector<std::size_t>& first,
                      std::size_t rows, std::size_t blocks, int head_count,
                      const Column& column) {
    for (std::size_t i = 0; i < rows; ++i) {
        auto& anchor = first[column(0, i)];
        if (anchor == rows) {
            anchor = i;
        } else {
            state.add_difference(i, anchor, blocks, head_count, column);
        }
    }
}

template <class Q, class Column>
std::vector<int> reduce_incidence(std::size_t rows, std::size_t blocks,
                                  int head_count, int width, const Column& column) {
    std::vector<std::size_t> first(head_count, rows);
    Echelon<Q> state(width - head_count);
    reduce_incidence(state, first, rows, blocks, head_count, column);
    return state.selected(head_count);
}

template <class Column>
std::vector<int> ordered_pivots(std::size_t rows, std::size_t blocks,
                                int heads, int width, const Column& column, bool& multiprecision) {
    try {
        return reduce_incidence<Small>(rows, blocks, heads, width, column);
    } catch (const ArithmeticRange&) {
        multiprecision = true;
        return reduce_incidence<Big>(rows, blocks, heads, width, column);
    }
}

}  // namespace powermean::k3

#endif
