// Exact witness action for the frozen FE quotient, using K3's sparse RREF.
#include "binary_round.h"
#include "k3_echelon.h"
#include "k3_layout.h"

#include <cmath>
#include <map>
#include <new>
#include <stdexcept>
#include <type_traits>

namespace powermean {
namespace {

using k3::Big;
using Integer = boost::multiprecision::cpp_int;
using Witness = std::map<int, Big>;
using TailTerm = std::pair<int, const Big*>;
using RevealedTail = std::vector<std::vector<TailTerm>>;
static_assert(std::numeric_limits<double>::is_iec559 &&
              std::numeric_limits<double>::digits == 53,
              "The exact basis action requires IEEE binary64");

struct ResourceLimit {};

class Budget {
public:
    explicit Budget(std::size_t bytes)
        : remaining(bytes), unlimited(bytes == std::numeric_limits<std::size_t>::max()) {}

    void take(std::size_t count, std::size_t bytes) {
        if (unlimited || !count) return;
        if (bytes > remaining / count) throw ResourceLimit{};
        remaining -= count * bytes;
    }

    void release(std::size_t count, std::size_t bytes) {
        if (!unlimited) remaining += count * bytes;
    }

    std::size_t available() const { return remaining; }

private:
    std::size_t remaining;
    bool unlimited;
};

// Reserve row capacities as they grow. The extra per-rational allowance
// uses the Hadamard bound for minors of the tail incidence matrix.
class SparseWorkspace {
public:
    SparseWorkspace(Budget& limit, int tail, std::size_t blocks)
        : budget(limit), rational_bytes(rational_allowance(tail, blocks)) {}

    template <class T>
    void grow(std::vector<T>& values, std::size_t size) {
        const std::size_t prior = values.capacity();
        if (size <= prior) return;
        const std::size_t target = std::max(size, prior + prior / 2);
        std::size_t bytes = sizeof(T);
        if constexpr (std::is_same_v<T, k3::Row<Big>::value_type>) bytes += rational_bytes;
        budget.take(target, bytes);
        values.reserve(target);
        budget.take(values.capacity() - target, bytes);
        budget.release(prior, bytes);
    }

    void add_holder(std::unordered_set<int>& holders, int holder) {
        const std::size_t required = holders.size() + 1;
        const auto old = holders.bucket_count();
        if (old <= 1 || required > old) {
            const auto requested = 2 * required;
            budget.take(requested, sizeof(void*));
            holders.reserve(requested);
            budget.take(holders.bucket_count() - requested, sizeof(void*));
            if (old > 1) budget.release(old, sizeof(void*));
        }
        budget.take(1, holder_bytes);
        holders.insert(holder);
    }

    void remove_holder(std::unordered_set<int>& holders, int holder) {
        budget.release(holders.erase(holder), holder_bytes);
    }

    void clear_holders(std::unordered_set<int>& holders) {
        budget.release(holders.size(), holder_bytes);
        holders.clear();
    }

    std::size_t witness_entry_bytes() const {
        return sizeof(Witness::value_type) + 4 * sizeof(void*) + rational_bytes;
    }

private:
    static std::size_t rational_allowance(int tail, std::size_t blocks) {
        const auto row_terms = std::min<std::size_t>(tail, 2 * (blocks - 1));
        std::size_t log_terms = 1;
        for (std::size_t terms = row_terms; terms > 1; terms >>= 1) ++log_terms;
        // One allowance per entry covers numerator and denominator storage.
        // Reserve extra space for products and transient normalization limbs.
        if (tail && log_terms > (std::numeric_limits<std::size_t>::max() - 64) /
                                  static_cast<std::size_t>(tail)) throw ResourceLimit{};
        return 64 + static_cast<std::size_t>(tail) * log_terms;
    }

    static constexpr auto holder_bytes = sizeof(int) + 3 * sizeof(void*);
    Budget& budget;
    std::size_t rational_bytes;
};

void update(double& value, double term) {
    value += term;
    if (!std::isfinite(value)) throw std::overflow_error("basis witness action");
}

class TemporaryWitness {
public:
    TemporaryWitness(Budget& workspace, std::size_t bytes)
        : budget(workspace), entry_bytes(bytes) {}

    ~TemporaryWitness() { budget.release(entries.size(), entry_bytes); }

    void add(int index, const Big& amount) {
        if (amount == 0) return;
        auto where = entries.lower_bound(index);
        if (where == entries.end() || where->first != index) {
            budget.take(1, entry_bytes);
            entries.emplace_hint(where, index, amount);
        } else {
            where->second += amount;
            if (where->second == 0) {
                entries.erase(where);
                budget.release(1, entry_bytes);
            }
        }
    }

    void apply(std::vector<double>& coefficients, double effect) const {
        for (const auto& term : entries) {
            // Materialize the multiplication at binary64 precision before
            // the addition, even when compiler FMA contraction is enabled.
            volatile double product = effect * nearest_binary64(term.second);
            update(coefficients[term.first], product);
        }
    }

private:
    Budget& budget;
    std::size_t entry_bytes;
    Witness entries;
};

std::vector<int> level_links(const BasisResult& basis, const k3::Grouping& groups) {
    std::vector<int> next(basis.component_by_level.size(), -1);
    std::vector<int> first(basis.n_components, -1);
    for (std::size_t b = 0; b < groups.blocks; ++b) {
        for (int v = groups.offset[b]; v < groups.offset[b + 1]; ++v) {
            first[basis.component_by_level[v]] = -1;
        }
        for (int v = groups.offset[b + 1]; v-- > groups.offset[b];) {
            int& root = first[basis.component_by_level[v]];
            next[v] = root;
            root = v;
        }
    }
    return next;
}

class WitnessGeometry {
public:
    WitnessGeometry(const k3::Grouping& g, const BasisResult& b,
                    const std::vector<std::size_t>& starts, const RevealedTail& tails)
        : groups(g), basis(b), anchors(starts), relations(tails),
          next(level_links(b, g)), head(b.n_components, -1) {
        for (int h = groups.levels[0]; h-- > 0;) head[basis.component_by_level[h]] = h;
    }

    void independent_reference(std::vector<double>& result, int v, double effect) const {
        for (int h = head[basis.component_by_level[v]]; h >= 0; h = next[h]) {
            update(result[basis.position_by_level[h]], effect);
        }
        for (int peer = next[v]; peer >= 0; peer = next[peer]) {
            update(result[basis.position_by_level[peer]], -effect);
        }
    }

    void reference(TemporaryWitness& witness, int v) const {
        for (int h = head[basis.component_by_level[v]]; h >= 0; h = next[h]) {
            witness.add(basis.position_by_level[h], Big{1});
        }
        for (int peer = next[v]; peer >= 0; peer = next[peer]) {
            const int slot = basis.position_by_level[peer];
            if (slot >= 0) witness.add(slot, Big{-1});
            else dependent(witness, peer, true);
        }
    }

    void dependent(TemporaryWitness& witness, int v, bool reversed = false) const {
        const auto& relation = relations[basis.column_by_level[v] - groups.levels[0]];
        for (const auto& term : relation) {
            witness.add(term.first, reversed ? k3::q_neg(*term.second) : *term.second);
        }
        for (int h = head[basis.component_by_level[v]]; h >= 0; h = next[h]) {
            Big value = 0;
            for (std::size_t b = 1; b < groups.blocks; ++b) {
                const int anchor_vertex = groups.vertex(b, anchors[h]);
                if (anchor_vertex == v) ++value;
                const int slot = basis.position_by_level[anchor_vertex];
                if (slot < 0) continue;
                const auto found = std::lower_bound(relation.begin(), relation.end(), slot,
                    [](const TailTerm& term, int target) { return term.first < target; });
                if (found != relation.end() && found->first == slot) value -= *found->second;
            }
            witness.add(basis.position_by_level[h], reversed ? k3::q_neg(value) : value);
        }
    }

private:
    const k3::Grouping& groups;
    const BasisResult& basis;
    const std::vector<std::size_t>& anchors;
    const RevealedTail& relations;
    std::vector<int> next;
    std::vector<int> head;
};

RevealedTail expose_relations(const std::vector<k3::Row<Big>>& tails,
                              const BasisResult& basis, int heads, Budget& budget,
                              SparseWorkspace& workspace) {
    budget.take(tails.size(), sizeof(RevealedTail::value_type));
    RevealedTail relations(tails.size());
    for (std::size_t slot = heads; slot < basis.selected.size(); ++slot) {
        const int pivot = basis.selected[slot] - heads;
        for (const auto& entry : tails[pivot]) {
            auto& relation = relations[entry.first];
            workspace.grow(relation, relation.size() + 1);
            relation.emplace_back(static_cast<int>(slot), &entry.second);
        }
    }
    return relations;
}

BasisActionResult apply_action(const k3::Grouping& groups, BasisResult& basis, int width,
                               const double* effects, Budget& budget) {
    const int heads = groups.levels[0];
    const auto column = [&](std::size_t b, std::size_t i) {
        return basis.column_by_level[groups.vertex(b, i)];
    };
    const bool full = groups.blocks <= 2 || width == heads ||
        k3::full_tail_rank_mod2(groups.rows, groups.blocks, heads, width, column, budget.available());
    std::vector<std::size_t> anchors;
    const int tail = full ? 0 : width - heads;
    SparseWorkspace workspace(budget, tail, groups.blocks);
    budget.take(tail, sizeof(k3::Row<Big>) + sizeof(std::unordered_set<int>) +
                      sizeof(int) + sizeof(unsigned char));
    k3::Echelon<Big, SparseWorkspace> reducer(tail, workspace);
    if (full) {
        basis.selected.resize(width);
        std::iota(basis.selected.begin(), basis.selected.end(), 0);
    } else {
        anchors.assign(heads, groups.rows);
        k3::reduce_incidence(reducer, anchors, groups.rows, groups.blocks, heads, column);
        basis.selected = reducer.selected(heads);
    }
    const bool aliased = static_cast<int>(basis.selected.size()) != width;
    k3::map_selected_positions(width, basis);
    const auto relations = aliased
        ? expose_relations(reducer.pivot_tails(), basis, heads, budget, workspace)
        : RevealedTail{};
    const WitnessGeometry geometry(groups, basis, anchors, relations);
    std::vector<double> coefficients(basis.selected.size());
    for (std::size_t b = 0; b < groups.blocks; ++b) {
        for (int v = groups.offset[b]; v < groups.offset[b + 1]; ++v) {
            const int slot = basis.position_by_level[v];
            if (slot >= 0) update(coefficients[slot], effects[v]);
        }
        for (int v = groups.offset[b]; v < groups.offset[b + 1]; ++v) {
            if (basis.position_by_level[v] >= 0 || effects[v] == 0) continue;
            if (!aliased) {
                geometry.independent_reference(coefficients, v, effects[v]);
                continue;
            }
            TemporaryWitness witness(budget, workspace.witness_entry_bytes());
            if (basis.column_by_level[v] < 0) geometry.reference(witness, v);
            else geometry.dependent(witness, v);
            witness.apply(coefficients, effects[v]);
        }
    }
    return {BasisStatus::ok, std::move(basis.position_by_level), std::move(coefficients)};
}

BasisActionResult action_failure(BasisStatus status) {
    return {status, {}, {}};
}

}  // namespace

BasisActionResult exact_fe_action(const int* codes, std::size_t n_rows,
                                 const int* level_counts, std::size_t n_blocks,
                                 const double* effects, std::size_t n_effects,
                                 std::size_t workspace_bytes) {
    const auto shape = k3::grouping_shape(codes, n_rows, level_counts, n_blocks);
    if (shape != BasisStatus::ok) return action_failure(shape);
    if (!effects) return action_failure(BasisStatus::invalid_shape);
    std::size_t total = 0;
    for (std::size_t b = 0; b < n_blocks; ++b) {
        const int count = level_counts[b];
        if (count <= 0) return action_failure(BasisStatus::invalid_shape);
        if (static_cast<std::size_t>(count) > n_rows) return action_failure(BasisStatus::missing_level);
        if (static_cast<std::size_t>(count) > std::numeric_limits<int>::max() - total) {
            return action_failure(BasisStatus::capacity);
        }
        total += count;
    }
    if (total != n_effects) return action_failure(BasisStatus::invalid_shape);
    for (std::size_t v = 0; v < total; ++v) {
        if (!std::isfinite(effects[v])) return action_failure(BasisStatus::not_representable);
    }
    try {
        Budget budget(workspace_bytes);
        budget.take(n_blocks + 1, sizeof(int));
        // Linear grouping, pivot maps, level links, anchors and coefficients.
        budget.take(total, 8 * sizeof(int) + sizeof(std::size_t) + sizeof(double));
        k3::Grouping groups{codes, n_rows, level_counts, n_blocks, std::vector<int>(n_blocks + 1)};
        BasisResult basis;
        int width;
        const auto mapped = k3::map_grouping(groups, basis, width);
        if (mapped != BasisStatus::ok) return action_failure(mapped);
        return apply_action(groups, basis, width, effects, budget);
    } catch (const ResourceLimit&) {
        return action_failure(BasisStatus::resource_limit);
    } catch (const std::overflow_error&) {
        return action_failure(BasisStatus::not_representable);
    } catch (const std::bad_alloc&) {
        return action_failure(BasisStatus::allocation);
    } catch (const std::length_error&) {
        return action_failure(BasisStatus::capacity);
    }
}

}  // namespace powermean
