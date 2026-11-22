#include "gravity_proof.h"

#include <cmath>
#include <limits>
#include <new>
#include <stdexcept>
#include <vector>

namespace powermean {
namespace {

using Status = GravityStatus;

constexpr signed char unobserved = -1;
constexpr signed char zero_outcome = -2;
constexpr double rank_threshold = 1e-9;

struct Pair {
    int exporter = -1;
    int importer = -1;
    int positives = 0;
    int next_exporter = -1;
    int next_importer = -1;
};

struct Vertex {
    int first_pair = -1;
    int parent = -1;
    int parent_pair = -1;
    int depth = 0;
    std::int64_t potential = 0;
};

template<class T>
bool aligned(GravityView<T> view, std::size_t size) {
    return view.data != nullptr && view.size == size;
}

bool in_range(int code, int count) {
    return code >= 0 && code < count;
}

double conservative_bound(std::int64_t contrast, std::int64_t length,
                          std::size_t ones) {
    // Each factor converts exactly to double. Round products and the quotient
    // outwards so a retained column has a strict lower bound above the cutoff.
    const double difference = static_cast<double>(contrast);
    const double numerator = std::nextafter(difference * difference, 0.0);
    const double denominator = std::nextafter(
        static_cast<double>(2 * length) * static_cast<double>(ones),
        std::numeric_limits<double>::infinity());
    return std::nextafter(numerator / denominator, 0.0);
}

GravityProof check_gravity(const GravityInput& input) {
    const std::size_t size = input.pair.size;
    if (size == 0 || size > std::numeric_limits<int>::max() ||
        !aligned(input.pair, size) || !aligned(input.exporter, size) ||
        !aligned(input.importer, size) || !aligned(input.year, size) ||
        !aligned(input.y, size) || !aligned(input.x, size) ||
        input.n_pairs < 1 || input.n_exporters < 1 || input.n_importers < 1 ||
        input.n_years < 1) {
        return {Status::invalid_input};
    }
    if (input.n_years < 2 || input.n_exporters > input.n_pairs ||
        input.n_importers > input.n_pairs ||
        input.n_exporters > std::numeric_limits<int>::max() - input.n_importers ||
        static_cast<std::size_t>(input.n_pairs) > size ||
        static_cast<std::size_t>(input.n_years) > size ||
        static_cast<std::size_t>(input.n_pairs) > 2 * size / input.n_years) {
        return {Status::unsupported_layout};
    }

    const std::size_t years = input.n_years;
    const int vertex_count = input.n_exporters + input.n_importers;
    std::vector<signed char> cells(static_cast<std::size_t>(input.n_pairs) * years,
                                   unobserved);
    std::vector<Pair> pairs(input.n_pairs);
    std::size_t ones = 0;
    for (std::size_t row = 0; row < size; ++row) {
        const int pair_id = input.pair[row];
        const int exporter = input.exporter[row];
        const int importer = input.importer[row];
        const int year = input.year[row];
        const double y = input.y[row];
        const double x = input.x[row];
        if (!in_range(pair_id, input.n_pairs) ||
            !in_range(exporter, input.n_exporters) ||
            !in_range(importer, input.n_importers) ||
            !in_range(year, input.n_years) || !std::isfinite(y) || y < 0) {
            return {Status::invalid_input};
        }
        if (x != 0 && x != 1) {
            return {Status::unsupported_slope};
        }
        ones += x == 1;
        Pair& pair = pairs[pair_id];
        const int importer_vertex = input.n_exporters + importer;
        if (pair.exporter == -1) {
            pair.exporter = exporter;
            pair.importer = importer_vertex;
        } else if (pair.exporter != exporter || pair.importer != importer_vertex) {
            return {Status::inconsistent_pair};
        }

        auto& cell = cells[static_cast<std::size_t>(pair_id) * years + year];
        if (cell != unobserved) {
            return {Status::duplicate_cell};
        }
        cell = y > 0 ? static_cast<signed char>(x) : zero_outcome;
        pair.positives += y > 0;
    }

    std::vector<Vertex> vertices(vertex_count);
    std::vector<int> complete_pairs;
    complete_pairs.reserve(input.n_pairs);
    for (int id = 0; id < input.n_pairs; ++id) {
        Pair& pair = pairs[id];
        if (pair.positives == 0) {
            return {Status::missing_positive_pair};
        }
        pair.next_exporter = vertices[pair.exporter].first_pair;
        vertices[pair.exporter].first_pair = id;
        pair.next_importer = vertices[pair.importer].first_pair;
        vertices[pair.importer].first_pair = id;
        if (static_cast<std::size_t>(pair.positives) == years) {
            complete_pairs.push_back(id);
        }
    }

    // A gravity pair label must identify one distinct exporter/importer pair.
    std::vector<int> last_exporter(input.n_importers, -1);
    for (int exporter = 0; exporter < input.n_exporters; ++exporter) {
        for (int id = vertices[exporter].first_pair; id != -1;
             id = pairs[id].next_exporter) {
            int& last = last_exporter[pairs[id].importer - input.n_exporters];
            if (last == exporter) {
                return {Status::inconsistent_pair};
            }
            last = exporter;
        }
    }

    // A spanning tree of pairs observed positively on every date suffices to
    // extend any positive-sample FE null relation to all observed rows.
    std::vector<int> order{0};
    order.reserve(vertex_count);
    vertices[0].parent = 0;
    for (std::size_t index = 0; index < order.size(); ++index) {
        const int vertex = order[index];
        const bool exporter_side = vertex < input.n_exporters;
        int id = vertices[vertex].first_pair;
        while (id != -1) {
            const Pair& pair = pairs[id];
            const int neighbor = exporter_side ? pair.importer : pair.exporter;
            if (static_cast<std::size_t>(pair.positives) == years &&
                vertices[neighbor].parent == -1) {
                Vertex& next = vertices[neighbor];
                next.parent = vertex;
                next.parent_pair = id;
                next.depth = vertices[vertex].depth + 1;
                order.push_back(neighbor);
            }
            id = exporter_side ? pair.next_exporter : pair.next_importer;
        }
    }
    if (order.size() != static_cast<std::size_t>(vertex_count)) {
        return {Status::disconnected_positive_graph};
    }
    if (ones == 0 || complete_pairs.size() == order.size() - 1) {
        return {Status::no_contrast};
    }

    GravityProof best{Status::no_contrast};
    for (int year = 1; year < input.n_years; ++year) {
        const auto date_difference = [&](int pair) {
            const std::size_t offset = static_cast<std::size_t>(pair) * years;
            return cells[offset + year] - cells[offset];
        };
        for (std::size_t index = 1; index < order.size(); ++index) {
            Vertex& vertex = vertices[order[index]];
            vertex.potential = date_difference(vertex.parent_pair) -
                               vertices[vertex.parent].potential;
        }
        for (int id : complete_pairs) {
            const Pair& pair = pairs[id];
            const Vertex& exporter = vertices[pair.exporter];
            const Vertex& importer = vertices[pair.importer];
            const std::int64_t contrast = date_difference(id) -
                                          exporter.potential - importer.potential;
            if (contrast == 0) {
                continue;
            }

            // The tree path and this edge form a cycle. Alternating signs on
            // two dates annihilate every FE and give exactly this x contrast.
            // Shared tree ancestors shorten the cycle, making this bound safe.
            const std::int64_t length =
                static_cast<std::int64_t>(exporter.depth) + importer.depth + 1;
            const double bound = conservative_bound(contrast, length, ones);
            if (bound > rank_threshold) {
                return {Status::proved, bound, contrast, length, id, year};
            }
            if (bound > best.lower_bound) {
                best = {Status::weak_contrast, bound, contrast, length, id, year};
            }
        }
    }
    return best;
}

} // namespace

GravityProof gravity_positive_proof(const GravityInput& input) {
    try {
        return check_gravity(input);
    } catch (const std::bad_alloc&) {
        return {Status::unavailable};
    } catch (const std::length_error&) {
        return {Status::unavailable};
    }
}

} // namespace powermean
