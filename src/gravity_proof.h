#ifndef POWERMEAN_GRAVITY_PROOF_HPP
#define POWERMEAN_GRAVITY_PROOF_HPP

#include <cstddef>
#include <cstdint>

namespace powermean {

template<class T>
struct GravityView {
    const T* data;
    std::size_t size;

    const T& operator[](std::size_t row) const { return data[row]; }
};

struct GravityInput {
    GravityView<int> pair;
    GravityView<int> exporter;
    GravityView<int> importer;
    GravityView<int> year;
    GravityView<double> y;
    GravityView<double> x;
    int n_pairs;
    int n_exporters;
    int n_importers;
    int n_years;
};

enum class GravityStatus {
    proved,
    invalid_input,
    unsupported_layout,
    duplicate_cell,
    inconsistent_pair,
    unsupported_slope,
    missing_positive_pair,
    disconnected_positive_graph,
    no_contrast,
    weak_contrast,
    unavailable
};

struct GravityProof {
    GravityStatus status;
    double lower_bound = 0;
    std::int64_t contrast = 0;
    std::int64_t cycle_length_bound = 0;
    int pair = -1;
    int year = -1;
};

GravityProof gravity_positive_proof(const GravityInput& input);

} // namespace powermean

#endif
