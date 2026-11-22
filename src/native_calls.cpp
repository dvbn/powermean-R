#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include "native_calls.h"
#include "gravity_proof.h"
#include "k2.h"
#include "k3.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace {

struct Dimensions {
    std::size_t rows;
    std::size_t columns;
};

Dimensions matrix_size(SEXP matrix, SEXPTYPE type) {
    SEXP dim = Rf_getAttrib(matrix, R_DimSymbol);
    if (TYPEOF(matrix) != type || TYPEOF(dim) != INTSXP || Rf_xlength(dim) != 2 ||
        INTEGER(dim)[0] <= 0 || INTEGER(dim)[1] <= 0) {
        Rf_error("Native kernels require a nonempty matrix of the expected storage type.");
    }
    return {static_cast<std::size_t>(INTEGER(dim)[0]),
            static_cast<std::size_t>(INTEGER(dim)[1])};
}

std::size_t level_total(SEXP levels, std::size_t blocks, std::size_t* offsets = nullptr) {
    if (TYPEOF(levels) != INTSXP || Rf_xlength(levels) != static_cast<R_xlen_t>(blocks)) {
        Rf_error("Every fixed-effect block must have an integer level count.");
    }
    std::size_t total = 0;
    if (offsets != nullptr) offsets[0] = 0;
    for (std::size_t block = 0; block < blocks; ++block) {
        const int count = INTEGER(levels)[block];
        if (count <= 0 || total > static_cast<std::size_t>(INT_MAX - count)) {
            Rf_error("Fixed-effect level counts exceed the native range.");
        }
        total += count;
        if (offsets != nullptr) offsets[block + 1] = total;
    }
    return total;
}

bool valid_budget(double value) {
    return std::isfinite(value) && value == std::floor(value) && value >= 0 &&
        value <= 0x1p53;
}

}  // namespace

SEXP pm_design_scale_status(SEXP design) {
    const Dimensions dim = matrix_size(design, REALSXP);
    const double* column = REAL(design);
    bool rounded = false;
    for (std::size_t j = 0; j < dim.columns; ++j, column += dim.rows) {
        double largest = 0;
        double smallest = std::numeric_limits<double>::infinity();
        for (std::size_t i = 0; i < dim.rows; ++i) {
            const double magnitude = std::abs(column[i]);
            largest = std::max(largest, magnitude);
            if (magnitude > 0) smallest = std::min(smallest, magnitude);
        }
        int exponent = 0;
        std::frexp(largest, &exponent);
        // Dividing by 2^(exponent - 1) rounds to zero through the
        // midpoint between zero and the smallest subnormal.
        if (exponent > 1 && smallest <= std::scalbn(
                std::numeric_limits<double>::denorm_min(), exponent - 2)) {
            return Rf_ScalarInteger(1);
        }
        // Normal quotients by a power of two are exactly reversible.
        const double first_normal = std::scalbn(std::numeric_limits<double>::min(), exponent - 1);
        if (!rounded && exponent > 1 && smallest < first_normal) {
            for (std::size_t i = 0; i < dim.rows; ++i) {
                const double value = column[i];
                if (std::abs(value) < first_normal &&
                    std::scalbn(std::scalbn(value, 1 - exponent), exponent - 1) != value) {
                    rounded = true;
                    break;
                }
            }
        }
    }
    return Rf_ScalarInteger(rounded ? 2 : 0);
}

SEXP project_call(SEXP codes, SEXP levels, SEXP weight, SEXP rhs, SEXP control) {
    const auto groups = matrix_size(codes, INTSXP);
    const auto design = matrix_size(rhs, REALSXP);
    if (groups.rows != design.rows || (weight != R_NilValue &&
        (TYPEOF(weight) != REALSXP || Rf_xlength(weight) != static_cast<R_xlen_t>(groups.rows)))) {
        Rf_error("Rows, group codes and projection weights must align.");
    }
    const auto controls = Rf_xlength(control);
    if (TYPEOF(control) != REALSXP || (controls != 3 && controls != 4) ||
        !valid_budget(REAL(control)[2]) ||
        (controls == 4 && REAL(control)[3] != 0 && REAL(control)[3] != 1)) {
        Rf_error("Projection controls require two tolerances and an integer iteration budget.");
    }
    const auto method = controls == 3 ? powermean::ProjectionMethod::native_iterations :
        REAL(control)[3] == 1 ? powermean::ProjectionMethod::accelerated_sweeps :
                               powermean::ProjectionMethod::forward_sweeps;
    const auto n_offsets = groups.columns + 1;
    auto* offsets = reinterpret_cast<std::size_t*>(R_alloc(n_offsets, sizeof(std::size_t)));
    level_total(levels, groups.columns, offsets);
    const powermean::ProjectionShape shape{
        groups.rows, groups.columns, design.columns, groups.rows, design.rows
    };
    powermean::ProjectionSizes sizes{};
    if (powermean::projection_sizes(shape, {offsets, n_offsets}, sizes) !=
        powermean::ProjectionStatus::converged ||
        sizes.workspace > static_cast<std::size_t>(R_XLEN_T_MAX)) {
        Rf_error("Projection workspace sizes exceed the native range.");
    }
    auto* work = reinterpret_cast<double*>(R_alloc(sizes.workspace, sizeof(double)));
    auto* diagnostics = reinterpret_cast<powermean::ProjectionDiagnostic*>(
        R_alloc(sizes.diagnostics, sizeof(powermean::ProjectionDiagnostic))
    );
    SEXP values = PROTECT(Rf_duplicate(rhs));
    const auto status = powermean::project_fixed_effects(
        shape, {offsets, n_offsets}, {INTEGER(codes), sizes.codes},
        {weight == R_NilValue ? nullptr : REAL(weight), weight == R_NilValue ? 0 : groups.rows},
        {REAL(values), sizes.values},
        {REAL(control)[0], REAL(control)[1], static_cast<std::uint64_t>(REAL(control)[2]), method},
        {work, sizes.workspace}, {diagnostics, sizes.diagnostics}
    );
    SEXP output = PROTECT(Rf_allocVector(VECSXP, 4));
    SET_VECTOR_ELT(output, 0, values);
    SET_VECTOR_ELT(output, 1, Rf_ScalarInteger(static_cast<int>(status)));
    if (status == powermean::ProjectionStatus::converged) {
        const bool wide_budget = REAL(control)[2] > static_cast<double>(INT_MAX);
        SEXP iterations = Rf_allocVector(wide_budget ? REALSXP : INTSXP, design.columns);
        SET_VECTOR_ELT(output, 2, iterations);
        SEXP residual = Rf_allocVector(REALSXP, design.columns);
        SET_VECTOR_ELT(output, 3, residual);
        for (std::size_t col = 0; col < design.columns; ++col) {
            if (wide_budget) REAL(iterations)[col] = static_cast<double>(diagnostics[col].iterations);
            else INTEGER(iterations)[col] = static_cast<int>(diagnostics[col].iterations);
            REAL(residual)[col] = diagnostics[col].max_abs_group_mean;
        }
    }
    UNPROTECT(2);
    return output;
}

SEXP project_effect_call(SEXP codes, SEXP levels, SEXP weight, SEXP rhs, SEXP control, SEXP initial) {
    const auto groups = matrix_size(codes, INTSXP);
    const auto design = matrix_size(rhs, REALSXP);
    if (groups.rows != design.rows || (weight != R_NilValue &&
        (TYPEOF(weight) != REALSXP || Rf_xlength(weight) != static_cast<R_xlen_t>(groups.rows)))) {
        Rf_error("Rows, group codes and projection weights must align.");
    }
    const auto controls = Rf_xlength(control);
    if (TYPEOF(control) != REALSXP || (controls != 3 && controls != 4) ||
        !valid_budget(REAL(control)[2]) ||
        (controls == 4 && REAL(control)[3] != 0 && REAL(control)[3] != 1)) {
        Rf_error("Projection controls require two tolerances and an integer iteration budget.");
    }
    const auto method = controls == 3 ? powermean::ProjectionMethod::native_iterations :
        REAL(control)[3] == 1 ? powermean::ProjectionMethod::accelerated_sweeps :
                               powermean::ProjectionMethod::forward_sweeps;
    const auto n_offsets = groups.columns + 1;
    auto* offsets = reinterpret_cast<std::size_t*>(R_alloc(n_offsets, sizeof(std::size_t)));
    const auto total = level_total(levels, groups.columns, offsets);
    if (initial != R_NilValue) {
        const auto seed = matrix_size(initial, REALSXP);
        if (seed.rows != total || seed.columns != design.columns) {
            Rf_error("Initial effects need one row per fixed-effect level and one column per RHS.");
        }
    }
    const powermean::ProjectionShape shape{
        groups.rows, groups.columns, design.columns, groups.rows, design.rows
    };
    powermean::ProjectionEffectSizes sizes{};
    if (powermean::projection_effect_sizes(shape, {offsets, n_offsets}, sizes) !=
        powermean::ProjectionStatus::converged ||
        sizes.workspace > static_cast<std::size_t>(R_XLEN_T_MAX) ||
        sizes.effects > static_cast<std::size_t>(R_XLEN_T_MAX)) {
        Rf_error("Effect projection workspace sizes exceed the native range.");
    }
    auto* work = reinterpret_cast<double*>(R_alloc(sizes.workspace, sizeof(double)));
    auto* diagnostics = reinterpret_cast<powermean::ProjectionEffectDiagnostic*>(
        R_alloc(design.columns, sizeof(powermean::ProjectionEffectDiagnostic))
    );

    SEXP values = PROTECT(Rf_allocMatrix(REALSXP, design.rows, design.columns));
    SEXP effects = PROTECT(Rf_allocMatrix(REALSXP, total, design.columns));
    Rf_setAttrib(values, R_DimNamesSymbol, Rf_getAttrib(rhs, R_DimNamesSymbol));
    const auto status = powermean::project_fixed_effect_candidates(
        shape, {offsets, n_offsets}, {INTEGER(codes), sizes.projection.codes},
        {weight == R_NilValue ? nullptr : REAL(weight), weight == R_NilValue ? 0 : groups.rows},
        {REAL(rhs), sizes.projection.values},
        {initial == R_NilValue ? nullptr : REAL(initial), initial == R_NilValue ? 0 : sizes.effects},
        {REAL(values), sizes.projection.values}, {REAL(effects), sizes.effects},
        {REAL(control)[0], REAL(control)[1], static_cast<std::uint64_t>(REAL(control)[2]), method},
        {work, sizes.workspace}, {diagnostics, design.columns}
    );

    SEXP output = PROTECT(Rf_allocVector(VECSXP, 7));
    SET_VECTOR_ELT(output, 0, values);
    SET_VECTOR_ELT(output, 1, effects);
    SET_VECTOR_ELT(output, 2, Rf_ScalarInteger(static_cast<int>(status)));
    SEXP iterations = Rf_allocVector(REALSXP, design.columns);
    SET_VECTOR_ELT(output, 3, iterations);
    SEXP maximum = Rf_allocVector(REALSXP, design.columns);
    SET_VECTOR_ELT(output, 4, maximum);
    SEXP endings = Rf_allocVector(INTSXP, design.columns);
    SET_VECTOR_ELT(output, 5, endings);
    SEXP formed = Rf_allocVector(LGLSXP, design.columns);
    SET_VECTOR_ELT(output, 6, formed);
    for (std::size_t col = 0; col < design.columns; ++col) {
        const auto& diagnostic = diagnostics[col];
        REAL(iterations)[col] = static_cast<double>(diagnostic.projection.iterations);
        const double group_mean = diagnostic.projection.max_abs_group_mean;
        REAL(maximum)[col] = std::isfinite(group_mean) ? group_mean : NA_REAL;
        INTEGER(endings)[col] = static_cast<int>(diagnostic.projection.status);
        LOGICAL(formed)[col] = diagnostic.formed;
        if (!diagnostic.formed) {
            std::fill_n(REAL(values) + col * design.rows, design.rows, NA_REAL);
            std::fill_n(REAL(effects) + col * total, total, NA_REAL);
        }
    }
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 7));
    const char* const fields[] = {
        "values", "effects", "status", "iterations", "group_residual", "column_status", "formed"
    };
    for (std::size_t item = 0; item < 7; ++item) {
        SET_STRING_ELT(names, item, Rf_mkChar(fields[item]));
    }
    Rf_setAttrib(output, R_NamesSymbol, names);
    UNPROTECT(4);
    return output;
}

SEXP basis_call(SEXP codes, SEXP levels) {
    const auto shape = matrix_size(codes, INTSXP);
    const auto total = level_total(levels, shape.columns);
    const int* coded = INTEGER(codes);
    const int* counts = INTEGER(levels);
    SEXP output = PROTECT(Rf_allocVector(VECSXP, 8));
    for (int item = 0; item < 4; ++item) {
        SET_VECTOR_ELT(output, item, Rf_allocVector(INTSXP, total));
    }
    int rank = 0;
    int components = 0;
    int status = 0;
    bool wide = false;
    {
        const auto basis = powermean::exact_fe_basis(coded, shape.rows, counts, shape.columns);
        status = static_cast<int>(basis.status);
        if (basis.status == powermean::BasisStatus::ok) {
            std::copy(basis.component_by_level.begin(), basis.component_by_level.end(),
                      INTEGER(VECTOR_ELT(output, 0)));
            std::copy(basis.column_by_level.begin(), basis.column_by_level.end(),
                      INTEGER(VECTOR_ELT(output, 1)));
            std::copy(basis.position_by_level.begin(), basis.position_by_level.end(),
                      INTEGER(VECTOR_ELT(output, 2)));
            auto selected = INTEGER(VECTOR_ELT(output, 3));
            std::fill_n(selected, total, NA_INTEGER);
            std::copy(basis.selected.begin(), basis.selected.end(), selected);
            rank = static_cast<int>(basis.selected.size());
            components = basis.n_components;
            wide = basis.multiprecision;
        }
    }
    SET_VECTOR_ELT(output, 4, Rf_ScalarInteger(rank));
    SET_VECTOR_ELT(output, 5, Rf_ScalarInteger(components));
    SET_VECTOR_ELT(output, 6, Rf_ScalarLogical(wide));
    SET_VECTOR_ELT(output, 7, Rf_ScalarInteger(status));
    UNPROTECT(1);
    return output;
}

SEXP gravity_proof_call(SEXP codes, SEXP roles, SEXP maps, SEXP levels, SEXP y, SEXP x) {
    const auto shape = matrix_size(codes, INTSXP);
    if (shape.columns != 3 || TYPEOF(roles) != INTSXP || Rf_xlength(roles) != 3 ||
        TYPEOF(levels) != INTSXP || Rf_xlength(levels) != 4 ||
        TYPEOF(maps) != VECSXP || Rf_xlength(maps) != 2 ||
        TYPEOF(y) != REALSXP || TYPEOF(x) != REALSXP ||
        Rf_xlength(y) != static_cast<R_xlen_t>(shape.rows) ||
        Rf_xlength(x) != static_cast<R_xlen_t>(shape.rows)) {
        Rf_error("Invalid gravity support proof inputs.");
    }
    const int* order = INTEGER(roles);
    if (order[0] < 0 || order[0] >= 3 || order[1] < 0 || order[1] >= 3 ||
        order[2] < 0 || order[2] >= 3 || order[0] == order[1] ||
        order[0] == order[2] || order[1] == order[2]) {
        Rf_error("Gravity support proof requires three distinct fixed-effect roles.");
    }
    const auto n = shape.rows;
    const int* ptr = INTEGER(codes);
    const int* exporters = ptr + n * order[1];
    const int* importers = ptr + n * order[2];
    SEXP emap = VECTOR_ELT(maps, 0);
    SEXP imap = VECTOR_ELT(maps, 1);
    const auto edim = matrix_size(emap, INTSXP);
    const auto idim = matrix_size(imap, INTSXP);
    if (edim.columns != 2 || idim.columns != 2) Rf_error("Invalid port and year maps.");
    const int* e = INTEGER(emap);
    const int* i = INTEGER(imap);
    auto* ex = reinterpret_cast<int*>(R_alloc(n, sizeof(int)));
    auto* im = reinterpret_cast<int*>(R_alloc(n, sizeof(int)));
    auto* year = reinterpret_cast<int*>(R_alloc(n, sizeof(int)));
    for (std::size_t row = 0; row < n; ++row) {
        const auto left = exporters[row];
        const auto right = importers[row];
        if (left < 0 || right < 0 || static_cast<std::size_t>(left) >= edim.rows ||
            static_cast<std::size_t>(right) >= idim.rows ||
            e[edim.rows + left] != i[idim.rows + right]) {
            return Rf_ScalarLogical(FALSE);
        }
        ex[row] = e[left];
        im[row] = i[right];
        year[row] = e[edim.rows + left];
    }
    const int* count = INTEGER(levels);
    const powermean::GravityInput input{
        {ptr + n * order[0], n}, {ex, n}, {im, n}, {year, n},
        {REAL(y), n}, {REAL(x), n}, count[0], count[1], count[2], count[3]
    };
    const auto outcome = powermean::gravity_positive_proof(input);
    return Rf_ScalarLogical(outcome.status == powermean::GravityStatus::proved);
}
