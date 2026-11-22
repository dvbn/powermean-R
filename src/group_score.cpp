#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>

#include "native_calls.h"
#include "groups.h"

#include <algorithm>
#include <cmath>
#include <cstddef>

SEXP group_weighted_sums_call(SEXP codes, SEXP levels, SEXP block, SEXP design, SEXP score) {
    SEXP code_dimensions = Rf_getAttrib(codes, R_DimSymbol);
    SEXP design_dimensions = Rf_getAttrib(design, R_DimSymbol);
    if (TYPEOF(codes) != INTSXP || TYPEOF(levels) != INTSXP || TYPEOF(block) != INTSXP ||
        TYPEOF(design) != REALSXP || TYPEOF(score) != REALSXP ||
        TYPEOF(code_dimensions) != INTSXP || TYPEOF(design_dimensions) != INTSXP ||
        Rf_xlength(code_dimensions) != 2 || Rf_xlength(design_dimensions) != 2 ||
        Rf_xlength(block) != 1) {
        Rf_error("Grouped score inputs have invalid storage.");
    }
    const int rows = INTEGER(code_dimensions)[0];
    const int blocks = INTEGER(code_dimensions)[1];
    const int selected = INTEGER(block)[0];
    const int columns = INTEGER(design_dimensions)[1];
    if (rows < 1 || columns < 0 || selected < 1 || selected > blocks || Rf_xlength(levels) != blocks ||
        INTEGER(design_dimensions)[0] != rows || Rf_xlength(score) != rows) {
        Rf_error("Grouped score inputs do not align.");
    }
    const int capacity = INTEGER(levels)[selected - 1];
    if (capacity < 1 || capacity > rows) Rf_error("Invalid group level count.");
    SEXP output = PROTECT(Rf_allocMatrix(REALSXP, capacity, columns));
    const int* group = INTEGER(codes) + static_cast<std::size_t>(rows) * (selected - 1);
    if (powermean::group_weighted_sums(group, rows, capacity, REAL(design), REAL(score),
                                      columns, REAL(output)) != powermean::GroupStatus::ok) {
        Rf_error("Invalid group code.");
    }
    UNPROTECT(1);
    return output;
}

SEXP group_score_call(SEXP codes, SEXP levels, SEXP score, SEXP scale) {
    const SEXP dimensions = Rf_getAttrib(codes, R_DimSymbol);
    if (TYPEOF(codes) != INTSXP || TYPEOF(dimensions) != INTSXP ||
        Rf_xlength(dimensions) != 2 || TYPEOF(levels) != INTSXP ||
        TYPEOF(score) != REALSXP || TYPEOF(scale) != REALSXP) {
        Rf_error("Fixed-effect score inputs have invalid storage.");
    }
    const int rows = INTEGER(dimensions)[0];
    const int blocks = INTEGER(dimensions)[1];
    if (rows < 1 || blocks < 1 || Rf_xlength(levels) != blocks ||
        Rf_xlength(score) != rows || Rf_xlength(scale) != rows) {
        Rf_error("Fixed-effect score inputs do not align.");
    }
    const int* block_levels = INTEGER(levels);
    const int capacity = *std::max_element(block_levels, block_levels + blocks);
    if (capacity < 1 || capacity > rows) Rf_error("Invalid fixed-effect level count.");

    auto* numerator = reinterpret_cast<double*>(R_alloc(capacity, sizeof(double)));
    auto* denominator = reinterpret_cast<double*>(R_alloc(capacity, sizeof(double)));
    const double* raw_score = REAL(score);
    const double* bound = REAL(scale);
    double largest = 0;
    for (int block = 0; block < blocks; ++block) {
        const int count = block_levels[block];
        if (count < 1) Rf_error("Invalid fixed-effect level count.");
        std::fill_n(numerator, count, 0.0);
        std::fill_n(denominator, count, 0.0);
        const int* group = INTEGER(codes) + static_cast<std::size_t>(rows) * block;
        for (int row = 0; row < rows; ++row) {
            const int level = group[row];
            if (level < 0 || level >= count) Rf_error("Invalid fixed-effect code.");
            numerator[level] += raw_score[row];
            denominator[level] += bound[row];
        }
        for (int level = 0; level < count; ++level) {
            const double ratio = std::abs(numerator[level]) / denominator[level];
            if (!std::isfinite(ratio) || !std::isfinite(denominator[level]) ||
                denominator[level] <= 0) {
                return Rf_ScalarReal(R_PosInf);
            }
            largest = std::max(largest, ratio);
        }
    }
    return Rf_ScalarReal(largest);
}

SEXP group_center_call(SEXP codes, SEXP levels, SEXP values) {
    const SEXP dimensions = Rf_getAttrib(values, R_DimSymbol);
    if (TYPEOF(codes) != INTSXP || TYPEOF(levels) != INTSXP ||
        Rf_xlength(levels) != 1 || TYPEOF(values) != REALSXP ||
        TYPEOF(dimensions) != INTSXP || Rf_xlength(dimensions) != 2) {
        Rf_error("Group centering inputs have invalid storage.");
    }
    const int rows = INTEGER(dimensions)[0];
    const int columns = INTEGER(dimensions)[1];
    const int capacity = INTEGER(levels)[0];
    if (rows < 1 || columns < 1 || capacity < 1 || capacity > rows || Rf_xlength(codes) != rows) {
        Rf_error("Group centering inputs do not align.");
    }
    auto* count = reinterpret_cast<int*>(R_alloc(capacity, sizeof(int)));
    auto* group_total = reinterpret_cast<double*>(R_alloc(capacity, sizeof(double)));
    std::fill_n(count, capacity, 0);
    const int* group = INTEGER(codes);
    for (int row = 0; row < rows; ++row) {
        const int level = group[row];
        if (level < 0 || level >= capacity) Rf_error("Invalid group code.");
        ++count[level];
    }
    SEXP result = PROTECT(Rf_duplicate(values));
    double* buffer = REAL(result);
    for (int column = 0; column < columns; ++column) {
        std::fill_n(group_total, capacity, 0.0);
        double* data = buffer + static_cast<std::size_t>(rows) * column;
        for (int row = 0; row < rows; ++row) group_total[group[row]] += data[row];
        for (int level = 0; level < capacity; ++level) {
            if (count[level] > 0) group_total[level] /= count[level];
        }
        for (int row = 0; row < rows; ++row) data[row] -= group_total[group[row]];
    }
    UNPROTECT(1);
    return result;
}
