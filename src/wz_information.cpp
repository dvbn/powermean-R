#define R_NO_REMAP
#include <R.h>
#include <R_ext/Lapack.h>
#include <R_ext/Utils.h>
#include <Rinternals.h>

#include "native_calls.h"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace {

#ifdef La_INT
using LapackInt = La_INT;
#else
using LapackInt = int;
#endif

// For each date pair, add the positive outer product of the difference
// between the two sparse nuisance rows. This is D' (diag(mu)-mu*pi') D.
void add_pair_information(const double* mean, const double* probability, int pair,
                          int pairs, int dates, const int* left, const int* right,
                          int dimension, double* information) {
    const std::size_t base = static_cast<std::size_t>(pair) * dates;
    const int sign[] = {1, 1, -1, -1};
    for (int first = 0; first < dates; ++first) {
        const std::size_t a = pair + static_cast<std::size_t>(first) * pairs;
        for (int second = first + 1; second < dates; ++second) {
            const std::size_t b = pair + static_cast<std::size_t>(second) * pairs;
            const double weight = mean[a] < mean[b] ? mean[a] * probability[b] : mean[b] * probability[a];
            const int columns[] = {left[base + first] - 1, right[base + first] - 1,
                                   left[base + second] - 1, right[base + second] - 1};
            for (int j = 0; j < 4; ++j) {
                if (columns[j] < 0) continue;
                double* column = information + static_cast<std::size_t>(columns[j]) * dimension;
                for (int i = 0; i < 4; ++i) {
                    if (columns[i] >= 0) column[columns[i]] += sign[i] * sign[j] * weight;
                }
            }
        }
    }
}

bool equilibrate(double* matrix, int size, double* scale) {
    std::fill_n(scale, size, 0.0);
    for (int col = 0; col < size; ++col) {
        const double* column = matrix + static_cast<std::size_t>(col) * size;
        for (int row = 0; row < size; ++row) {
            const double value = std::abs(column[row]);
            if (!std::isfinite(value)) return false;
            scale[row] = std::max(scale[row], value);
        }
    }
    for (int row = 0; row < size; ++row) {
        if (scale[row] <= 0) return false;
        scale[row] = std::sqrt(scale[row]);
    }
    for (int col = 0; col < size; ++col) {
        double* column = matrix + static_cast<std::size_t>(col) * size;
        for (int row = 0; row < size; ++row) column[row] = column[row] / scale[row] / scale[col];
    }
    return true;
}

bool unscale_nuisance(double* matrix, const double* scale, int size) {
    for (int col = 0; col < size; ++col) {
        for (int row = 0; row <= col; ++row) {
            const auto upper = row + static_cast<std::size_t>(col) * size;
            const auto lower = col + static_cast<std::size_t>(row) * size;
            const double value = matrix[upper];
            matrix[lower] = (value / scale[col]) / scale[row];
            matrix[upper] = (value / scale[row]) / scale[col];
            if (!std::isfinite(matrix[upper]) || !std::isfinite(matrix[lower])) {
                return false;
            }
        }
    }
    return true;
}

SEXP invert_nuisance(SEXP information, const double* scale, int size) {
    auto* matrix = REAL(information);
    auto* work = reinterpret_cast<double*>(R_alloc(size, sizeof(double)));
    const LapackInt n = static_cast<LapackInt>(size);
    const double norm = F77_CALL(dlange)("1", &n, &n, matrix, &n, work FCONE);
    LapackInt info = 0;
    F77_CALL(dpotrf)("U", &n, matrix, &n, &info FCONE);
    if (info != 0) return R_NilValue;
    F77_CALL(dpotri)("U", &n, matrix, &n, &info FCONE);
    if (info != 0) return R_NilValue;
    const double inverse_norm = F77_CALL(dlansy)("1", "U", &n, matrix, &n, work FCONE FCONE);
    const double reciprocal = 1 / (norm * inverse_norm);
    // Request the original R condition estimate throughout the boundary band.
    if (reciprocal >= 1e-14 && reciprocal <= 1e-12) return Rf_ScalarLogical(TRUE);
    if (!std::isfinite(reciprocal) || reciprocal <= 1e-13) return R_NilValue;
    return unscale_nuisance(matrix, scale, size) ? information : Rf_ScalarLogical(FALSE);
}

// Verify that distinct exporter countries have zero cross terms in the
// upper triangle used by both Cholesky paths.
bool has_block_diagonal(const double* information, int size, const int* widths, int blocks) {
    int start = 0;
    for (int block = 0; block < blocks; ++block) {
        const int end = start + widths[block];
        for (int col = start; col < end; ++col) {
            const double* column = information + static_cast<std::size_t>(col) * size;
            for (int row = 0; row < start; ++row) {
                if (column[row] != 0) return false;
            }
        }
        start = end;
    }
    return true;
}

bool form_block_inverse(double* inverse, int size, const int* widths, int blocks, int split) {
    const LapackInt n = size, left = split, right = size - split;
    double* cross = inverse + static_cast<std::size_t>(left) * n;
    double* schur = cross + left;
    const double one = 1, minus_one = -1;
    LapackInt info = 0;
    int start = 0;
    // For A = R_A' R_A, compute W = R_A^{-T} B and S = D - W'W.
    for (int block = 0; block < blocks; ++block) {
        const LapackInt width = widths[block];
        double* diagonal = inverse + static_cast<std::size_t>(start) * (n + 1);
        F77_CALL(dpotrf)("U", &width, diagonal, &n, &info FCONE);
        if (info != 0) return false;
        F77_CALL(dtrsm)("L", "U", "T", "N", &width, &right, &one,
                        diagonal, &n, cross + start, &n FCONE FCONE FCONE FCONE);
        start += width;
    }
    F77_CALL(dsyrk)("U", "T", &right, &left, &minus_one, cross, &n, &one,
                    schur, &n FCONE FCONE);
    F77_CALL(dpotrf)("U", &right, schur, &n, &info FCONE);
    if (info != 0) return false;
    start = 0;
    // C = A^{-1} B = R_A^{-1} W. Use the same factors for A^{-1}.
    for (int block = 0; block < blocks; ++block) {
        const LapackInt width = widths[block];
        double* diagonal = inverse + static_cast<std::size_t>(start) * (n + 1);
        F77_CALL(dtrsm)("L", "U", "N", "N", &width, &right, &one,
                        diagonal, &n, cross + start, &n FCONE FCONE FCONE FCONE);
        F77_CALL(dpotri)("U", &width, diagonal, &n, &info FCONE);
        if (info != 0) return false;
        start += width;
    }
    // With S = R_S' R_S and Q = C R_S^{-1}, the inverse blocks are
    // A^{-1} + QQ', -Q R_S^{-T}, and S^{-1}.
    F77_CALL(dtrsm)("R", "U", "N", "N", &left, &right, &one,
                    schur, &n, cross, &n FCONE FCONE FCONE FCONE);
    F77_CALL(dsyrk)("U", "N", &left, &right, &one, cross, &n, &one,
                    inverse, &n FCONE FCONE);
    F77_CALL(dtrsm)("R", "U", "T", "N", &left, &right, &minus_one,
                    schur, &n, cross, &n FCONE FCONE FCONE FCONE);
    F77_CALL(dpotri)("U", &right, schur, &n, &info FCONE);
    return info == 0;
}

bool check_block_inverse(const double* information, const double* inverse, int size, int split) {
    const LapackInt n = size, stride = 1;
    auto* work = reinterpret_cast<double*>(R_alloc(3 * size, sizeof(double)));
    const double norm = F77_CALL(dlansy)("1", "U", &n, information, &n, work FCONE FCONE);
    const double inverse_norm = F77_CALL(dlansy)("1", "U", &n, inverse, &n, work FCONE FCONE);
    const double reciprocal = 1 / (norm * inverse_norm);
    if (!std::isfinite(reciprocal) || reciprocal <= 1e-8) return false;
    double* probe = work;
    double* solution = work + size;
    double* residual = solution + size;
    const double one = 1, zero = 0, minus_one = -1;
    for (int check = 0; check < 4; ++check) {
        std::fill_n(probe, size, check == 0 ? 1.0 : 0.0);
        if (check == 1) {
            for (int row = 0; row < size; ++row) probe[row] = row % 2 ? -1.0 : 1.0;
        } else if (check == 2) {
            probe[split - 1] = 1.0;
        } else if (check == 3) {
            probe[size - 1] = 1.0;
        }
        F77_CALL(dsymv)("U", &n, &one, inverse, &n, probe, &stride, &zero,
                        solution, &stride FCONE);
        std::copy_n(probe, size, residual);
        F77_CALL(dsymv)("U", &n, &one, information, &n, solution, &stride, &minus_one,
                        residual, &stride FCONE);
        double error = 0, magnitude = 0;
        for (int row = 0; row < size; ++row) {
            if (!std::isfinite(residual[row]) || !std::isfinite(solution[row])) return false;
            error = std::max(error, std::abs(residual[row]));
            magnitude = std::max(magnitude, std::abs(solution[row]));
        }
        if (error > 1e-8 || error > 1e-10 * (norm * magnitude + 1)) return false;
    }
    return true;
}

}  // namespace

SEXP wz_nuisance_matrix_call(SEXP mean, SEXP probability, SEXP left, SEXP right, SEXP dimensions) {
    const SEXP shape = Rf_getAttrib(mean, R_DimSymbol);
    if (TYPEOF(mean) != REALSXP || TYPEOF(probability) != REALSXP || TYPEOF(shape) != INTSXP ||
        Rf_xlength(shape) != 2 || TYPEOF(left) != INTSXP || TYPEOF(right) != INTSXP ||
        TYPEOF(dimensions) != INTSXP || Rf_xlength(dimensions) < 1) {
        Rf_error("Invalid WZ nuisance information inputs.");
    }
    const int pairs = INTEGER(shape)[0];
    const int dates = INTEGER(shape)[1];
    const int size = INTEGER(dimensions)[0];
    // A scalar dimension requests the equilibrated matrix and its scales.
    const R_xlen_t dimension_args = Rf_xlength(dimensions);
    const bool need_inverse = dimension_args > 1 && INTEGER(dimensions)[1] == 1;
    if (dimension_args > 1 && !need_inverse) Rf_error("Invalid WZ nuisance request.");
    const R_xlen_t rows = Rf_xlength(mean);
    if (pairs < 1 || dates < 2 || size < 1 || Rf_xlength(probability) != rows ||
        Rf_xlength(left) != rows || Rf_xlength(right) != rows) {
        Rf_error("WZ nuisance information inputs do not align.");
    }
    const double* m = REAL(mean);
    const double* pi = REAL(probability);
    const int* l = INTEGER(left);
    const int* r = INTEGER(right);
    const int* widths = dimension_args > 2 ? INTEGER(dimensions) + 2 : nullptr;
    const int blocks = dimension_args > 2 && dimension_args - 2 < size ?
        static_cast<int>(dimension_args - 2) : 0;
    int split = 0;
    for (int block = 0; block < blocks; ++block) {
        if (widths[block] < 1 || widths[block] > dates || widths[block] >= size - split) {
            split = 0;
            break;
        }
        split += widths[block];
    }
    for (R_xlen_t row = 0; row < rows; ++row) {
        if (l[row] < 0 || l[row] > size || r[row] < 0 || r[row] > size) {
            Rf_error("Invalid WZ nuisance coordinate.");
        }
        if (!std::isfinite(m[row]) || !std::isfinite(pi[row]) || m[row] <= 0 || pi[row] < 0 || pi[row] > 1) {
            return R_NilValue;
        }
        if (split && (l[row] > split || (r[row] && r[row] <= split))) split = 0;
    }
    SEXP information = PROTECT(Rf_allocMatrix(REALSXP, size, size));
    SEXP scale = PROTECT(Rf_allocVector(REALSXP, size));
    std::fill_n(REAL(information), Rf_xlength(information), 0.0);
    for (int pair = 0; pair < pairs; ++pair) {
        if ((pair & 8191) == 0) R_CheckUserInterrupt();
        add_pair_information(m, pi, pair, pairs, dates, l, r, size, REAL(information));
    }
    if (!equilibrate(REAL(information), size, REAL(scale))) {
        UNPROTECT(2);
        return R_NilValue;
    }
    if (need_inverse) {
        if (split && has_block_diagonal(REAL(information), size, widths, blocks)) {
            SEXP candidate = PROTECT(Rf_allocMatrix(REALSXP, size, size));
            double* values = REAL(candidate);
            const double* original = REAL(information);
            for (int col = 0; col < size; ++col) {
                const auto offset = static_cast<std::size_t>(col) * size;
                std::copy_n(original + offset, col + 1, values + offset);
            }
            if (form_block_inverse(values, size, widths, blocks, split) &&
                check_block_inverse(original, values, size, split) &&
                unscale_nuisance(values, REAL(scale), size)) {
                UNPROTECT(3);
                return candidate;
            }
            UNPROTECT(1);
        }
        SEXP output = invert_nuisance(information, REAL(scale), size);
        UNPROTECT(2);
        return output;
    }
    SEXP output = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(output, 0, information);
    SET_VECTOR_ELT(output, 1, scale);
    UNPROTECT(3);
    return output;
}
