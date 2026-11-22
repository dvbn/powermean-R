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

struct MatrixView {
    const double* value;
    std::size_t rows;
    std::size_t columns;

    double get(std::size_t row, std::size_t column) const {
        return value[row + rows * column];
    }
};

MatrixView matrix_view(SEXP value) {
    SEXP dimensions = Rf_getAttrib(value, R_DimSymbol);
    if (TYPEOF(value) != REALSXP || TYPEOF(dimensions) != INTSXP ||
        Rf_xlength(dimensions) != 2 || INTEGER(dimensions)[0] < 0 || INTEGER(dimensions)[1] < 0) {
        Rf_error("WZ arrays require double matrices.");
    }
    return {REAL(value), static_cast<std::size_t>(INTEGER(dimensions)[0]),
            static_cast<std::size_t>(INTEGER(dimensions)[1])};
}

template <class T>
T* workspace(std::size_t size) {
    return reinterpret_cast<T*>(R_alloc(size, sizeof(T)));
}

double nuisance_map(const MatrixView& inverse, int a, int b) {
    return a && b ? inverse.get(a - 1, b - 1) : 0;
}

struct PairSolveResult {
    const char* failure = nullptr;
    double rcond = NA_REAL;
    int pivot = 0;
};

PairSolveResult pair_solve(std::size_t dates, double* matrix, double* score, double* scale,
                           LapackInt* pivot, double* work, LapackInt* integer_work, double cutoff) {
    std::fill_n(scale, dates, 0.0);
    for (std::size_t j = 0; j < dates; ++j) {
        for (std::size_t i = 0; i < dates; ++i) {
            const double magnitude = std::abs(matrix[i + dates * j]);
            if (!std::isfinite(magnitude)) return {"nonfinite"};
            scale[i] = std::max(scale[i], magnitude);
            scale[j] = std::max(scale[j], magnitude);
        }
    }
    for (std::size_t i = 0; i < dates; ++i) {
        if (!(scale[i] > 0)) return {"zero_scale"};
        scale[i] = std::sqrt(scale[i]);
        score[i] /= scale[i];
        if (!std::isfinite(score[i])) return {"nonfinite"};
    }
    double norm = 0;
    double distance_from_identity = 0;
    for (std::size_t j = 0; j < dates; ++j) {
        double column = 0;
        double deviation = 0;
        for (std::size_t i = 0; i < dates; ++i) {
            double& cell = matrix[i + dates * j];
            cell = (cell / scale[i]) / scale[j];
            const double magnitude = std::abs(cell);
            column += magnitude;
            deviation += i == j ? std::abs(cell - 1) : magnitude;
        }
        norm = std::max(norm, column);
        distance_from_identity = std::max(distance_from_identity, deviation);
    }
    const LapackInt n = static_cast<LapackInt>(dates);
    const LapackInt one = 1;
    LapackInt info = 0;
    F77_CALL(dgetrf)(&n, &n, matrix, &n, pivot, &info);
    if (info > 0) return {"zero_pivot", NA_REAL, static_cast<int>(info)};
    if (info < 0) return {"factorization"};
    // With ||M - I||_1 <= 1/2 the Neumann bound gives rcond(M) >= 1/3.
    // For every other pair keep the release condition estimator and cutoff.
    if (distance_from_identity > 0.5 || cutoff >= 1.0 / 3.0) {
        double reciprocal = 0;
        F77_CALL(dgecon)("1", &n, matrix, &n, &norm, &reciprocal, work, integer_work, &info FCONE);
        if (info != 0 || !std::isfinite(reciprocal)) return {"condition_estimate"};
        if (reciprocal <= cutoff) return {"ill_conditioned", reciprocal};
    }
    F77_CALL(dgetrs)("N", &n, &one, matrix, &n, pivot, score, &n, &info FCONE);
    if (info != 0) return {"factorization"};
    for (std::size_t i = 0; i < dates; ++i) {
        score[i] /= scale[i];
        if (!std::isfinite(score[i])) return {"nonfinite"};
    }
    return {};
}

SEXP pair_failure(const PairSolveResult& result) {
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 3));
    SEXP names = PROTECT(Rf_allocVector(STRSXP, 3));
    SET_VECTOR_ELT(out, 0, Rf_mkString(result.failure));
    SET_VECTOR_ELT(out, 1, Rf_ScalarReal(result.rcond));
    SET_VECTOR_ELT(out, 2, Rf_ScalarInteger(result.pivot));
    SET_STRING_ELT(names, 0, Rf_mkChar("failure"));
    SET_STRING_ELT(names, 1, Rf_mkChar("rcond"));
    SET_STRING_ELT(names, 2, Rf_mkChar("pivot"));
    Rf_setAttrib(out, R_NamesSymbol, names);
    UNPROTECT(2);
    return out;
}

SEXP wz_pair_meats(SEXP mean_value, SEXP probability_value, SEXP score_value,
                   SEXP slope_value, SEXP mapped_value, SEXP nuisance_value,
                   SEXP from_codes, SEXP to_codes, double cutoff) {
    const MatrixView mean = matrix_view(mean_value);
    const MatrixView probability = matrix_view(probability_value);
    const MatrixView score = matrix_view(score_value);
    const MatrixView slope = matrix_view(slope_value);
    const MatrixView mapped = matrix_view(mapped_value);
    const MatrixView inverse = matrix_view(nuisance_value);
    const std::size_t dates = mean.columns;
    const std::size_t groups = mean.rows;
    const std::size_t p = slope.columns;
    if (groups < 2 || dates < 2 || !p || slope.rows != groups * dates ||
        mapped.rows != slope.rows || mapped.columns != p ||
        probability.rows != groups || probability.columns != dates ||
        score.rows != groups || score.columns != dates || inverse.rows != inverse.columns ||
        TYPEOF(from_codes) != INTSXP || TYPEOF(to_codes) != INTSXP ||
        static_cast<std::size_t>(Rf_xlength(from_codes)) != slope.rows ||
        static_cast<std::size_t>(Rf_xlength(to_codes)) != slope.rows) {
        Rf_error("WZ pair arrays do not align.");
    }
    const int* from = INTEGER(from_codes);
    const int* to = INTEGER(to_codes);
    for (std::size_t row = 0; row < slope.rows; ++row) {
        if (from[row] < 0 || to[row] < 0 || static_cast<std::size_t>(from[row]) > inverse.rows ||
            static_cast<std::size_t>(to[row]) > inverse.rows) {
            Rf_error("WZ nuisance coordinates are outside the identified basis.");
        }
    }
    auto* left = workspace<double>(dates * dates);
    auto* solved = workspace<double>(dates);
    auto* scale = workspace<double>(dates);
    auto* pivot = workspace<LapackInt>(dates);
    auto* work = workspace<double>(4 * dates);
    auto* integer_work = workspace<LapackInt>(dates);
    auto* ordinary = workspace<double>(p);
    auto* adjusted = workspace<double>(p);

    SEXP raw = PROTECT(Rf_allocMatrix(REALSXP, p, p));
    SEXP corrected = PROTECT(Rf_allocMatrix(REALSXP, p, p));
    std::fill_n(REAL(raw), p * p, 0.0);
    std::fill_n(REAL(corrected), p * p, 0.0);
    for (std::size_t group = 0; group < groups; ++group) {
        if ((group & 1023U) == 0) R_CheckUserInterrupt();
        const std::size_t begin = group * dates;
        for (std::size_t j = 0; j < dates; ++j) {
            double average = 0;
            const std::size_t col = begin + j;
            for (std::size_t i = 0; i < dates; ++i) {
                const std::size_t row = begin + i;
                double value = nuisance_map(inverse, from[row], from[col]) +
                    nuisance_map(inverse, to[row], from[col]) +
                    nuisance_map(inverse, from[row], to[col]) +
                    nuisance_map(inverse, to[row], to[col]);
                for (std::size_t k = 0; k < p; ++k) value += mapped.get(row, k) * slope.get(col, k);
                left[i + dates * j] = value;
                average += probability.get(group, i) * value;
            }
            for (std::size_t i = 0; i < dates; ++i) {
                double& value = left[i + dates * j];
                value = static_cast<double>(i == j) - mean.get(group, i) * (value - average);
            }
            solved[j] = score.get(group, j);
        }
        const auto result = pair_solve(dates, left, solved, scale, pivot, work, integer_work, cutoff);
        if (result.failure != nullptr) {
            UNPROTECT(2);
            return pair_failure(result);
        }
        for (std::size_t k = 0; k < p; ++k) {
            ordinary[k] = adjusted[k] = 0;
            for (std::size_t i = 0; i < dates; ++i) {
                const double x = slope.get(begin + i, k);
                ordinary[k] += x * score.get(group, i);
                adjusted[k] += x * solved[i];
            }
        }
        for (std::size_t j = 0; j < p; ++j) {
            for (std::size_t i = 0; i < p; ++i) {
                REAL(raw)[i + p * j] += ordinary[i] * ordinary[j];
                REAL(corrected)[i + p * j] += adjusted[i] * ordinary[j];
            }
        }
    }
    SEXP out = PROTECT(Rf_allocVector(VECSXP, 2));
    SET_VECTOR_ELT(out, 0, raw);
    SET_VECTOR_ELT(out, 1, corrected);
    UNPROTECT(3);
    return out;
}

}  // namespace

SEXP wz_pair_meats_call(SEXP mean, SEXP probability, SEXP score, SEXP slope,
                        SEXP mapped, SEXP nuisance, SEXP origin, SEXP destination) {
    return wz_pair_meats(mean, probability, score, slope, mapped, nuisance,
                         origin, destination, 1e-13);
}

SEXP wz_pair_meats_controlled_call(SEXP mean, SEXP probability, SEXP score,
                                   SEXP slope, SEXP mapped, SEXP nuisance,
                                   SEXP origin, SEXP destination, SEXP threshold) {
    const double cutoff = Rf_asReal(threshold);
    if (!std::isfinite(cutoff) || cutoff <= 0 || cutoff >= 1) {
        Rf_error("Invalid WZ reciprocal condition threshold.");
    }
    return wz_pair_meats(mean, probability, score, slope, mapped, nuisance,
                         origin, destination, cutoff);
}
