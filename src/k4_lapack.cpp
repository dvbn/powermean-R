// Adapted from Powermean 0.1.0, BSD-3-Clause. See inst/LICENSE.K4.
// Frozen source: 7b22dcf01c6fea74176823e24a4b18df559cd386
// src/powermean/_native/iols/bind.cpp
// SHA-256 7d12b4308bbff9f4ee6f3dc3f72f77c06874462c7e762d6bc7c8ad707060e4e5.

#include "k4_core.h"
#include <R_ext/Lapack.h>

namespace powermean::k4 {

// R versions before La_INT used ordinary LAPACK int.
#ifdef La_INT
using K4LapackInt = La_INT;
#else
using K4LapackInt = int;
#endif

double stable_l2_norm(
    const Eigen::Ref<const Eigen::VectorXd>& values) {
    double norm = 0.0;
    for (Eigen::Index idx = 0; idx < values.size(); ++idx) {
        norm = std::hypot(norm, std::abs(values(idx)));
    }
    return norm;
}

Eigen::MatrixXd min_norm_lstsq(
    const Eigen::Ref<const Eigen::MatrixXd>& a,
    const Eigen::Ref<const Eigen::MatrixXd>& b) {
    if (a.rows() == 0 || a.cols() == 0) {
        throw std::invalid_argument("least-squares matrix must be non-empty");
    }
    if (a.rows() != b.rows()) {
        throw std::invalid_argument("least-squares row dimensions do not match");
    }
    if (a.rows() >
            static_cast<Eigen::Index>(std::numeric_limits<K4LapackInt>::max()) ||
        a.cols() >
            static_cast<Eigen::Index>(std::numeric_limits<K4LapackInt>::max()) ||
        b.cols() >
            static_cast<Eigen::Index>(std::numeric_limits<K4LapackInt>::max())) {
        throw std::overflow_error("least-squares dimensions exceed LAPACK int range");
    }
    K4LapackInt m = static_cast<K4LapackInt>(a.rows());
    K4LapackInt n = static_cast<K4LapackInt>(a.cols());
    K4LapackInt nrhs = static_cast<K4LapackInt>(b.cols());
    K4LapackInt lda = std::max<K4LapackInt>(1, m);
    K4LapackInt ldb = std::max<K4LapackInt>(1, std::max(m, n));
    K4LapackInt minmn = std::min(m, n);

    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> a_lapack =
        a;
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor> b_lapack =
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>::Zero(
            ldb,
            nrhs);
    b_lapack.topRows(m) = b;
    Eigen::VectorXd singular(minmn);

    double rcond =
        std::numeric_limits<double>::epsilon() *
        static_cast<double>(std::max(a.rows(), a.cols()));
    K4LapackInt rank = 0;
    K4LapackInt info = 0;

    K4LapackInt smlsiz = 25;
    K4LapackInt nlvl = 0;
    if (minmn > smlsiz + 1) {
        nlvl = static_cast<K4LapackInt>(
            std::log(static_cast<double>(minmn) / static_cast<double>(smlsiz + 1)) /
            std::log(2.0)) +
               1;
    }
    K4LapackInt liwork =
        std::max<K4LapackInt>(
            1,
            3 * minmn * std::max<K4LapackInt>(0, nlvl) + 11 * minmn);
    std::vector<K4LapackInt> iwork(static_cast<std::size_t>(liwork), 0);

    double work_query = 0.0;
    K4LapackInt lwork = -1;
    F77_CALL(dgelsd)(
        &m,
        &n,
        &nrhs,
        a_lapack.data(),
        &lda,
        b_lapack.data(),
        &ldb,
        singular.data(),
        &rcond,
        &rank,
        &work_query,
        &lwork,
        iwork.data(),
        &info);
    if (info != 0 || !std::isfinite(work_query) || work_query <= 0.0 ||
        work_query > static_cast<double>(std::numeric_limits<K4LapackInt>::max())) {
        throw std::runtime_error("LAPACK DGELSD workspace query failed");
    }
    lwork = static_cast<K4LapackInt>(std::ceil(work_query));
    std::vector<double> work(
        static_cast<std::size_t>(std::max<K4LapackInt>(1, lwork)),
        0.0);
    F77_CALL(dgelsd)(
        &m,
        &n,
        &nrhs,
        a_lapack.data(),
        &lda,
        b_lapack.data(),
        &ldb,
        singular.data(),
        &rcond,
        &rank,
        work.data(),
        &lwork,
        iwork.data(),
        &info);
    if (info != 0) {
        throw std::runtime_error("LAPACK DGELSD least-squares solve failed");
    }
    return Eigen::MatrixXd(b_lapack.topRows(n));
}

Eigen::MatrixXd min_norm_lstsq_anderson(
    const Eigen::Ref<const Eigen::MatrixXd>& a,
    const Eigen::Ref<const Eigen::MatrixXd>& b) {
    return min_norm_lstsq(a, b);
}

} // namespace powermean::k4
