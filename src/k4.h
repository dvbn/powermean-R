#ifndef POWERMEAN_K4_H
#define POWERMEAN_K4_H

#include <string>
#include <vector>

namespace powermean {

struct DenseIolsOptions {
    int anderson_memory = 5;
    double eps = 1e-12;
    int stopcrit = 1000000;
    int max_phase1_iterations = 100000;
    int max_phase2_evaluations = 5000000;
    std::vector<double> pilot_deltas;
    double rho_log_min = -7;
    double rho_log_max = 14;
    void (*interrupt)() = nullptr;
};

struct DenseIolsResult {
    bool converged;
    int iterations;
    std::string phase;
    int restarts;
    int phase1_attempts;
    int phase1_accepted;
    int phase1_rejected;
    int phase1_evaluations;
    int phase2_evaluations;
    int response_unit_exponent;
    double intercept;
    double residual;
    double score_residual;
    double intercept_score_residual;
    double final_rho;
    double final_criterion;
    std::string failure_reason;
    std::string score_validation_status;
    std::vector<double> slopes;
    std::vector<double> index;
};

// Column-major slopes without an intercept, in the caller's normalized units.
DenseIolsResult dense_iols(const double* y, const double* slopes, int n, int p,
                          double tau, const DenseIolsOptions& options = {});

} // namespace powermean

#endif
