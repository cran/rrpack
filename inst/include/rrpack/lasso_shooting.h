#ifndef RRPACK_LASSO_SHOOTING_H
#define RRPACK_LASSO_SHOOTING_H

#include <RcppArmadillo.h>

namespace rrpack
{
    // convert arma vector type to Rcpp vector type
    template <typename T>
    inline Rcpp::NumericVector arma2rvec(const T& x) {
        return Rcpp::NumericVector(x.begin(), x.end());
    }
    // sign function
    inline double sign(const double x)
    {
        if (x < 0) {
            return - 1.0;
        }
        if (x > 0) {
            return 1.0;
        }
        return 0.0;
    }
    // soft-thresholding
    inline double soft_threshold(const double z,
                                 const double lambda,
                                 arma::uvec& is_active,
                                 const size_t idx)
    {
        double tmp { std::abs(z) - lambda };
        if (tmp < 0) {
            is_active(idx) = 0;
            return 0;
        }
        return tmp * sign(z);
    }
    // one active cycle
    inline void one_active_cycle(const arma::mat& xtx,
                                 const arma::vec& xty,
                                 const double lambda,
                                 arma::vec& beta,
                                 arma::uvec& is_active
        )
    {
        for (size_t j { 0 }; j < xtx.n_rows; ++j) {
            if (is_active(j) == 0) {
                continue;
            }
            double tmp { 0.0 };
            for (size_t k { 0 }; k < xtx.n_rows; ++k) {
                if (k == j) {
                    continue;
                }
                tmp += xtx(k, j) * beta(k);
            }
            double z { xty(j) - tmp };
            beta(j) = rrpack::soft_threshold(
                z, lambda, is_active, j) / xtx(j, j);
        }
    }
    inline void run_active_cycles(const arma::mat& xtx,
                                  const arma::mat& xty,
                                  const double lambda,
                                  arma::vec& beta,
                                  arma::uvec& is_active,
                                  const double epsilon,
                                  const unsigned int max_iter)
    {
        arma::vec beta0 { beta };
        for (size_t i { 0 }; i < max_iter; ++i) {
            one_active_cycle(xtx, xty, lambda, beta, is_active);
            if (arma::max(arma::abs(beta - beta0)) < epsilon) {
                break;
            }
            beta0 = beta;
        }
    }
}

#endif /* RRPACK_LASSO_SHOOTING_H */
