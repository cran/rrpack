#include <RcppArmadillo.h>

namespace lasso
{
    // convert arma vector type to Rcpp vector type
    template <typename T>
    inline Rcpp::NumericVector arma2rvec(const T& x) {
        return Rcpp::NumericVector(x.begin(), x.end());
    }
    // sign function
    double sign(const double x)
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
    double soft_threshold(const double z,
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
    void one_active_cycle(const arma::mat& xtx,
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
            beta(j) = lasso::soft_threshold(
                z, lambda, is_active, j) / xtx(j, j);
        }
    }
    void run_active_cycles(const arma::mat& xtx,
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

// lasso shooting for a fix tuning parameter
// [[Rcpp::export]]
Rcpp::NumericVector lasso_shooting(const arma::mat& xtx,
                                   const arma::vec& xty,
                                   const double lambda,
                                   const double epsilon = 1e-6,
                                   const unsigned int max_iter = 1e5)
{
    const size_t p { xtx.n_rows };
    arma::vec beta { arma::zeros(p) };
    const arma::uvec all_active { arma::ones<arma::uvec>(p) };
    arma::uvec is_active { all_active }, is_active0 { is_active };
    for (size_t i {0}; i < max_iter; ++i) {
        lasso::run_active_cycles(xtx, xty, lambda, beta, is_active,
                                 epsilon, max_iter);
        is_active0 = is_active;
        is_active = all_active;
        lasso::one_active_cycle(xtx, xty, lambda, beta, is_active);
        if (arma::accu(arma::abs(is_active0 - is_active)) == 0) {
            break;
        }
    }
    return lasso::arma2rvec(beta);
}
