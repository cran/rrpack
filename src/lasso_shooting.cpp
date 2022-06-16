#include <RcppArmadillo.h>

#include <rrpack.h>

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
        rrpack::run_active_cycles(xtx, xty, lambda, beta, is_active,
                                  epsilon, max_iter);
        is_active0 = is_active;
        is_active = all_active;
        rrpack::one_active_cycle(xtx, xty, lambda, beta, is_active);
        if (arma::accu(arma::abs(is_active0 - is_active)) == 0) {
            break;
        }
    }
    return rrpack::arma2rvec(beta);
}
