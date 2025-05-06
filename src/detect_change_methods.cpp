#include <RcppArmadillo.h>
#include <iostream>
#include <iomanip>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat C_cusum(arma::mat& ts,
                  const double& threshold,
                  const arma::uword& nb_samples,
                  const arma::uword& n_times) {

    // Aux variables
    arma::mat ci(ts.n_rows, n_times, arma::fill::zeros);
    arma::mat change(ts.n_rows, n_times, arma::fill::zeros);

    double change, mean;


    // For each pixel
    for (arma::uword i = 0; i < ts.n_rows; i++) {
        // Probability to be a Forest
        arma::colvec p_for(n_times, arma::fill::zeros);
        // Probability to be a Non-Forest
        arma::colvec p_nfor(n_times, arma::fill::zeros);
        // Probability to be a Non-Forest in the past
        arma::colvec p_nfor_past(n_times, arma::fill::zeros);

        // Aux variables
        arma::uword col_idx = 0;
        bool update_res = false;

        // For each band
        for (arma::uword c = 0; c < ts.n_cols; c = c + n_times) {
            // Deseasonlize time series
            if (quantile_values.size() > 1) {
                ts.submat(i, c, i, c + n_times - 1) = C_bayts_calc_sub(
                    ts.submat(i, c, i, c + n_times - 1),
                    quantile_values.submat(0, c, 0, c + n_times - 1)
                );
            }
            // Estimate a normal distribution based on Forest stats
            p_for = C_dnorm(
                ts.submat(i, c, i, c + n_times - 1).t(),
                mean(0, col_idx),
                sd(0, col_idx)
            );
            // Estimate a normal distribution based on Non-Forest stats
            p_nfor = C_dnorm(
                ts.submat(i, c, i, c + n_times - 1).t(),
                mean(1, col_idx),
                sd(1, col_idx)
            );
            // Clean values lower than 0.00001
            p_nfor.elem(arma::find(p_nfor < 0.00001)).zeros();
            // Estimate a conditional prob for each positive distribution value
            p_nfor.elem(arma::find(p_nfor > 0)) = C_bayts_calc_pcond(
                p_nfor.elem(arma::find(p_nfor > 0)),
                p_for.elem(arma::find(p_nfor > 0))
            );
            // Fix the range of prob values
            p_nfor.elem(arma::find(p_nfor < bwf(0))).fill(bwf(0));
            p_nfor.elem(arma::find(p_nfor > bwf(1))).fill(bwf(1));

            // Update NF probabilities with a Bayesian approach
            if (update_res) {
                arma::uvec p1 = arma::find_finite(p_nfor);
                arma::uvec p2 = arma::find_finite(p_nfor_past);

                arma::uvec non_na_idxs = arma::intersect(p1, p2);

                p_nfor(non_na_idxs) = C_bayts_calc_pbayes(
                    p_nfor(non_na_idxs), p_nfor_past(non_na_idxs)
                );

                arma::uvec p1_na = arma::find_nonfinite(p_nfor);
                p_nfor(p1_na) = p_nfor_past(p1_na);
            }
            // Update Non-Forest probs
            p_nfor_past = p_nfor;
            update_res = true;
            col_idx++;
        }
        // Get the probs for NF values
        p_res.submat(i, 1, i, n_times) = p_nfor.t();
    }
    // Return the probs results
    return p_res;
}
