#include <RcppArmadillo.h>
#include <iostream>
#include <iomanip>
#include <cmath>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;

arma::vec C_calc_pbayes(const arma::vec& prior, const arma::vec& post) {
    return (prior % post) / ((prior % post) + ((1 - prior) % (1 - post)));
}

double C_calc_pbayes(const double& prior, const double& post) {
    double res = (prior * post) / ((prior * post) + ((1 - prior) * (1 - post)));
    return (std::floor(res * 1000000000000000.0) / 1000000000000000.0);
}


// [[Rcpp::export]]
arma::rowvec C_calc_probcusum(const arma::rowvec& ts,
                              const arma::uword& warmup_period,
                              const double& threshold) {
    arma::uword current_t = warmup_period;
    arma::rowvec res(ts.n_cols, arma::fill::zeros);
    double mean = arma::mean(ts.subvec(0, warmup_period - 1));
    double std = arma::stddev(ts.subvec(0, warmup_period - 1));
    double sum, stand_sum, p_obs, probability = 0;
    for (arma::uword c = warmup_period; c < ts.n_cols; c++) {
        sum += ts.at(c) - mean;
        stand_sum = sum  / (std * pow(c, 0.5));
        p_obs = arma::normcdf(abs(stand_sum));
        probability = 2 * (1 - p_obs);
        res(c) = 1 - probability;
        // Find a change
        if (probability < threshold) {
            c++;
            if ((c + warmup_period) >= ts.n_cols) {
                break;
            }
            mean = arma::mean(ts.subvec(c, c + warmup_period - 1));
            std = arma::stddev(ts.subvec(c, c + warmup_period - 1));
            sum = 0;
            c += warmup_period;
        }
    }
    return(res);
}

// [[Rcpp::export]]
arma::mat C_cusum(const arma::mat& ts,
                  const arma::uword& warmup_period,
                  const double& threshold,
                  const arma::uword& n_times,
                  const arma::uword& n_bands) {

    // Aux variables
    arma::mat res(ts.n_rows, n_times, arma::fill::zeros);

    // For each pixel
    for (arma::uword i = 0; i < ts.n_rows; i++) {
        // Aux variables
        arma::uword col_idx = 0;
        bool update_res = false;
        // Currently probability of a change
        arma::rowvec p_change(n_times - 1, arma::fill::zeros);
        // Past probability of a change
        arma::rowvec p_change_past(n_times - 1, arma::fill::zeros);
        // For each band
        for (arma::uword c = 0; c < ts.n_cols; c = c + n_times) {
            // ...
            p_change = C_calc_probcusum(
                ts.submat(i, c, i, c + n_times - 1), warmup_period,
                threshold
            );

            // Update NF probabilities with a Bayesian approach
            if (update_res) {
                arma::uvec p1 = arma::find_finite(p_change);
                arma::uvec p2 = arma::find_finite(p_change_past);

                arma::uvec non_na_idxs = arma::intersect(p1, p2);
                p_change(non_na_idxs) = C_calc_pbayes(
                    p_change(non_na_idxs), p_change_past(non_na_idxs)
                );

                arma::uvec p1_na = arma::find_nonfinite(p_change);
                p_change(p1_na) = p_change_past(p1_na);
            }
            // ...
            p_change_past = p_change;
            update_res = true;
        }
        res.submat(i, 0, i, n_times - 1) = p_change;
    }
    // Return the results
    return res;
}

// [[Rcpp::export]]
arma::mat C_bocd(const arma::mat& ts,
                 const arma::uword& n_times,
                 const arma::uword& n_bands) {

    // Aux variables
    arma::mat res(ts.n_rows, n_times, arma::fill::zeros);

    // For each pixel
    for (arma::uword i = 0; i < ts.n_rows; i++) {
        // For each band
        for (arma::uword c = 0; c < ts.n_cols; c = c + n_times) {


        }
        res.submat(i, 0, i, n_times - 1) = 1;
    }
    // Return the results
    return res;
}


arma::vec bocd_cal(const arma::vec& ts) {
    arma::uword T = 0;

}

arma::vec hazard(const arma::uword& lambda,
                 const arma::vec& r) {
    arma::vec v(r.size(), arma::fill::ones);

    return (v / lambda);
}

arma::vec hazard(const arma::uword& lambda,
                 const arma::uword& r) {
    arma::vec v(r, arma::fill::ones);

    return (v / lambda);
}


