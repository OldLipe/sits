#include <RcppArmadillo.h>
#include <iostream>
#include <iomanip>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

// [[Rcpp::export]]
arma::mat C_cesbio_calc_rcr(const arma::mat& values,
                            const arma::uword xa,
                            const arma::uword xb) {
    // Create the result matrix
    arma::mat res(
            values.n_rows, values.n_cols, arma::fill::value(arma::datum::nan)
    );
    arma::mat values_after(1, xb, arma::fill::zeros);
    arma::mat values_before(1, xa, arma::fill::zeros);

    // Get the days before and after
    arma::uword db = xb - 1;
    arma::uword da = xa - 1;
    // Aux variables
    double mean_a, mean_b = 0;
    // Index date
    arma::uword i_idx = 0;
    arma::uvec valid_idxa;
    arma::uvec valid_idxb;

    // // For each pixel
    for (arma::uword c = 0; c < values.n_rows; ++c) {
        // For each date
        // We subtract by 2 because the date is inclusive
        for (arma::uword i = 0; (i + db) <= (values.n_cols - da - 2); ++i) {
            i_idx = i + db;

            values_before = values.submat(c, i, c, i_idx);
            valid_idxb = arma::find_finite(values_before);

            mean_b = arma::mean(
                arma::mean( values_before.cols(valid_idxb), 1)
            );

            values_after = values.submat(c, i_idx + 1, c, i_idx + 1 + da);
            valid_idxa = arma::find_finite(values_after);

            mean_a = arma::mean(
                arma::mean( values_after.cols(valid_idxa), 1)
            );
            res(c, i_idx) = mean_a / mean_b;
        }
    }

    return res;
}

// [[Rcpp::export]]
arma::mat C_cesbio_detect_shadow(const arma::mat& rcr,
                                 const double shadow_value) {
    // Create the result matrix
    arma::mat res(
            rcr.n_rows, 1, arma::fill::value(arma::datum::nan)
    );
    arma::uword min_idx = 0;
    for (arma::uword c = 0; c < rcr.n_rows; ++c) {
        min_idx = rcr.row(c).index_min();
        if (rcr.at(c, min_idx) <= shadow_value) {
            res(c, 0) = rcr.row(c).index_min() + 1;
        }
    }

    return res;
}

// [[Rcpp::export]]
arma::mat C_cesbio_detect_neigh(const arma::mat& rcr,
                                 const double neigh_value) {
    // Create the result matrix
    arma::mat res(
            rcr.n_rows, 1, arma::fill::value(arma::datum::nan)
    );
    arma::uword min_idx = 0;
    for (arma::uword c = 0; c < rcr.n_rows; ++c) {
        min_idx = rcr.row(c).index_min();
        if (rcr.at(c, min_idx) <= neigh_value) {
            res(c, 0) = rcr.row(c).index_min() + 1;
        }
    }

    return res;
}
