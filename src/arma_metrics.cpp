#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::mat min_ts_2(arma::mat& x) {

    int nrows = x.n_rows;
    //int ncols = mtx.ncol();

    for (int i = 0; i < nrows; i++) {
        arma::vec v = x.col(i);

        x(i, 0) = arma::min(v);
    }



    return x;

}

// [[Rcpp::export]]
arma::vec max_mat(const arma::mat& mtx) {


    return  arma::max(mtx, 1);;
}

// [[Rcpp::export]]
arma::vec min_mat(const arma::mat& mtx) {

    return  arma::min(mtx, 1);;
}

// [[Rcpp::export]]
arma::vec mean_mat(const arma::mat& mtx) {

    return  arma::mean(mtx, 1);;
}

// [[Rcpp::export]]
arma::vec std_mat(const arma::mat& mtx) {

    return  arma::stddev(mtx, 0, 1);
}

// [[Rcpp::export]]
arma::vec amplitude_ts(const arma::mat& mtx) {

    return  arma::max(mtx, 1) - arma::min(mtx, 1);
}

// [[Rcpp::export]]
arma::vec fslope_ts(const arma::mat& mtx) {

    return  arma::max(arma::abs(arma::diff(mtx, 1, 1)), 1);
}

// [[Rcpp::export]]
arma::vec abs_sum_ts(const arma::mat& mtx) {

    return  arma::sum(arma::abs(mtx), 1);
}

// [[Rcpp::export]]
arma::vec amd_ts(const arma::mat& mtx) {

    return  arma::mean(arma::abs(arma::diff(mtx, 1, 1)), 1);
}

// [[Rcpp::export]]
arma::mat max_ts_2(arma::mat& x) {

    int nrows = x.n_rows;
    //int ncols = mtx.ncol();

    for (int i = 0; i < nrows; i++) {
        arma::rowvec v = x.row(i);
        Rcpp::Rcout << "The value is " <<v<< std::endl;

        Rcpp::Rcout << "max " <<arma::max(v)<< std::endl;

        Rcpp::Rcout << "min " <<arma::min(v)<< std::endl;

        Rcpp::Rcout << "mean " <<arma::min(v)<< std::endl;

        Rcpp::Rcout << "std " <<arma::stddev(v)<< std::endl;

        Rcpp::Rcout << "amplitude_ts " << arma::max(v) - arma::min(v)<< std::endl;

        Rcpp::Rcout << "fslope_ts " << arma::max(arma::abs(arma::diff(v))) << std::endl;

        Rcpp::Rcout << "abs_sum_ts " << arma::sum(arma::abs(v)) << std::endl;

        Rcpp::Rcout << "amd_ts " << arma::mean(arma::abs(arma::diff(v))) << std::endl;

        Rcpp::Rcout << "mse_ts " << arma::mean(arma::square(arma::abs(arma::fft(v)))) << std::endl;


        arma::vec P_025 = {0.25};
        arma::vec P_050 = {0.50};
        arma::vec P_075 = {0.75};
        Rcpp::Rcout << "fqr_ts " << arma::quantile(v, P_025) << std::endl;

        Rcpp::Rcout << "tqr_ts " << arma::quantile(v, P_075) << std::endl;

        Rcpp::Rcout << "sqr_ts " << arma::quantile(v, P_050) << std::endl;

        Rcpp::Rcout << "iqr_ts " << arma::quantile(v, P_075) - arma::quantile(v, P_025) << std::endl;
    }

    return x;

}
