#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
arma::mat normalize_scale_data(const arma::mat& data) {

    // this functions implements scale(x, center = TRUE, scale = TRUE)

    // int nrows = data.n_rows;
    // int ncols = data.n_cols;

    // arma::mat new_data(nrows, ncols);

    // mean column
    arma::mat center = arma::mean(data, 0);
    Rcpp::Rcout << "1 " << center << std::endl;
    arma::mat Y = data.each_row() - center;

    arma::vec scale = Y.each_col([](const arma::vec& b){
        sqrt(arma::sum(arma::square(b)) / (arma::max(b.size() - 1)));
        //arma::mean(b);
        });

    Rcpp::Rcout << "2 " << scale << std::endl;

    arma::mat x = Y.each_row() / scale;

    return Y;
}

// // [[Rcpp::export]]
// double scale_norm(const arma::vec& Y) {
//
//     double dd =
//     return dd;
// }
