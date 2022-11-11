// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// average_probs
NumericMatrix average_probs(const List& data_lst);
RcppExport SEXP _sits_average_probs(SEXP data_lstSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type data_lst(data_lstSEXP);
    rcpp_result_gen = Rcpp::wrap(average_probs(data_lst));
    return rcpp_result_gen;
END_RCPP
}
// weighted_probs
NumericMatrix weighted_probs(const List& data_lst, const NumericVector& weights);
RcppExport SEXP _sits_weighted_probs(SEXP data_lstSEXP, SEXP weightsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type data_lst(data_lstSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type weights(weightsSEXP);
    rcpp_result_gen = Rcpp::wrap(weighted_probs(data_lst, weights));
    return rcpp_result_gen;
END_RCPP
}
// weighted_uncert_probs
NumericMatrix weighted_uncert_probs(const List& data_lst, const List& unc_lst);
RcppExport SEXP _sits_weighted_uncert_probs(SEXP data_lstSEXP, SEXP unc_lstSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const List& >::type data_lst(data_lstSEXP);
    Rcpp::traits::input_parameter< const List& >::type unc_lst(unc_lstSEXP);
    rcpp_result_gen = Rcpp::wrap(weighted_uncert_probs(data_lst, unc_lst));
    return rcpp_result_gen;
END_RCPP
}
// C_kernel_median
NumericVector C_kernel_median(const NumericMatrix& x, int ncols, int nrows, int band, int window_size);
RcppExport SEXP _sits_C_kernel_median(SEXP xSEXP, SEXP ncolsSEXP, SEXP nrowsSEXP, SEXP bandSEXP, SEXP window_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type ncols(ncolsSEXP);
    Rcpp::traits::input_parameter< int >::type nrows(nrowsSEXP);
    Rcpp::traits::input_parameter< int >::type band(bandSEXP);
    Rcpp::traits::input_parameter< int >::type window_size(window_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(C_kernel_median(x, ncols, nrows, band, window_size));
    return rcpp_result_gen;
END_RCPP
}
// C_kernel_sum
NumericVector C_kernel_sum(const NumericMatrix& x, int ncols, int nrows, int band, int window_size);
RcppExport SEXP _sits_C_kernel_sum(SEXP xSEXP, SEXP ncolsSEXP, SEXP nrowsSEXP, SEXP bandSEXP, SEXP window_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type ncols(ncolsSEXP);
    Rcpp::traits::input_parameter< int >::type nrows(nrowsSEXP);
    Rcpp::traits::input_parameter< int >::type band(bandSEXP);
    Rcpp::traits::input_parameter< int >::type window_size(window_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(C_kernel_sum(x, ncols, nrows, band, window_size));
    return rcpp_result_gen;
END_RCPP
}
// C_kernel_mean
NumericVector C_kernel_mean(const NumericMatrix& x, int ncols, int nrows, int band, int window_size);
RcppExport SEXP _sits_C_kernel_mean(SEXP xSEXP, SEXP ncolsSEXP, SEXP nrowsSEXP, SEXP bandSEXP, SEXP window_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type ncols(ncolsSEXP);
    Rcpp::traits::input_parameter< int >::type nrows(nrowsSEXP);
    Rcpp::traits::input_parameter< int >::type band(bandSEXP);
    Rcpp::traits::input_parameter< int >::type window_size(window_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(C_kernel_mean(x, ncols, nrows, band, window_size));
    return rcpp_result_gen;
END_RCPP
}
// C_kernel_sd
NumericVector C_kernel_sd(const NumericMatrix& x, int ncols, int nrows, int band, int window_size);
RcppExport SEXP _sits_C_kernel_sd(SEXP xSEXP, SEXP ncolsSEXP, SEXP nrowsSEXP, SEXP bandSEXP, SEXP window_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type ncols(ncolsSEXP);
    Rcpp::traits::input_parameter< int >::type nrows(nrowsSEXP);
    Rcpp::traits::input_parameter< int >::type band(bandSEXP);
    Rcpp::traits::input_parameter< int >::type window_size(window_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(C_kernel_sd(x, ncols, nrows, band, window_size));
    return rcpp_result_gen;
END_RCPP
}
// C_kernel_var
NumericVector C_kernel_var(const NumericMatrix& x, int ncols, int nrows, int band, int window_size);
RcppExport SEXP _sits_C_kernel_var(SEXP xSEXP, SEXP ncolsSEXP, SEXP nrowsSEXP, SEXP bandSEXP, SEXP window_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type ncols(ncolsSEXP);
    Rcpp::traits::input_parameter< int >::type nrows(nrowsSEXP);
    Rcpp::traits::input_parameter< int >::type band(bandSEXP);
    Rcpp::traits::input_parameter< int >::type window_size(window_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(C_kernel_var(x, ncols, nrows, band, window_size));
    return rcpp_result_gen;
END_RCPP
}
// C_kernel_min
NumericVector C_kernel_min(const NumericMatrix& x, int ncols, int nrows, int band, int window_size);
RcppExport SEXP _sits_C_kernel_min(SEXP xSEXP, SEXP ncolsSEXP, SEXP nrowsSEXP, SEXP bandSEXP, SEXP window_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type ncols(ncolsSEXP);
    Rcpp::traits::input_parameter< int >::type nrows(nrowsSEXP);
    Rcpp::traits::input_parameter< int >::type band(bandSEXP);
    Rcpp::traits::input_parameter< int >::type window_size(window_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(C_kernel_min(x, ncols, nrows, band, window_size));
    return rcpp_result_gen;
END_RCPP
}
// C_kernel_max
NumericVector C_kernel_max(const NumericMatrix& x, int ncols, int nrows, int band, int window_size);
RcppExport SEXP _sits_C_kernel_max(SEXP xSEXP, SEXP ncolsSEXP, SEXP nrowsSEXP, SEXP bandSEXP, SEXP window_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type ncols(ncolsSEXP);
    Rcpp::traits::input_parameter< int >::type nrows(nrowsSEXP);
    Rcpp::traits::input_parameter< int >::type band(bandSEXP);
    Rcpp::traits::input_parameter< int >::type window_size(window_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(C_kernel_max(x, ncols, nrows, band, window_size));
    return rcpp_result_gen;
END_RCPP
}
// C_kernel_bayes_mean
NumericMatrix C_kernel_bayes_mean(const NumericMatrix& x, int ncols, int nrows, int window_size);
RcppExport SEXP _sits_C_kernel_bayes_mean(SEXP xSEXP, SEXP ncolsSEXP, SEXP nrowsSEXP, SEXP window_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type ncols(ncolsSEXP);
    Rcpp::traits::input_parameter< int >::type nrows(nrowsSEXP);
    Rcpp::traits::input_parameter< int >::type window_size(window_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(C_kernel_bayes_mean(x, ncols, nrows, window_size));
    return rcpp_result_gen;
END_RCPP
}
// C_kernel_bayes_var
NumericMatrix C_kernel_bayes_var(const NumericMatrix& x, int ncols, int nrows, int window_size);
RcppExport SEXP _sits_C_kernel_bayes_var(SEXP xSEXP, SEXP ncolsSEXP, SEXP nrowsSEXP, SEXP window_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type ncols(ncolsSEXP);
    Rcpp::traits::input_parameter< int >::type nrows(nrowsSEXP);
    Rcpp::traits::input_parameter< int >::type window_size(window_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(C_kernel_bayes_var(x, ncols, nrows, window_size));
    return rcpp_result_gen;
END_RCPP
}
// C_bayes_posterior
NumericMatrix C_bayes_posterior(const NumericMatrix& x, const NumericVector& s, const NumericMatrix& m, const NumericMatrix& v);
RcppExport SEXP _sits_C_bayes_posterior(SEXP xSEXP, SEXP sSEXP, SEXP mSEXP, SEXP vSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type s(sSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const NumericMatrix& >::type v(vSEXP);
    rcpp_result_gen = Rcpp::wrap(C_bayes_posterior(x, s, m, v));
    return rcpp_result_gen;
END_RCPP
}
// C_label_max_prob
arma::colvec C_label_max_prob(const arma::mat& x);
RcppExport SEXP _sits_C_label_max_prob(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(C_label_max_prob(x));
    return rcpp_result_gen;
END_RCPP
}
// linear_interp
NumericMatrix linear_interp(NumericMatrix& mtx);
RcppExport SEXP _sits_linear_interp(SEXP mtxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix& >::type mtx(mtxSEXP);
    rcpp_result_gen = Rcpp::wrap(linear_interp(mtx));
    return rcpp_result_gen;
END_RCPP
}
// linear_interp_vec
NumericVector linear_interp_vec(NumericVector& vec);
RcppExport SEXP _sits_linear_interp_vec(SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector& >::type vec(vecSEXP);
    rcpp_result_gen = Rcpp::wrap(linear_interp_vec(vec));
    return rcpp_result_gen;
END_RCPP
}
// batch_calc
arma::mat batch_calc(const int& n_pixels, const int& max_lines_per_batch);
RcppExport SEXP _sits_batch_calc(SEXP n_pixelsSEXP, SEXP max_lines_per_batchSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const int& >::type n_pixels(n_pixelsSEXP);
    Rcpp::traits::input_parameter< const int& >::type max_lines_per_batch(max_lines_per_batchSEXP);
    rcpp_result_gen = Rcpp::wrap(batch_calc(n_pixels, max_lines_per_batch));
    return rcpp_result_gen;
END_RCPP
}
// C_nnls_solver_batch
arma::mat C_nnls_solver_batch(const arma::mat& x, const arma::mat& em, const bool rmse, const int max_it, const float tol);
RcppExport SEXP _sits_C_nnls_solver_batch(SEXP xSEXP, SEXP emSEXP, SEXP rmseSEXP, SEXP max_itSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type em(emSEXP);
    Rcpp::traits::input_parameter< const bool >::type rmse(rmseSEXP);
    Rcpp::traits::input_parameter< const int >::type max_it(max_itSEXP);
    Rcpp::traits::input_parameter< const float >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(C_nnls_solver_batch(x, em, rmse, max_it, tol));
    return rcpp_result_gen;
END_RCPP
}
// C_nnls_solver
arma::mat C_nnls_solver(const arma::mat& x, const arma::mat& em, const bool rmse, const int max_it, const float tol);
RcppExport SEXP _sits_C_nnls_solver(SEXP xSEXP, SEXP emSEXP, SEXP rmseSEXP, SEXP max_itSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type em(emSEXP);
    Rcpp::traits::input_parameter< const bool >::type rmse(rmseSEXP);
    Rcpp::traits::input_parameter< const int >::type max_it(max_itSEXP);
    Rcpp::traits::input_parameter< const float >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(C_nnls_solver(x, em, rmse, max_it, tol));
    return rcpp_result_gen;
END_RCPP
}
// C_normalize_data
arma::mat C_normalize_data(const arma::mat& data, const arma::rowvec& min, const arma::rowvec& max);
RcppExport SEXP _sits_C_normalize_data(SEXP dataSEXP, SEXP minSEXP, SEXP maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type min(minSEXP);
    Rcpp::traits::input_parameter< const arma::rowvec& >::type max(maxSEXP);
    rcpp_result_gen = Rcpp::wrap(C_normalize_data(data, min, max));
    return rcpp_result_gen;
END_RCPP
}
// C_normalize_data_0
NumericMatrix C_normalize_data_0(const NumericMatrix& data, const double& min, const double& max);
RcppExport SEXP _sits_C_normalize_data_0(SEXP dataSEXP, SEXP minSEXP, SEXP maxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const double& >::type min(minSEXP);
    Rcpp::traits::input_parameter< const double& >::type max(maxSEXP);
    rcpp_result_gen = Rcpp::wrap(C_normalize_data_0(data, min, max));
    return rcpp_result_gen;
END_RCPP
}
// max_sampling
DataFrame max_sampling(const IntegerMatrix& data, const int band, const int img_nrow, const int img_ncol, const int window_size);
RcppExport SEXP _sits_max_sampling(SEXP dataSEXP, SEXP bandSEXP, SEXP img_nrowSEXP, SEXP img_ncolSEXP, SEXP window_sizeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const IntegerMatrix& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const int >::type band(bandSEXP);
    Rcpp::traits::input_parameter< const int >::type img_nrow(img_nrowSEXP);
    Rcpp::traits::input_parameter< const int >::type img_ncol(img_ncolSEXP);
    Rcpp::traits::input_parameter< const int >::type window_size(window_sizeSEXP);
    rcpp_result_gen = Rcpp::wrap(max_sampling(data, band, img_nrow, img_ncol, window_size));
    return rcpp_result_gen;
END_RCPP
}
// bayes_smoother
arma::mat bayes_smoother(const arma::mat& m, const arma::uword m_nrow, const arma::uword m_ncol, const arma::mat& w, const arma::mat& sigma, bool covar_sigma0, const double neigh_fraction);
RcppExport SEXP _sits_bayes_smoother(SEXP mSEXP, SEXP m_nrowSEXP, SEXP m_ncolSEXP, SEXP wSEXP, SEXP sigmaSEXP, SEXP covar_sigma0SEXP, SEXP neigh_fractionSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type m_nrow(m_nrowSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type m_ncol(m_ncolSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type w(wSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigma(sigmaSEXP);
    Rcpp::traits::input_parameter< bool >::type covar_sigma0(covar_sigma0SEXP);
    Rcpp::traits::input_parameter< const double >::type neigh_fraction(neigh_fractionSEXP);
    rcpp_result_gen = Rcpp::wrap(bayes_smoother(m, m_nrow, m_ncol, w, sigma, covar_sigma0, neigh_fraction));
    return rcpp_result_gen;
END_RCPP
}
// bilateral_smoother
arma::mat bilateral_smoother(const arma::mat& m, const arma::uword m_nrow, const arma::uword m_ncol, const arma::mat& w, double tau);
RcppExport SEXP _sits_bilateral_smoother(SEXP mSEXP, SEXP m_nrowSEXP, SEXP m_ncolSEXP, SEXP wSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type m(mSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type m_nrow(m_nrowSEXP);
    Rcpp::traits::input_parameter< const arma::uword >::type m_ncol(m_ncolSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(bilateral_smoother(m, m_nrow, m_ncol, w, tau));
    return rcpp_result_gen;
END_RCPP
}
// smooth_sg
arma::vec smooth_sg(const arma::vec& data, const arma::mat& f_res, const int& p, const int& n);
RcppExport SEXP _sits_smooth_sg(SEXP dataSEXP, SEXP f_resSEXP, SEXP pSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type f_res(f_resSEXP);
    Rcpp::traits::input_parameter< const int& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(smooth_sg(data, f_res, p, n));
    return rcpp_result_gen;
END_RCPP
}
// smooth_sg_mtx
arma::mat smooth_sg_mtx(const arma::mat& data, const arma::mat& f_res, const int& p, const int& n);
RcppExport SEXP _sits_smooth_sg_mtx(SEXP dataSEXP, SEXP f_resSEXP, SEXP pSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type f_res(f_resSEXP);
    Rcpp::traits::input_parameter< const int& >::type p(pSEXP);
    Rcpp::traits::input_parameter< const int& >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(smooth_sg_mtx(data, f_res, p, n));
    return rcpp_result_gen;
END_RCPP
}
// smooth_whit
NumericVector smooth_whit(const NumericVector& data, const double& lambda, const int& length);
RcppExport SEXP _sits_smooth_whit(SEXP dataSEXP, SEXP lambdaSEXP, SEXP lengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const double& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const int& >::type length(lengthSEXP);
    rcpp_result_gen = Rcpp::wrap(smooth_whit(data, lambda, length));
    return rcpp_result_gen;
END_RCPP
}
// smooth_whit_mtx
NumericMatrix smooth_whit_mtx(const NumericMatrix& data, const double& lambda, const int& length);
RcppExport SEXP _sits_smooth_whit_mtx(SEXP dataSEXP, SEXP lambdaSEXP, SEXP lengthSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericMatrix& >::type data(dataSEXP);
    Rcpp::traits::input_parameter< const double& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const int& >::type length(lengthSEXP);
    rcpp_result_gen = Rcpp::wrap(smooth_whit_mtx(data, lambda, length));
    return rcpp_result_gen;
END_RCPP
}
// C_entropy_probs
arma::mat C_entropy_probs(const arma::mat& x);
RcppExport SEXP _sits_C_entropy_probs(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(C_entropy_probs(x));
    return rcpp_result_gen;
END_RCPP
}
// C_margin_probs
arma::mat C_margin_probs(const arma::mat& x);
RcppExport SEXP _sits_C_margin_probs(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(C_margin_probs(x));
    return rcpp_result_gen;
END_RCPP
}
// C_least_probs
arma::mat C_least_probs(const arma::mat& x);
RcppExport SEXP _sits_C_least_probs(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(C_least_probs(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_sits_average_probs", (DL_FUNC) &_sits_average_probs, 1},
    {"_sits_weighted_probs", (DL_FUNC) &_sits_weighted_probs, 2},
    {"_sits_weighted_uncert_probs", (DL_FUNC) &_sits_weighted_uncert_probs, 2},
    {"_sits_C_kernel_median", (DL_FUNC) &_sits_C_kernel_median, 5},
    {"_sits_C_kernel_sum", (DL_FUNC) &_sits_C_kernel_sum, 5},
    {"_sits_C_kernel_mean", (DL_FUNC) &_sits_C_kernel_mean, 5},
    {"_sits_C_kernel_sd", (DL_FUNC) &_sits_C_kernel_sd, 5},
    {"_sits_C_kernel_var", (DL_FUNC) &_sits_C_kernel_var, 5},
    {"_sits_C_kernel_min", (DL_FUNC) &_sits_C_kernel_min, 5},
    {"_sits_C_kernel_max", (DL_FUNC) &_sits_C_kernel_max, 5},
    {"_sits_C_kernel_bayes_mean", (DL_FUNC) &_sits_C_kernel_bayes_mean, 4},
    {"_sits_C_kernel_bayes_var", (DL_FUNC) &_sits_C_kernel_bayes_var, 4},
    {"_sits_C_bayes_posterior", (DL_FUNC) &_sits_C_bayes_posterior, 4},
    {"_sits_C_label_max_prob", (DL_FUNC) &_sits_C_label_max_prob, 1},
    {"_sits_linear_interp", (DL_FUNC) &_sits_linear_interp, 1},
    {"_sits_linear_interp_vec", (DL_FUNC) &_sits_linear_interp_vec, 1},
    {"_sits_batch_calc", (DL_FUNC) &_sits_batch_calc, 2},
    {"_sits_C_nnls_solver_batch", (DL_FUNC) &_sits_C_nnls_solver_batch, 5},
    {"_sits_C_nnls_solver", (DL_FUNC) &_sits_C_nnls_solver, 5},
    {"_sits_C_normalize_data", (DL_FUNC) &_sits_C_normalize_data, 3},
    {"_sits_C_normalize_data_0", (DL_FUNC) &_sits_C_normalize_data_0, 3},
    {"_sits_max_sampling", (DL_FUNC) &_sits_max_sampling, 5},
    {"_sits_bayes_smoother", (DL_FUNC) &_sits_bayes_smoother, 7},
    {"_sits_bilateral_smoother", (DL_FUNC) &_sits_bilateral_smoother, 5},
    {"_sits_smooth_sg", (DL_FUNC) &_sits_smooth_sg, 4},
    {"_sits_smooth_sg_mtx", (DL_FUNC) &_sits_smooth_sg_mtx, 4},
    {"_sits_smooth_whit", (DL_FUNC) &_sits_smooth_whit, 3},
    {"_sits_smooth_whit_mtx", (DL_FUNC) &_sits_smooth_whit_mtx, 3},
    {"_sits_C_entropy_probs", (DL_FUNC) &_sits_C_entropy_probs, 1},
    {"_sits_C_margin_probs", (DL_FUNC) &_sits_C_margin_probs, 1},
    {"_sits_C_least_probs", (DL_FUNC) &_sits_C_least_probs, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_sits(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
