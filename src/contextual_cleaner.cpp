#include <Rcpp.h>

using namespace Rcpp;

// compute outside indices of a vector as a mirror
IntegerVector locus_mirror3(int size, int leg) {
    IntegerVector res(size + 2 * leg);
    for (int i = 0; i < res.length(); ++i) {
        if (i < leg)
            res(i) = leg - i - 1;
        else if (i < size + leg)
            res(i) = i - leg;
        else
            res(i) = 2 * size + leg - i - 1;
    }
    return res;
}

// This function was implemented by Robert Hijmans in the
// terra  package (GPL>=3) and adapted by Felipe Carvalho.
// The source code can be found in
// https://github.com/rspatial/terra/blob/bcce14dd1778a36a43e2a211704feb8128f2c953/src/vecmath.h
double modal(const NumericVector& neigh) {

    std::map<double, size_t> count;
    for(size_t i=0; i<neigh.size(); i++) {
        if (std::isnan(neigh[i])) {
            return NAN;
        } else {
            count[neigh[i]]++;
        }
    }

    std::map<double, size_t>::iterator mode =
        std::max_element(count.begin(), count.end(),[] (const std::pair<double, size_t>& a,
                                     const std::pair<double, size_t>& b)->bool{ return a.second < b.second; } );

    return mode->first;
}

// [[Rcpp::export]]
NumericVector context_cleaner(const NumericMatrix& x, int ncols, int nrows,
                              int band, int window_size, int target_class,
                              int mode_class) {
    // initialize result vectors
    NumericVector res(x.nrow());
    NumericVector neigh(window_size * window_size);
    double modal_value;
    if (window_size < 1) {
        res = x(_, band);
        return res;
    }
    // compute window leg
    int leg = window_size / 2;
    // compute locus mirror
    IntegerVector loci = locus_mirror3(nrows, leg);
    IntegerVector locj = locus_mirror3(ncols, leg);
    // compute values for each pixel
    for (int i = 0; i < nrows; ++i) {
        for (int j = 0; j < ncols; ++j) {
            // window
            for (int wi = 0; wi < window_size; ++wi)
                for (int wj = 0; wj < window_size; ++wj)
                    neigh(wi * window_size + wj) =
                        x(loci(wi + i) * ncols + locj(wj + j), band);

            if (x(i * ncols + j, band) == target_class) {
                modal_value = modal(neigh);
                if (modal_value == mode_class) {
                    res(i * ncols + j) = mode_class;
                } else {
                    res(i * ncols + j) = target_class;
                }
            } else {
                res(i * ncols + j) = x(i * ncols + j, band);
            }
        }
    }
    return res;
}
