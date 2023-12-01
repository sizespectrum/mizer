#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix inner_project_loop(int no_sp, 
                                 NumericMatrix n, 
                                 NumericMatrix A, 
                                 NumericMatrix B, 
                                 NumericMatrix S,
                                 NumericVector w_min_idx,
                                 NumericVector w_max_idx) {
    
    for (int i = 0; i < no_sp; i++) {
        for (int j = w_min_idx[i]; j < w_max_idx[i]; j++) {
            n(i,j) = (S(i,j) - A(i,j)*n(i,j-1)) / B(i,j);
        }
    }
    return n;
}
