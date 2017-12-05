#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix inner_project_loop(int no_sp, int no_w, 
                NumericMatrix n, NumericMatrix A, NumericMatrix B,
                NumericMatrix S, NumericVector w_min_idx) {
    
    for (int i = 0; i < no_sp; i++) {
        for (int j = w_min_idx[i]; j < no_w; j++) {
            n(i,j) = (S(i,j) - A(i,j)*n(i,j-1)) / B(i,j);
        }
    }
    return n;
}