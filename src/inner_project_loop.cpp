#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
NumericMatrix inner_project_loop(int no_sp, int no_w, 
                NumericMatrix n, NumericMatrix A, NumericMatrix B,
                NumericMatrix S, IntegerVector w_min_idx) {
    n = Rcpp::clone(n);
    
    for (int i = 0; i < no_sp; i++) {
        int start = w_min_idx[i];
        if (start < 1) start = 1;
        for (int j = start; j < no_w; j++) {
            n(i,j) = (S(i,j) - A(i,j)*n(i,j-1)) / B(i,j);
        }
    }
    return n;
}
