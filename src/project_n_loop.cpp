#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix project_n_loop(NumericMatrix n, NumericMatrix a, NumericMatrix b, NumericMatrix c, 
                             NumericMatrix S, NumericVector w_min_idx) {
    int no_sp = n.nrow();
    int no_w = n.ncol();
    
    // Temporary vectors for Thomas algorithm
    // Allocated once to be reused across species
    NumericVector c_prime(no_w);
    NumericVector d_prime(no_w);
    
    for (int i = 0; i < no_sp; i++) {
        // R uses 1-based indexing for w_min_idx, so subtract 1
        int j_start = w_min_idx[i] - 1; 
        
        if (j_start >= no_w) continue;
        
        // Thomas Algorithm
        // Solve A * n = S for the species range [j_start, no_w-1]
        
        // Forward elimination
        double b_val = b(i, j_start);
        if (b_val == 0) b_val = 1e-10; // Avoid division by zero
        
        c_prime[j_start] = c(i, j_start) / b_val;
        d_prime[j_start] = S(i, j_start) / b_val;
        
        for (int j = j_start + 1; j < no_w; j++) {
            double a_val = a(i, j);
            double temp = b(i, j) - a_val * c_prime[j - 1];
            if (temp == 0) temp = 1e-10; // Avoid division by zero
            
            if (j < no_w - 1) {
                c_prime[j] = c(i, j) / temp;
            }
            d_prime[j] = (S(i, j) - a_val * d_prime[j - 1]) / temp;
        }
        
        // Backward substitution
        n(i, no_w - 1) = d_prime[no_w - 1];
        for (int j = no_w - 2; j >= j_start; j--) {
            n(i, j) = d_prime[j] - c_prime[j] * n(i, j + 1);
        }
        
        // Note: values of n for j < j_start remain unchanged as desired
    }
    
    return n;
}
