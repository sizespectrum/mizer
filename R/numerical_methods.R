#' Thomas algorithm for solving tridiagonal system
#'
#' Solves a tridiagonal system of linear equations A * x = d
#' where A is a tridiagonal matrix defined by diagonals a, b, c.
#'
#' @param a Lower diagonal (length N). a_1 is ignored/0.
#' @param b Main diagonal (length N).
#' @param c Upper diagonal (length N). c_N is ignored/0.
#' @param d Right hand side vector (length N).
#'
#' @return Solution vector x (length N).
#' @noRd
thomas_solve <- function(a, b, c, d) {
    n <- length(d)
    c_prime <- numeric(n)
    d_prime <- numeric(n)
    
    # Forward elimination
    c_prime[1] <- c[1] / b[1]
    d_prime[1] <- d[1] / b[1]
    
    if (n > 1) {
        for (i in 2:n) {
            temp <- b[i] - a[i] * c_prime[i - 1]
            if (i < n) {
                c_prime[i] <- c[i] / temp
            }
            d_prime[i] <- (d[i] - a[i] * d_prime[i - 1]) / temp
        }
    }
    
    # Backward substitution
    x <- numeric(n)
    x[n] <- d_prime[n]
    if (n > 1) {
        for (i in (n - 1):1) {
            x[i] <- d_prime[i] - c_prime[i] * x[i + 1]
        }
    }
    
    return(x)
}
