test_that("thomas_solve solves a tridiagonal system", {
    # A well-conditioned tridiagonal system solved against base R's solve()
    a <- c(0, 1, 2, 1)      # lower diagonal, a[1] unused
    b <- c(4, 5, 6, 4)      # main diagonal
    c <- c(1, 2, 1, 0)      # upper diagonal, c[n] unused
    d <- c(7, 8, 9, 10)

    A <- diag(b)
    for (i in 2:4) A[i, i - 1] <- a[i]
    for (i in 1:3) A[i, i + 1] <- c[i]

    expect_equal(thomas_solve(a, b, c, d), as.vector(solve(A, d)))
})

test_that("thomas_solve handles a system of size one", {
    expect_equal(thomas_solve(0, 2, 0, 6), 3)
})

test_that("thomas_solve agrees with solve() on a random diagonally dominant system", {
    set.seed(42)
    n <- 20
    a <- c(0, runif(n - 1, 0.1, 1))
    c <- c(runif(n - 1, 0.1, 1), 0)
    # make it diagonally dominant so the algorithm is stable without pivoting
    b <- abs(a) + abs(c) + runif(n, 1, 2)
    d <- rnorm(n)

    A <- diag(b)
    for (i in 2:n) A[i, i - 1] <- a[i]
    for (i in 1:(n - 1)) A[i, i + 1] <- c[i]

    expect_equal(thomas_solve(a, b, c, d), as.vector(solve(A, d)))
})
