loc <- cbind(runif(10), runif(10))
dist <- as.matrix(dist(loc))

# Check that results is okay
# Should return matrix
# Result should be positive definite
# Output should be symmetric (and square)
c1 <- try(covCauchyFn@fun(
  pars = c(rexp(2, rate = 1/10),
           runif(1, 1e-10, 2),
           rexp(1, rate = 1)),
  dist = dist))
c2 <- try(covPowExpFn@fun(
  pars = c(rexp(2, rate = 1/10),
           runif(1, 1e-10, 2)),
  dist = dist))
c3 <- try(covExpFn@fun(
  pars = rexp(2, rate = 1/10),
  dist = dist))
c4 <- try(covPowExpFn@fun(
  pars = c(rexp(2, rate = 1), runif(1,0,2)),
  dist = dist))
# Boundary cases
c5 <- try(covMaternFn@fun(
  pars = c(0.1,1, Inf),
  dist = dist))
c6 <- try(covPowExpFn@fun(
  pars = c(1,0.1,2),
  dist = dist))
c7 <- try(covCauchyFn@fun(
  pars = c(1,0.1,2, 0.1),
  dist = dist))
c8 <- try(covCauchyFn@fun(
  pars = c(1,0.1,2, Inf),
  dist = dist))

# Check output matrix is symmetric
#  and positive definite
tinytest::expect_true(
  isTRUE(all(
    !inherits(c1, "try-error"),
    is.matrix(c1),
    isSymmetric(c1),
    all(eigen(c1, only.values = TRUE)$values > 0)
    )))
tinytest::expect_true(
  isTRUE(all(
    !inherits(c2, "try-error"),
    is.matrix(c2),
    isSymmetric(c2),
    all(eigen(c2, only.values = TRUE)$values > 0)
  )))
tinytest::expect_true(
  isTRUE(all(
    !inherits(c3, "try-error"),
    is.matrix(c3),
    isSymmetric(c3),
    all(eigen(c3, only.values = TRUE)$values > 0)
  )))
tinytest::expect_true(
  isTRUE(all(
    !inherits(c4, "try-error"),
    is.matrix(c4),
    isSymmetric(c4),
    all(eigen(c4, only.values = TRUE)$values > 0)
  )))
tinytest::expect_true(
  isTRUE(all(
    !inherits(c5, "try-error"),
    is.matrix(c5),
    isSymmetric(c5),
    all(eigen(c5, only.values = TRUE)$values > 0)
  )))
tinytest::expect_true(
  isTRUE(all(
    !inherits(c6, "try-error"),
    is.matrix(c6),
    isSymmetric(c6),
    all(eigen(c6, only.values = TRUE)$values > 0)
  )))
tinytest::expect_true(
  isTRUE(all(
    !inherits(c7, "try-error"),
    is.matrix(c7),
    isSymmetric(c7),
    all(eigen(c7, only.values = TRUE)$values > 0)
  )))
tinytest::expect_true(
  isTRUE(all(
    !inherits(c8, "try-error"),
    is.matrix(c8),
    isSymmetric(c8),
    all(eigen(c8, only.values = TRUE)$values > 0)
  )))

# Check number of arguments
tinytest::expect_error(
  covMaternFn@fun(rep(1,2), dist = dist),
  )
tinytest::expect_error(
  covMaternFn@fun(rep(1,4), dist = dist),
)
tinytest::expect_error(
  covCauchyFn@fun(rep(1,3), dist = dist),
)
tinytest::expect_error(
  covPowExpFn@fun(rep(1,2), dist = dist),
)

# Check the case of vector is handled
#  correctly for Matern
dout <- covMaternFn@fun(pars = c(0.1, 1, 0.5),
                        dist = c(dist[1,]))
tinytest::expect_true(isTRUE(
  all(length(dout) == length(dist[1,]),
      dout > 0)))

# Check that returns variance for distance 0
tinytest::expect_true(
  covMaternFn@fun(pars = c(0.1, 1, rexp(1, rate = 1/100)),
                  dist = 0) == 0.1)
tinytest::expect_true(
  covExpFn@fun(pars = c(0.1, 1),
                  dist = 0) == 0.1)
tinytest::expect_true(
  covCauchyFn@fun(pars = c(0.1,
                           rexp(1),
                           runif(1,0,2),
                           rexp(1, rate = 1/10)),
               dist = 0) == 0.1)
tinytest::expect_true(
  covPowExpFn@fun(pars = c(0.1,
                           rexp(1),
                           runif(1,0,2)),
                  dist = 0) == 0.1)
