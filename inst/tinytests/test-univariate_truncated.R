# Check density
xp <- seq(0, 10, by = 1)
tinytest::expect_equal(dtnorm(xp, a = 0,
                              log = TRUE) - log(2),
                       dnorm(xp, log = TRUE))
tinytest::expect_equal(dtnorm(-xp, b = 0,
                              log = TRUE) - log(2),
                       dnorm(-xp, log = TRUE))
tinytest::expect_equal(dtnorm(xp, log = TRUE),
                       dnorm(xp, log = TRUE))


# Tests for univariate truncated Normal distribution
set.seed(1234)
# Check unconditional
kst <- ks.test(
  x = rtnorm(
    n = 1000L,
    mean = 2,
    sd = 1,
    a = -Inf,
    b = Inf
  ),
  y = pnorm,
  mean = 2,
  sd = 1
)
tinytest::expect_true(kst$p.value > 0.01)

# Check restricted
n <- 1e5L
a = 1
b = 2
mu = 20
sigma = 1
alpha <- (a - mu) / sigma
beta <- (b - mu) / sigma
samp <- rtnorm(
  n = n,
  mean = mu,
  sd = sigma,
  a = a,
  b = b
)
mean_tnorm <-
  mu + (dnorm(alpha) - dnorm(beta)) /
  (pnorm(beta) - pnorm(alpha)) * sigma
sd_tnorm <-
  sqrt(sigma ^ 2 * (
    1 + (alpha * dnorm(alpha) - beta * dnorm(beta)) / (pnorm(beta) - pnorm(alpha)) - (dnorm(alpha) -
 dnorm(beta)) ^ 2 / (pnorm(beta) - pnorm(alpha)) ^ 2
  ))
tinytest::expect_true(a < min(samp) & max(samp) < b)
tinytest::expect_equal(mean(samp),
                       mean_tnorm,
                       tolerance = qnorm(0.99) *
                         sd_tnorm / sqrt(n))
tinytest::expect_error(rtnorm(n = 10, a = 38))
