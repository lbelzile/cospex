# Test that gradients are correctly implemented
# also with vector tau argument -
# grad returns a vector

x <- rnorm(100)
censoring <- x < 0
mu <- runif(100)
tau <- 1.2
tinytest::expect_equal(
  numDeriv::grad(function(sc){
  cospec_cond_loglik(x = x, censoring = censoring, mu = mu, tau = sc)}, x = tau),
  cospec_cond_loglik_grad_tau(
  x = x,
  censoring = censoring,
  mu = mu,
  tau = tau)
)

tinytest::expect_equal(
  numDeriv::grad(function(sc){
    cospec_cond_loglik_grad_tau(x = x, censoring = censoring, mu = mu, tau = sc)}, x = tau),
  cospec_cond_loglik_hessian_tau(
    x = x,
    censoring = censoring,
    mu = mu,
    tau = tau)
)


tinytest::expect_equal(
  numDeriv::hessian(function(sc){
    cospec_cond_loglik(x = x, censoring = censoring, mu = mu, tau = sc)}, x = tau)[1,1],
  cospec_cond_loglik_hessian_tau(
    x = x,
    censoring = censoring,
    mu = mu,
    tau = tau)
)

# Test that gradients are correctly implemented
# also with vector tau argument -
# grad returns a vector

x <- rnorm(100)
censoring <- x < 0
mu <- runif(100)
tau <- rexp(100)

tinytest::expect_equal(
  numDeriv::grad(function(sc){
    cospec_cond_loglik(x = x, censoring = censoring, mu = mu, tau = sc)}, x = tau),
  cospec_cond_loglik_grad_tau(
    x = x,
    censoring = censoring,
    mu = mu,
    tau = tau)
)

tinytest::expect_equal(
  numDeriv::grad(function(sc){
    cospec_cond_loglik_grad_tau(x = x, censoring = censoring, mu = mu, tau = sc)}, x = tau),
  cospec_cond_loglik_hessian_tau(
    x = x,
    censoring = censoring,
    mu = mu,
    tau = tau)
)


tinytest::expect_equal(
  numDeriv::hessian(function(sc){
    cospec_cond_loglik(x = x, censoring = censoring, mu = mu, tau = sc)}, x = tau)[1,1],
  cospec_cond_loglik_hessian_tau(
    x = x,
    censoring = censoring,
    mu = mu,
    tau = tau)
)

tinytest::expect_equal(
  numDeriv::grad(function(muv){
    cospec_cond_loglik(x = x, censoring = censoring, mu = muv, tau = tau)}, x = mu),
  cospec_cond_loglik_grad_mu(
    x = x,
    censoring = censoring,
    mu = mu,
    tau = tau)
)

# The hessian is a diagonal matrix
# the output from the function is awkward, because it is meant to be used in conjunction with the chain rule
tinytest::expect_equal(
  diag(numDeriv::hessian(function(muv){
    cospec_cond_loglik(x = x, censoring = censoring, mu = muv, tau = tau)}, x = mu)),
  cospec_cond_loglik_hessian_mu(
    x = x,
    censoring = censoring,
    mu = mu,
    tau = tau)$pmixed
)
