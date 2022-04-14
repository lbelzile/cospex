# Testing univariate Laplace update
#
# Run a Markov chain with a conjugate model
# to compare with the answer and check if the
# model converges to the posterior

# Beta-binomial
loglik <- function(x, n, p){
  dbinom(x = x, size = n, prob = p, log = TRUE)
  }
loglik_grad <- function(x, n, p){
  x/p - (n-x)/(1-p)
  }
loglik_hessian <- function(x, n, p){
  -x/p^2 - (n-x)/(1-p)^2
}

logprior <- function(p,  alpha = 1, beta = 1){
  dbeta(x = p, shape1 = alpha, shape2 = beta, log = TRUE)
}
logprior_grad <- function(p,  alpha = 1, beta = 1){
  (alpha-1)/p - (beta-1)/(1-p)
}
logprior_hessian <- function(p, alpha = 1, beta = 1){
  -(alpha-1)/p^2 - (beta-1)/(1-p)^2
}

# Check derivatives
x0 <- 10
n <- 40
p0 <- 0.4
alpha = 1
beta = 2
isTRUE(all.equal(numDeriv::grad(function(p){
  loglik(x = x0, n = n, p = p)}, x = p0),
  loglik_grad(x = x0, n = n, p = p0)))

isTRUE(all.equal(numDeriv::hessian(function(p){
  loglik(x = x0, n = n, p = p)}, x = p0)[1,1],
  loglik_hessian(x = x0, n = n, p = p0)))

isTRUE(all.equal(numDeriv::grad(function(p){
  logprior(alpha = alpha, beta = beta, p = p)}, x = p0),
  logprior_grad(alpha = alpha, beta = beta, p = p0)))

isTRUE(all.equal(numDeriv::hessian(function(p){
  logprior(alpha = alpha, beta = beta, p = p)}, x = p0)[1,1],
  logprior_hessian(alpha = alpha, beta = beta, p = p0)))


chain <- numeric(length = 1e4L)
chain[1] <- 0.8
for(i in 2:length(chain)){
  chain[i] <-
    univariate_laplace_update(
      par_curr = chain[i-1],
      par_name = "p",
      loglik = loglik,
      loglik_grad = loglik_grad,
      loglik_hessian = loglik_hessian,
      logprior = logprior,
      logprior_grad = logprior_grad,
      logprior_hessian = logprior_hessian,
      lb = 0,
      ub = 1,
      damping = 0.8,
      x = x0,
      n = n)
}
chain0 <- chain[-(1:100)]
acf(chain0)[1]
mean(chain0) - (1+x0)/(2 + n)
var(chain0)

hist(chain0,
     probability = TRUE,
     breaks = 100)
curve(dbeta(x,
            shape1 = 1 + x0,
            shape2 = 1 + n - x0),
      from = 0,
      to = 1,
      add = TRUE,
      col = 2)

# Run a vanilla random walk Metropolis Sampler
chain2 <- numeric(1e4L)
curr <- 0.5
sd_prop <- 0.05
for(i in seq_along(chain2)){
  prop <- rtnorm(n = 1,
                 a = 0,
                 b = 1,
                 mean = curr,
                 sd = sd_prop)
 ratio <-
   loglik(x = x0, n = n, p = prop) -
    loglik(x = x0, n = n, p = curr) +
  logprior(p = prop) - logprior(p = curr) +
   dtnorm(x = curr, a = 0, b = 1,
          mean = prop, sd = sd_prop, log = TRUE) -
   dtnorm(x = prop, a = 0, b = 1,
          mean = curr, sd = sd_prop, log = TRUE)
  if(ratio > log(runif(1))){
    chain2[i] <- prop
    curr <- prop
  } else{
    chain2[i] <- curr
  }
}
hist(chain2, probability = TRUE, breaks = 100)
curve(dbeta(x,
            shape1 = 1 + x0,
            shape2 = 1 + n - x0),
      from = 0,
      to = 1,
      add = TRUE,
      col = 2)
