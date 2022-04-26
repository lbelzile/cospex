if(false){

# We simulate data from a conjugate model to check the algorithm
# Use a simple linear regression model
# Y | beta ~ No(Xbeta, sigmasq*I)
# beta ~ No(0, Omega)
set.seed(1234)
n <- 100
p <- 5
sigma <- sqrt(2)
Omega <- rWishart(n = 1,
                  df = p + 5,
                  Sigma = diag(p) +
                    matrix(1, nrow = p, ncol = p))[,,1]
X <- matrix(data = runif(n*p), nrow = n, ncol = p)
beta0 <- rnorm(p)
y <- X %*% beta0 + rnorm(n = n, sd = sigma)

logLikFn <- function(betaf, y, Xmat, sigma){
  sum(dnorm(x = y,
        mean = as.vector(Xmat %*% betaf),
        sd = sigma,
        log = TRUE))
}
#Check function
logLikFn(beta0, y = y, Xmat = X, sigma = sigma)
# Compute mean and variance of posterior
postmean <- c(Omega %*% t(X) %*% solve(sigma^2*diag(n) + X %*% Omega %*% t(X)) %*% y)
postcov <- solve(solve(Omega) + crossprod(X)/ sigma^2)


nrep <- 1e5L + 250
beta_samp <- matrix(NA, nrow = nrep, ncol = p)
beta_samp[1,] <- rep(0, p)
for(i in 2:nrep){
  beta_samp[i,] <-
    elliptical_slice_sampling(
    curState = beta_samp[i-1,],
    logLikFn = logLikFn,
    normSamp = c(mvtnorm::rmvnorm(n = 1, sigma = Omega)),
    y = y,
    Xmat = X,
    sigma = sigma)$x
}
beta_samp <- beta_samp[-(1:250),]
# rho <- acf(beta_samp[,1], plot = FALSE)$acf[2]
# var(beta_samp[,1])*(1+rho)/(1-rho)
plot(beta_samp[,1], type = "l")
abline(h = postmean[1] +
         c(0, qnorm(c(0.025,0.975))*sqrt(postcov[1,1])),
       col = 2)
confint <- postmean[1] +
  qnorm(c(0.025,0.975))*
  sqrt(postcov[1,1])
mean(beta_samp[,1] > confint[2] | beta_samp[,1] < confint[1])
confint <- postmean[1] +
  qnorm(c(0.05,0.95))*
  sqrt(postcov[1,1])
mean(beta_samp[,1] > confint[2] | beta_samp[,1] < confint[1])

(colMeans(beta_samp) - postmean)/postmean

acf(beta_samp[,1])

# Ellipse or prior to posterior
plot(ellipse::ellipse(Omega[1:2,1:2]),
     xlim = c(-10, 10),
     ylim = c(-10, 10),
     type = "l")
lines(ellipse::ellipse(postcov[1:2,1:2],
                       centre = postmean[1:2]),
      col = 2)
# Sample from cospex model
nsites <- 100
pts <- cbind(runif(nsites + 1), runif(nsites + 1))
distm <- mev::distg(pts, scale = 1, rho = 0)
ds0 <- distm[1,-1]
x0 <- rexp(1)
a <- x0*exp(-ds0/2)
b <- x0^0.5
Sigma <- fields::Matern(d = distm[-1,-1], phi = 0.5)
X <- a + b*mvtnorm::rmvnorm(n = 1, sigma = Sigma)[1,] +
  rnorm(nrow(Sigma), sd = 0.1)
cens <- X < 0
X <- pmax(0, X)
logLikFn <- function(w){
  cospec_cond_loglik(
    x = X,
    censoring = cens,
    mu = a+b*w,
    tau = 0.1)
  }
# Plot log-likelihood surface
# val_x1 <- seq(-3,2, length.out = 1e3)
# val_x2 <- seq(-2,3, length.out = 1e3)
# logLikSurf <- matrix(NA,
#                      nrow = length(val_x1),
#                      ncol = length(val_x2))
# for(i in seq_along(val_x1)){
#   for(j in seq_along(val_x2)){
#     logLikSurf[i,j] <- logLikFn(c(val_x1[i], val_x2[j]))
#   }
# }
# fields::image.plot(val_x1, val_x2, logLikSurf)

p <- nrow(Sigma)
nrep <- 1e4L + 250
beta_samp <- matrix(NA, nrow = nrep, ncol = p)
beta_samp[1,] <- rep(0, p)
niter <- rep(0, nrep)
for(i in 2:nrep){
  slice <- elliptical_slice_sampling(
    curState = beta_samp[i-1,],
    logLikFn = logLikFn,
    normSamp = c(mvtnorm::rmvnorm(n = 1, sigma = Sigma)))
  beta_samp[i,] <- slice$x
  niter[i] <- slice$niter
}
beta_samp <- beta_samp[-(1:250),]
niter <- niter[-(1:250)]

plot(beta_samp[,which(cens)[1]],type = "l")
plot(beta_samp[,which(!cens)[1]],type = "l")
acf(beta_samp[,which(cens)[1]])
acf(beta_samp[,which(!cens)[1]])

}
