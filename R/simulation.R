#' Simulate from conditional spatial extremes model
#'
#' @param n [integer] number of replicates
#' @param dist [matrix] \code{d} by \code{d} matrix of pairwise distance between sites
#' @param dist0 [numeric] \code{d} vector of distance to conditioning site
#' @param covFn [cospexFn] covariance function for the residual Gaussian process
#' @param locFn [cospexFn] location function
#' @param locxFn [cospexFn] location function
#' @param scaleFn [cospexFn] scale function
#' @param nuggetVar [numeric] nugget variance; default to zero
#' @param return an \code{n} by \code{d+1} matrix of replications whose first column is the value at the conditioning site
#'
#' @examples
#' n <- 100L
#' loc <- cbind(runif(11), runif(11))
#' d <- as.matrix(dist(loc))
#' dist0 <- d[1, -1]
#' dist <-  d[-1, -1]
#' rcospex(n = n,
#'     dist = dist,
#'     dist0 = dist0,
#'     covFn = covExpFn,
#'     locxFn = locxExpFn,
#'     scaleFn = scale2Fn,
#'     locxPars = c(0.5, 1),
#'     scalePars = 0.5,
#'     covPars = c(2, 0.5),
#'     nuggetVar = sqrt(0.05)
#'     )
rcospex <- function(n,
                    dist,
                    dist0,
                    locFn = NULL,
                    locxFn,
                    scaleFn,
                    covFn,
                    locPars,
                    locxPars,
                    scalePars,
                    covPars,
                    nuggetVar = 0,
                    threshold = 0) {
  n <- as.integer(n)
  stopifnot("Missing \"dist\" argument." = !missing(dist))
  dist <- try(as.matrix(dist))
  stopifnot("Invalid \"dist\" argument." = !inherits(dist, "try-catch"),
    "Number of sites must be positive." = n > 0,
    "\"dist\" must be a matrix." = is.matrix(dist),
    "\"dist\" must be symmetric" = isSymmetric(dist),
    "\"dist0\" must be positive." = isTRUE(all(dist0 > 0)),
    "Invalid location function" = inherits(locxFn, "cospexFn"),
    "Invalid scale function" = inherits(scaleFn, "cospexFn"),
    "Invalid covariance function" = inherits(covFn, "cospexFn"),
    "\"nuggetVar\" should be a scalar" = length(nuggetVar) == 1L,
    "\"nuggetVar\" should be non-negative." = isTRUE(nuggetVar >= 0),
    "\"threshold\" should be a scalar" = length(threshold) == 1L,
    "\"threshold\" should be non-negative." = isTRUE(threshold >= 0)
  )

stopifnot(
  "Invalid type for location function" =  type(locxFn) == "locx",
  "Invalid type for scale function." = type(scaleFn) == "scale",
    "Invalid type for covariance function" = type(covFn) == "cov")
  if(!is.null(locFn)){
    stopifnot("Invalid constant location function" =  inherits(locFn, "cospexFn"))
    stopifnot("Invalid type for constant location function" =  type(locFn) == "loc")
  }

  d <- nrow(dist)
  stopifnot("Distance \"dist0\" argument missing" = !missing(dist0),
            "Length of \"dist0\" and size of \"dist\" do not match." = isTRUE(length(dist0) == d))
  # Simulate exceedance at conditioning site
  x0 <- rexp(n) + threshold
 # Conditional covariance matrix given Z(s0)=0
  Sigma <- covFn@fun(par = covPars,
           dist = dist) -
    tcrossprod(covFn@fun(par = covPars,
             dist = as.numeric(dist0)))/
    covFn@fun(par = covPars,
             dist = 0)
  # Simulate Gaussian residual process
  Z0 <- matrix(rnorm(n = n * nrow(dist)),
               ncol = d,
               nrow = n) %*%
    chol(Sigma)
  loc <- outer(x0, locxFn@fun(par = locxPars,
                         dist = dist0),
               FUN = "*")
  if (!is.null(locFn)) {
    loc <- loc + locFn@fun(par = locPars,
                           dist = dist)
  }
  sim <- loc + scaleFn@fun(par = scalePars, x = x0) * Z0
  if (nuggetVar > 0) {
   # add column by column
    sim <- sim + rnorm(n, sd = sqrt(nuggetVar))
  }
  return(cbind(x0, data))
}
