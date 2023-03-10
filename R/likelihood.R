#' Negative log likelihood of conditional spatial extremes
#'
#' This is just a wrapper around the Gaussian density
#' but it dispatches arguments of the location,
#' scale, variogram/covariance and marginals, in this order
#' @param pars [numeric] vector of parameters for \code{loc}, \code{locx}, \code{scale}, \code{cov}
#' @param data [matrix] \code{n} by \code{d} matrix of observations
#' @param x0 [numeric] \code{n} vector of values at the conditioning site
#' @param dist [matrix] \code{d} by \code{d} distance matrix
#' @param dist0 [numeric] \code{d} vector of distance to the conditioning site
#' @param locFn [\code{cospexFn}] type 1 location function. Default to \code{NULL}
#' @param locxFn [\code{cospexFn}] type 2 location function
#' @param scaleFn [\code{cospexFn}] scaling function
#' @param covFn [\code{covFn}] covariance function
#' @param ... additional arguments, currently ignored
#' @return the value of the negative log likelihood
#' @author Leo Belzile
#' @export
#'
#' @examples
#' n <- 100L
#' loc <- cbind(runif(11), runif(11))
#' d <- as.matrix(dist(loc))
#' dist0 <- d[1, -1]
#' dist <-  d[-1, -1]
#' data <-
#' rcospex(n = n,
#'     dist = dist,
#'     dist0 = dist0,
#'     covFn = covExpFn,
#'     locxFn = locxExpFn,
#'     scaleFn = scale2Fn,
#'     locxPars = c(0.5, 1),
#'     scalePars = 0.5,
#'     covPars = c(2, 0.5)
#'     )
#'  cospex_nll(
#'   pars = c(0.5, 1, 0.5, 2, 0.5),
#'   data = data[,-1],
#'   x0 = data[,1],
#'   dist = dist,
#'   dist0 = dist0,,
#'   covFn = covExpFn,
#'   locxFn = locxExpFn,
#'   scaleFn = scale2Fn
#'  )
cospex_nll <- function(pars,
                       data,
                       x0,
                       dist,
                       dist0,
                       locFn = NULL,
                       locxFn,
                       scaleFn,
                       covFn,
                       ...){
  stopifnot(isTRUE(is.matrix(data)))
  D <- ncol(data)
  n <- nrow(data)
  # Dispatch parameters to function
  # Keep in mind some can be fixed in the optimization
  # so rather than @npar, use is.fixed
  # Non-Gaussian marginals require Jacobian of transformation
  # Need distance matrix for covariance and
  # distance to conditioning site
  npar_locFn <- ifelse(is.null(locFn),
                       0,
                       sum(!is.fixed(locFn)))
  npar_locxFn <- sum(!is.fixed(locxFn))
  npar_scaleFn <- sum(!is.fixed(scaleFn))
  npar_covFn <- sum(!is.fixed(covFn))
  npars <- c(npar_locFn, npar_locxFn,
            npar_scaleFn, npar_covFn)
  cpars <- c(0, cumsum(npars))
  npar_total <- sum(npars)
  stopifnot(length(pars) == npar_total,
            is.numeric(pars))
    pars_loc <- pars[cpars[1] +
      seq_len(length.out = npar_locFn)]
  pars_locx <- pars[cpars[2] +
    seq_len(length.out = npar_locxFn)]
  pars_scale <- pars[cpars[3] +
    seq_len(length.out = npar_scaleFn)]
  pars_cov <- pars[cpars[4] +
                       seq_len(length.out = npar_covFn)]

  # The function a() - locxFn, returns typically a vector but we have outer x0 * loc2
  loc <- outer(x0,
                locxFn@fun(par = pars_locx,
                           dist = dist0),
                FUN = "*")
  stopifnot(isTRUE(is.matrix(loc)),
            isTRUE(all.equal(
              target = dim(loc),
              current = c(n, D),
              check.attributes = FALSE)))
  # The above computes the outer product
  # so is not memory friendly...
  if(!is.null(locFn)){
    # The function gamma (constant location field)
    # should return a d-vector
    loc_cst <- locFn@fun(par = pars_loc, dist = dist)
    stopifnot(length(loc_cst) == ncol(loc))
    loc <- loc +
      matrix(loc_cst,
             nrow = nrow(loc),
             ncol = ncol(loc))
  }
  # This returns either a vector
  # (if the scale is the same for all sites)
  # or a n x d matrix if the scale is
  # spatially varying
  scale <- scaleFn@fun(par = pars_scale,
                     x = x0,
                     dist = dist0)
  if(is.matrix(scale)){
    scale_sv <- TRUE
    stopifnot(isTRUE(is.matrix(scale)),
              all.equal(dim(scale), c(n, D),
                        check.attributes = FALSE))
  } else{
    scale_sv <- FALSE
    stopifnot(length(scale) == n)
  }
  covmat <- covFn@fun(par = pars_cov,
                      dist = dist) -
    tcrossprod(covFn@fun(par = pars_cov,
                         dist = as.numeric(dist0)))/
    covFn@fun(par = pars_cov,
              dist = 0)

  stopifnot(nrow(covmat) == D,
            ncol(covmat) == D)
 # Benchmarking shows this is less efficient in
 # small cases
  cholesky <- chol(covmat)
  # cholinv <- backsolve(cholesky, diag(nrow(cholinv)))
  # apply(data, 1, function(v){
  # sum(crossprod(v, cholinv)^2)})}

  ## Negative log likelihood
   0.5*n*D*log(2*pi) +
    ifelse(scale_sv, 1, D)*sum(log(scale)) +
    n*sum(log(diag(cholesky))) +
     # determinant is 2*diag(R), but factor 1/2
    0.5*sum(apply((data-loc)/scale, 1,
                  function(v){
                    crossprod(backsolve(cholesky, v, transpose = TRUE))}))
}
