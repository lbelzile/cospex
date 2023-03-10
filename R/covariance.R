# Covariance models
# The exact forms are extracted from
# the Handbook of spatial statistics,
# except for the Matérn which is parametrized
# in terms

#' Matérn covariance matrix
#'
#' @param pars [numeric] parameter vector, parametrized in terms of variance, scale and shape
#' @param dist [matrix] square matrix of distance
#' @param dmat [logical] if \code{TRUE}, a matrix of pairwise distance with zero diagonal entries
#' @param check_dist [logical] if \code{TRUE}; check that entries are positives and matrix symmetrical
#' @param ... additional arguments, currently ignored
#' @export
cov.Matern <-
  function (pars,
            dist,
            dmat = TRUE,
            check_dist = FALSE,
            ...)
{
    stopifnot(isTRUE(all(pars > 0)),
              length(pars) == 3L)
    if(isTRUE(dmat)){
      if(is.vector(dist)){
        dmat <- FALSE
      }
      if(is.matrix(dist)){
        if(nrow(dist) != ncol(dist)){
          dmat <- FALSE
        }
      }
    }
    if(isTRUE(check_dist) & dmat){
      stopifnot(isSymmetric.matrix(dist),
                nrow(dist) == ncol(dist),
                isTRUE(all(diag(dist) == 0))
      )
    }
    scale <- pars[1]
    rrange <- pars[2] #kappa
    shape <- pars[3] #nu
    if(isTRUE(dmat)){
      d <- dist[lower.tri(dist, diag = FALSE)] / rrange
    } else{
      d <- dist / rrange
    }
  if (isTRUE(all.equal(shape, 0.5,
                       check.attributes = FALSE))) {
    covE <- scale * exp(-d)
  } else if (isTRUE(all.equal(shape, 1.5,
                              check.attributes = FALSE))) {
    covE <- scale * (1 + d) * exp(-d)
  } else if(shape > 1e10){
    covE <- scale*exp(-0.5*d^2/rrange^2)
    } else{
    matern <- function(d, sigma, nu){
      ifelse(d == 0, sigma, exp((1-nu)*log(2) - lgamma(nu) + nu*log(d) + log(besselK(d, nu = nu))))
    }
   covE <- matern(d = d,
                  sigma = scale,
                  nu = shape)
    }
  if(dmat){
    covM <- matrix(scale,
                   nrow = nrow(dist),
                   ncol = ncol(dist))
    covM[lower.tri(covM, diag = FALSE)] <- covE
    covM[upper.tri(covM, diag = FALSE)] <- t(covM)[upper.tri(covM, diag = FALSE)]
    return(covM)
  } else{
    return(covE)
  }
  }

#' Exponential covariance function
#'
#' @pars [numeric] vector of parameters, for variance and scale
#' @pars dist [numeric] distance vector or matrix
#' @pars ... additional arguments, currently ignored
#' @export
cov.exp <- function(pars,
                    dist,
                    ...){
  stopifnot(length(pars) == 2L)
  cov.Matern(pars = c(pars, 0.5),
             dist = dist, ...)
}

#' Powered exponential covariance
#'
#' @param pars [numeric] vector of variance, scale and shape
#' @param dist matrix or vector of distance
#' @export
cov.powerexp <- function(pars, dist, ...){
  stopifnot(length(pars) == 3L,
            isTRUE(all(pars > 0)),
            pars[3] <= 2)
  pars[1]*exp(-(dist/pars[2])^pars[3])
}

#' Cauchy covariance model
#'
#' The Cauchy covariance model is of the form
#' \deqn{c(d) = \sigma^2 \left\{1+\left(\frac{d}{\theta}\right)^\alpha\right\}^{-\beta/\alpha}}
#' where the vector of parameters consists of the variance \eqn{\sigma^2}, scale \eqn{\theta}, shape1 \eqn{\alpha} and shape2 \eqn{beta}, in this order
#'
#' @param pars [numeric] vector of parameters of length 4
#' @export
cov.Cauchy <- function(pars, dist, ...){
  stopifnot(length(pars) == 4L,
            isTRUE(all(pars > 0)),
            pars[3] <= 2)
  pars[1]*(1+(dist/pars[2])^pars[3])^(-pars[4]/pars[3])
}
