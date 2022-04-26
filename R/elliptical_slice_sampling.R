#' Elliptical slice sampling
#'
#' This function implements the elliptical slice sampler of
#' Murray, Adams and MacKay (2012) given a random sample from a Gaussian distribution with mean zero and covariance matrix. Sampling is performed outside of the function to allow users to leverage existing functions for sparse matrices.
#'
#' @param curState current value of the random state vector
#' @param logLikFn log likelihood function for the model, whose first arguemnt is \code{x}, corresponding to \code{curState}
#' @param normSamp Gaussian random vector
#' @param ... additional arguments passed to \code{logLikFn}
#' @return a list with elements
#' \itemize{
#' \item{x}{vector of \code{length(curState)}}
#' \item{niter}{number of iterations}
#' }
#' @export
elliptical_slice_sampling <-
  function(curState, logLikFn, normSamp, ...){
    if(length(curState) != length(normSamp)){
      stop("The sample \"normSamp\" and the current value \"curState\" are not of the same length.")
    }
    # Check that it's a function and signature
    stopifnot(is.function(logLikFn))
    fnArgs <- formalArgs(logLikFn)
    ellipsis <- list(...)
    if(!isTRUE(all(names(ellipsis) %in% fnArgs[-1]))){
      stop("Some arguments passed to the function not called by \"logLikFn\".")
    }
  logy <- do.call(what = logLikFn,
                  args = c(ellipsis,
                           `names<-`(list(curState), fnArgs[1]))) + log(runif(1))
  theta <- runif(1, 0, 2*pi)
  thetaBounds <- c(theta-2*pi, theta)
  output <- FALSE
  niter <- 0L
  while(!output){
    niter <- niter + 1L
    # Proposal vector
    propState <- curState*cos(theta) + normSamp*sin(theta)
    logLikProp <-
      do.call(what = logLikFn,
              args = c(ellipsis,
              `names<-`(list(propState), fnArgs[1])))
    if(isTRUE(logLikProp > logy)){
      output <- TRUE
    } else{
      if(theta < 0){
        thetaBounds[1] <- theta
      } else{
        thetaBounds[2] <- theta
      }
      theta <- runif(n = 1,
                     min = thetaBounds[1],
                     max = thetaBounds[2])
    }
  }
  return(list(x = propState, niter = niter))
}
