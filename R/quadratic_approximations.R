

#' Univariate updates based on Laplace approximation
#'
#' Following Rue and Held (2004), we perform
#' a random walk based on a normal approximation
#' at the curr parameter value.
#'
#' The algorithm uses a random walk proposal
#' and a Metropolis acceptance step.
#'
#' @param cur [double] curr value of the parameter
#' @param par_name [string] the name of the argument for the function
#' @param loglik [function] log likelihood function
#' @param loglik_gradient [function] derivative of the log likelihood with respect to parameter of interest
#' @param loglik_hessian [function] second derivative of the log likelihood with respect to the parameter of interest
#' @param logprior [function] log prior of parameter
#' @param logprior_gradient [function] derivative of the log prior with respect to parameter of interest
#' @param logprior_hessian [function] second derivative of the log prior with respect to the parameter of interest
#' @param lb [scalar] lower bound for the parameter
#' @param ub [scalar] upper bound for the parameter
#' @param damping [scalar] contraction factor for the Newton update
#' @param ... additional arguments passed to the log likelihood function and its derivatives
#' @return a new value for the parameter of interest
#' @export
#' @references Rue, H. and L. Held (2005), Gaussian Markov random fields, CRC press, section 4.4.1
univariate_laplace_update <-
  function(par_curr,
           par_name,
           loglik,
           loglik_grad,
           loglik_hessian,
           logprior,
           logprior_grad,
           logprior_hessian,
           lb = -Inf,
           ub = Inf,
           damping = 1,
           ...){
    ellipsis <- list(...)
    args <- formalArgs(loglik)
    stopifnot(length(par_curr) == 1L,
              is.finite(par_curr),
              length(lb) == 1,
              length(ub) == 1,
              lb < ub,
              length(damping) == 1L,
              damping > 0,
              damping <= 1,
              isTRUE(all(args == formalArgs(loglik_grad))),
              isTRUE(all(args == formalArgs(loglik_hessian)))
    )
    # Parameter value to be returned
    par_new <- par_curr # if all things fail
    # Copy arguments into a list to call function
    loglik_args <- ellipsis[names(ellipsis) %in% args]
    # Override value of parameter
    loglik_args[[par_name]] <- par_curr
    # TODO determine whether we allow users to pass
    # other arguments to logprior
    logpost_grad_curr <-
      do.call(what = loglik_grad,
              args = loglik_args) +
      logprior_grad(par_curr)
    logpost_hessian_curr <-
      do.call(what = loglik_hessian,
              args = loglik_args) +
      logprior_hessian(par_curr)
    mean_curr <- par_curr -
      damping * logpost_grad_curr/logpost_hessian_curr
    precision_curr <- -logpost_hessian_curr
    if(!isTRUE(precision_curr > 0)){
      return(par_new)
    }
    par_prop <- rtnorm(n = 1,
                       a = lb,
                       b = ub,
                       mean = mean_curr,
                       sd = sqrt(1/precision_curr))
    logpost_curr <-
      do.call(what = loglik,
              args = loglik_args) +
      logprior(par_curr)
    # Reverse move
    loglik_args[[par_name]] <- par_prop
    # TODO determine whether we allow users to pass
    # other arguments to logprior
    logpost_grad_prop <-
      do.call(what = loglik_grad,
              args = loglik_args) +
      logprior_grad(par_prop)
    logpost_hessian_prop <-
      do.call(what = loglik_hessian,
              args = loglik_args) +
      logprior_hessian(par_prop)
    mean_prop <- par_prop -
      damping * logpost_grad_prop/logpost_hessian_prop
    precision_prop <- -logpost_hessian_prop
    if(!isTRUE(precision_prop > 0)){
      return(par_new)
    }
    logpost_prop <-
      do.call(what = loglik,
              args = loglik_args) +
      logprior(par_prop)
    log_MH_ratio <-
      logpost_prop - logpost_curr +
      dtnorm(par_curr,
             a = lb,
             b = ub,
             mean = mean_prop,
             sd = sqrt(1/precision_prop),
             log = TRUE) -
      dtnorm(par_prop,
             a = lb,
             b = ub,
             mean = mean_curr,
             sd = sqrt(1/precision_curr),
             log = TRUE)
    if(log_MH_ratio > log(runif(1))){
      par_new <- par_prop
    }
    return(par_new)
  }


#' Multivariate updates based on Laplace approximation
#'
#' Following Rue and Held (2004), we perform
#' a random walk based on a multivariate normal approximation
#' at the \code{cur} parameter value.
#'
#' @param cur [double] current value of the vector parameter
#' @param par_name [string] the name of the argument for the function
#' @param loglik [function] log likelihood function
#' @param loglik_gradient [function] derivative of the log likelihood with respect to parameter of interest
#' @param loglik_hessian [function] hessian matrix, the second derivative of the log likelihood with respect to the parameter of interest
#' @param logprior [function] log prior of parameter
#' @param logprior_gradient [function] derivative of the log prior with respect to parameter of interest
#' @param logprior_hessian [function] second derivative of the log prior with respect to the parameter of interest
#' @param lb [double] vector of lower bounds for the parameters
#' @param ub [double] vector of upper bounds for the parameters
#' @param damping [double] double or vector of contraction factor for the Newton update
#' @param ... additional arguments passed to the log likelihood function and its derivatives
#' @return a new vector for the parameter of interest
#' @export
#' @references Rue, H. and L. Held (2005), Gaussian Markov random fields, CRC press, section 4.4.2
multivariate_laplace_update <-
  function(par_curr,
           par_name,
           loglik,
           loglik_grad,
           loglik_hessian,
           logprior,
           logprior_grad,
           logprior_hessian,
           lb = -Inf,
           ub = Inf,
           damping = 1,
           ...){
    ellipsis <- list(...)
    args <- formalArgs(loglik)
    stopifnot(length(par_curr) == 1L,
              is.finite(par_curr),
              length(lb) == 1,
              length(ub) == 1,
              lb < ub,
              length(damping) == 1L,
              damping > 0,
              damping <= 1,
              isTRUE(all(args == formalArgs(loglik_grad))),
              isTRUE(all(args == formalArgs(loglik_hessian)))
    )
    # Parameter value to be returned
    par_new <- par_curr # if all things fail
    # Copy arguments into a list to call function
    loglik_args <- ellipsis[names(ellipsis) %in% args]
    # Override value of parameter
    loglik_args[[par_name]] <- par_curr
    # TODO determine whether we allow users to pass
    # other arguments to logprior
    logpost_grad_curr <-
      do.call(what = loglik_grad,
              args = loglik_args) +
      logprior_grad(par_curr)
    logpost_hessian_curr <-
      do.call(what = loglik_hessian,
              args = loglik_args) +
      logprior_hessian(par_curr)
    eigen_precision <- eigen(-logpost_hessian_curr)
    if(!isTRUE(all(eigen_precision$values > 0))){
      return(par_new)
    }
    covar_curr <- tcrossprod(eigen_precision$vectors %*% diag(1/sqrt(eigen_precision$values)))
    mean_curr <- par_curr -
      damping * logpost_grad_curr %*% covar_curr
    precision_curr <- -logpost_hessian_curr
    par_prop <- TruncatedNormal::rtmvnorm(
      n = 1,
      a = lb,
      b = ub,
      mu = mean_curr,
      sigma = covar_curr)
    logpost_curr <-
      do.call(what = loglik,
              args = loglik_args) +
      logprior(par_curr)
    # Reverse move
    loglik_args[[par_name]] <- par_prop
    # TODO determine whether we allow users to pass
    # other arguments to logprior
    logpost_grad_prop <-
      do.call(what = loglik_grad,
              args = loglik_args) +
      logprior_grad(par_prop)
    logpost_hessian_prop <-
      do.call(what = loglik_hessian,
              args = loglik_args) +
      logprior_hessian(par_prop)
    eigen_precision <- eigen(-logpost_hessian_prop)
    if(!isTRUE(all(eigen_precision$values > 0))){
      return(par_new)
    }
    covar_prop <- tcrossprod(eigen_precision$vectors %*% diag(1/sqrt(eigen_precision$values)))
    mean_curr <- as.numeric(par_prop -
      damping * logpost_grad_prop %*% covar_prop)
    logpost_prop <-
      do.call(what = loglik,
              args = loglik_args) +
      logprior(par_prop)
    log_MH_ratio <-
      logpost_prop - logpost_curr +
      TruncatedNormal::dtmvnorm(par_curr,
             a = lb,
             b = ub,
             mu = mean_prop,
             sigma = covar_prop,
             log = TRUE) -
      TruncatedNormal::dtmvnorm(par_prop,
             a = lb,
             b = ub,
             mu = mean_curr,
             sigma = covar_curr,
             log = TRUE)
    if(log_MH_ratio > log(runif(1))){
      par_new <- par_prop
    }
    return(par_new)
  }
