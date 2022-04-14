#' Conditional censored likelihood
#'
#' Censored likelihood for observations from the conditional
#' spatial extremes model, given the imputed random effects
#' Each spatio-temporal observation is conditionally independent
#' of one another, so the likelihood is formed from the
#' product of univariate normal density and distribution functions
#'
#' @param x [numeric] vector of observations
#' @param censoring [logical] if \code{TRUE}, observation is left-censored, observed otherwise
#' @param mu [numeric] mean vector
#' @param tau [numeric] standard deviation vector (or scalar)
#' @export
cospec_cond_loglik <-
  function(x,
           censoring,
           mu,
           tau){
  stopifnot(length(x) == length(censoring),
            length(mu) == length(x),
            length(tau) %in% c(1L, length(mu)),
            is.logical(censoring),
            isTRUE(all(tau > 0))
            )
  sum(ifelse(censoring,
         pnorm(q = x,
               mean = mu,
               sd = tau,
               log.p = TRUE),
         dnorm(x = x,
               mean = mu,
               sd = tau,
               log = TRUE)
         ))
}

#' Gradient of conditional censored log likelihood
#'
#' Gradient with respect to scale \eqn{\tau}.
#'
#' @inheritParams cospec_cond_loglik
#' @export
cospec_cond_loglik_grad_tau <-
  function(x,
           censoring,
           mu,
           tau){
    stopifnot(length(x) == length(censoring),
              length(mu) == length(x),
              length(tau) %in% c(1L, length(mu)),
              is.logical(censoring)
    )
   # Centered inputs
   xc <- (x - mu) / tau
   sum(ifelse(censoring,
           -xc/tau*
             exp(dnorm(x = xc, log = TRUE) -
                 pnorm(q = xc, log.p = TRUE)),
           -1/tau + xc^2/tau
    ))
  }

#' Hessian of conditional censored log likelihood
#'
#' Hessian of conditional censored log likelihood
#'  with respect to scale \eqn{\tau}.
#'
#' @inheritParams cospec_cond_loglik
#' @export
cospec_cond_loglik_hessian_tau <-
  function(x,
           censoring,
           mu,
           tau){
    stopifnot(length(x) == length(censoring),
              length(mu) == length(x),
              length(tau) %in% c(1L, length(mu)),
              is.logical(censoring)
    )
    # Centered inputs
    xc <- (x - mu) / tau
    rdc <- exp(dnorm(x = xc, log = TRUE) -
                pnorm(q = xc, log.p = TRUE))
    sum(ifelse(censoring,
           2*xc/tau^2*rdc - xc^2/tau^2*(xc*rdc + rdc^2),
           1/tau^2 -3*xc^2/tau^2
    ))
  }


#' Univariate truncated normal sampler
#'
#' This implementations works for a < 37
#' or b > 37
#' @param n [integer] sample size
#' @param mean [numeric] scalar mean
#' @param sd [numeric] standard deviation
#' @param a [numeric] lower bound
#' @param b [numeric] upper bound
#' @export
rtnorm <- function(n = 1,
                   mean = 0,
                   sd = 1,
                   a = -Inf,
                   b = Inf
                   ){
  stopifnot(length(a) == 1L,
            length(b) == 1L,
            length(mean) == 1L,
            length(sd) == 1L,
            isTRUE(a < b),
            isTRUE(sd > 0))
  if(a < 0){
    Fa <- pnorm((a - mean) / sd)
    Fb <- pnorm((b - mean) / sd)
    mean + sd * qnorm(Fa + runif(n) * (Fb - Fa))
  } else {
    Fa <- pnorm((a - mean) / sd, lower.tail = FALSE)
    Fb <- pnorm((b - mean) / sd, lower.tail = FALSE)
    mean + sd * (-qnorm(Fa - (Fa - Fb) * runif(n)))
  }
}

#' Density of a univariate truncated Normal
#' @inheritParams rtnorm
#' @param x [numeric] observation vector
#' @param log [logical] if \code{TRUE}, return the log density
#' @export
dtnorm <- function(x,
                   a = -Inf,
                   b = Inf,
                   mean = 0,
                   sd = 1,
                   log = FALSE){
  stopifnot(length(a) == 1L,
            length(b) == 1L,
            length(mean) == 1L,
            length(sd) == 1L,
            isTRUE(a < b),
            isTRUE(sd > 0))
  dens <- dnorm(x, mean = mean, sd = sd, log = TRUE)
  if(is.finite(a) & is.finite(b)){
  dens <- dens  - log(
    pnorm(q = b, mean = mean, sd = sd) -
    pnorm(q = a, mean = mean, sd = sd))
  } else if(is.infinite(a) & is.finite(b)){
  dens <- dens - pnorm(q = b, mean = mean, sd = sd, log.p = TRUE)
  } else if(is.finite(a) & is.infinite(b)){
    dens <- dens - pnorm(q = a, mean = mean, sd = sd,
                         log.p = TRUE, lower.tail = FALSE)
  }
    return(dens)
}


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
