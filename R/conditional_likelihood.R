#' Conditional censored likelihood
#'
#' Censored likelihood for observations from the conditional
#' spatial extremes model, given the imputed random effects
#' Each spatio-temporal observation is conditionally independent
#' of one another, so the likelihood is formed from the
#' product of univariate normal density and distribution functions
#'
#' @param x [numeric] vector of observations
#' @param censoring [logical] vector; if \code{TRUE}, observation is left-censored, observed otherwise
#' @param mu [numeric] mean vector
#' @param tau [numeric] standard deviation vector (or scalar)
#' @export
cospex_cond_loglik <-
  function(x,
           censoring,
           mu,
           tau){
  stopifnot(is.vector(x),
            is.vector(censoring),
            is.vector(mu),
            is.vector(tau),
            length(x) == length(censoring),
            length(mu) == length(x),
            length(tau) == 1L,
            is.logical(censoring),
            isTRUE(tau>0)
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
#' @inheritParams cospex_cond_loglik
#' @export
cospex_cond_loglik_grad_tau <-
  function(x,
           censoring,
           mu,
           tau){
    stopifnot(length(x) == length(censoring),
              length(mu) == length(x),
              length(tau) == 1L,
              isTRUE(tau > 0),
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
#' @inheritParams cospex_cond_loglik
#' @export
cospex_cond_loglik_hessian_tau <-
  function(x,
           censoring,
           mu,
           tau){
    stopifnot(length(x) == length(censoring),
              length(mu) == length(x),
              length(tau) == 1L,
              isTRUE(tau > 0),
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


#' Gradient of conditional censored log likelihood
#'
#' Gradient with respect to location parameter.
#'
#' @inheritParams cospex_cond_loglik
#' @export
cospex_cond_loglik_grad_mu <-
  function(x,
           censoring,
           mu,
           tau){
    stopifnot(length(x) == length(censoring),
              length(mu) == length(x),
              length(tau) == 1L,
              isTRUE(tau > 0),
              is.logical(censoring)
    )
    # Centered inputs
    xc <- (x - mu) / tau
    ifelse(censoring,
               -1/tau*
                 exp(dnorm(x = xc, log = TRUE) -
                       pnorm(q = xc, log.p = TRUE)),
               xc/tau
    )
  }


#' Hessian of conditional censored log likelihood
#'
#' Hessian with respect to location parameter.
#' @return a list with two vectors, to be used with the chain rule
#' @inheritParams cospex_cond_loglik
#' @export
cospex_cond_loglik_hessian_mu <-
  function(x,
           censoring,
           mu,
           tau){
    stopifnot(length(x) == length(censoring),
              length(mu) == length(x),
              length(tau) == 1L,
              isTRUE(tau > 0),
              is.logical(censoring)
    )
    # Centered inputs
    xc <- (x - mu) / tau
    mfact <- ifelse(censoring,
           exp(dnorm(x = xc, log = TRUE) -
                   pnorm(q = xc, log.p = TRUE)), 0)
    mat1 <- ifelse(censoring,
                   mfact/tau^2*(xc+mfact),
                   1/tau^2)
    mat2 <- ifelse(censoring,
                   mfact/tau,
                   -xc/tau)
    list(pmixed = -mat1,
         psecond = -mat2)
  }
