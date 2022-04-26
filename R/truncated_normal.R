#' Univariate truncated normal sampler
#'
#' This implementations works for a < 37
#' or b > 37
#' @param n [integer] sample size
#' @param mean [numeric] scalar location parameter
#' @param sd [numeric] scalar scale parameter
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
  a_std <- (a - mean) / sd
  b_std <- (b - mean) / sd
  if(b_std > -37 | a_std < 37){
    stop("Interval requested is beyond numerical tolerance.\nUse \"TruncatedNormal\" package \"rtnorm\" instead for rare events simulation.")
  }
  if(a_std < 0){
    Fa <- pnorm(a_std)
    Fb <- pnorm(b_std)
    mean + sd * qnorm(Fa + runif(n) * (Fb - Fa))
  } else {
    Fa <- pnorm(a_std, lower.tail = FALSE)
    Fb <- pnorm(b_std, lower.tail = FALSE)
    mean + sd * (-qnorm(Fa - (Fa - Fb) * runif(n)))
  }
}

#' Density of a univariate truncated Normal
#' @inheritParams rtnorm
#' @param x [numeric] observation vector
#' @param log [logical] if \code{TRUE}, return the log density
#' @export
dtnorm <- function(x,
                   mean = 0,
                   sd = 1,
                   a = -Inf,
                   b = Inf,
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
