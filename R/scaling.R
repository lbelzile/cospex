# Scaling functions for the conditional extremes model
# Two place function: distance versus value at site s0
# Spatially varying (include potential for anisotropy)
# with
# 1) default values of parameters
# 2) parameter constraints
# 3) gradient with respect to each parameter
# 4) spatially varying (logical)
# 5) spatial anisotropy ? If outside, then problems with gradients...
#  need logical switch to deal with fixed distance (not recomputed)

# Location-scale family with
# a(s,x) + gamma(s) + b(s,x)Z0(s)
# where a(s, x) = alpha(s)*x such that alpha(s0)=1
# gamma(s0) = 0.
# For ergodicity, need b=1 as h = d(s, s0) -> oo
# Allow for fixed parameters for each submodel

setClass("cospexFn",
         slots = c(npar = "integer",   # number of parameters for optimization
                   default = "numeric",# default values for parameters / starting values for optimization
                   # pars = "ANY",
                   fun = "function",  # affine normalization function
                   grad = "ANY",  # gradient function (same signature as func)
                   type = "character", # type of function, location or scale
                   lower = "numeric",  # lower bound
                   upper = "numeric",  # upper bound
                   cst = "logical",    # FALSE if spatially varying
                   lag = "logical", # whether there is lag-asymptotic dependence, only for locx functions
                   fixed = "logical"), # TRUE if arguments are held fixed for optimization
         prototype = prototype(npar = 1L,
                               # pars = NULL,
                               default = 1,
                               fun = function(x){1},
                               grad = function(x){0},
                               type = "scale",
                               lower = -Inf,
                               upper = Inf,
                               cst = TRUE,
                               lag = FALSE,
                               fixed = TRUE),
         validity = function(object) {
           length_args <- c(
             length(object), #number of parameters
             length(default(object)),
             length(lowerBound(object)),
             length(upperBound(object)),
             length(is.fixed(object))
           )
           if(length(unique(length_args)) != 1L){
             return("Bounds and constraints must have the same length as the parameter vector.")
           }
           if(length(object@type) != 1L || !object@type %in% c("scale","loc","locx","cov")){
             return("\"type\" must be a character, either \"cov\", \"scale\",\"loc\" or \"locx\".")
           }
           if(length(object@cst) != 1L){
             return("\"cst\" must be a logical of length one.")
           }
           if(! (is.null(object@grad) | is.function(object@grad))){
             return("Optional argument \"grad\" must be a function or NULL.")
           }
           # if(!is.null(object@pars) | length(object@pars) == object@npar){
           #   return("\"pars\" must be NULL or a vector of length \"npar\".")
           # }
           if(!isTRUE(all(object@lower < object@upper))){
             return("Lower bound \"lower\" and upper bound \"upper\" are not ordered: must provide \"lower\" < \"upper\".")
           }
           if(!isTRUE(all(object@lower <= object@default, object@upper >= object@default))){
             return("Some default values for parameters are out of bound: either \"default\" > \"upper\" or \"default\" < \"lower\".")
           }
           if(!is.null(grad(object))){
           if(!isTRUE(all.equal(formalArgs(fun(object)), formalArgs(grad(object))))){
             return("Arguments of the scaling function and its gradient do not match.")
           }
           }
           TRUE
         }
)

# Initialization

#' Affine renormalization function
#'
#' This class is used to define a location or scaling function.
#'
#' @param npar number of parameters passed as first argument to \code{func} and \code{grad}
#' @param fun function of \code{npar} argument and distance
#' @param grad gradient of \code{func} with respect to parameters
#' @param cst logical; is the function constant. If \code{FALSE}, this indicates that \code{func} is a two-place function that returns
#' @param type integer indicating the type of function, either \code{scale}, \code{loc} or \code{locx} for locations (\code{locx}) gets multiplied by \code{x}, the value at the conditioning site
#' @param lower \code{npar} vector of lower bounds for the parameter vector
#' @param upper \code{npar} vector of upper bounds for the parameter vector
#' @param lag logical scalar indicating whether the location model is lag-asymptotic (\code{locx}) or else if the function is compatible with lag-asymptotic dependence (\code{scale}). Argument ignored for pure location components (\code{loc}).
#' @param fixed \code{npar} logical vector indicating whether parameters are held fixed to their \code{default} values
cospexFn <-
  function(npar,
           default,
           fun,
           grad = NULL,
           type = "scale",
           lower = rep(-Inf, npar),
           upper = rep(Inf, npar),
           cst = TRUE,
           lag = FALSE,
           fixed = rep(FALSE, npar)
           )
  {
  new("cospexFn",
      npar = npar,
      default = default,
      fun = fun,
      grad = grad,
      lower = lower,
      upper = upper,
      type = type,
      cst = cst,
      lag = lag,
      fixed = fixed)
}

# Define methods and 'getters'
setMethod("length", "cospexFn", function(x) x@npar)

setGeneric("default", function(x) standardGeneric("default"))
setMethod("default", "cospexFn", function(x) x@default)

setGeneric("grad", function(x) standardGeneric("grad"))
setMethod("grad", "cospexFn", function(x) x@grad)

setGeneric("fun", function(x) standardGeneric("fun"))
setMethod("fun", "cospexFn", function(x) x@fun)

setGeneric("lowerBound", function(x) standardGeneric("lowerBound"))
setMethod("lowerBound", "cospexFn", function(x) x@lower)

setGeneric("upperBound", function(x) standardGeneric("upperBound"))
setMethod("upperBound", "cospexFn", function(x) x@upper)

setGeneric("is.fixed", function(x) standardGeneric("is.fixed"))
setMethod("is.fixed", "cospexFn", function(x) x@fixed)

setGeneric("type", function(x) standardGeneric("type"))
setMethod("type", "cospexFn", function(x) x@type)


setGeneric("is.lagDep", function(x) standardGeneric("is.lagDep"))
setMethod("is.lagDep", "cospexFn", function(x){
  ifelse(x@type == "locx", x@lag, NULL) })

setMethod("show", "cospexFn",
          function(object){
            cat(class(object), "instance with", length(object), "parameters.", "\n")
            }
          )

# setters are bound to fail unless we replace bounds with other of the same length
# this mean we cannot modify the length of the parameter vector
setGeneric("lowerBound<-", function(x, value) standardGeneric("lowerBound<-"))
setReplaceMethod("lowerBound", "cospexFn", function(x, value) {x@lower <- value; validObject(x); x})

setGeneric("upperBound<-", function(x, value) standardGeneric("upperBound<-"))
setReplaceMethod("upperBound", "cospexFn", function(x, value) {x@upper <- value; validObject(x); x})

setGeneric("is.fixed<-", function(x, value) standardGeneric("is.fixed<-"))
setReplaceMethod("is.fixed", "cospexFn", function(x, value) {x@fixed <- value; validObject(x); x})

setGeneric("is.lagDep<-", function(x, value) standardGeneric("is.lagDep<-"))
setReplaceMethod("is.lagDep", "cospexFn", function(x, value) {x@lag <- value; validObject(x); x})


setGeneric("default<-", function(x, value) standardGeneric("default<-"))
setReplaceMethod("default", "cospexFn", function(x, value) {x@default <- value; validObject(x); x})

#' Exponential correlation function
#'
#' @param par [numeric] parameter vector, consisting of positive scale \eqn{\sigma>0} and smoothness \eqn{\nu \in (0,2)}.
#' @param dist [numeric] vector of \code{d} distance
#' @param ... additional parameters, currently ignored
#' @return a positive vector of length \code{d}
powExpCorFn <- function(par, dist,...){
  stopifnot(isTRUE(all(dist >= 0)),
            # length(x) == 1L,
            # is.numeric(x),
            # is.finite(x),
            # isTRUE(x > 0),
            par[2] <= 2,
            par[1] > 0,
            par[2] > 0
  )
 exp(-(dist/par[1])^par[2])
}

#' Gradient of the powered exponential correlation location function
#' @inheritParams powExpCorFunc
powExpCorGrad <- function(par, dist,...){
  stopifnot(isTRUE(all(dist >= 0)),
            # length(x) == 1L,
            # is.numeric(x),
            # x > 0,
            par[2] <= 2,
            par[1] > 0,
            par[2] > 0
  )
  exp(-(dist/par[1])^par[2])*(dist/par[1])^par[2] *
    c(par[2]/par[1], log(par[2]))
}




# Define scaling functions used in the various papers

#' Powered exponential correlation function
#'
#' This model defines asymptotically independent processes.
#' @export
locxExpFn <- new("cospexFn",
                 npar = 2L,
                 default = c(1, 1),
                 fun = powExpCorFn,
                 grad = powExpCorGrad,
                 type = "locx",
                 lower = c(0, 0),
                 upper = c(Inf, 2),
                 cst = FALSE,
                 lag = FALSE,
                 fixed = rep(FALSE, 2))

#' @export
locxCstFn <- new("cospexFn",
                 npar = 0L,
                 default = numeric(),
                 fun = identity,
                 grad = function(x){0},
                 type = "locx",
                 lower = numeric(),
                 upper = numeric(),
                 cst = TRUE,
                 lag = FALSE,
                 fixed = logical())



#' Scaling function for lag-dependent data
#'
#' @param par [numeric] vector of length two with scale parameter \eqn{lambda \ge 0} and shape parameter \eqn{beta < 0}
#' @param x [numeric] scalar; positive value at the conditioning site
#' @param ... additional arguments, currently ignored
s1Fn <- function(par, x, ...){
  stopifnot(par[1] >= 0,
            par[2] < 0,
            length(x) == 1L,
            is.numeric(x),
            is.finite(x),
            isTRUE(x > 0)
            )
1/(1+par[1]*x^par[2])
}

#' Gradient of first scaling function
#' @inheritParams s1Func
s1Grad <- function(par, x, ...){
  stopifnot(length(par) == 2L,
            par[1] >= 0,
            par[2] < 0,
            length(x) == 1L,
            is.numeric(x),
            is.finite(x),
            isTRUE(x > 0)
  )
  c(-par[1]*x^par[2]*log(x)/(par[1]*x^par[2] + 1)^2,
    -x^par[2]/(par[1]*x^par[2] + 1)^2)
}

#' @export
scale1Fn <- new("cospexFn",
                    npar = 2L,
                    default = c(1, 0),
                    fun = s1Fn,
                    grad = s1Grad,
                    type = "scale",
                    lower = c(0, -Inf),
                    upper = c(Inf, 0),
                    cst = FALSE,
                    lag = TRUE,
                    fixed = rep(FALSE, 2))


#' Second scaling function
#'
#' @param par [numeric] scalar, shape parameter \eqn{0 \le \beta < 0}
#' @param x [numeric] scalar; positive value at the conditioning site
#' @param ... additional arguments, currently ignored
s2Fn <- function(par, x, ...){
  stopifnot(length(par) == 1L,
            par[1] >= 0,
            par[1] < 1,
            is.numeric(x),
            is.finite(x),
            isTRUE(all(x > 0))
  )
  x^par[1]
}

#' Gradient of second scaling function
#' @inheritParams s2Func
s2Grad <- function(par, x, ...){
  stopifnot(length(par) == 1L,
            par[1] >= 0,
            par[1] < 1,
            is.numeric(x),
            is.finite(x),
            isTRUE(all(x > 0))
  )
  x^par*log(x)
}

#' @export
scale2Fn <- new("cospexFn",
                 npar = 1L,
                 default = 0,
                 fun = s2Fn,
                 grad = s2Grad,
                 type = "scale",
                 lower = 0,
                 upper = 1,
                 cst = FALSE,
                 lag = FALSE,
                 fixed = FALSE)

# Define covariance models
#' @export
covMaternFn <- new(
  "cospexFn",
    npar = 3L,
    default = c(1, 1, 0.5),
    fun = cov.Matern,
    grad = NULL,
    type = "cov",
    lower = rep(1e-10, 3L),
    upper = rep(Inf, 3L),
    cst = FALSE,
    lag = FALSE,
    fixed = rep(FALSE, 3))

#' @export
covExpFn <- new(
  "cospexFn",
  npar = 2L,
  default = c(1, 1),
  fun = cov.exp,
  grad = NULL,
  type = "cov",
  lower = rep(1e-10, 2L),
  upper = rep(Inf, 2L),
  cst = FALSE,
  lag = FALSE,
  fixed = rep(FALSE, 2L))

#' @export
covCauchyFn <- new(
  "cospexFn",
  npar = 4L,
  default = rep(1, 4L),
  fun = cov.Cauchy,
  grad = NULL,
  type = "cov",
  lower = rep(1e-10, 4),
  upper = c(Inf, Inf, 2, Inf),
  cst = FALSE,
  lag = FALSE,
  fixed = rep(FALSE, 4L))

#' @export
covPowExpFn <- new(
  "cospexFn",
  npar = 3L,
  default = rep(1, 3L),
  fun = cov.powerexp,
  grad = NULL,
  type = "cov",
  lower = rep(1e-10, 3L),
  upper = c(Inf, Inf, 2),
  cst = FALSE,
  lag = FALSE,
  fixed = rep(FALSE, 3L))
