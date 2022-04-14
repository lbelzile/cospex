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

# Location-scale family with a(s,x) + gamma(s) + b(s,x)Z0(s)
# where a(s, x) = alpha(s)*x such that alpha(s0)=1
# gamma(s0) = 0.
# For ergodicity, need b=1 as h = d(s, s0) -> oo
# Allow for fixed parameters for each submodel

setClass("ScalingFunction",
         slots = c(npar = "integer",   # number of parameters for optimization
                   default = "numeric",# default values for parameters / starting values for optimization
                   func = "function",  # affine normalization function
                   grad = "function",  # gradient function (same signature as func)
                   type = "character", # type of function, location or scale
                   lower = "numeric",  # lower bound
                   upper = "numeric",  # upper bound
                   cst = "logical",    # FALSE if spatially varying
                   fixed = "logical"), # TRUE if arguments are held fixed for optimization
         prototype = prototype(npar = 1L,
                               default = 1,
                               func = function(x){1},
                               grad = function(x){0},
                               type = "scale",
                               lower = -Inf,
                               upper = Inf,
                               cst = TRUE,
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
           if(length(object@type) != 1L || !object@type %in% c("scale","loc","locx")){
             return("\"type\" must be a character, either \"scale\",\"loc\" or \"locx\".")
           }
           if(length(object@cst) != 1L){
             return("\"cst\" must be a logical of length one.")
           }
           if(!isTRUE(all(object@lower < object@upper))){
             return("Lower bound \"lower\" and upper bound \"upper\" are not ordered: must provide \"lower\" < \"upper\".")
           }
           if(!isTRUE(all(object@lower <= object@default, object@upper >= object@default))){
             return("Some default values for parameters are out of bound: either \"default\" > \"upper\" or \"default\" < \"lower\".")
           }
           if(!isTRUE(all.equal(formalArgs(func(object)), formalArgs(grad(object))))){
             return("Arguments of the scaling function and its gradient do not match.")
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
#' @param func function of \code{npar} argument and distance
#' @param grad gradient of \code{func} with respect to parameters
#' @param cst logical; is the function constant. If \code{FALSE}, this indicates that \code{func} is a two-place function that returns
#' @param type integer indicating the type of function, either \code{scale}, \code{loc} or \code{locx} for locations (\code{locx}) gets multiplied by \code{x}, the value at the conditioning site
#' @param lower \code{npar} vector of lower bounds for the parameter vector
#' @param upper \code{npar} vector of upper bounds for the parameter vector
#' @param fixed \code{npar} logical vector indicating whether parameters are held fixed to their \code{default} values
ScalingFunction <-
  function(npar,
           default,
           func,
           grad,
           type = "scale",
           lower = rep(-Inf, npar),
           upper = rep(Inf, npar),
           cst = TRUE,
           fixed = rep(FALSE, npar)
           )
  {
  new("ScalingFunction",
      npar = npar,
      default = default,
      func = func,
      grad = grad,
      lower = lower,
      upper = upper,
      type = type,
      cst = TRUE,
      fixed = fixed)
}

# Define methods and 'getters'
setMethod("length", "ScalingFunction", function(x) x@npar)

setGeneric("default", function(x) standardGeneric("default"))
setMethod("default", "ScalingFunction", function(x) x@default)

setGeneric("grad", function(x) standardGeneric("grad"))
setMethod("grad", "ScalingFunction", function(x) x@grad)

setGeneric("func", function(x) standardGeneric("func"))
setMethod("func", "ScalingFunction", function(x) x@func)

setGeneric("lowerBound", function(x) standardGeneric("lowerBound"))
setMethod("lowerBound", "ScalingFunction", function(x) x@lower)

setGeneric("upperBound", function(x) standardGeneric("upperBound"))
setMethod("upperBound", "ScalingFunction", function(x) x@upper)

setGeneric("is.fixed", function(x) standardGeneric("is.fixed"))
setMethod("is.fixed", "ScalingFunction", function(x) x@fixed)

setMethod("show", "ScalingFunction",
          function(object){
            cat(class(object), "instance with", length(object), "parameters.", "\n")
            }
          )

# setters are bound to fail unless we replace bounds with other of the same length
# this mean we cannot modify the length of the parameter vector
setGeneric("lowerBound<-", function(x, value) standardGeneric("lowerBound<-"))
setReplaceMethod("lowerBound", "ScalingFunction", function(x, value) {x@lower <- value; validObject(x); x})

setGeneric("upperBound<-", function(x, value) standardGeneric("upperBound<-"))
setReplaceMethod("upperBound", "ScalingFunction", function(x, value) {x@upper <- value; validObject(x); x})

setGeneric("is.fixed<-", function(x, value) standardGeneric("is.fixed<-"))
setReplaceMethod("is.fixed", "ScalingFunction", function(x, value) {x@fixed <- value; validObject(x); x})

setGeneric("default<-", function(x, value) standardGeneric("default<-"))
setReplaceMethod("default", "ScalingFunction", function(x, value) {x@default <- value; validObject(x); x})


expCorFunc <- function(par, x, dist,...){
  stopifnot(isTRUE(all(dist > 0)),
            length(x) == 1L,
            is.numeric(x),
            x > 0,
            par[2] < 2,
            par[1] > 0,
            par[2] > 0
  )
 x*exp(-(dist/par[1])^par[2])
}
expCorGrad <- function(par, x, dist,...){
  stopifnot(isTRUE(all(dist > 0)),
            length(x) == 1L,
            is.numeric(x),
            x > 0,
            par[2] < 2,
            par[1] > 0,
            par[2] > 0
  )
  x*exp(-(dist/par[1])^par[2])*(dist/par[1])^par[2]
    c(par[2]/par[1], log(par[2]))
}
# Define scaling functions used in the various papers

locFunExp <- new("ScalingFunction",
                 npar = 2L,
                 default = c(1, 1),
                 func = expCorFunc,
                 grad = expCorGrad,
                 type = "locx",
                 lower = c(0, 0),
                 upper = c(Inf, 2),
                 cst = FALSE,
                 fixed = rep(FALSE, 2))
