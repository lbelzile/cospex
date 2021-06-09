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

setClass("ScalingFunction", contains = "VIRTUAL",
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
             length(isFixed(object))
           )
           if(length(unique(length_args)) != 1L){
             return("Bounds and constraints must have the same length as the parameter vector.")
           }
           if(length(type) != 1L || !type %in% c("scale","loc","locx")){
             return("`type` must be a character, either `scale`,`loc` or `locx`.")
           }
           if(length(cst) != 1L){
             return("`cst` must be a logical of length one.")
           }
           if(!isTRUE(all(lower < upper))){
             return("Lower bound `lower` and upper bound `upper` are not ordered: must provide `lower` < `upper`.")
           }
           if(!isTRUE(all(lower <= default, upper >= default))){
             return("Some default values for parameters are out of bound: either `default` > `upper` or `default` < `lower`.")
           }
           if(!isTRUE(all.equal(functionBody(func(object)), functionBody(grad(object))))){
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
setMethod("is.fixed", "ScalingFunction", function(x) x@isFixed)

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

expCor <- function(par, x, dist,...){
  stopifnot(isTRUE(all(dist > 0)),
            length(x) == 1L,
            is.numeric(x)
            x > 0,
            par[2] < 2,
            par[1] > 0,
            par[2] > 0
  )
 x*exp(-(dist/par[1])^par[2])
}

# Define scaling functions used in the various papers
#
locFunExp <- ScalingFunction(npar = 2,
                             func = expCor,
                             default = c(1, 1),
                             lower = c(0, 0),
                             upper = c(Inf, 2),
                             cst = FALSE,
                             fixed = c(FALSE, FALSE))

             )
