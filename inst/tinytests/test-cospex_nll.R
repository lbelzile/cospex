# Check that the full likelihood
loc <- cbind(runif(10), runif(10))
dist <- as.matrix(dist(loc))
x0 <- rexp(100)


cospex_nll(pars = ,
           dist = ,
           covFn = covMaternFn,
           scaleFn = scale2Fn,
           locxFn = locFunExp,
           )
