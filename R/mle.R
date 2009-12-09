## I should use stats4:::mle, but this does not work well for me with
## computing the variance-covariance matrix.
find.mle <- function(func, x.init, ...) {
  UseMethod("find.mle")
}

find.mle.default <- function(func, x.init, ...) {
  ans <- optim(x.init, func, ...)
  names(ans)[names(ans) == "value"] = "lnLik"
  if ( ans$convergence != 0 )
    warning("Convergence problems in find.mle: ",
            tolower(ans$message))
  class(ans) <- "mle"
  ans
}

find.mle.bisse <- function(func, x.init, control=list(),
                           lower=NULL, fail.value=NULL, ...) {
  if ( is.null(lower) )
    lower <- c(1e-4, 1e-4, 0, 0, 1e-4, 1e-4)
  if ( is.null(fail.value) || is.na(fail.value) )
    fail.value <- func(x.init) - 1000

  dx <- 1e-5
  par.names <- c("lambda0", "lambda1", "mu0", "mu1", "q01", "q10")

  if ( inherits(func, c("constrained", "fixed")) ) {
    free <- environment(func)$free
    control <- modifyList(list(fnscale=-1, ndeps=rep(dx, sum(free))),
                          control)
    if ( length(lower)  == length(free) )
      lower  <- lower[free]
    if ( length(x.init) == length(free) )
      x.init <- x.init[free]
    if ( is.null(names(x.init)) )
      names(x.init) <- par.names[free]
  } else {
    if ( is.null(names(x.init)) )
      names(x.init) <- par.names
    control <- modifyList(list(fnscale=-1, ndeps=rep(dx, 6)), control)
  }

  ans <- find.mle.default(func, x.init, control=control, lower=lower,
                          method="L-BFGS-B", fail.value=fail.value,
                          ...)
  class(ans) <- c(class(ans), "mle.bisse")

  at.edge <- ans$par - lower == 0 & lower > 0
  if ( any(at.edge) ) {
    at.edge <- at.edge[at.edge]
    warning(sprintf("Parameter(s) %s at edge of parameter space",
                    paste(which(at.edge), collapse=", ")))
  }

  ans
}





## find.mle <- function(func, x.init, lower, upper=Inf,
##                      constraint=NULL, fixed=NULL, ...) {
##   names.full <- c("lambda0", "lambda1", "mu0", "mu1", "q01", "q10")
##   full <- is.null(constraint) && is.null(fixed)
  
##   if ( missing(lower) )
##     lower <- c(1e-6, 1e-6, 0, 0, 1e-6, 1e-6)

##   if ( full ) {
##     f <- func
##     npar <- 6
##     ok <- 1:6
##   } else if ( !is.null(fixed) ) {
##     f <- function(x)
##       bisse(obj, spread(x, fixed), ...)
##     ok <- which(is.na(fixed))
##     npar <- length(ok)
##   } else if ( !is.null(constraint) ) {
##     f <- function(x)
##       bisse(obj, add.constrained(x, constraint), ...)
##     ok <- which(is.na(constraint))
##     npar <- length(ok)
##   } else {
##     stop("Not yet implemented!")
##   }

##   names(x.init) <- names.full[ok]

##   if ( !full ) {
##     if ( length(lower) == 6 ) lower <- lower[ok]
##     if ( length(upper) == 6 ) upper <- upper[ok]
##     if ( length(x.init) == 6 ) x.init <- x.init[ok]
##   }
  
##   est <- optim(x.init, f, method="L-BFGS-B", lower=lower, upper=upper,
##                control=list(fnscale=-1))

##   if ( est$convergence != 0 )
##     warning("Convergence was not complete; check the output carefully")
##   else
##     est <- est[c("par", "value")]
  
##   names(est)[2] <- "logLik"
##   class(est) <- "mle"

##   if ( !is.null(fixed) )
##     est$allPars <- add.fixed(est$par, fixed)
##   else if ( !is.null(constraint) )
##     est$allPars <- add.constrained(est$par, constraint)

##   if (!full)
##     names(est$allPars) <- names.full

##   est
## }

## mle.phylo <- mle.phylo4 <- function(obj, states, unresolved=NULL,
##                                     ...) {
##   mle(bissecache(obj, states, unresolved), ...)
## }
