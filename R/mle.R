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

# used by bisse and gse2 methods
find.mle.core <- function(func, x.init, control, lower, fail.value, 
                          par.names, ...)
{
    if ( is.null(lower) )
        lower <- c(1e-4, 1e-4, 0, 0, 1e-4, 1e-4)  # 6 vs 7 param danger here
    if ( is.null(fail.value) || is.na(fail.value) )
        fail.value <- func(x.init) - 1000

    dx <- 1e-5

    if ( inherits(func, c("constrained", "fixed")) )
    {
        free <- environment(func)$free
        control <- modifyList(list(fnscale=-1, ndeps=rep(dx, sum(free))),
                              control)
        if ( length(lower)  == length(free) )
            lower  <- lower[free]
        if ( length(x.init) == length(free) )
            x.init <- x.init[free]
        if ( is.null(names(x.init)) )
            names(x.init) <- par.names[free]
    } else    # FIXME: optim() error about x.init when unconstrained?
    {
        control <- modifyList(list(fnscale=-1, 
                              ndeps=rep(dx, length(par.names))), control)
        if ( is.null(names(x.init)) )
            names(x.init) <- par.names
    }

    ans <- find.mle.default(func, x.init, control=control, lower=lower,
                            method="L-BFGS-B", fail.value=fail.value,
                            ...)
    at.edge <- ans$par - lower == 0 & lower > 0
    if ( any(at.edge) )
    {
        at.edge <- at.edge[at.edge]
        warning(sprintf("Parameter(s) %s at edge of parameter space",
                        paste(which(at.edge), collapse=", ")))
    }

    return(ans)
}

find.mle.bisse <- function(func, x.init, control=list(),
                           lower=NULL, fail.value=NULL, ...)
{
    par.names <- c("lambda0", "lambda1", "mu0", "mu1", "q01", "q10")
    ans <- find.mle.core(func, x.init, control, lower, fail.value, 
                         par.names, ...)
    class(ans) <- c(class(ans), "mle.bisse")
    return(ans)
}

find.mle.gse2 <- function(func, x.init, control=list(),
                           lower=NULL, fail.value=NULL, ...)
{
    par.names <- c("sA", "sB", "sAB", "xA", "xB", "dA", "dB")
    ans <- find.mle.core(func, x.init, control, lower, fail.value, 
                         par.names, ...)
    class(ans) <- c(class(ans), "mle.gse2")
    return(ans)
}
