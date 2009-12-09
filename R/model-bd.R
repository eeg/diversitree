## The B-D model is different to the others in that I am not using
## most of the infrastructure - instead the entire calculations are
## done at once.
bd.ll <- function(cache, pars, prior=NULL) {
  N <- cache$N
  x <- cache$x

  r <- pars[1] - pars[2]
  a <- pars[2] / pars[1]
  ## This allows r < 0, a > 1.  a < 0 is not allowed though
  ## Also, if r < 0, then a must > 1 (and vv.).  XOR captures this.
  ## if ( a < 0 || sign(1-a) != sign(r) ) return(-Inf)
  if ( a < 0 || xor(a > 1, r < 0) )
    return(-Inf)

  loglik <- lfactorial(N - 1) + (N - 2) * log(abs(r)) + r * sum(x[3:N]) +
    N * log(abs(1 - a)) - 2*sum(log(abs(exp(r * x[2:N]) - a)))

  if ( !is.null(prior) )
    loglik <- loglik + prior.bd(pars)

  loglik
}

make.bd <- function(tree, times=branching.times(tree),
                            sampling.f=NULL, unresolved=NULL) {
  if ( !is.null(sampling.f) ) .NotYetUsed("sampling.f")
  if ( !is.null(unresolved) ) .NotYetUsed("unresolved")

  if ( !missing(times) && !missing(tree) )
    stop("times cannot be specified if tree given")
  x <- c(NA, times)
  N <- length(x)
  cache <- list(N=N, x=x)
  f <- function(pars, ..., fail.value=NULL) {
    if (!is.null(fail.value)) 
      protect(bd.ll, fail.value)(cache, pars, ...)
    else bd.ll(cache, pars, ...)
  }
  class(f) <- c("bd", "function")
  f
}

argnames.bd <- function(x, ...) {
  ret <- attr(x, "argnames")
  if ( is.null(ret) )
    c("lambda", "mu")
  else
    ret
}
`argnames<-.bd` <- function(x, value) {
  if ( length(value) != 2 )
    stop("Invalid names length")
  attr(x, "argnames") <- value
  x  
}

find.mle.bd <- function(func, x.init,
                        method=c("nlm", "L-BFGS-B", "Nelder-Mead",
                          "subplex"),
                        control=list(),
                        fail.value=NULL, hessian=FALSE, ...) {
  method <- match.arg(method)
  names <- argnames(func)

  ## I really should use parameters estimated from the Yule model
  ## here.  Currently using parameters from nowhere.
  if ( missing(x.init) )
    x.init <- structure(c(.2, .1), names=argnames(func))
  if ( is.null(names(x.init)) )
    names(x.init) <- argnames(func)
  
  if ( inherits(func, "constrained") ) {
    ## Identify the parameters we do have:
    arg.idx <- match(names, argnames(environment(func)$f))
    if ( length(x.init) == 2 )
      x.init <- x.init[arg.idx]
  }

  if ( method=="nlm" ) {
    if ( is.null(fail.value) || is.na(fail.value) )
      fail.value <- func(x.init, ...) - 1000
    fit <- nlm(invert(func), x.init, fail.value=fail.value, ...)
    names(fit)[names(fit) == "estimate"] <- "par"
    names(fit)[names(fit) == "minimum"] <- "lnLik"
    names(fit)[names(fit) == "iterations"] <- "counts"
    fit$lnLik <- -fit$lnLik 
    names(fit$par) <- names(x.init)
    fit <- fit[c("par", "lnLik", "counts", "code", "gradient")]
  } else if ( method == "subplex" ) {
    if ( !require(subplex) )
      stop("The subplex package is required")
    control <- modifyList(list(reltol=.Machine$double.eps^0.25),
                          control, ...)
    if ( is.null(fail.value) || is.na(fail.value) )
      fail.value <- -Inf
    fit <- subplex(x.init, invert(func), control=control, fail.value)
    fit$value <- -fit$value
    fit$hessian <- NULL
  } else {
    if ( is.null(fail.value) || is.na(fail.value) )
      fail.value <- func(x.init, ...) - 1000
    control <- modifyList(list(fnscale=-1), control)
    fit <- optim(x.init, func, control=control, method=method,
                 fail.value=fail.value, ...)
  }

  names(fit)[names(fit) == "value"] <- "lnLik"
  if ( method == "nlm" ) {
    if ( fit$code > 2 )
      warning("Convergence problems in find.mle: code = ",
              fit$code, " (see ?nlm for details)")
  } else {
    if ( fit$convergence != 0 )
      warning("Convergence problems in find.mle: ",
              tolower(fit$message))
  }

  if ( hessian ) {
    if ( !require(numDeriv) )
      stop("The package numDeriv is required to compute the hessian")
    fit$hessian <- hessian(func, fit$par, ...)
  }

  class(fit) <- "mle.bd"
  fit
}

logLik.mle.bd <- function(object, ...) {
  ll <- object$lnLik
  attr(ll, "df") <- length(object$par)
  class(ll) <- "logLik"
  ll
}

anova.mle.bd <- anova.mle.bisse
prior.bd <- function(pars, r)
  - sum(pars * r)

print.bd <- function(x, ...) {
  cat("Constant rate birth-death likelihood function:\n")
  print(unclass(x))
}
