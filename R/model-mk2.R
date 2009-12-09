## The Mk2 model of character evolution.  This is likely to end up
## being the Mkn model, really.  However, a simple 'mk2' interface
## will remain that does not require specifying the entire transition
## matrix.

## Special for k=2
branches.mk2 <- function(y, len, pars, t0) {
  q01 <- pars[1]
  q10 <- pars[2]
  if ( q01 + q10 > 0 ) {
    x <- exp(-(q01+q10)*len) * (y[1] - y[2])
    z <- q10 * y[1] + q01 * y[2]
    ret <- cbind(z + x * q01, z - x * q10) / (q01 + q10)
  } else {
    ret <- matrix(rep(y, length(len)), length(len), 2, TRUE)
  }

  q <- apply(ret, 1, min)
  i <- q > 0
  ret[i,] <- ret[i,] / q[i]
  lq <- q
  lq[i] <- log(q[i])
  cbind(lq, ret, deparse.level=0)
}

## The n-state version is not much different:
make.branches.mkn <- function(k) {
  if ( k == 2 )
    warning("Two states is faster with Mk2")
  make.ode <- diversitree:::make.ode
  mkn.ode <- make.ode("derivs_mkn", "diversitree", "initmod_mkn", k, FALSE)
  RTOL <- ATOL <- 1e-8
  eps <- 0
  
  function(y, len, pars, t0) {
    ret <- t(mkn.ode(y, c(t0, t0+len), pars, rtol=RTOL, atol=ATOL)[-1,-1])
    if ( all(ret >= eps) ) {
      q <- apply(ret, 1, min)
      i <- q > 0
      ret[i,] <- ret[i,] / q[i]
      lq <- q
      lq[i] <- log(q[i])
      cbind(lq, ret, deparse.level=0)
    } else {
      ti <- len[length(len)]/2
      len1 <- c(len[len <= ti], ti)
      len2 <- len[len > ti] - ti
      n1 <- length(len1)
      ret1 <- Recall(y, len1, pars, t0)
      ret2 <- Recall(ret1[n1,-1], len2, pars, t0 + ti)
      ret2[,1] <- ret2[,1] + ret1[n1,1]
      rbind(ret1[-n1,], ret2)
    }
  }
}

## Most everything else applies to both 2-state and k-state:
initial.conditions.mkn <- function(init, pars, t, is.root=FALSE)
  init[1,] * init[2,]
branches.unresolved.mkn <- function(...)
  stop("Cannot use unresolved clades with Mk2/Mkn")

make.cache.mkn <- function(tree, states, k=2) {
  if ( is.null(names(states)) )
    stop("The states vector must contain names")
  if ( !all(tree$tip.label %in% names(states)) )
    stop("Not all species have state information")
  states <- states[tree$tip.label]
  names(states) <- tree$tip.label
  cache <- make.cache(tree)
  cache$k <- k
  cache$tip.state  <- states
  cache$y <- initial.tip.mkn(cache)
  cache
}
initial.tip.mkn <- function(cache) {
  k <- cache$k
  tip.state <- cache$tip.state
  if ( any(tip.state < 1 | tip.state > k, na.rm=TRUE) )
    stop(sprintf("tip states must be in the range [1, %d]", k))
  
  y <- matrix(rep(rep(0, k), k + 1), k+1, k, TRUE)
  y[k+1,] <- diag(y[1:k,1:k]) <- 1
  i <- cache$tip.state
  i[is.na(i)] <- k + 1
  list(y=y, i=i, types=sort(unique(i)))
}
root.mkn <- function(vals, lq, pars, root, root.p=NULL) {
  k <- length(vals)
  if ( !is.null(root.p) ) {
    if ( root != ROOT.GIVEN )
      warning("Ignoring specified root state")
    else if ( length(root.p) != k )
      stop("root.p of wrong length")
  }

  if ( root == ROOT.FLAT )
    p <- rep(1/k, k)
  else if ( root == ROOT.EQUI )
    p <- stationary.freq.mkn(pars)
  else if ( root == ROOT.OBS )
    p <- vals/sum(vals)
  else if ( root == ROOT.GIVEN )
    p <- root.p
  else if ( root != ROOT.BOTH )
    stop("Invalid root mode")

  logcomp <- sum(lq)
  if ( root == ROOT.BOTH )
    loglik <- log(vals) + logcomp
  else
    loglik <- log(sum(p * vals)) + logcomp
  loglik
}
stationary.freq.mkn <- function(pars) {
  .NotYetImplemented()
}
prior.mkn <- function(pars, r)
  - sum(pars * r)


## And the likelihood functions:
make.mk2 <- function(tree, states) {
  cache <- make.cache.mkn(tree, states + 1, 2)
  f <- function(pars, ..., fail.value = NULL) {
    if (!is.null(fail.value)) 
      protect(mk2.ll, fail.value)(cache, pars, ...)
    else mk2.ll(cache, pars, ...)
  }
  class(f) <- c("mkn", "function")
  f
}

mk2.ll <- function(cache, pars, prior=NULL, root=ROOT.OBS,
                   root.p=NULL) {
  if ( length(pars) != 2 )
    stop("Invalid length parameters")
  if ( any(pars < 0) || any(!is.finite(pars)) )
    return(-Inf)
  ans <- all.branches(pars, cache, initial.conditions.mkn, branches.mk2,
                      branches.unresolved.mkn)
  loglik <- root.mkn(ans$init[cache$root,], ans$lq, pars, root, root.p)
  if ( !is.null(prior) )
    loglik <- loglik + prior.mkn(pars, prior)
  loglik
}

make.mkn <- function(tree, states, k) {
  cache <- make.cache.mkn(tree, states)
  branches <- make.branches.mkn(k)
  qmat <- matrix(0, k, k)
  idx <- cbind(rep(1:k, each=k-1),
               unlist(lapply(1:k, function(i) (1:k)[-i])))

  mkn.ll <- function(cache, pars, prior=NULL, root=ROOT.OBS,
                     root.p=NULL) {
    if ( length(pars) != k*(k-1) )
      stop("Invalid length parameters")
    if ( any(!is.finite(pars)) )
      return(-Inf)
    qmat[idx] <- pars
    diag(qmat) <- -rowSums(qmat)
    
    ans <- all.branches(qmat, cache, initial.conditions.mkn, branches,
                        branches.unresolved.mkn)
    loglik <- root.mkn(ans$init[cache$root,], ans$lq, pars, root, root.p)
    if ( !is.null(prior) )
      loglik <- loglik + prior.mkn(pars, prior)
    loglik
  }

  f <- function(pars, ..., fail.value = NULL) {
    if (!is.null(fail.value)) 
      protect(mkn.ll, fail.value)(cache, pars, ...)
    else mkn.ll(cache, pars, ...)
  }
  class(f) <- c("mkn", "function")
  f
}

## ## For reference, as this will be useful for doing the
## ## reconstructions?
## Pij <- function(t, q01, q10) {
##   x <- exp(-(q01+q10)*t)
##   cbind((x*q01 + q10), (1 - x)*q01,
##         (1 - x)*q10, (x*q10 + q01)) / (q01 + q10)
## }
find.mle.mkn <- function(func, x.init,
                        method=c("nlminb", "L-BFGS-B", "Nelder-Mead",
                          "subplex"),
                        control=list(),
                        fail.value=NULL, hessian=FALSE, ...) {
  method <- match.arg(method)
  names <- argnames(func)

  ## Parameters pulled from nowhere - should at least get these from
  ## the tree?  That said, no other function gets starting help.
  if ( missing(x.init) )
    x.init <- structure(c(.1, .1), names=argnames(func))
  if ( is.null(names(x.init)) )
    names(x.init) <- argnames(func)

  if ( inherits(func, "constrained") ) {
    ## Identify the parameters we do have:
    arg.idx <- match(names, argnames(environment(func)$f))
    if ( length(x.init) == 2 ) x.init <- x.init[arg.idx]
  }

  if ( method=="nlminb" ) {
    if ( is.null(fail.value) || is.na(fail.value) )
      fail.value <- func(x.init, ...) - 1000
    fit <- nlminb(x.init, invert(func), lower=c(0, 0),
                  upper=c(Inf, Inf), ...)
    names(fit)[names(fit) == "objective"] <- "lnLik"
    names(fit)[names(fit) == "evaluations"] <- "counts"
    fit$lnLik <- -fit$lnLik
    fit <- fit[c("par", "lnLik", "counts", "convergence", "message",
                 "iterations")]
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
  if ( fit$convergence != 0 )
    warning("Convergence problems in find.mle: ",
            tolower(fit$message))

  if ( hessian ) {
    if ( !require(numDeriv) )
      stop("The package numDeriv is required to compute the hessian")
    fit$hessian <- hessian(func, fit$par, ...)
  }

  class(fit) <- "mle.mkn"
  fit
}

argnames.mkn <- function(x, ...) {
  ret <- attr(x, "argnames")
  if ( is.null(ret) ) {
    k <- environment(x)$cache$k
    if ( k == 2 )
      c("q01", "q10")
    else
      sprintf("q%d%d", rep(1:k, each=k), rep(1:k, k))
  } else {
    ret
  }
}

`argnames<-.mkn` <- function(x, value) {
  k <- environment(x)$cache$k  
  if ( length(value) != k*k )
    stop("Invalid names length")
  attr(x, "argnames") <- value
  x  
}

logLik.mle.mkn <- function(object, ...) {
  ll <- object$lnLik
  attr(ll, "df") <- length(object$par)
  class(ll) <- "logLik"
  ll
}

anova.mle.mkn <- anova.mle.bisse
prior.mkn <- function(pars, r)
  - sum(pars * r)

print.mkn <- function(x, ...) {
  cat("Mk-n likelihood function:\n")
  print(unclass(x))
}
