# used to be solve.gse2.C.guts
RTOL <- ATOL <- 1e-8
eps <- 0
branches.gse2 <- function(y, len, pars, t0) {
  ret <- t(gse2.ode(y, c(t0, t0+len), pars, rtol=RTOL, atol=ATOL)[-1,-1])
  if ( all(ret[,4:6] >= eps) ) {
    # TODO: should compare all of 4, 5, and 6 on next line
    q <- ret[cbind(seq_along(len), as.integer(ret[,5] > ret[,6]) + 5)]
    i <- q > 0
    ret[i,4:6] <- ret[i,4:6] / q[i]
    lq <- q
    lq[i] <- log(q[i])
    cbind(lq, ret, deparse.level=0)
  } else {
    ti <- len[length(len)]/2
    len1 <- c(len[len <= ti], ti)
    len2 <- len[len > ti] - ti
    n1 <- length(len1)
    ret1 <- Recall(y, len1, pars, t0)
    ret2 <- Recall(ret1[n1,2:7], len2, pars, t0 + ti)
    ret2[,1] <- ret2[,1] + ret1[n1,1]
    rbind(ret1[-n1,], ret2)
  }
}

branches.unresolved.gse2 <- function(...)
  stop("Cannot use unresolved clades with GSE2")

initial.conditions.gse2 <- function(init, pars, t, is.root=FALSE) {
  # E.0, E.1, E.2
  e <- init[1, c(1,2,3)]

  # D.1, D.2  (Eq. 6bc)
  d12 <- init[1, c(5,6)] * init[2, c(5,6)]

  if ( !is.root )
  {
    # pars[1:3] = sA, sB, sAB
    d12 <- d12 * pars[c(1,2)]

    # D.0 (Eq. 6a)
    d0 <- 0.5 * sum(init[1, c(4,5)] * init[2, c(5,4)] * pars[1] + 
                    init[1, c(4,6)] * init[2, c(6,4)] * pars[2] +
                    init[1, c(5,6)] * init[2, c(6,5)] * pars[3])
  } else
  {
    d0 <- 0.5 * sum(init[1, c(4,5)] * init[2, c(5,4)] + 
                    init[1, c(4,6)] * init[2, c(6,4)] +
                    init[1, c(5,6)] * init[2, c(6,5)])
  }
  d <- c(d0, d12)

  c(e, d)
}

ROOT.FLAT  <- 1
ROOT.EQUI  <- 2
ROOT.OBS   <- 3
ROOT.GIVEN <- 4
ROOT.ALL   <- 5   # rather than bisse's ROOT.BOTH

gse2.ll <- function(cache, pars, prior=NULL, root=ROOT.OBS,
                     condition.surv=TRUE, root.p=NA,
                     intermediates=FALSE) {
  if ( any(pars < 0) || any(!is.finite(pars)) )
    return(-Inf)
  ans <- all.branches(pars, cache, initial.conditions.gse2,
                      branches.gse2, branches.unresolved.gse2)
  loglik <- root.gse2(ans$init[cache$root,], ans$lq, pars, root,
                       condition.surv, root.p)
  if ( !is.null(prior) )
    loglik <- loglik + prior.gse2(pars, prior)
  if ( intermediates ) {
    attr(loglik, "cache") <- cache
    attr(loglik, "intermediates") <- ans
    attr(loglik, "vals") <- ans$init[cache$root,]
    attr(loglik, "logComp") <- sum(ans$lq)
  }

  loglik
}


prior.gse2 <- function(pars, r)
  - sum(pars * r)

make.gse2 <- function(tree, states, unresolved=NULL, sampling.f=NULL,
                       nt.extra=10) {
  cache <- make.cache.gse2(tree, states, unresolved=unresolved,
                            sampling.f=sampling.f, nt.extra=10)
  f <- function(pars, ..., fail.value = NULL) {
    if (!is.null(fail.value)) 
      protect(gse2.ll, fail.value)(cache, pars, ...)
    else gse2.ll(cache, pars, ...)
  }
  class(f) <- c("gse2", "function")
  f
}

## GSE2 root calculations
# used to be in gse2.ll
root.gse2 <- function(vars, lq, pars, root, condition.surv, root.p=NA) {
  if ( !is.na(root.p[1]) )
  {
    if ( root != ROOT.GIVEN )
      warning("Ignoring specified root state")
    if ( length(root.p) != 3 )
      stop("root.p must be a vector of length 3")
  }
  
  e.root <- vars[1:3]
  d.root <- vars[4:6]

  # p = vector of weights for root state (not a scalar like bisse's root.p0)
  if ( root == ROOT.FLAT )
      p <- rep(1/3., 3)
  else if ( root == ROOT.EQUI)
      p <- gse2.stationary.freq(pars)
  else if ( root == ROOT.OBS )
      p <- d.root / sum(d.root)
  else if ( root == ROOT.GIVEN )
      p <- root.p
  else if ( root == ROOT.ALL )
      p <- rep(1., 3)
  else
      stop(paste("Invalid root mode", root))

  if ( condition.surv )
    d.root <- d.root / (1-e.root)^2
  logcomp <- sum(lq)
  if ( root == ROOT.ALL )
    loglik <- log(d.root)# + logcomp
  else
    loglik <- log(sum(p * d.root)) + logcomp  # okay for ROOT.ALL?
  loglik
}

# returns a vector, unlike bisse.stationary.freq()
gse2.stationary.freq <- function(pars)
{
  sA  <- pars[1]
  sB  <- pars[2]
  sAB <- pars[3]
  xA  <- pars[4]
  xB  <- pars[5]
  dA  <- pars[6]
  dB  <- pars[7]

  A <- matrix(c(
               -xA - xB - sAB,   dA,             dB,
                sA + xB + sAB,   sA - dA - xA,   0,
                sB + xA + sAB,   0,              sB - dB - xB
               ), nrow=3, byrow=T)

  # continuous time, so the dominant eigenvalue is the largest one
  # find its index and get that eigenvector
  # return it, normalized so the elements sum to 1
  evA <- eigen(A)
  i <- which(evA$values == max(evA$values))
  freqs <- evA$vectors[,i] / sum(evA$vectors[,i])
  return(freqs)
}

# dunno if this is much better than nothing...
# TODO: come up with something more sensible
starting.point.gse2 <- function(tree, q.div=5) {
  fit <- suppressWarnings(birthdeath(tree))
  r <- fit$para[2]
  e <- fit$para[1]
  p <- rep(c(r/(1-e), r*e/(1-e), r/q.div), each=2)
  names(p) <- c("sA", "sB", "sAB", "xA", "xB", "dA", "dB")
  p
}

#almost identical to make.cache.bisse(), but uses three states
# used to be gse2cache()
make.cache.gse2 <- function(tree, states, unresolved=NULL,
                             sampling.f=NULL, nt.extra=10) {
  if ( is.null(names(states)) )
    stop("The states vector must contain names")

  ## Check 'sampling.f'
  if ( !is.null(sampling.f) && !is.null(unresolved) )
    stop("Cannot specify both sampling.f and unresolved")
  else if ( is.null(sampling.f) )
    sampling.f <- c(1, 1, 1)
  else if ( length(sampling.f) != 3 )
    stop("sampling.f must be of length 3 (or NULL)")
  else if ( max(sampling.f) > 1 || min(sampling.f) < 0 )
    stop("sampling.f must be on range [0,1]")
  
  if ( !is.null(unresolved) ) {
      stop("Unresolved clades not yet available for GSE2")
  }
  
  ## Check that we know about all required species (this requires
  ## processing the unresolved clade information).
  known <- names(states)
  if ( !is.null(unresolved) )
    known <- unique(c(known, as.character(unresolved$tip.label)))
  if ( !all(tree$tip.label %in% known) )
    stop("Not all species have state information")
  states <- states[tree$tip.label]
  names(states) <- tree$tip.label

  cache <- make.cache(tree)
  cache$tip.state  <- states
  cache$unresolved <- unresolved
  cache$sampling.f <- sampling.f
  cache$nt.extra   <- nt.extra
  cache$y <- initial.tip.gse2(cache)
  cache
}

## Initial conditions at the tips are given by their tip states:
## There are four types of initial condition in gse2:
##   state0: c(f_0, 0,   0,    1-f_0, 1-f_1, 1-f_2)
##   state1: c(0,   f_1, 0,    1-f_0, 1-f_1, 1-f_2)
##   state3: c(0,   0,   f_2,  1-f_0, 1-f_1, 1-f_2)
##   state?: c(f_0, f_1, f_2,  1-f_0, 1-f_1, 1-f_2)
## Build this small matrix of possible initial conditions (y), then
## work out which of the three types things are (i).  Note that y[i,]
## gives the full initial conditions.  The element 'types' contains
## the different possible conditions.
# used to be in gse2.branches()
initial.tip.gse2 <- function(cache) {
  f <- cache$sampling.f
  y <- rbind(c(1-f, f[1], 0, 0),
             c(1-f, 0, f[2], 0),
             c(1-f, 0, 0, f[3]),
             c(1-f, f))
  i <- cache$tip.state + 1
  if ( !is.null(cache$unresolved) )
    i <- i[-cache$unresolved$i]
  i[is.na(i)] <- 4
  list(y=y, i=i, types=sort(unique(i)))
}

# For GSE2, think about what a sensible set of default models would be.
# all.models.bisse <- function(f, p, ...) { ... }


# copied from constrain.R:

argnames.gse2 <- function(x, ...) {
  ret <- attr(x, "argnames")
  if ( is.null(ret) )
    c("sA", "sB", "sAB", "xA", "xB", "dA", "dB")
  else
    ret
}

`argnames<-.gse2` <- function(x, value) {
  if ( length(value) != 7 )
    stop("Invalid names length")
  attr(x, "argnames") <- value
  x  
}

print.gse2 <- function(x, ...) {
  cat("GSE2 likelihood function:\n")
  print(unclass(x))
}

# copied from find.mle.R:

find.mle.gse2 <- function(func, x.init,
                           method=c("L-BFGS-B", "Nelder-Mead", "subplex"),
                           control=list(),
                           fail.value=NULL, hessian=FALSE, ...) {
  method <- match.arg(method)
  names <- argnames(func)
  npar <- length(names)
  
  if ( is.null(names(x.init)) )
    names(x.init) <- names

  if ( inherits(func, c("constrained", "fixed")) ) {
    ## Identify the parameters we do have:
    arg.idx <- match(names, argnames(environment(func)$f))
    if ( length(x.init) == 7 ) x.init <- x.init[arg.idx]
  }

  if ( method == "subplex" ) {
    if ( !require(subplex) )
      stop("The subplex package is required")
    if ( is.null(fail.value) || is.na(fail.value) )
      fail.value <- -Inf
    control <- modifyList(list(reltol=.Machine$double.eps^0.25),
                          control)
    ans <- subplex(x.init, invert(func), control=control)
    ans$value <- -ans$value
    ans$hessian <- NULL
  } else {
    if ( is.null(fail.value) || is.na(fail.value) )
      fail.value <- func(x.init, ...) - 1000
    dx <- 1e-5
    control <- modifyList(list(fnscale=-1, ndeps=rep(dx, npar)),
                          control)
    ans <- optim(x.init, func, control=control, method="L-BFGS-B",
                 fail.value=fail.value, lower=0, ...)
    ans$hessian <- NULL    
  }
  
  names(ans)[names(ans) == "value"] <- "lnLik"
  if ( ans$convergence != 0 )
    warning("Convergence problems in find.mle: ",
            tolower(ans$message))
  class(ans) <- "mle.gse2"

  if ( hessian ) {
    if ( !require(numDeriv) )
      stop("The package numDeriv is required to compute the hessian")
    ans$hessian <- hessian(func, ans$par, ...)
  }

  ans
}

logLik.mle.gse2 <- function(object, ...) {
  ll <- object$lnLik
  attr(ll, "df") <- length(object$par)
  class(ll) <- "logLik"
  ll
}

## Code based on MASS:::anova.negbin and ape:::anova.ace
anova.mle.gse2 <- function(object, ...) {
  mlist <- c(list(object), list(...))
  if ( length(mlist) == 1L )
    stop("Need to specify more than one model")
  if ( is.null(names(mlist)) )
    names(mlist) <-
      c("full", model=sprintf("model %d", seq_len(length(mlist)-1)))
  else
    names(mlist)[1] <- "full"

  ll <- lapply(mlist, logLik)
  ll.val <- sapply(ll, as.numeric)
  chisq <- c(NA, abs(2*(ll.val[1] - ll.val[-1])))
  df <- sapply(ll, attr, "df")
  ddf <- c(NA, abs(df[1] - df[-1]))
  
  out <- data.frame(Df=df,
                    lnLik=sapply(ll, as.numeric),
                    AIC=sapply(mlist, AIC),
                    ChiSq=chisq,
                    "Pr(>|Chi|)"=1 - pchisq(chisq, ddf),
                    check.names=FALSE)
  rownames(out) <- names(mlist)
    
  class(out) <- c("anova", "data.frame")
  out
}
