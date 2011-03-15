## Models should provide:
##   1. make
##   2. print
##   3. argnames / argnames<-
##   4. find.mle
## Generally, make will require:
##   5. make.cache (also initial.tip, root)
##   6. ll
##   7. initial.conditions
##   8. branches

# The full PunctSSE parameter structure is a vector with, for n states:
#     n * n * (n+1) / 2  speciation rates (lambda_ijk, j <= k)
#     n                  extinction rates (mu_i)
#     n * n - n          transition rates (q_ij, i != j)
#   = (n + 3) * n^2 / 2  elements

## 1: make
make.punctsse <- function(tree, states, k, sampling.f=NULL, strict=TRUE,
                       safe=FALSE) {
  cache <- make.cache.musse(tree, states, k, sampling.f, strict)
  branches <- make.branches.punctsse(k, safe)

  ll.punctsse <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                       root.p=NULL, intermediates=FALSE) {
    npars <- (k + 3) * k * k / 2
    if ( length(pars) != npars )
      stop(sprintf("Invalid length parameters (expected %d)",
                   npars))
    if ( any(!is.finite(pars)) || any(pars < 0) )
      return(-Inf)
    if ( !is.null(root.p) &&  root != ROOT.GIVEN )
      warning("Ignoring specified root state")

    ll.xxsse.punctsse(pars, cache, initial.conditions.punctsse, branches,
             condition.surv, root, root.p, intermediates)
  }

  ll <- function(pars, ...) ll.punctsse(pars, ...)
  class(ll) <- c("punctsse", "function")
  attr(ll, "k") <- k
  ll
}

## 2: print
print.punctsse <- function(x, ...) {
  cat("PunctSSE likelihood function:\n")
  print(unclass(x))
}

## 3: argnames / argnames <-
argnames.punctsse <- function(x, k=attr(x, "k"), ...) {
  ret <- attr(x, "argnames")
  if ( is.null(ret) ) {
    fmt <- sprintf("%%0%dd", ceiling(log10(k + .5)))
    sstr <- sprintf(fmt, 1:k)
    lambda.names <- sprintf("lambda%s%s%s", rep(sstr, each=k*(k+1)/2),
                            rep(rep(sstr, times=seq(k,1,-1)), k), 
                            unlist(lapply(1:k, function(i) sstr[i:k])))
    mu.names <- sprintf("mu%s", sstr)
    q.names <- sprintf("q%s%s", rep(sstr, each=k-1), 
                       unlist(lapply(1:k, function(i) sstr[-i])))
    c(lambda.names, mu.names, q.names)
  } else {
    ret
  }
}
`argnames<-.punctsse` <- function(x, value) {
  k <- attr(x, "k")
  if ( length(value) != (k+3)*k*k/2 )
    stop("Invalid names length")
  if ( any(duplicated(value)) )
    stop("Duplicate argument names")
  attr(x, "argnames") <- value
  x
}

# Note: Might eventually want utilities for converting from/to list-of-arrays
# parameter specification: flatten.pars.punctsse() and listify.pars.punctsse().

## 4: find.mle
find.mle.punctsse <- function(func, x.init, method, fail.value=NA, ...) {
  if ( missing(method) )
    method <- "subplex"
  NextMethod("find.mle", method=method, class.append="fit.mle.punctsse")
}

mcmc.punctsse <- mcmc.lowerzero

## Make requires the usual functions:
## 5: make.cache (initial.tip)
# Note: punctsse uses the same functions as musse here:
#       initial.tip.punctsse() = initial.tip.musse()
#       make.cache.punctsse() = make.cache.musse()

## 6: ll.punctsse is done within make.punctsse

## 7: initial.conditions:
initial.conditions.punctsse <- function(init, pars, t, is.root=FALSE) {
  n <- length(init[[1]])/2 # called k elsewhere, but k is used as an index below
  nseq <- seq_len(n)
  nlam <- n*(n+1)/2

  # E_i(t), same for N and M
  e <- init[[1]][nseq]

  # D_i(t), formed from N and M
  DM <- init[[1]][(n+1):(2*n)]
  DN <- init[[2]][(n+1):(2*n)]
  j <- rep(nseq, times=seq(n,1,-1))
  k <- unlist(lapply(1:n, function(i) nseq[i:n]))
  get.di <- function(i)
  {
    x <- seq((i-1)*nlam+1, i*nlam)
    sum(pars[x] * 0.5 * (DM[j] * DN[k] + DM[k] * DN[j]))
  }
  d <- unlist(lapply(nseq, get.di))

  c(e, d)
}

## 8: branches
make.branches.punctsse <- function(k, safe=FALSE) {
  RTOL <- ATOL <- 1e-8

  qmat <- matrix(0, k, k)
  idx.qmat <- cbind(rep(1:k, each=k-1),
               unlist(lapply(1:k, function(i) (1:k)[-i])))
  x <- k * k * (k + 1) / 2 + k
  idx.lm <- seq_len(x)
  idx.q <- seq(x+1, (k+3)*k*k/2)

  punctsse.ode <- make.ode("derivs_punctsse", "diversitreeGP",
                        "initmod_punctsse", 2*k, FALSE)

  branches.punctsse <- function(y, len, pars, t0) {
    qmat[idx.qmat] <- pars[idx.q]
    diag(qmat) <- -rowSums(qmat)
    pars <- c(pars[idx.lm], qmat)
    t(punctsse.ode(y, c(t0, t0+len), pars, rtol=RTOL, atol=ATOL)[-1,-1])
  }
  make.branches(branches.punctsse, (k+1):(2*k))
}

# don't see a need for punctsse.Q()

# based on starting.point.geosse()
starting.point.punctsse <- function(tree, k, eps=0.5) {
  if (eps == 0) {
    s <- (log(Ntip(tree)) - log(2)) / max(branching.times(tree))
    x <- 0
    d <- s/10
  } else {
    n <- Ntip(tree)
    r <- ( log( (n/2) * (1 - eps*eps) + 2*eps + (1 - eps)/2 *
           sqrt( n * (n*eps*eps - 8*eps + 2*n*eps + n))) - log(2)
         ) / max(branching.times(tree))
    s <- r / (1 - eps)
    x <- s * eps
    q <- s - x
  }
  p <- c( rep(s / (k*(k+1)/2), k*k*(k+1)/2 ), rep(x, k), rep(q, k*(k-1)) )
  names(p) <- argnames.punctsse(NULL, k=k)
  p
}

## NOTE:
## Functions below are taken from diversitree-branches.R.  Internal
## modifications were necessary, but the parent functions could likely be
## generalized by passing some more functions as arguments.  At the moment,
## though, separate seems better than integrated.

# modified from diversitree-branches.R: root.p.xxsse()
#   returned p is always a vector of length k (or NULL)
#   equilibrium freqs not available even for k = 2
root.p.punctsse <- function(vals, pars, root, root.p=NULL) {
  k <- length(vals) / 2
  d.root <- vals[(k+1):(2*k)]

  if ( root == ROOT.FLAT )
    p <- rep(1/k, k)
  else if ( root == ROOT.EQUI )
    stop("Equilibrium root freqs not available for PunctSSE")
  else if ( root == ROOT.OBS )
    p <- d.root / sum(d.root)
  else if ( root == ROOT.GIVEN ) {
    if ( length(root.p) != length(d.root) )
      stop("Invalid length for root.p")
    p <- root.p
  } else if ( root == ROOT.ALL )
    p <- NULL
  else
    stop("Invalid root mode")
  p
}

## modified from diversitree-branches.R: root.xxsse()
##   the only difference is lambda
root.punctsse <- function(vals, pars, lq, condition.surv, root.p) {
  logcomp <- sum(lq)

  k <- length(vals) / 2
  i <- seq_len(k)

  e.root <- vals[i]
  d.root <- vals[-i]

  if ( condition.surv )
  {
    # species in state i are subject to all lambda_ijk speciation rates
    nsum <- k*(k+1)/2
    lambda <- colSums(matrix(pars[1:(nsum*k)], nrow=nsum))
    d.root <- d.root / (lambda * (1-e.root)^2)
  }

  if ( is.null(root.p) ) # ROOT.BOTH
    loglik <- log(d.root) + logcomp
  else
    loglik <- log(sum(root.p * d.root)) + logcomp
  loglik
}

## modified from diversitree-branches.R: ll.xxsse()
##   only difference is names of root function calls (the above functions)
ll.xxsse.punctsse <- function(pars, cache, initial.conditions,
                     branches, condition.surv, root, root.p,
                     intermediates) {
  ans <- all.branches(pars, cache, initial.conditions, branches)
  vals <- ans$init[[cache$root]]
  root.p <- root.p.punctsse(vals, pars, root, root.p)
  loglik <- root.punctsse(vals, pars, ans$lq, condition.surv, root.p)
  ans$root.p <- root.p
  cleanup(loglik, pars, intermediates, cache, ans)
}
