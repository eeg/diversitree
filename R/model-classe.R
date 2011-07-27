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

# The full ClaSSE parameter structure is a vector with, for n states:
#     n * n * (n+1) / 2  speciation rates (lambda_ijk, j <= k)
#     n                  extinction rates (mu_i)
#     n * n - n          transition rates (q_ij, i != j)
#   = (n + 3) * n^2 / 2  elements

## 1: make
make.classe <- function(tree, states, k, sampling.f=NULL, strict=TRUE,
                       control=list()) {
  control <- check.control.ode(control)
  backend <- control$backend

  unresolved <- NULL
  nt.extra <- 10
  cache <- make.cache.musse(tree, states, k, sampling.f, strict)
  
  # some cvodes stuff will come here
  branches <- make.branches.classe(cache, control)

  initial.conditions.classe <- make.initial.conditions.classe(k)

  ll.classe <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                       root.p=NULL, intermediates=FALSE) {
    npars <- (k + 3) * k * k / 2
    if ( length(pars) != npars )
      stop(sprintf("Invalid length parameters (expected %d)",
                   npars))
    if ( any(!is.finite(pars)) || any(pars < 0) )
      return(-Inf)
    if ( !is.null(root.p) &&  root != ROOT.GIVEN )
      warning("Ignoring specified root state")

    # some cvodes stuff will come here
    ll.xxsse.classe(pars, cache, initial.conditions.classe, branches,
             condition.surv, root, root.p, intermediates)
  }

  # ll <- function(pars, ...) ll.classe(pars, ...)
  class(ll.classe) <- c("classe", "function")
  attr(ll.classe, "k") <- k
  ll.classe
}

## 2: print
print.classe <- function(x, ...) {
  cat("ClaSSE likelihood function:\n")
  print(unclass(x))
}

## 3: argnames / argnames <-
argnames.classe <- function(x, k=attr(x, "k"), ...) {
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
`argnames<-.classe` <- function(x, value) {
  k <- attr(x, "k")
  if ( length(value) != (k+3)*k*k/2 )
    stop("Invalid names length")
  if ( any(duplicated(value)) )
    stop("Duplicate argument names")
  attr(x, "argnames") <- value
  x
}

# Note: Might eventually want utilities for converting from/to list-of-arrays
# parameter specification: flatten.pars.classe() and listify.pars.classe().

## 4: find.mle
find.mle.classe <- function(func, x.init, method, fail.value=NA, ...) {
  if ( missing(method) )
    method <- "subplex"
  NextMethod("find.mle", method=method, class.append="fit.mle.classe")
}

mcmc.classe <- mcmc.lowerzero

## Make requires the usual functions:
## 5: make.cache (initial.tip)
# Note: classe uses the same functions as musse here:
#       initial.tip.classe() = initial.tip.musse()
#       make.cache.classe() = make.cache.musse()

## 6: ll.classe is done within make.classe

## 7: initial.conditions:
# save on index computations by wrappping initial.conditions.classe
make.initial.conditions.classe <- function(n)
{
  # n = number of states; called k elsewhere but k is used as an index below
  nseq <- seq_len(n)
  lam.idx <- matrix(seq_len(n*n*(n+1)/2), byrow=T, nrow=n)

  idxD <- (n+1):(2*n)
  j <- rep(nseq, times=seq(n,1,-1))
  k <- unlist(lapply(1:n, function(i) nseq[i:n]))
  d <- rep(NA, n)

  initial.conditions.classe <- function(init, pars, t, is.root=FALSE) {
    # E_i(t), same for N and M
    e <- init[nseq,1]

    # D_i(t), formed from N and M
    DM <- init[idxD,1]
    DN <- init[idxD,2]
    DM.DN <- 0.5 * (DM[j] * DN[k] + DM[k] * DN[j])
    for (i in nseq) d[i] <- sum(pars[lam.idx[i,]] * DM.DN) # slower with apply

    # a touch slower but cleaner:
    #   idxlam = seq_len(n*n*(n+1)/2)
    #     d = colSums(matrix(pars[idxlam], ncol=n) * DM.DN)
    # or slightly better (but still not faster than for):
    #   lamseq = seq_len(n*n*(n+1)/2)
    #   lam.mat = matrix(lamseq, ncol=n)
    #     lam.mat[lamseq] = pars[lamseq]
    #     d = colSums(lam.mat * DM.DN)

    c(e, d)
  }
}

## 8: branches
make.branches.classe <- function(cache, control) {
  k <- cache$k
  neq <- as.integer(2*k)
  np <- as.integer((k+3)*k*k/2)
  comp.idx <- as.integer((k+1):(2*k))

  # These indices are used for flattening the parameter structure.
  qmat <- matrix(0, k, k)
  idx.qmat <- cbind(rep(1:k, each=k-1),
               unlist(lapply(1:k, function(i) (1:k)[-i])))
  x <- k * k * (k + 1) / 2 + k
  idx.lm <- seq_len(x)
  idx.q <- seq(x+1, np)
  idx.pars <- list(qmat=qmat, idx.qmat=idx.qmat, idx.q=idx.q, idx.lm=idx.lm)

  # Consequently, can't use the generic make.ode.branches().
  make.ode.branches.classe("classe", "diversitreeGP", neq, np, comp.idx,
                             control, idx.pars)
}

## This is instead of make.ode.branches() in diversitree-branches.R.
## It's separate because parameter re-organization happens within branches().
make.ode.branches.classe <- function(model, dll, neq, np, comp.idx, control,
                                       idx.pars) {
  backend <- control$backend
  safe <- control$safe
  tol <- control$tol
  eps <- control$eps

  qmat <- idx.pars$qmat
  idx.qmat <- idx.pars$idx.qmat
  idx.q <- idx.pars$idx.q
  idx.lm <- idx.pars$idx.lm

  # some cvodes stuff will go here
  initfunc <- sprintf("initmod_%s", model)
  derivs <- sprintf("derivs_%s", model)

  RTOL <- ATOL <- tol
  ode <- make.ode(derivs, dll, initfunc, neq, safe)
  branches <- function(y, len, pars, t0) {
    qmat[idx.qmat] <- pars[idx.q]
    diag(qmat) <- -rowSums(qmat)
    pars <- c(pars[idx.lm], qmat)
    ode(y, c(t0, t0+len), pars, rtol=RTOL, atol=ATOL)[-1,-1,drop=FALSE]
    # old: t(classe.ode(y, c(t0, t0+len), pars, rtol=RTOL, atol=ATOL)[-1,-1])
  }

  make.branches(branches, comp.idx, eps)
}

# don't see a need for classe.Q()

# based on starting.point.geosse()
starting.point.classe <- function(tree, k, eps=0.5) {
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
  names(p) <- argnames.classe(NULL, k=k)
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
root.p.classe <- function(vals, pars, root, root.p=NULL) {
  k <- length(vals) / 2
  d.root <- vals[(k+1):(2*k)]

  if ( root == ROOT.FLAT )
    p <- rep(1/k, k)
  else if ( root == ROOT.EQUI )
    stop("Equilibrium root freqs not available for ClaSSE")
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
root.classe <- function(vals, pars, lq, condition.surv, root.p) {
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
    # d.root <- d.root / (lambda * (1-e.root)^2) # old
    d.root <- d.root / sum(root.p * lambda * (1 - e.root)^2)
  }

  if ( is.null(root.p) ) # ROOT.BOTH
    loglik <- log(d.root) + logcomp
  else
    loglik <- log(sum(root.p * d.root)) + logcomp
  loglik
}

## modified from diversitree-branches.R: ll.xxsse()
##   only difference is names of root function calls (the above functions)
ll.xxsse.classe <- function(pars, cache, initial.conditions,
                     branches, condition.surv, root, root.p,
                     intermediates) {
  ans <- all.branches.matrix(pars, cache, initial.conditions, branches)
  # vals <- ans$init[[cache$root]]
  vals <- ans$init[,cache$root]
  root.p <- root.p.classe(vals, pars, root, root.p)
  loglik <- root.classe(vals, pars, ans$lq, condition.surv, root.p)

  if ( intermediates ) {
    ans$root.p <- root.p
    attr(loglik, "intermediates") <- ans
    attr(loglik, "vals") <- vals
  }

  loglik
}
