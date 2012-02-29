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
  cache <- make.cache.musse(tree, states, k, sampling.f, strict)
  
  if ( backend == "CVODES" )
    all.branches <- make.all.branches.C.classe(cache, control)
  else
    branches <- make.branches.classe(cache, control)

  initial.conditions <- make.initial.conditions.classe(k)

  f.pars <- make.classe.pars(k)

  ll.classe <- function(pars, condition.surv=TRUE, root=ROOT.OBS,
                       root.p=NULL, intermediates=FALSE) {
    check.pars.classe(pars, k)
    if ( !is.null(root.p) &&  root != ROOT.GIVEN )
      warning("Ignoring specified root state")
    if ( root == ROOT.EQUI )
      stop("Equilibrium root freqs not available for ClaSSE")

    if ( backend == "CVODES" ) {
      ans <- all.branches(f.pars(pars), intermediates)
      lq <- ans[[1]]
      vals <- ans[[2]]

      if ( intermediates ) {
        ans <- .Call("r_get_vals", ptr, PACKAGE="diversitree")
        names(ans) <- c("init", "base", "lq")
      }
    } else {
      ans <- all.branches.matrix(f.pars(pars), cache,
                                 initial.conditions, branches)
      lq <- ans$lq
      vals <- ans$init[,cache$root]
    }
    root.p <- root.p.xxsse(vals, pars, root, root.p)
    loglik <- root.classe(vals, pars, lq, condition.surv, root.p)
  
    if ( intermediates ) {
      ans$root.p <- root.p
      attr(loglik, "intermediates") <- ans
      attr(loglik, "vals") <- vals
    }
    loglik
  }

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
## save on index computations by wrappping initial.conditions.classe
make.initial.conditions.classe <- function(n) {
  ## n = number of states; called k elsewhere but k is used as an index below
  nseq <- seq_len(n)
  lam.idx <- matrix(seq_len(n*n*(n+1)/2), byrow=TRUE, nrow=n)

  idxD <- (n+1):(2*n)
  j <- rep(nseq, times=seq(n,1,-1))
  k <- unlist(lapply(1:n, function(i) nseq[i:n]))
  d <- rep(NA, n)

  initial.conditions.classe <- function(init, pars, t, is.root=FALSE) {
    ## E_i(t), same for N and M
    e <- init[nseq,1]

    ## D_i(t), formed from N and M
    DM <- init[idxD,1]
    DN <- init[idxD,2]
    DM.DN <- 0.5 * (DM[j] * DN[k] + DM[k] * DN[j])
    for (i in nseq)
      d[i] <- sum(pars[lam.idx[i,]] * DM.DN) # slower with apply

    ## a touch slower but cleaner:
    ##   idxlam = seq_len(n*n*(n+1)/2)
    ##     d = colSums(matrix(pars[idxlam], ncol=n) * DM.DN)
    ## or slightly better (but still not faster than for):
    ##   lamseq = seq_len(n*n*(n+1)/2)
    ##   lam.mat = matrix(lamseq, ncol=n)
    ##     lam.mat[lamseq] = pars[lamseq]
    ##     d = colSums(lam.mat * DM.DN)

    c(e, d)
  }
}

## 8: branches
make.branches.classe <- function(cache, control) {
  k <- cache$k
  neq <- as.integer(2*k)
  np0 <- as.integer((k+3)*k*k/2)  # number of actual params
  np <- np0 + k                   # number, including Q's diagonal elements
  comp.idx <- as.integer((k+1):(2*k))
  make.ode.branches("classe", "diversitree", neq, np, comp.idx,
                    control)
}

make.all.branches.C.classe <- function(cache, control) {
  k <- cache$k
  neq <- as.integer(2*k)
  np0 <- as.integer((k+3)*k*k/2)  # number of actual params
  np <- np0 + k                   # number, including Q's diagonal elements
  comp.idx <- as.integer((k+1):(2*k))
 
  make.all.branches.C(cache, "classe", "diversitree",
                      neq, np, comp.idx, control)
}

## don't see a need for classe.Q()

## todo, at least for k = 2: stationary.freq.classe

## based on starting.point.geosse()
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
    ## species in state i are subject to all lambda_ijk speciation rates
    nsum <- k*(k+1)/2
    lambda <- colSums(matrix(pars[1:(nsum*k)], nrow=nsum))
    ## d.root <- d.root / (lambda * (1-e.root)^2) # old
    d.root <- d.root / sum(root.p * lambda * (1 - e.root)^2)
  }

  if ( is.null(root.p) ) # ROOT.BOTH
    loglik <- log(d.root) + logcomp
  else
    loglik <- log(sum(root.p * d.root)) + logcomp
  loglik
}

check.pars.classe <- function(pars, k) {
  npars <- (k + 3) * k * k / 2
  if ( length(pars) != npars )
    stop(sprintf("Invalid length parameters (expected %d)",
                 npars))
  if ( any(!is.finite(pars)) || any(pars < 0) )
    stop("Parameters must be non-negative and finite")
  TRUE
}

make.classe.pars <- function(k) {
  np0 <- as.integer((k+3)*k*k/2)  # number of actual params
  np <- np0 + k                   # number, including Q's diagonal elements
  
  qmat <- matrix(0, k, k)
  idx.qmat <- cbind(rep(1:k, each=k-1),
               unlist(lapply(1:k, function(i) (1:k)[-i])))
  x <- k * k * (k + 1) / 2 + k
  idx.lm <- seq_len(x)
  idx.q <- seq(x+1, np0)
  
  function(pars) {
    qmat[idx.qmat] <- pars[idx.q]
    diag(qmat) <- -rowSums(qmat)
    c(pars[idx.lm], qmat)
  }
}
