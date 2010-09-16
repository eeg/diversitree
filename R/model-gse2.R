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
##   (no longer branches.unresolved)

## 1: make
make.gse2 <- function(tree, states, unresolved=NULL, sampling.f=NULL,
                       nt.extra=10, safe=FALSE, strict=TRUE) {
  cache <- make.cache.gse2(tree, states, unresolved=unresolved,
                            sampling.f=sampling.f, nt.extra=nt.extra,
                            strict=strict)
  branches <- make.branches.gse2(safe)
  ll <- function(pars, ...) ll.gse2(cache, pars, branches, ...)
  class(ll) <- c("gse2", "function")
  ll
}

## 2: print
print.gse2 <- function(x, ...) {
  cat("GSE2 likelihood function:\n")
  print(unclass(x))
}

## 3: argnames / argnames<-
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

## 4: find.mle
find.mle.gse2 <- function(func, x.init, method,
                           fail.value=NA, ...) {
  if ( missing(method) )
    method <- "subplex"
  NextMethod("find.mle", method=method, class.append="fit.mle.gse2")
}

## Make requires the usual functions:
## 5: make.cache (initial.tip, root)
# almost identical to make.cache.bisse(), but uses three states
make.cache.gse2 <- function(tree, states, unresolved=NULL,
                             sampling.f=NULL, nt.extra=10,
                             strict=TRUE) {
  ## RGF: There is a potential issue here with states, as
  ## 'unresolved' may contain one of the states.  For now I am
  ## disabling the check, but this is not great.
  if ( strict && !is.null(unresolved) ) {
    strict <- FALSE
  }

  tree <- check.tree(tree)
  states <- check.states(tree, states, strict=strict, strict.vals=0:2)

  # check unresolved
  if ( !is.null(unresolved) ) {
      stop("Unresolved clades not yet available for GSE2")
  }
  if ( inherits(tree, "clade.tree") ) {
    if ( !is.null(unresolved) )
      stop("'unresolved' cannot be specified where 'tree' is a clade.tree")
    #unresolved <- make.unresolved(tree$clades, states)
  }

  ## Check 'sampling.f'
  if ( !is.null(sampling.f) && !is.null(unresolved) )
    stop("Cannot specify both sampling.f and unresolved")
  else
    sampling.f <- check.sampling.f(sampling.f, 3)

  # would also need: unresolved <- check.unresolved(cache, unresolved, nt.extra)

  cache <- make.cache(tree)
  cache$tip.state  <- states
  cache$sampling.f <- sampling.f
  cache$unresolved <- unresolved # would need more here

  cache$y <- initial.tip.gse2(cache)
  cache
}

## By the time this hits, unresolved clades and any other non-standard
## tips have been removed.  We have an index "tips" (equal to 1:n.tip
## for plain bisse/gse2) that is the "index" (in phy$edge numbering) of the
## tips, and a state vector cache$tip.state, both of the same length.
## The length of the terminal branches is cache$len[cache$tips].
##
## Initial conditions at the tips are given by their tip states:
## There are four types of initial condition in gse2:
##   state0: c(f_0, 0,   0,    1-f_0, 1-f_1, 1-f_2)
##   state1: c(0,   f_1, 0,    1-f_0, 1-f_1, 1-f_2)
##   state2: c(0,   0,   f_2,  1-f_0, 1-f_1, 1-f_2)
##   state?: c(f_0, f_1, f_2,  1-f_0, 1-f_1, 1-f_2)
initial.tip.gse2 <- function(cache) {
  f <- cache$sampling.f
  y <- list(c(1-f, f[1], 0, 0),
            c(1-f, 0, f[2], 0),
            c(1-f, 0, 0, f[3]),
            c(1-f, f))
  y.i <- cache$tip.state + 1
  y.i[is.na(y.i)] <- 4

  tips <- cache$tips

  ## This would return a data structure appropriate for the more basic
  ## tip treatment:
  ##   dt.tips.ordered(y[y.i], tips, cache$len[tips])
  dt.tips.grouped(y, y.i, tips, cache$len[tips])
}

## 6: ll (note: condition.surv=TRUE in bisse)
ll.gse2 <- function(cache, pars, branches, prior=NULL,
                     condition.surv=FALSE, root=ROOT.OBS, root.p=NULL,
                     intermediates=FALSE) {
if ( !is.null(prior) )
    stop("'prior' argument to likelihood function no longer accepted")
 if ( length(pars) != 7 )
    stop("Invalid parameter length (expected 7)")
  if ( any(pars < 0) || any(!is.finite(pars)))
    return(-Inf)

  if ( !is.null(root.p) && root != ROOT.GIVEN )
    warning("Ignoring specified root state")

  # would need something here for unresolved

  ll.xxsse.gse2(pars, cache, initial.conditions.gse2, branches,
           condition.surv, root, root.p, intermediates)
}

## 7: initial.conditions:
initial.conditions.gse2 <- function(init, pars, t, is.root=FALSE) {
  # E.0, E.1, E.2
  e <- init[[1]][c(1,2,3)]

  # D.1, D.2  (Eq. 6bc)
  d12 <- init[[1]][c(5,6)] * init[[2]][c(5,6)] * pars[c(1,2)]

  # D.0 (Eq. 6a)
  d0 <- 0.5 * sum(init[[1]][c(4,5)] * init[[2]][c(5,4)] * pars[1] + 
                  init[[1]][c(4,6)] * init[[2]][c(6,4)] * pars[2] +
                  init[[1]][c(5,6)] * init[[2]][c(6,5)] * pars[3])
  d <- c(d0, d12)

  c(e, d)
}

## 8: branches
make.branches.gse2 <- function(safe=FALSE) {
  RTOL <- ATOL <- 1e-8

  gse2.ode <- make.ode("gse2_derivs", "diversitreeGSE", "gse2_initmod", 6, safe)
  branches <- function(y, len, pars, t0)
    t(gse2.ode(y, c(t0, t0+len), pars, rtol=RTOL, atol=ATOL)[-1,-1])
  
  make.branches(branches, 4:6)
}

## 9: branches.unresolved
branches.unresolved.gse2 <- function(...)
  stop("Cannot yet use unresolved clades with GSE2")

## Additional functions
stationary.freq.gse2 <- function(pars) {
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

# dunno if this is much better than nothing... should come up with something more sensible
starting.point.gse2 <- function(tree, q.div=5, yule=FALSE) {
  ## RGF: Use qs estimated from Mk2?  Can be slow is the only reason
  ## I have not set this up by default.
  ## find.mle(constrain(make.mk2(phy, phy$tip.state), q10 ~ q01), .1)$par
  pars.bd <- suppressWarnings(starting.point.bd(tree, yule))
  if  ( pars.bd[1] > pars.bd[2] )
    p <- rep(c(pars.bd, (pars.bd[1] - pars.bd[2]) / q.div), each=2)
  else
    p <- rep(c(pars.bd, pars.bd[1] / q.div), each=2)
  p <- c(p[1], p)
  names(p) <- argnames.gse2(NULL)
  p
}

# For GSE2, think about what a sensible set of default models would be.
# all.models.gse2 <- function(f, p, ...) { ... }

# check.unresolved would go here
# mle and anova stuff is covered generically in mle.R

# NOTE:
# Functions below are taken from diversitree-branches.R.  Internal
# modifications were necessary, but the parent functions could likely be
# generalized with the aid of classes.

# modified from diversitree-branches.R: root.p.xxsse()
#   allows ROOT.EQUI for gse2
#   returned p is always a vector of length k (or NULL)
root.p.gse2 <- function(vals, pars, root, root.p=NULL) {
  k <- length(vals) / 2
  d.root <- vals[-seq_len(k)]

  if ( root == ROOT.FLAT )
    p <- rep(1/k, k)
  else if ( root == ROOT.EQUI )
    p <- stationary.freq.gse2(pars)
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

# modified from diversitree-branches.R: root.xxsse()
#   the only difference is lambda
root.gse2 <- function(vals, pars, lq, condition.surv, root.p) {
  logcomp <- sum(lq)

  k <- length(vals) / 2
  i <- seq_len(k)

  # note: AB species are subject to all three speciation rates
  lambda <- c(sum(pars[1:3]), pars[1:2])
  e.root <- vals[i]
  d.root <- vals[-i]

  if ( condition.surv )
    d.root <- d.root / (lambda * (1-e.root)^2)

  if ( is.null(root.p) ) # ROOT.BOTH
    loglik <- log(d.root) + logcomp
  else
    loglik <- log(sum(root.p * d.root)) + logcomp
  loglik
}

# modified from diversitree-branches.R: ll.xxsse()
#   only difference is names of root function calls (the above functions)
ll.xxsse.gse2 <- function(pars, cache, initial.conditions,
                     branches, condition.surv, root, root.p,
                     intermediates) {
  ans <- all.branches(pars, cache, initial.conditions, branches)
  vals <- ans$init[[cache$root]]
  root.p <- root.p.gse2(vals, pars, root, root.p)
  loglik <- root.gse2(vals, pars, ans$lq, condition.surv, root.p)
  ans$root.p <- root.p
  cleanup(loglik, pars, intermediates, cache, ans)
}
