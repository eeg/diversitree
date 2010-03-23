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
##   9. branches.unresolved

## 1: make
make.gse2 <- function(tree, states, unresolved=NULL, sampling.f=NULL,
                       nt.extra=10, safe=FALSE) {
  cache <- make.cache.gse2(tree, states, unresolved=unresolved,
                            sampling.f=sampling.f, nt.extra=10)
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
    method <- "optim"
  NextMethod("find.mle", method=method, class.append="fit.mle.gse2")
}

## Make requires the usual functions:
## 5: make.cache (initial.tip, root)
# almost identical to make.cache.bisse(), but uses three states
make.cache.gse2 <- function(tree, states, unresolved=NULL,
                             sampling.f=NULL, nt.extra=10) {
  if ( !inherits(tree, "phylo") )
    stop("'tree' must be a valid phylo tree")
  if ( is.null(names(states)) )
    stop("The states vector must contain names")

  ## Check 'sampling.f'
  if ( !is.null(sampling.f) && !is.null(unresolved) )
    stop("Cannot specify both sampling.f and unresolved")
  else if ( is.null(sampling.f) )
    sampling.f <- c(1, 1, 1)
  else if ( length(sampling.f) != 3 )
    stop("sampling.f must be of length 3 (or NULL)")
  else if ( max(sampling.f) > 1 || min(sampling.f) <= 0 )
    stop("sampling.f must be on range (0,1]")

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

## 6: ll
ll.gse2 <- function(cache, pars, branches, prior=NULL,
                     condition.surv=FALSE, root=ROOT.OBS, root.p=NULL,
                     intermediates=FALSE) {
 if ( length(pars) != 7 )
    stop("Invalid parameter length (expected 7)")
  if ( any(pars < 0) || any(!is.finite(pars)))
    return(-Inf)

  if ( !is.null(root.p) && root != ROOT.GIVEN )
    warning("Ignoring specified root state")

  gse2.ll(pars, cache, initial.conditions.gse2,
           branches, branches.unresolved.gse2,
           condition.surv, root, root.p,
           prior, intermediates)
}
# gse2.ll() was here

## 7: initial.conditions:
initial.conditions.gse2 <- function(init, pars, t, is.root=FALSE) {
  # E.0, E.1, E.2
  e <- init[1, c(1,2,3)]

  # D.1, D.2  (Eq. 6bc)
  d12 <- init[1, c(5,6)] * init[2, c(5,6)] * pars[c(1,2)]

  # D.0 (Eq. 6a)
  d0 <- 0.5 * sum(init[1, c(4,5)] * init[2, c(5,4)] * pars[1] + 
                  init[1, c(4,6)] * init[2, c(6,4)] * pars[2] +
                  init[1, c(5,6)] * init[2, c(6,5)] * pars[3])
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

# dunno if this is much better than nothing...
# TODO: come up with something more sensible
starting.point.gse2 <- function(tree, q.div=5, yule=FALSE) {
  ## TODO: Use qs estimated from Mk2?  Can be slow is the only reason
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

# mle and anova stuff that used to be here is now covered generically in mle.R

# modified from diversitree-branches.R: xxsse.ll()
gse2.ll <- function(pars, cache, initial.conditions,
                     branches, branches.unresolved, 
                     condition.surv, root, root.p,
                     prior, intermediates) {
  ans <- all.branches(pars, cache, initial.conditions,
                      branches, branches.unresolved)

  vals <- ans$init[cache$root,]
  root.p <- root.p.gse2(vals, pars, root, root.p)
  loglik <- root.gse2(vals, pars, ans$lq, condition.surv, root.p)
  cleanup(loglik, pars, prior, intermediates, cache, ans)
}

# modified from diversitree-branches.R: root.xxsse()
root.gse2 <- function(vals, pars, lq, condition.surv, root.p) {
  logcomp <- sum(lq)

  k <- 3 # number of states
  i <- seq_len(k)
  lambda <- c(pars[1:2], sum(pars[1:3]))
  e.root <- vals[i]
  d.root <- vals[-i]

  # root.mode stuff moved to root.p.gse2()
  
  if ( condition.surv )
    d.root <- d.root / (lambda * (1-e.root)^2)

  if ( is.null(root.p) ) # ROOT.BOTH
    loglik <- log(d.root) + logcomp
  else
    loglik <- log(sum(root.p * d.root)) + logcomp
  loglik
}

# unlike root.p.xxsse(), returned p is always a vector of length k
root.p.gse2 <- function(vals, pars, root, root.p=NULL) {
  k <- 3 # number of states
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
