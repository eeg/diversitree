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

# The parameter structure used here is no longer a vector, but a 
# list with these elements (k states):
#     lambda = k x k x k array of speciation rates
#     mu = length k vector of extinction rates
#     q = k x k array of transition rates
#     nstates = as.integer(k)
# Functions to convert between param list and vector are included.

## 1: make
make.gpsdd <- function(tree, states, k, unresolved=NULL, sampling.f=NULL,
                       nt.extra=10, safe=FALSE) {
  cache <- make.cache.gpsdd(tree, states, k, unresolved=unresolved,
                            sampling.f=sampling.f, nt.extra=10)
  branches <- make.branches.gpsdd(k, safe)
  ll <- function(pars, ...) ll.gpsdd(cache, pars, branches, ...)
  class(ll) <- c("gpsdd", "function")
  attr(ll, "k") <- k
  ll
}

## 2: print
print.gpsdd <- function(x, ...) {
  cat("GP-SDD likelihood function:\n")
  print(unclass(x))
}

## 3: argnames / argnames<-
# only for when the parameter structure must be flattened
# may want to change for double-digit states
argnames.gpsdd <- function(x, k, ...) {
  ret <- attr(x, "argnames")
  if ( is.null(ret) ) {
    if ( missing(k) )
      k <- attr(x, "k")
    else
      if ( !is.null(x) )
        stop("k can only be given if x is null")
    lambda.names <- sprintf("lam%d%d%d", rep(1:k, each=k*(k+1)/2), 
            rep(rep(1:k, times=seq(k,1,-1)), k), 
            unlist(lapply(1:k, function(i) i:k)))
    mu.names <- sprintf("mu%d", seq_len(k))
    q.names <- sprintf("q%d%d", rep(1:k, each=k-1),
            unlist(lapply(1:k, function(i) (1:k)[-i])))
    c(lambda.names, mu.names, q.names)
  } else {
    ret
  }
}
`argnames<-.gpsdd` <- function(x, value) {
  k <- environment(x)$cache$k
  if ( length(value) != k*k*(k+1)/2 + k + k*(k-1))
    stop("Invalid names length")
  attr(x, "argnames") <- value
  x  
}

## convert between list of parameter arrays and single vector

# recycle logic from argnames
# tried omitting elements from unlist(parlist), but that's no cleaner
flatten.pars.gpsdd <- function(parlist)
{
  k <- parlist$nstates

  i.lam <- cbind(rep(1:k, each=k*(k+1)/2), rep(rep(1:k, times=seq(k,1,-1)), k),
                 unlist(lapply(1:k, function(i) i:k)))

  i.q <- cbind(rep(1:k, each=k-1), unlist(lapply(1:k, function(i) (1:k)[-i])))

  parvec <- c(parlist$lambda[i.lam], parlist$mu, parlist$q[i.q])
  names(parvec) <- diversitreeGP:::argnames.gpsdd(NULL, k)
  return(parvec)
}

# reverse the logic in flatten.pars.gpsdd
# unused elements in lambda and q are assigned 0 (would be nice if 
#   they could be NA, but the C code needs help for that)
listify.pars.gpsdd <- function(parvec, k)
{
  if ( length(parvec) != k*k*(k+1)/2 + k + k*(k-1))
    stop("Invalid length of parameter vector.")

  Lam <- array(0, dim=rep(k, 3))  # 3 = parent + 2 daughters
  idx <- cbind(rep(1:k, each=k*(k+1)/2), rep(rep(1:k, times=seq(k,1,-1)), k),
               unlist(lapply(1:k, function(i) i:k)))
  j <- length(idx[,1])
  Lam[idx] <- parvec[seq(j)]

  Mu <- parvec[seq(j+1, j+k)]

  Q <- array(0, dim=rep(k, 2))
  idx <- cbind(rep(1:k, each=k-1), unlist(lapply(1:k, function(i) (1:k)[-i])))
  Q[idx] <- parvec[seq(j+k+1, length(parvec))]

  list(lambda=Lam, mu=Mu, q=Q, nstates=as.integer(k))
}

## 4: find.mle
find.mle.gpsdd <- function(func, x.init, method,
                           fail.value=NA, ...) {
  if ( missing(method) )
    method <- "optim"
  NextMethod("find.mle", method=method, class.append="fit.mle.gpsdd")
}

## Make requires the usual functions:
## 5: make.cache (initial.tip, root)
# same as make.cache.bisse() but disallows unresolved clades
make.cache.gpsdd <- function(tree, states, k, unresolved=NULL,
                             sampling.f=NULL, nt.extra=10) {
  if ( !inherits(tree, "phylo") )
    stop("'tree' must be a valid phylo tree")
  if ( is.null(names(states)) )
    stop("The states vector must contain names")
  if ( inherits(tree, "clade.tree") ) {
    if ( !is.null(unresolved) )
      stop("'unresolved' cannot be specified where 'tree' is a clade.tree")
    unresolved <- make.unresolved(tree$clades, states)
  }

  ## Check 'sampling.f'
  if ( !is.null(sampling.f) && !is.null(unresolved) )
    stop("Cannot specify both sampling.f and unresolved")
  else if ( is.null(sampling.f) )
    sampling.f <- rep(1, k)
  else if ( length(sampling.f) != k )
    stop("sampling.f must be of length nstates (or NULL)")
  else if ( max(sampling.f) > 1 || min(sampling.f) <= 0 )
    stop("sampling.f must be on range (0,1]")

  if ( !is.null(unresolved) ) {
      stop("Unresolved clades not yet available for GP-SDD")
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
  cache$k <- k
  cache$tip.state  <- states
  cache$unresolved <- unresolved
  cache$sampling.f <- sampling.f
  cache$nt.extra   <- nt.extra
  cache$y <- initial.tip.gpsdd(cache)
  cache
}

## Initial conditions at the tips are given by their tip states:
## There are k + 1 types of initial condition:
##
##              E.1    E.2         E.k      D.1  D.2    D.k
##   state1:  c(1-f_1, 1-f_2, ..., 1-f_k,   f_1, 0, ..., 0)
##   state2:  c(1-f_1, 1-f_2, ..., 1-f_k,   0, f_2, ..., 0)
##   etc.
##   stateNA: c(1-f_1, 1-f_2, ..., 1-f_k,   f_1, f_2, ...f_k)
##
## Build this small-ish matrix of possible initial conditions (y), then
## work out which of the three types things are (i).  Note that y[i,]
## gives the full initial conditions.  The element 'types' contains
## the different possible conditions.
initial.tip.gpsdd <- function(cache) {
  k <- cache$k
  tip.state <- cache$tip.state
  if ( any(tip.state < 1 | tip.state > k, na.rm=TRUE) )
    stop(sprintf("tip states must be in the range [1, %d]", k))

  f <- cache$sampling.f
  y <- matrix(rep(c(1-f, rep(0,k)), k), nrow=k, byrow=T)
  idx <- cbind(seq(k), seq(k)+k)
  y[idx] <- f
  y <- rbind(y, c(1-f, f))

  i <- cache$tip.state
  if ( !is.null(cache$unresolved) )
    i <- i[-cache$unresolved$i]
  i[is.na(i)] <- k + 1
  list(y=y, i=i, types=sort(unique(i)))
}

## 6: ll
ll.gpsdd <- function(cache, pars, branches, prior=NULL,
                     condition.surv=TRUE, root=ROOT.OBS, root.p=NULL,
                     intermediates=FALSE) {
  k = cache$k
  if (class(pars) != "list")
    if (length(pars) == k*k*(k+1)/2 + k + k*(k-1))
      pars <- listify.pars.gpsdd(pars, k)
    else
      stop("Invalid parameter vector")
  # todo: lots more error checking somewhere

  if (any(unlist(lapply(pars, function(x) any(x < 0)))) || 
    any(unlist(lapply(pars, function(x) !is.finite(x)))))
    return(-Inf)

  if ( !is.null(root.p) && root != ROOT.GIVEN )
    warning("Ignoring specified root state")

  gpsdd.ll(pars, cache, initial.conditions.gpsdd,
           branches, branches.unresolved.gpsdd,
           condition.surv, root, root.p,
           prior, intermediates)
}

## 7: initial.conditions:
initial.conditions.gpsdd <- function(init, pars, t, is.root=FALSE)
{
  n <- pars$nstates     # called k elsewhere, but k is used as an index below
  # note: init[1,] = clade N, init[2,] = clade M; init[,1:n]] = E, init[,(n+1):(2*n)] = D

  # E.1, E.2
  e <- init[1, seq(n)]

  # D.1, D.2
  d <- rep(0, n)
  for (i in seq(n)) {
    lambda <- pars$lambda[i,,]
    for (j in seq(n))
      for (k in seq(j))
        d[i] <- d[i] + 0.5 * lambda[j,k] * 
                       sum(init[1, c(n+j,n+k)] * init[2, c(n+k,n+j)])
  }
  # todo: be more R-like, without loops

  c(e, d)
}

## 8: branches
make.branches.gpsdd <- function(k, safe=FALSE) {
  RTOL <- ATOL <- 1e-8
  gpsdd.ode <- make.ode("derivs_gp", "diversitreeGP", "initmod_gp", 2*k, safe)
  branches <- function(y, len, pars, t0)
    t(gpsdd.ode(y, c(t0, t0+len), pars, rtol=RTOL, atol=ATOL)[-1,-1])
  
  make.branches(branches, seq(k+1, 2*k))
}

## 9: branches.unresolved
branches.unresolved.gpsdd <- function(...)
  stop("Cannot yet use unresolved clades with GP-SDD")

## Additional functions

# TODO: will need to collapse lambda_ijk into transition matrix
stationary.freq.gpsdd <- function(pars) {
  stop("Cannot yet use stationary frequency with GP-SDD")
}

# TODO: could use bisse's, and pad with zeroes
starting.point.gpsdd <- function(tree, q.div=5, yule=FALSE) {
  stop("Cannot yet use starting.point with GP-SDD")
}

# modified from diversitree-branches.R: xxsse.ll()
gpsdd.ll <- function(pars, cache, initial.conditions,
                     branches, branches.unresolved, 
                     condition.surv, root, root.p,
                     prior, intermediates) {
  ans <- all.branches(pars, cache, initial.conditions,
                      branches, branches.unresolved)

  vals <- ans$init[cache$root,]
  root.p <- root.p.gpsdd(vals, pars, root, root.p)
  loglik <- root.gpsdd(vals, pars, ans$lq, condition.surv, root.p)
  cleanup(loglik, pars, prior, intermediates, cache, ans)
}

# modified from diversitree-branches.R: root.xxsse()
root.gpsdd <- function(vals, pars, lq, condition.surv, root.p) {
  logcomp <- sum(lq)

  i <- seq_len(pars$nstates)
  e.root <- vals[i]
  d.root <- vals[-i]

  # sum the lambdas for each parent state; any type of speciation 
  #   could happen at the root
  lambda <- apply(pars$lambda, 1, sum)
  if ( condition.surv )
    d.root <- d.root / (lambda * (1-e.root)^2)

  if ( is.null(root.p) ) # ROOT.BOTH
    loglik <- log(d.root) + logcomp
  else
    loglik <- log(sum(root.p * d.root)) + logcomp
  loglik
}

# only re-defined for stationary.freq.gpsdd(), which doesn't even exist yet
root.p.gpsdd <- function(vals, pars, root, root.p=NULL) {
  k <- pars$nstates
  d.root <- vals[-seq_len(k)]

  if ( root == ROOT.FLAT )
    p <- rep(1/k, k)
  else if ( root == ROOT.EQUI )
    p <- stationary.freq.gpsdd(pars)
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
