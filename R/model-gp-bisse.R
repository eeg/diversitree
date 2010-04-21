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
# list with these elements:
#     lambda = 2x2x2 array of speciation rates
#     mu = length 2 vector of extinction rates
#     q = 2x2 array of transition rates
#     nstates = as.integer(2)

## 1: make
make.gpbisse <- function(tree, states, unresolved=NULL, sampling.f=NULL,
                       nt.extra=10, safe=FALSE) {
  cache <- make.cache.gpbisse(tree, states, unresolved=unresolved,
                            sampling.f=sampling.f, nt.extra=10)
  branches <- make.branches.gpbisse(safe)
  ll <- function(pars, ...) ll.gpbisse(cache, pars, branches, ...)
  class(ll) <- c("gpbisse", "function")
  ll
}

## 2: print
print.gpbisse <- function(x, ...) {
  cat("GP-BiSSE likelihood function:\n")
  print(unclass(x))
}

## 3: argnames / argnames<-
# only for when the parameter structure must be flattened
argnames.gpbisse <- function(x, ...) {
  ret <- attr(x, "argnames")
  if ( is.null(ret) )
  {
    k <- 2
    lambda.names <- sprintf("lam%d%d%d", rep(1:k, each=k*(k+1)/2), 
            rep(rep(1:k, times=seq(k,1,-1)), k), 
            unlist(lapply(1:k, function(i) i:k)))
    mu.names <- sprintf("mu%d", seq_len(k))
    q.names <- sprintf("q%d%d", rep(1:k, each=k-1),
            unlist(lapply(1:k, function(i) (1:k)[-i])))
    c(lambda.names, mu.names, q.names)
  }
  else
    ret
}
`argnames<-.gpbisse` <- function(x, value) {
  k <- 2
  if ( length(value) != k*k*(k+1)/2 + k + k*(k-1))
    stop("Invalid names length")
  attr(x, "argnames") <- value
  x  
}

## 4: find.mle
find.mle.gpbisse <- function(func, x.init, method,
                           fail.value=NA, ...) {
  if ( missing(method) )
    method <- "optim"
  NextMethod("find.mle", method=method, class.append="fit.mle.gpbisse")
}

## Make requires the usual functions:
## 5: make.cache (initial.tip, root)
# same as make.cache.bisse() but disallows unresolved clades
make.cache.gpbisse <- function(tree, states, unresolved=NULL,
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
    sampling.f <- c(1, 1)
  else if ( length(sampling.f) != 2 )
    stop("sampling.f must be of length 2 (or NULL)")
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
  cache$tip.state  <- states
  cache$unresolved <- unresolved
  cache$sampling.f <- sampling.f
  cache$nt.extra   <- nt.extra
  cache$y <- initial.tip.bisse(cache)
  cache
}

## Initial conditions at the tips are given by their tip states:
## There are three types of initial condition in bisse:
##   state0: c(f_0, 0,   1-f_0, 1-f_1)
##   state1: c(0,   f_1, 1-f_0, 1-f_1)
##   state?: c(f_0, f_1, 1-f_0, 1-f_1)
## Build this small matrix of possible initial conditions (y), then
## work out which of the three types things are (i).  Note that y[i,]
## gives the full initial conditions.  The element 'types' contains
## the different possible conditions.

# initial.tip.bisse() can be used for gpbisse

## 6: ll
ll.gpbisse <- function(cache, pars, branches, prior=NULL,
                     condition.surv=TRUE, root=ROOT.OBS, root.p=NULL,
                     intermediates=FALSE,
                     root.p0=NA, root.p1=NA) {
    if (class(pars) != "list")
        # todo: lots more error checking somewhere
        stop("Invalid parameter structure (expecting a list)")
    if (any(unlist(lapply(pars, function(x) any(x < 0)))) || 
        any(unlist(lapply(pars, function(x) !is.finite(x)))))
        return(-Inf)
#--------------------------------------------------
#   # see argnames above for generalizing to k states
#   if ( length(pars) != 10 )
#     stop("Invalid parameter length (expected 10)")
#   if ( any(pars < 0) || any(!is.finite(pars)) )
#     return(-Inf)
#-------------------------------------------------- 

  if ( !is.na(root.p0) ) {
    warning("root.p0 is deprecated: please use root.p instead")
    root.p <- c(root.p0, 1-root.p0)
  } else if ( !is.na(root.p1) ) {
    warning("root.p1 is deprecated: please use root.p instead")
    root.p <- c(1-root.p1, root.p1)
  }
  if ( !is.null(root.p) &&  root != ROOT.GIVEN )
    warning("Ignoring specified root state")

  gpbisse.ll(pars, cache, initial.conditions.gpbisse,
           branches, branches.unresolved.gpbisse,
           condition.surv, root, root.p,
           prior, intermediates)
}

## 7: initial.conditions:
initial.conditions.gpbisse <- function(init, pars, t, is.root=FALSE)
{
  n <- pars$nstates
  # note: init[1,] = clade N, init[2,] = clade M; init[,1:2]] = E, init[,3:4] = D

  # E.1, E.2
  e <- init[1, c(1,2)]

  # D.1, D.2
  d <- rep(0, n)
  for (i in seq(n))
  {
    lambda <- pars$lambda[i,,]
    for (j in seq(n))
    {
      for (k in seq(j))
      {
        d[i] <- d[i] + 0.5 * lambda[j,k] * 
                       sum(init[1, c(n+j,n+k)] * init[2, c(n+k,n+j)])
      }
    }
  }
  # todo: be more R-like, without loops
#--------------------------------------------------
#   # pars starts with: "lam111" "lam112" "lam122" "lam211" "lam212" "lam222"
#   # todo: can probably be much cleaner for k states
#   d1 <- 0.5 * sum(init[1, c(3,3)] * init[2, c(3,3)] * pars[1] + 
#                   init[1, c(3,4)] * init[2, c(4,3)] * pars[2] +
#                   init[1, c(4,4)] * init[2, c(4,4)] * pars[3])
#   d2 <- 0.5 * sum(init[1, c(3,3)] * init[2, c(3,3)] * pars[4] + 
#                   init[1, c(3,4)] * init[2, c(4,3)] * pars[5] +
#                   init[1, c(4,4)] * init[2, c(4,4)] * pars[6])
#   d <- c(d1, d2)
#-------------------------------------------------- 
  c(e, d)
}

## 8: branches
make.branches.gpbisse <- function(safe=FALSE) {
  RTOL <- ATOL <- 1e-8
  gpbisse.ode <- make.ode("derivs_gp", "diversitreeGP", "initmod_gp", 4, safe)
  branches <- function(y, len, pars, t0)
    t(gpbisse.ode(y, c(t0, t0+len), pars, rtol=RTOL, atol=ATOL)[-1,-1])
  
  make.branches(branches, 3:4)
}

## 9: branches.unresolved
branches.unresolved.gpbisse <- function(...)
  stop("Cannot yet use unresolved clades with GP-SDD")

## Additional functions

# TODO: will need to collapse lambda_ijk into transition matrix
stationary.freq.gpbisse <- function(pars) {
  stop("Cannot yet use stationary frequency with GP-SDD")
}

# TODO: could use bisse's, and pad with zeroes
starting.point.gpbisse <- function(tree, q.div=5, yule=FALSE) {
  stop("Cannot yet use starting.point with GP-SDD")
}

# modified from diversitree-branches.R: xxsse.ll()
gpbisse.ll <- function(pars, cache, initial.conditions,
                     branches, branches.unresolved, 
                     condition.surv, root, root.p,
                     prior, intermediates) {
  ans <- all.branches(pars, cache, initial.conditions,
                      branches, branches.unresolved)

  vals <- ans$init[cache$root,]
  root.p <- root.p.gpbisse(vals, pars, root, root.p)
  loglik <- root.gpbisse(vals, pars, ans$lq, condition.surv, root.p)
  cleanup(loglik, pars, prior, intermediates, cache, ans)
}

# modified from diversitree-branches.R: root.xxsse()
root.gpbisse <- function(vals, pars, lq, condition.surv, root.p) {
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

# only re-defined for stationary.freq.gpbisse(), which doesn't even exist yet
root.p.gpbisse <- function(vals, pars, root, root.p=NULL) {
  k <- pars$nstates
  d.root <- vals[-seq_len(k)]

  if ( root == ROOT.FLAT )
    p <- rep(1/k, k)
  else if ( root == ROOT.EQUI )
    p <- stationary.freq.gpbisse(pars)
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
