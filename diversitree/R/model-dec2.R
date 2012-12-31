## Simple case of the Ree & Smith DEC model.
## Two regions.  Focus on rate estimation, not ancestral state reconstruction.  
## Like GeoSSE but without speciation and global extinction rates to be estimated.
## Like Mkn (k=3) but allowing state changes at nodes.

## States:
##   0 AB
##   1 A
##   2 B
## Parameters:
##   xA 0 -> 2
##   xB 0 -> 1
##   dA 1 -> 0
##   dB 2 -> 0

## 1: make
make.dec2 <- function(tree, states, strict=TRUE, control=list()) {
  control <- check.control.dec2(control)
  cache <- make.cache.dec2(tree, states, strict, control)
  all.branches <- make.all.branches.mkn(cache, control)
  rootfunc <- rootfunc.dec2
  f.pars <- make.pars.dec2()

  ll <- function(pars, root=ROOT.OBS, root.p=NULL, intermediates=FALSE) {
    qmat <- f.pars(pars)
    ans <- all.branches(qmat, intermediates)
    rootfunc(ans, qmat, root, root.p, intermediates)
  }
  class(ll) <- c("dec2", "dtlik", "function")
  ll
}

## 2: info
make.info.dec2 <- function(phy) {
  k <- 3
  list(name="dec2",
       name.pretty="DEC(2)",
       ## Parameters:
       np=as.integer(k * k),
       argnames=default.argnames.dec2(),
       ## Variables:
       ny=as.integer(k), # TODO/NEW: only for ode version...
       k=as.integer(k),
       idx.e=integer(0),
       idx.d=seq_len(k),
       ## Phylogeny:
       phy=phy,
       ## Inference:
       ml.default="subplex",
       mcmc.lowerzero=TRUE,
       ## These are optional
       doc=NULL,
       reference=c(
         "Ree et al. (2005)",
         "Ree & Smith (2008)"))
}
default.argnames.dec2 <- function() {
  c("xA", "xB", "dA", "dB")
}

## 3: make.cache (& initial.tip)
make.cache.dec2 <- function(tree, states, strict, control) {
  method <- control$method

  tree <- check.tree(tree)
  if ( !is.null(states) ) # for multitrait
    states <- check.states(tree, states, strict=strict,
                           strict.vals=0:2)
  cache <- make.cache(tree)
  cache$info <- make.info.dec2(tree)
  cache$states  <- states
  if ( method == "ode" ) {
    cache$y <- initial.tip.dec2.ode(cache)
    cache$info$name.ode <- "mknode"
  }

  cache
}

# mostly initial.tip.mkn.ode()
initial.tip.dec2.ode <- function(cache) {
  k <- cache$info$k
  y <- matrix(0, k+1, k, TRUE)
  y[k+1,] <- diag(y[1:k,]) <- 1
  y <- matrix.to.list(y)
  y.i <- cache$states
  y.i <- y.i + 1L    # new: see base.zero in initial.tip.xxsse()
  dt.tips.grouped(y, y.i, cache)
}

## 4: initial conditions
initial.conditions.dec2 <- function(init, pars, t, idx) {
  # in init, row = state, col = daughter

  # 0 -> 0 + 1
  # 0 -> 0 + 2
  # 0 -> 1 + 2
  d0 <- 0.5 * sum(
          init[c(1,2), 1] * init[c(2,1), 2] +
          init[c(1,3), 1] * init[c(3,1), 2] +
          init[c(2,3), 1] * init[c(3,2), 2]
  )

  # 1 -> 1 + 1
  # 2 -> 2 + 2
  d12 <- init[2:3, 1] * init[2:3, 2]

  c(d0, d12)
}

## 5: rootfunc
rootfunc.dec2 <- function(res, pars, root, root.p, intermediates) {
  d.root <- res$vals
  lq <- res$lq
  k <- length(d.root)

  # root.p <- root.p.calc(d.root, pars, root, root.p,
  #                       stationary.freq.dec2)
  root.p <- root.p.calc(d.root, pars, root, root.p)
  if ( root == ROOT.ALL )
    loglik <- log(d.root) + sum(lq)
  else
    loglik <- log(sum(root.p * d.root)) + sum(lq)

  if ( intermediates ) {
    res$root.p <- root.p
    attr(loglik, "intermediates") <- res
    attr(loglik, "vals") <- d.root
  }

  loglik
}

make.all.branches.dec2 <- function(cache, control) {
  if ( control$method == "ode" ) {
    if ( !is.null(control$backend) && control$backend == "expokit" )
      make.all.branches.mkn.expokit(cache, control)
    else
      make.all.branches.dtlik(cache, control, initial.conditions.mkn) # currently only allow this option
  } else {
    make.all.branches.mkn.exp(cache, control)
  }
}

######################################################################
## Additional functions:

## Parameter manipulation:
## Makes a function that converts the parameter vector into a k^2 Q
## matrix.
make.pars.dec2 <- function() {
  k <- 3
  qmat <- matrix(0, k, k)
  #             xA,      xB,      dA,      dB
  idx <- rbind( c(1, 3), c(1, 2), c(2, 1), c(3, 1) )
  npar <- 4

  function(pars) {
    check.pars.nonnegative(pars, npar)
    qmat[idx] <- pars
    diag(qmat) <- -rowSums(qmat)
    qmat
  }
}

## Checking:
check.control.dec2 <- function(control) {
  # only allowing "ode" for now; uses initial.conditions.dec2() above
  control <- modifyList(list(method="ode"), control)
  # methods <- c("exp", "ode")
  methods <- c("ode")
  if ( !(control$method %in% methods) )
    stop(sprintf("control$method must be in %s",
                 paste(methods, collapse=", ")))
  control
}

##   xA 0 -> 2
##   xB 0 -> 1
##   dA 1 -> 0
##   dB 2 -> 0
