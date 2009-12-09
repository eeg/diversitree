RTOL <- ATOL <- 1e-8
eps <- 0
branches.bisse <- function(y, len, pars, t0) {
  ret <- t(bisse.ode(y, c(t0, t0+len), pars, rtol=RTOL, atol=ATOL)[-1,-1])
  if ( all(ret[,3:4] >= eps) ) {
    q <- ret[cbind(seq_along(len), as.integer(ret[,3] > ret[,4]) + 3)]
    i <- q > 0
    ret[i,3:4] <- ret[i,3:4] / q[i]
    lq <- q
    lq[i] <- log(q[i])
    cbind(lq, ret, deparse.level=0)
  } else {
    ti <- len[length(len)]/2
    len1 <- c(len[len <= ti], ti)
    len2 <- len[len > ti] - ti
    n1 <- length(len1)
    ret1 <- Recall(y, len1, pars, t0)
    ret2 <- Recall(ret1[n1,2:5], len2, pars, t0 + ti)
    ret2[,1] <- ret2[,1] + ret1[n1,1]
    rbind(ret1[-n1,], ret2)
  }
}

branches.unresolved.bisse <- function(pars, len, unresolved) {
  Nc <- unresolved$Nc
  k <- unresolved$k
  nsc <- unresolved$nsc
  t <- len
  nt <- max(Nc) + unresolved$nt.extra
  
  lambda0 <- pars[1]
  lambda1 <- pars[2]
  mu0 <- pars[3]
  mu1 <- pars[4]
  q01 <- pars[5]
  q10 <- pars[6]
  ret <- bucexpl(nt, mu0, mu1, lambda0, lambda1, q01, q10, t,
          Nc, nsc, k)[,c(3,4,1,2)]
  q <- ret[cbind(seq_along(len), as.integer(ret[,3] > ret[,4]) + 3)]
  ret[,3:4] <- ret[,3:4] / q
  cbind(log(q), ret, deparse.level=0)
}

initial.conditions.bisse <- function(init, pars, t, is.root=FALSE) {
  e <- init[1,c(1,2)]
  d <- init[1,c(3,4)] * init[2,c(3,4)]
  if ( !is.root )
    d <- d * pars[c(1,2)]
  c(e, d)
}

ROOT.FLAT  <- 1
ROOT.EQUI  <- 2
ROOT.OBS   <- 3
ROOT.GIVEN <- 4
ROOT.BOTH  <- 5

bisse.ll <- function(cache, pars, prior=NULL, root=ROOT.OBS,
                     condition.surv=TRUE, root.p0=NA, root.p1=NA,
                     intermediates=FALSE) {
  if ( any(pars < 0) || any(!is.finite(pars)) )
    return(-Inf)
  ans <- all.branches(pars, cache, initial.conditions.bisse,
                      branches.bisse, branches.unresolved.bisse)
  loglik <- root.bisse(ans$init[cache$root,], ans$lq, pars, root,
                       condition.surv, root.p0, root.p1)
  if ( !is.null(prior) )
    loglik <- loglik + prior.bisse(pars, prior)
  if ( intermediates ) {
    attr(loglik, "cache") <- cache
    attr(loglik, "intermediates") <- ans
    attr(loglik, "vals") <- ans$init[cache$root,]
    attr(loglik, "logComp") <- sum(ans$lq)
  }

  loglik
}


prior.bisse <- function(pars, r)
  - sum(pars * r)

make.bisse <- function(tree, states, unresolved=NULL, sampling.f=NULL,
                       nt.extra=10) {
  cache <- make.cache.bisse(tree, states, unresolved=unresolved,
                            sampling.f=sampling.f, nt.extra=10)
  f <- function(pars, ..., fail.value = NULL) {
    if (!is.null(fail.value)) 
      protect(bisse.ll, fail.value)(cache, pars, ...)
    else bisse.ll(cache, pars, ...)
  }
  class(f) <- c("bisse", "function")
  f
}

## BiSSE root calculations; this gets reused by a few things.
root.bisse <- function(vars, lq, pars, root, condition.surv, root.p0=NA,
                       root.p1=NA) {
  if ( !is.na(root.p0) || !is.na(root.p1) )
    if ( root != ROOT.GIVEN )
      warning("Ignoring specified root state")
  if ( !is.na(root.p1) )
    root.p0 <- 1 - root.p1
  
  e.root <- vars[c(1,2)]
  d.root <- vars[c(3,4)]

  if ( root == ROOT.FLAT )
    p <- 0.5
  else if ( root == ROOT.EQUI )
    p <- bisse.stationary.freq(pars)
  else if ( root == ROOT.OBS )
    p <- d.root[1]/sum(d.root)
  else if ( root == ROOT.GIVEN )
    p <- root.p0
  else if ( root != ROOT.BOTH )
    stop("Invalid root mode")

  if ( condition.surv )
    d.root <- d.root / (1-e.root)^2
  logcomp <- sum(lq)
  if ( root == ROOT.BOTH )
    loglik <- log(d.root)# + logcomp
  else
    loglik <- log(sum(c(p, 1-p) * d.root)) + logcomp
  loglik
}

bisse.stationary.freq <- function(pars) {
  lambda0 <- pars[1]
  lambda1 <- pars[2]
  mu0 <- pars[3]
  mu1 <- pars[4]
  q01 <- pars[5]
  q10 <- pars[6]

  g <- (lambda0 - mu0) - (lambda1 - mu1)
  eps <- (lambda0 + mu0 + lambda1 + mu1)*1e-14
  if ( abs(g) < eps ) {
    if ( q01 + q10 == 0 )
      0.5
    else
      q10/(q01 + q10)
  } else {
    roots <- quadratic.roots(g, q10+q01-g, -q10)
    roots <- roots[roots >= 0 & roots <= 1]
    if ( length(roots) > 1 )
      NA
    else
      roots
  }
}

quadratic.roots <- function(a, b, c)
  (-b + c(-1, 1) * sqrt(b*b - 4*a*c))/(2 * a)

starting.point.bisse <- function(tree, q.div=5) {
 fit <- suppressWarnings(birthdeath(tree))
 r <- fit$para[2]
 e <- fit$para[1]
 p <- rep(c(r/(1-e), r*e/(1-e), r/q.div), each=2)
 names(p) <- c("lambda0", "lambda1", "mu0", "mu1", "q01", "q10")
 p
}
bisse.starting.point <- function(tree, q.div=5) {
  .Deprecated("starting.point.bisse")
  starting.point.bisse(tree, q.div)
}
starting.point <- bisse.starting.point <- function(tree, q.div=5) {
  .Deprecated("starting.point.bisse")
  starting.point.bisse(tree, q.div)
}

make.cache.bisse <- function(tree, states, unresolved=NULL,
                             sampling.f=NULL, nt.extra=10) {
  if ( is.null(names(states)) )
    stop("The states vector must contain names")

  ## Check 'sampling.f'
  if ( !is.null(sampling.f) && !is.null(unresolved) )
    stop("Cannot specify both sampling.f and unresolved")
  else if ( is.null(sampling.f) )
    sampling.f <- c(1, 1)
  else if ( length(sampling.f) != 2 )
    stop("sampling.f must be of length 2 (or NULL)")
  else if ( max(sampling.f) > 1 || min(sampling.f) < 0 )
    stop("sampling.f must be on range [0,1]")
  
  ## Check 'unresolved' (there is certainly room to streamline this in
  ## the future).
  if ( !is.null(unresolved) && nrow(unresolved) == 0 ) {
    unresolved <- NULL
    warning("Ignoring empty 'unresolved' argument")
  }
  if ( inherits(tree, "clade.tree") ) {
    if ( !is.null(unresolved) )
      stop("'unresolved' cannot be specified where 'tree' is a clade.tree")
    unresolved <- make.unresolved(tree$clades, states)
  }
  if ( !is.null(unresolved) ) {
    required <- c("tip.label", "Nc", "n0", "n1")
    if ( !all(required %in% names(unresolved)) )
      stop("Required columns missing from unresolved clades")
    if ( !all(unresolved$tip.label %in% tree$tip.label) )
      stop("Unknown tip species in 'unresolved'")
    unresolved$k   <- unresolved$n1
    unresolved$nsc <- unresolved$n0 + unresolved$n1
    unresolved$i   <- match(unresolved$tip.label, tree$tip.label)
    unresolved <- as.list(unresolved)
    unresolved$nt.extra <- nt.extra
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
initial.tip.bisse <- function(cache) {
  f <- cache$sampling.f
  y <- rbind(c(1-f, f[1], 0),
             c(1-f, 0, f[2]),
             c(1-f, f))
  i <- cache$tip.state + 1
  if ( !is.null(cache$unresolved) )
    i <- i[-cache$unresolved$i]
  i[is.na(i)] <- 3
  list(y=y, i=i, types=sort(unique(i)))
}

## This is here for reference, but not exported yet.  It should be
## tweaked in several ways
##   1. Starting parameter guessing should be done internally, at
##      least as an option.
##   2. Better listing of arguments
##   3. Automatic parsing of results into some sort of table; this
##      proabably requires classing this.
all.models.bisse <- function(f, p, ...) {
  f3 <- constrain(f, lambda1 ~ lambda0, mu1 ~ mu0, q01 ~ q10)
  f4.lm <- constrain(f, lambda1 ~ lambda0, mu1 ~ mu0)
  f4.lq <- constrain(f, lambda1 ~ lambda0, q01 ~ q10)
  f4.mq <- constrain(f, mu1 ~ mu0, q01 ~ q10)
  f5.l <- constrain(f, lambda1 ~ lambda0)
  f5.m <- constrain(f, mu1 ~ mu0)
  f5.q <- constrain(f, q01 ~ q10)

  ## Fit six and three parameter models
  if ( length(p) != 3 )
    stop("Starting point must be of length 3 (lambda, mu, q)")
  ans3 <- find.mle(f3, p, ...)

  ## Using the values from the 3p model, fit the 4p and 5p models:
  l <- ans3$par[1]
  m <- ans3$par[2]
  q <- ans3$par[3]

  ## Start the searches from the best model of the previous type.
  ans4.lm <- find.mle(f4.lm, c(l, m, q, q), ...)
  ans4.lq <- find.mle(f4.lq, c(l, m, m, q), ...)
  ans4.mq <- find.mle(f4.mq, c(l, l, m, q), ...)

  p.l <- if ( ans4.lm$lnLik > ans4.lq$lnLik )
    ans4.lm$par[c(1:2,2:4)] else ans4.lq$par[c(1:4,4)]
  p.m <- if ( ans4.lm$lnLik > ans4.mq$lnLik )
    ans4.lm$par[c(1,1:4)] else ans4.mq$par[c(1:4,4)]
  p.q <- if ( ans4.lq$lnLik > ans4.mq$lnLik )
    ans4.lq$par[c(1,1:4)] else ans4.mq$par[c(1:3,3:4)]
  ans5.l  <- find.mle(f5.l, p.l, ...)
  ans5.m  <- find.mle(f5.m, p.m, ...)
  ans5.q  <- find.mle(f5.q, p.q, ...)

  tmp <- list(ans5.l, ans5.m, ans5.q)
  i <- which.max(sapply(tmp, "[[", "lnLik"))
  p6 <- tmp[[i]]$par
  j <- list(c(1, 1:5), c(1:2, 2:5), c(1:5, 5))
  ans6 <- find.mle(f, p6[j[[i]]], ...)

  list(ans6=ans6,
       ans5.l =ans5.l,  ans5.m =ans5.m,  ans5.q =ans5.q,
       ans4.lm=ans4.lm, ans4.lq=ans4.lq, ans4.mq=ans4.mq,
       ans3=ans3)
}
