make.all.branches.mkn.pij <- function(cache, control) {
  k <- cache$info$k
  
  f.pij <- make.pij.mkn(cache$info, control)

  n.tip <- cache$n.tip
  n <- length(cache$len)

  ## TODO: deal with unknown states; tip calculations for unknown tip
  ## states just return the sum of nonzero ys from the matrix
  ## multiplication, so they can just be done separately.
  map <- t(sapply(1:k, function(i) (1:k) + (i - 1) * k))
  if (!all(is.na(cache$states))) {
    idx.tip <- cbind(c(map[cache$states,]),
                     rep(seq_len(n.tip), k))
  } else {
    idx.tip <- cbind(NA,
                     rep(seq_len(n.tip), k))
  }
  len.uniq <- sort(unique(cache$len))
  len.idx <- match(cache$len, len.uniq)

  ## Alter things to make it more speedy.  The '.C' denotes C-style
  ## base-0 indices.
  children.C <- toC.int(t(cache$children))
  order.C    <- toC.int(cache$order)
  
  ## At this point, the parameters are assumed to be a Q matrix
  function(pars, intermediates, preset=NULL) {
    if ( !is.null(preset) )
      stop("Preset values not allowed")
    pij <- f.pij(len.uniq, pars)[,len.idx]

    lq <- numeric(n)
    branch.init <- branch.base <- matrix(NA, k, n)
    storage.mode(branch.init) <- "numeric"

    ## tips
    ans <- matrix(pij[idx.tip], n.tip, k)
    if (any(is.na(cache$states)))
    {
      ## now fix up ans for the multistate/unknown tips
      ## (this is unlikely to be the most elegant solution)
      wts <- do.call("rbind", attr(cache$states, "multistate")$states)[names(cache$states)[is.na(cache$states)],]
      if (!(all(rowSums(wts) == 1)))
          stop("Expecting weights for uncertain tips to sum to 1.")
      i.na <- which(is.na(cache$states))
      ans[i.na,] <- 0
      st <- cache$states
      for (i.k in seq(k))
      {
        st[i.na] <- i.k
        idx.tip1 <- cbind(c(map[st,]), rep(seq_len(n.tip), k))
        ans1 <- matrix(pij[idx.tip1], n.tip, k)
        ans[i.na,] <- ans[i.na,] + 
                      matrix(rep(wts[,i.k], k), ncol=k) * ans1[i.na,]
      }
    }

    q <- rowSums(ans)
    branch.base[,seq_len(n.tip)] <- t.default(ans/q)
    lq[seq_len(n.tip)] <- log(q)

    ans <- .C("r_mkn_core",
              k        = as.integer(k),
              n        = length(order.C) - 1L,
              order    = order.C,
              children = children.C,
              pij      = pij,
              init     = branch.init,
              base     = branch.base,
              lq       = lq,
              NAOK=TRUE, package="diversitree")
    list(init=ans$init,
         base=ans$base,
         lq=ans$lq,
         vals=ans$init[,cache$root],
         pij=pij)
  }
}

######################################################################
## Mkn-special stuff.
## The calculations here are quite different to the rest of the
## package.

## Compute Pij matrices:
pij.mk2 <- function(len, pars) {
  ## The 3,2 indices here are because pars is a Q matrix by the time
  ## this gets called.
  q01 <- pars[3]
  q10 <- pars[2]
  x <- exp(-(q01+q10)*len)
  rbind((x*q01 + q10),
        (1 - x)*q10,
        (1 - x)*q01,
        (x*q10 + q01)) / (q01 + q10)
}

make.pij.mkn <- function(info, control) {
  if ( control$method == "mk2" )
    return(pij.mk2)
  control <- check.control.ode(control)

  k <- info$k

  ## Make a brand new info list here (a bit ugly)
  info <- list(name="mknpij",
               ny=k*k, np=k*k,
               idx.d=integer(0),
               derivs=make.derivs.mkn.pij(k))
  pij.ode <- make.ode(info, control)
  
  yi <- diag(k) # initial conditions always same.

  function(len, pars)
    pij.ode(yi, c(0, len), pars)
}

make.derivs.mkn.pij <- function(k) {
  force(k)
  function(t, y, pars) {
    Q <- matrix(pars, k, k)
    y <- matrix(y, k, k)
    ret <- Q %*% y
    dim(ret) <- NULL
    ret
  }
}
