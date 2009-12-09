### order of states:
# 0 = in both regions
# 1 = endemic to region A
# 2 = endemic to region B

### order of parameters:
# (sA, sB, sAB, xA, xB, dA, dB)


### root state assumptions
ROOT.FLAT  <- 1    # 1/3 for each state
ROOT.EQUI  <- 2    # equilibrium proportions
ROOT.OBS   <- 3    # Rich's D proportions; the default
ROOT.GIVEN <- 4    # user specified constants
ROOT.ALL   <- 5    # 1.0 for each state


### construct the likelihood function object
make.gse2 <- function(tree, states, unresolved=NULL, sampling.f=NULL, ...)
{
    cache <- gse2cache(tree, states, unresolved, sampling.f, ...)

    f <- function(pars, ..., fail.value=NULL)
    {
        if ( !is.null(fail.value) )
            protect(gse2.ll, fail.value)(cache, pars, ...)
        else
            gse2.ll(cache, pars, ...)
    }

    class(f) <- c(class(f), "gse2")   # need a new class for find.mle()
    return(f)
}


gse2.ll <- function(cache, pars, prior=NULL, root=ROOT.OBS,
                     condition.surv=TRUE, intermediates=FALSE, root.p=NA)
{
    if ( any(pars < 0) )
        return(-Inf)

    # most of the likelihood calculation
    ans <- gse2.branches(pars, cache)
    e.root <- ans$init[cache$root,1:3]
    d.root <- ans$init[cache$root,4:6]

    # p = vector of weights for root state (not a scalar like bisse's root.p0)
    if ( root == ROOT.FLAT )
        p <- rep(1/3., 3)
    else if ( root == ROOT.EQUI)
        p <- stationary.freq.gse2(pars)
    else if ( root == ROOT.OBS )
        p <- d.root / sum(d.root)
    else if ( root == ROOT.GIVEN )
        p <- root.p
    else if ( root == ROOT.ALL )
        p <- rep(1., 3)
    else
        stop(paste("Invalid root mode", root))

    # condition on survival
    if ( condition.surv )
        d.root <- d.root / (1-e.root)^2

    # restore the likelihood factors removed during the calculation
    logcomp <- sum(ans$lq)
    loglik <- log(sum(p * d.root)) + logcomp  # okay for ROOT.ALL?

    # apply exponential prior
    if ( !is.null(prior) )
        loglik <- loglik + bisse.prior(pars, prior)

    # report extra info
    if ( intermediates )
    {
        attr(loglik, "cache") <- cache
        attr(loglik, "intermediates") <- ans
        attr(loglik, "vals") <- ans$init[cache$root,] # root E and D
        attr(loglik, "logComp") <- logcomp
    }

    return(loglik)
}


### cache stuff unchanged by parameter values; mostly topology
#      almost identical to bissecache(), but uses three states
#      (could easily be generalized to N states and shared with N-bisse)
gse2cache <- function(tree, states, unresolved=NULL, sampling.f=NULL,
                      nt.extra=10, ...)
{
    if ( !is.null(unresolved) && nrow(unresolved) == 0 )
    {
        unresolved <- NULL
        warning("Ignoring empty 'unresolved' argument")
    }

    if ( inherits(tree, "clade.tree") )
    {
        if ( !is.null(unresolved) )
            stop("'unresolved' cannot be specified where 'tree' is a clade.tree")
        unresolved <- make.unresolved(tree$clades, states)
        states <- states[tree$tip.label]
        names(states) <- tree$tip.label
    }

    if ( !is.null(sampling.f) && !is.null(unresolved) )
        stop("Cannot specify both sampling.f and unresolved")
    else if ( is.null(sampling.f) )
        sampling.f <- c(1, 1, 1)
    else if ( length(sampling.f) != 3 )
        stop("sampling.f must be of length 3 (or NULL)")
    else if ( max(sampling.f) > 1 || min(sampling.f) < 0 )
        stop("sampling.f must be in the range [0,1]")

    if ( !is.null(unresolved) )
    {
        stop("Unresolved clades not yet available for GSE2")
        ### below must be fixed
        # required <- c("tip.label", "Nc", "n0", "n1")
        # if ( !all(required %in% names(unresolved)) )
        #     stop("Required columns missing from unresolved clades")
        # if ( !all(unresolved$tip.label %in% tree$tip.label) )
        #     stop("Unknown tip species in 'unresolved'")
        # unresolved$k   <- unresolved$n1
        # unresolved$nsc <- unresolved$n0 + unresolved$n1
        # unresolved$i   <- match(unresolved$tip.label, tree$tip.label)
    }

    if ( is.null(names(states)) ) {
        stop("The states vector must contain names")
    } else 
    {
        known <- names(states)
        if ( !is.null(unresolved) )
            known <- unique(c(known, as.character(unresolved$tip.label)))
        if ( !all(tree$tip.label %in% known) )
            stop("Not all species have state information")
    }

    edge <- tree$edge
    idx <- seq_len(max(edge))
    ntip <- length(tree$tip.label)
    root <- ntip + 1
  
    is.tip <- idx <= ntip

    tip.state <- rep(NA, ntip)
    tip.state[] <- states[match(tree$tip.label, names(states))]

    children <- lapply(idx[!is.tip], function(x) edge[edge[,1] == x,2])
    if ( !all(sapply(children, length)==2) )
        stop("Multifircations/unbranched nodes in tree - must get rid of them")
    children <- rbind(matrix(NA, ntip, 2), t(matrix(unlist(children), 2)))

    ans <- list(len=tree$edge.len[match(idx, edge[,2])],
                is.tip=is.tip,
                tip.state=tip.state,
                tip.label=tree$tip.label,
                children=children,
                order=get.ordering(children, is.tip, root),
                root=root,
                unresolved=unresolved,
                sampling.f=sampling.f,
                nt.extra=nt.extra)
    class(ans) <- "bissecache"    # can just borrow this class
    return(ans)
}


### get.ordering() provided in bisse.R


### most of the likelihood calculation
gse2.branches <- function(pars, cache)
{
    # why defined within gse2.branches()?
    tips <- function(state)
    {
        # i = tips in this state that are resolved
        if ( is.na(state) )
            i <- setdiff(which(is.na(tip.state)), unresolved$i)
        else
            i <- setdiff(which(tip.state == state), unresolved$i)

        if ( length(i) > 0 )
        {
            t <- len[i]                 # branch lengths
            j <- tapply(t, t)           # get indices, repeated as necessary
            t.uniq <- sort(unique(t))   # allow reusing integrals when possible
            
            # yi = (E.0, E.1, E.2, D.0, D.1, D.2) for the tip
            #       E_i = 1 - f_i
            #       D_i = f_i if the tip is state i, 0 otherwise
            #       D_i = f_i for all tips if the state is unknown
            yi <- c(1 - sampling.f, 0, 0, 0)
            if ( is.na(state) )
                yi[4:6] <- sampling.f
            else
                yi[state + 4] <- sampling.f[state + 1]
            
            # ans has one row per time, with E's and D's as columns
            ans <- t(solve.gse2(yi, t.uniq, pars))
            n <- length(i)
            branch.init[i,] <<- matrix(rep(yi, n), n, 6, TRUE)
            branch.base[i,] <<- ans[j,-1]     # drop column of times
        }
    }

    ### needs any changes?
    # unresolved.tips <- function()
    # {
    #     i <- unresolved$i
    #     Nc <- unresolved$Nc
    #     k <- unresolved$k
    #     nsc <- unresolved$nsc
    #     t <- len[i]
    #     nt.extra <- cache$nt.extra
    #     branch.base[i,] <<- solve.unresolved(Nc, k, t, pars, nsc, nt.extra)
    # }

    len <- cache$len
    tip.state <- cache$tip.state
    children <- cache$children
    order <- cache$order
    root <- cache$root
    unresolved <- cache$unresolved
    nnode <- length(order)
    sampling.f <- cache$sampling.f

    # set up required space
    n <- length(len)
    branch.init <- branch.base <- matrix(NA, n, 6)
    lq <- rep(0, n)  # for underflow compensation

    ### compute transition probabilities for the tip branches
    if ( !is.null(unresolved) )
        unresolved.tips()
    tips(0)
    tips(1)
    tips(2)
    tips(NA)

    ### integrate down the rest of the branches
    for ( i in order[-nnode] )   # all internal nodes except the root
    {
        # combine at this node
        y.in <- initial.conditions.gse2(branch.base[children[i,],], pars)

        # extract underflow factor
        lq[i] <- y.in[7]
        y.in <- y.in[-7]

        if ( any(is.na(y.in) | y.in < 0) )
            stop("Invalid conditions")

        # integrate down the branch
        branch.init[i,] <- y.in
        branch.base[i,] <- solve.gse2(y.in, len[i], pars)[-1]
    }

    ### the final node join does not include the speciation event
    y.in <- initial.conditions.gse2(branch.base[children[root,],], pars, TRUE)
    lq[root] <- y.in[7]
    branch.init[root,] <- y.in[-7]
    return( list(init=branch.init, base=branch.base, lq=lq) )
}


### do the node computations, including underflow compensation
# named to avoid conflict with initial.conditions() in bisse.R; could instead be internal to gse2.branches()
initial.conditions.gse2 <- function(init, pars, is.root=FALSE)
{
    # E.0, E.1, E.2
    e <- init[1, c(1,2,3)]

    # D.1, D.2  (Eq. 6bc)
    d12 <- init[1, c(5,6)] * init[2, c(5,6)]
    if ( !is.root )
    {
        # pars[1:3] = sA, sB, sAB
        d12 <- d12 * pars[c(1,2)]

        # D.0 (Eq. 6a)
        d0 <- 0.5 * sum(init[1, c(4,5)] * init[2, c(5,4)] * pars[1] + 
                        init[1, c(4,6)] * init[2, c(6,4)] * pars[2] +
                        init[1, c(5,6)] * init[2, c(6,5)] * pars[3])
    } else
    {
        d0 <- 0.5 * sum(init[1, c(4,5)] * init[2, c(5,4)] + 
                        init[1, c(4,6)] * init[2, c(6,4)] +
                        init[1, c(5,6)] * init[2, c(6,5)])
    }
    d <- c(d0, d12)

    # avoid underflow by yanking out a factor 
    # (gets added in with the final calculation of the likelihood)
    q <- min(d)
    if ( q != 0 )
        return( c(e, d/q, log(q)) )
    else
        return( c(e, d, 0) )
}


### equilibrium frequencies
stationary.freq.gse2 <- function(pars)
{
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


### use bisse.prior() in bisse.R


### needs changes
#--------------------------------------------------
# solve.unresolved <- function(Nc, k, t, pars, nsc=Nc, nt.extra=10,
#                              nt=max(Nc)+nt.extra) {
#   lambda0 <- pars[1]
#   lambda1 <- pars[2]
#   mu0 <- pars[3]
#   mu1 <- pars[4]
#   q01 <- pars[5]
#   q10 <- pars[6]
#   bucexpl(nt, mu0, mu1, lambda0, lambda1, q01, q10, t,
#           Nc, nsc, k)[,c(3:4,1:2)]
# }
#-------------------------------------------------- 
