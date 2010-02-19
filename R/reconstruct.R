# EEG: This is my old reconstruction code, using node.fixing.

### Reconstruct node states using the "global, marginal" method.
### This is not necessarily the most appropriate method for many questions 
###   (Pagel 1999), but it is the quickest and most commonly used.  
###     global: use a single set of rates (the ML estimates) rather than re-
###           estimating after each node fixing
###     marginal: consider each node separately, rather than using a joint set 
###           of node states
### Currently only works for a binary character.

### Returns a data.frame with node number and probabilities of being in 
###   states 0 or 1.
reconstruct.bisse <- function(tree, states, pars, ...)
{
    if ( any(pars < 0) || any(!is.finite(pars)) )
        stop("Invalid parameters.")

    ntips <- Ntip(tree)
    nnodes <- tree$Nnode
    ntotal <- ntips + nnodes

    # will hold the desired probabilities of each state
    ans <- data.frame(num = seq(nnodes) + ntips, p0 = rep(NA, nnodes), 
                      p1 = rep(NA, nnodes))

    # for each internal node
    for (node in seq(ntips + 1, ntotal))
    {
        node.fixing <- rep(NA, nnodes)
        likes <- rep(NA, 2)            # will hold (not-log-) likelihoods
        for (s in c(0, 1))             # fix node to state 0, then 1
        {
            node.fixing[node-ntips] <- s
            f <- make.bisse(tree, states, node.fixing = node.fixing)
            loglike <- f(pars = pars, ...)
            likes[s+1] <- exp(loglike)
        }
        ans[node-ntips, -1] <- likes / sum(likes)
    }
    
    return(ans)
}
