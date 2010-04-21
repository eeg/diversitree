source("/home/emma/src/miscR/ttn.R")
mycol <- c(rgb(0.832, 0.367, 0), rgb(0, 0.445, 0.695))
ttn <- read.ttn("example-clean.ttn", nodes=F)

plot(ttn$tree, show.tip.label=F, no.margin=T)
tiplabels(text=as.character(seq(ttn$tree$Nnode+1)), bg=mycol[ttn$states+1], adj=0)
dev.off()

#--------------------------------------------------
# straight bisse
#-------------------------------------------------- 

library(diversitree)

lnL.1 <- make.bisse(ttn$tree, ttn$states)
params <- c(1.4, 0.4, 0.2, 0.1, 0.6, 0.35)
names(params) <- argnames(lnL.1)
lnL.1(params)       # -121.3362

#--------------------------------------------------
# gp-bisse
#-------------------------------------------------- 

library(diversitreeGP)

### construct parameter list

nstates <- 2

Lam <- array(0, dim=rep(nstates, 3))
dimnames(Lam) <- list(paste("p", seq(nstates), sep="."), paste("d1", seq(nstates), sep="."), paste("d2", seq(nstates), sep="."))
Lam[1,1,1] <- 1.4
Lam[2,2,2] <- 0.4
Lam[1,,]

Mu <- c(0.2, 0.1)

Q <- array(0, dim=rep(nstates, 2))
Q[1,2] <- 0.6
Q[2,1] <- 0.35
dimnames(Q) <- list(paste("from", seq(nstates), sep="."), paste("to", seq(nstates), sep="."))

params <- list(lambda=Lam, mu=Mu, q=Q, nstates=as.integer(nstates))

### single likelihood calculation

lnL.2 <- make.gpbisse(ttn$tree, ttn$states)
lnL.2(params)       # -121.3362
# agrees with make.bisse!

### flatten and re-list-ify the parameter structure (needed for find.mle)

# recycle logic from argnames
flatten.pars.gpbisse <- function(parlist)
{
  k <- parlist$nstates

  i.lam <- cbind(rep(1:k, each=k*(k+1)/2), rep(rep(1:k, times=seq(k,1,-1)), k),
                 unlist(lapply(1:k, function(i) i:k)))

  i.q <- cbind(rep(1:k, each=k-1), unlist(lapply(1:k, function(i) (1:k)[-i])))

  parvec <- c(parlist$lambda[i.lam], parlist$mu, parlist$q[i.q])
  names(parvec) <- diversitreeGP:::argnames.gpbisse(NULL, k)
  return(parvec)
}
parvec1 <- flatten.pars.gpbisse(params)

# omitting elements from unlist(parlist) isn't any cleaner
# wrong for k > 2, and slower than flatten anyway
unlist.pars.gpbisse <- function(parlist)
{
    k <- parlist$nstates

    i <- unlist(sapply(seq(by=2, length.out=k), 
                       function(n) k * (n-1) + seq(1, n, by=2)))
    idx.lam <- rep(i, k) + rep(seq(k), each=length(i))

    idx.mu <- k^3 + seq(k)

    i <- unlist(lapply(seq(k), function(x) seq(x, by=k, length.out=k)))
    idx.q <- k^3 + k + i[-seq(1, k^2, by=k+1)]

    parvec <- unlist(parlist)[c(idx.lam, idx.mu, idx.q)]
    names(parvec) <- diversitreeGP:::argnames.gpbisse(NULL, k)
    return(parvec)
}
parvec2 <- unlist.pars.gpbisse(params)

# test with more states
nstates <- 3
Lam <- array(seq(from=0.01, by=0.01, length.out=nstates^3), dim=rep(nstates, 3))
Mu <- c(0.2, 0.3, 0.4)
Q <- array(seq(from=0.15, by=0.1, length.out=nstates^2), dim=rep(nstates, 2))
params <- list(lambda=Lam, mu=Mu, q=Q, nstates=as.integer(nstates))

p1 <- flatten.pars.gpbisse(params)
p2 <- unlist.pars.gpbisse(params)
p1 == p2

system.time(flatten.pars.gpbisse(params))
system.time(unlist.pars.gpbisse(params))


# undo logic in flatten.pars.gpbisse
listify.pars.gpbisse <- function(parvec, k)
{
    if ( length(parvec) != k*k*(k+1)/2 + k + k*(k-1))
        stop("Invalid length of parameter vector.")

    Lam <- array(0, dim=rep(k, 3))
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

params
p1 <- flatten.pars.gpbisse(params)
p2 <- listify.pars.gpbisse(p1, nstates)
lnL.2(p2)


