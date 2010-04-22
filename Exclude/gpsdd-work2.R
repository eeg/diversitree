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


#--------------------------------------------------
# single likelihood calculation
#-------------------------------------------------- 

lnL.2 <- make.gpbisse(ttn$tree, ttn$states)
lnL.2(params)       # -121.3362
# agrees with make.bisse!


#--------------------------------------------------
# flatten and re-list-ify the parameter structure (needed for find.mle)
#-------------------------------------------------- 

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
lnL.2(p1)
lnL.2(p2)


#--------------------------------------------------
# ML estimates
#-------------------------------------------------- 

### bisse
lnL.1 <- make.bisse(ttn$tree, ttn$states)
par.1 <- c(1.4, 0.4, 0.2, 0.1, 0.6, 0.35)
names(par.1) <- argnames(lnL.1)
find.mle(lnL.1, par.1, method="subplex")
# $par
#      lambda0      lambda1          mu0          mu1          q01          q10 
# 1.109809e+00 3.053652e-01 3.430396e-06 2.187329e-06 5.952360e-01 3.878092e-01 
# $lnLik
# [1] -120.3201


### gpbisse

lnL.2 <- make.gpbisse(ttn$tree, ttn$states)
par.2 <- params     # list defined above

par.3 <- flatten.pars.gpbisse(par.2)
find.mle(lnL.2, par.3, method="subplex")
#       lam111       lam112       lam122       lam211       lam212       lam222 
# 1.109786e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 3.053679e-01 
#          mu1          mu2          q12          q21 
# 2.306045e-07 2.381490e-07 5.951741e-01 3.877046e-01 
# $lnLik
# [1] -120.3201

par.4 <- par.3
par.4[2:5] <- 0.1
find.mle(lnL.2, par.4, method="subplex")
#       lam111       lam112       lam122       lam211       lam212       lam222 
# 1.109740e+00 3.092353e-08 6.019280e-07 1.560118e-06 1.046539e-10 3.053645e-01 
#          mu1          mu2          q12          q21 
# 1.712336e-06 5.788996e-07 5.951114e-01 3.876538e-01 
# $lnLik
# [1] -120.3201

# estimation of punctuational lambdas is suspiciously good!

attr(lnL.2, "k") <- nstates
lnL.3 <- constrain(lnL.2, lam112 ~ 0, lam122 ~ 0, lam211 ~ 0, lam212 ~ 0)
find.mle(lnL.3, par.3[-seq(2,5)], method="subplex")
# $par
#       lam111       lam222          mu1          mu2          q12          q21 
# 1.109809e+00 3.053652e-01 3.430396e-06 2.187329e-06 5.952360e-01 3.878092e-01 
# $lnLik
# [1] -120.3201

# constraining works!

lnL.4 <- constrain(lnL.2, lam122 ~ lam112, lam211 ~ lam112, lam212 ~ lam112)
find.mle(lnL.4, par.3[-seq(3,5)], method="subplex")
#       lam111       lam112       lam222          mu1          mu2          q12 
# 1.109843e+00 2.902625e-07 3.053895e-01 1.480125e-05 9.099486e-05 5.952365e-01 
#          q21 
# 3.877907e-01 
# $lnLik
# [1] -120.3202
find.mle(lnL.4, par.3, method="subplex")
# same answer, but "guessing" warning

