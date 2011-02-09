### starts off like gpsdd-work2.R, then moves on to punct sim tree
### Jan 17, 2011

library(diversitreeGP)

source("/home/emma/src/miscR/ttn.R")
mycol <- c(rgb(0.832, 0.367, 0), rgb(0, 0.445, 0.695))

#--------------------------------------------------
# re-doing from gpsdd-work2.R
#-------------------------------------------------- 

ttn <- read.ttn("example-clean.ttn", nodes=F)

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

lnL.2 <- make.gpsdd(ttn$tree, ttn$states+1, 2)
lnL.2(params)       # -121.3362
# agrees with value from gpsdd-work2.R

### ML estimation

parvec <- flatten.pars.gpsdd(params)
find.mle(lnL.2, parvec, method="subplex")
#       lam111       lam112       lam122       lam211       lam212       lam222 
# 1.109783e+00 0.000000e+00 0.000000e+00 0.000000e+00 0.000000e+00 3.053672e-01 
#          mu1          mu2          q12          q21 
# 2.591077e-07 2.234287e-09 5.951669e-01 3.876896e-01 
# $lnLik
# [1] -120.3201
# nearly identical to gpsdd-work2.R

parvec2 <- parvec
parvec2[1:6] <- 1.0
find.mle(lnL.2, parvec2, method="subplex")
#       lam111       lam112       lam122       lam211       lam212       lam222 
# 1.109797e+00 4.478385e-07 5.553669e-06 9.589924e-07 1.231699e-07 3.053644e-01 
#          mu1          mu2          q12          q21 
# 6.027390e-07 5.057120e-07 5.952078e-01 3.878014e-01 
# $lnLik
# [1] -120.3201
# same conclusion!


#--------------------------------------------------
# new: simulated tree with punct-only character change
#--------------------------------------------------

ttn <- read.ttn("/data/UICwork/nsf_proposal/prelim_results/punct/set4/tree-100.ttn", nodes=T)
ttn$tree$node.label <- NULL

# true values (lam = 6.2, p = 0.2, q = 0; symmetric)
# lam111 lam112 lam122 lam211 lam212 lam222    mu1    mu2    q12    q21 
#   4.96   1.24      0      0   1.24   4.96      3      0      0      0

lnL.3 <- make.gpsdd(ttn$tree, ttn$states+1, 2)
lnL.3a <- constrain(lnL.3, lam122 ~ 0, lam211 ~ 0)

parvec3 <- rep(1, 8)
names(parvec3) <- c("lam111", "lam112", "lam212", "lam222", "mu1", "mu2", "q12", "q21")
find.mle(lnL.3a, parvec3, method="subplex")
#       lam111       lam112       lam212       lam222          mu1          mu2 
# 6.318710e+00 2.425054e-06 1.952252e-06 7.394043e+00 7.315748e-01 3.601709e+00 
#          q12          q21 
# 1.926215e+00 1.459083e+00 
# $lnLik
# [1] 54.04424
# way wrong :(

parvec3 <- c(5, 1, 1, 5, 3, 0.1, 0.1, 0.1)
names(parvec3) <- c("lam111", "lam112", "lam212", "lam222", "mu1", "mu2", "q12", "q21")
find.mle(lnL.3a, parvec3, method="subplex")
#       lam111       lam112       lam212       lam222          mu1          mu2 
# 6.318678e+00 2.014823e-07 3.537784e-07 7.394119e+00 7.314644e-01 3.601857e+00 
#          q12          q21 
# 1.926247e+00 1.459049e+00 
# $lnLik
# [1] 54.04431

ttn <- read.ttn("/data/UICwork/nsf_proposal/prelim_results/punct/set4/tree-21.ttn", nodes=T)
ttn$tree$node.label <- NULL
#       lam111       lam112       lam212       lam222          mu1          mu2 
# 8.197010e+00 4.821762e-08 7.645407e-11 7.719576e+00 2.317347e+00 2.946245e+00 
#          q12          q21 
# 3.187476e+00 2.283578e+00 
# $lnLik
# [1] 84.277

ttn <- read.ttn("/data/UICwork/nsf_proposal/prelim_results/punct/set2/tree-100.ttn", nodes=T)
ttn$tree$node.label <- NULL
# true values (lam = 1,5, p = 0,0.2, q = 0)
# lam111 lam112 lam122 lam211 lam212 lam222    mu1    mu2    q12    q21 
#      1      0      0      0      1      4      0      0      0      0
#
#       lam111       lam112       lam212       lam222          mu1          mu2 
# 5.480797e+00 3.702360e-11 8.183324e-10 3.836273e+00 2.889527e+00 1.140193e-05 
#          q12          q21 
# 6.823255e-06 1.035233e+00 
# $lnLik
# [1] -17.96415

# not looking good...
