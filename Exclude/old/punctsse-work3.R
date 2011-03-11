### Testing PunctSSE on a tree ###

library(diversitreeGP)

source("/home/emma/src/miscR/ttn.R")
ttn <- read.ttn("example-clean.ttn", nodes=F)

#--------------------------------------------------
# A: bisse-like
#--------------------------------------------------

nstates <- 2
nlambda <- nstates*nstates*(nstates+1)/2

lnL.Ab <- make.bisse(ttn$tree, ttn$states)
parsAb <- c(1.0, 0.5, 0.2, 0.1, 0.6, 0.4)
names(parsAb) <- argnames(lnL.Ab)

lnL.Ap <- make.punctsse(ttn$tree, ttn$states+1, nstates)
parsAp <- c(rep(0, nlambda), parsAb[-seq(nstates)])
names(parsAp) <- argnames(lnL.Ap)
parsAp[c('lambda111', 'lambda222')] <- parsAb[seq(nstates)]

lnL.Ab(parsAb)  # -122.2663
lnL.Ap(parsAp)  # -122.2663 :)

#--------------------------------------------------
# B: musse-like
#--------------------------------------------------

nstates <- 4
states4 <- c(ttn$states[seq(10)]+2, ttn$states[-seq(10)]) + 1

lnL.Bm <- make.musse(ttn$tree, states4, nstates)
parsBm <- c(1.11, 2.22, 3.33, 4.44, 0.8, 0.9, 1.0, 1.1, 3, 3.15, 3.20, 3.26, 3.31, 3.36, 3.42, 3.47, 3.52, 3.58, 3.63, 3.68)
names(parsBm) <- argnames(lnL.Bm)

lnL.Bp <- make.punctsse(ttn$tree, states4, nstates)
nlambda <- nstates*nstates*(nstates+1)/2
parsBp <- c(rep(0, nlambda), parsBm[-seq(nstates)])
names(parsBp) <- argnames(lnL.Bp)
parsBp[c('lambda111', 'lambda222', 'lambda333', 'lambda444')] <- parsBm[1:4]

lnL.Bm(parsBm)  # -216.7327
lnL.Bp(parsBp)  # -216.7327

#--------------------------------------------------
# C: geosse-like
#--------------------------------------------------

nstates <- 3
states3 <- c(ttn$states[seq(10)]+1, ttn$states[-seq(10)]) + 1

lnL.Cg <- make.geosse(ttn$tree, states3-1)
lnL.Cp <- make.punctsse(ttn$tree, states3, nstates)

parsCg <- rep(NA, 7)
names(parsCg) <- argnames(lnL.Cg)
parsCp <- 0 * starting.point.punctsse(ttn$tree, nstates)

parsCg['sA'] <- parsCp['lambda222'] <- parsCp['lambda112'] <- 2.22
parsCg['sB'] <- parsCp['lambda333'] <- parsCp['lambda113'] <- 3.33
parsCg['sAB'] <- parsCp['lambda123'] <- 1.23
parsCg['xA'] <- parsCp['mu2'] <- parsCp['q13'] <- 0.8
parsCg['xB'] <- parsCp['mu3'] <- parsCp['q12'] <- 0.9
parsCg['dA'] <- parsCp['q21'] <- 3.0
parsCg['dB'] <- parsCp['q31'] <- 3.5

lnL.Cg(parsCg, condition.surv=T)  # -309.1453
lnL.Cp(parsCp, condition.surv=T)  # -309.1453

lnL.Cg(parsCg, condition.surv=F)  # -307.7621
lnL.Cp(parsCp, condition.surv=F)  # -307.7621



#--------------------------------------------------
# debugging work
#--------------------------------------------------

source("aaaGP.R")

# make.punctsse()
tree <- ttn$tree
states <- ttn$states+1
k <- nstates
sampling.f <- NULL
strict <- TRUE
safe <- FALSE
cache <- make.cache.musse(tree, states, k, sampling.f, strict)
branches <- make.branches.punctsse(k, safe)

ll.xxsse(parsAp, cache, initial.conditions.punctsse, branches, TRUE, ROOT.OBS, NULL, FALSE) # Inf

# ll.xxsse()
initial.conditions <- initial.conditions.punctsse
pars <- parsAp
condition.surv <- TRUE
root <- ROOT.OBS
root.p <- NULL
intermediates <- FALSE
ans <- all.branches(pars, cache, initial.conditions, branches)
vals <- ans$init[[cache$root]]
root.p <- root.p.xxsse(vals, pars, root, root.p)

loglik <- root.xxsse(vals, pars, ans$lq, condition.surv, root.p) # Inf
