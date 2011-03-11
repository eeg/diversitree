### Testing punctsse-eqs.c ###
# (see http://math.acadiau.ca/ACMMaC/howtos/C_R.html)

#--------------------------------------------------
# Preparation
#--------------------------------------------------

library(deSolve)
source("../R/util.R")   # now have make.ode()

RTOL <- ATOL <- 1e-8
# paste in argnames.XXX() functions from model-XXX.R

# turn the Q pars into a flattened matrix; for musse, punctsse
get.ode.pars.punctsse <- function(pars, k)
{
  qmat <- matrix(0, k, k)
  idx.qmat <- cbind(rep(1:k, each=k-1),
               unlist(lapply(1:k, function(i) (1:k)[-i])))
  x <- k * k * (k + 1) / 2 + k
  idx.lm <- seq(x)
  idx.q <- seq(x+1, (k+3)*k*k/2)
  qmat[idx.qmat] <- pars[idx.q]
  diag(qmat) <- -rowSums(qmat)
  pars <- c(pars[idx.lm], qmat)
}
get.ode.pars.musse <- function(pars, k)
{
  qmat <- matrix(0, k, k)
  idx.qmat <- cbind(rep(1:k, each=k-1),
               unlist(lapply(1:k, function(i) (1:k)[-i])))
  idx.lm <- 1:(2*k)
  idx.q <- (2*k+1):(k*(1+k))
  qmat[idx.qmat] <- pars[idx.q]
  diag(qmat) <- -rowSums(qmat)
  pars <- c(pars[idx.lm], qmat)
}

# integration time points
t0 <- 1
len <- c(0, 0.01, 0.1, 0.5, 1)

nstates <- 3

### bisse
# gcc -c -I/usr/share/R/include -fPIC bisse-eqs.c
# gcc -shared -o bisse-eqs.so bisse-eqs.o
dyn.unload("../src/bisse-eqs.so")
dyn.load("../src/bisse-eqs.so")
assign("bisse.ode", make.ode("derivs_bisse", "bisse-eqs", "initmod_bisse", 4))

### musse
# gcc -c -I/usr/share/R/include -fPIC musse-eqs.c
# gcc -shared -o musse-eqs.so util.o musse-eqs.o
dyn.unload("../src/musse-eqs.so")
dyn.load("../src/musse-eqs.so")
assign("musse.ode", make.ode("derivs_musse", "musse-eqs", "initmod_musse", nstates*2))

### geosse
# gcc -c -I/usr/share/R/include -fPIC geosse-eqs.c
# gcc -shared -o geosse-eqs.so geosse-eqs.o
dyn.unload("../src/geosse-eqs.so")
dyn.load("../src/geosse-eqs.so")
assign("geosse.ode", make.ode("derivs_geosse", "geosse-eqs", "initmod_geosse", 6))

### punctsse
# gcc -c -I/usr/share/R/include -fPIC punctsse-eqs.c
# gcc -shared -o punctsse-eqs.so util.o punctsse-eqs.o
dyn.unload("../src/punctsse-eqs.so")  # after re-compliation
dyn.load("../src/punctsse-eqs.so")
assign("punctsse.ode", make.ode("derivs_punctsse", "punctsse-eqs", "initmod_punctsse", nstates*2))

is.loaded("initmod_punctsse")

#--------------------------------------------------
# A: bisse-like
#--------------------------------------------------

# set nstates <- 2 above

parsAp <- c(1.11, 0, 0, 0, 0, 2.22, 0.8, 0.9, 3, 3.5)
names(parsAp) <- argnames.punctsse(parsAp, k=nstates)
parsAp2 <- get.ode.pars.punctsse(parsAp, nstates)

parsAb <- parsAp[c(1, 6, 7:10)]
#names(parsAb) <- argnames.bisse(parsAb)

# E_1, E_2, D_1, D_2
y <- c((seq(nstates)+5)/40, seq(nstates)/10 + 0.11)

ans.p <- punctsse.ode(y, c(t0, t0+len), parsAp2, rtol=RTOL, atol=ATOL)[-1,-1]
ans.b <- bisse.ode(y, c(t0, t0+len), parsAb, rtol=RTOL, atol=ATOL)[-1,-1]
ans.p - ans.b < ATOL

#       [,1]      [,2]      [,3]      [,4]       [,5]
# [1,] 0.150 0.1560461 0.2036809 0.3341365 0.41688459
# [2,] 0.175 0.1783534 0.2085161 0.3150286 0.39086064
# [3,] 0.210 0.2095535 0.1971129 0.1158479 0.06103487
# [4,] 0.310 0.2995603 0.2285053 0.1074512 0.05707963

#--------------------------------------------------
# B: musse-like
#--------------------------------------------------

# set nstates <- 4 above

parsBm <- c(1.11, 2.22, 3.33, 4.44, 0.8, 0.9, 1.0, 1.1, 3, 3.15, 3.20, 3.26, 3.31, 3.36, 3.42, 3.47, 3.52, 3.58, 3.63, 3.68)
#names(parsBm) <- argnames.musse(parsBm, nstates)
parsBm2 <- get.ode.pars.musse(parsBm, nstates)

parsBp <- c(rep(0, 40), parsBm[5:20])
names(parsBp) <- argnames.punctsse(parsBp, k=nstates)
parsBp[c('lambda111', 'lambda222', 'lambda333', 'lambda444')] <- parsBm[1:4]
parsBp2 <- get.ode.pars.punctsse(parsBp, nstates)

y <- c((seq(nstates)+5)/40, seq(nstates)/10 + 0.11)

ans.p <- punctsse.ode(y, c(t0, t0+len), parsBp2, rtol=RTOL, atol=ATOL)[-1,-1]
ans.m <- musse.ode(y, c(t0, t0+len), parsBm2, rtol=RTOL, atol=ATOL)[-1,-1]
ans.p - ans.m < ATOL

#       [,1]      [,2]      [,3]      [,4]       [,5]
# [1,] 0.150 0.1596113 0.2165531 0.3058746 0.34520916
# [2,] 0.175 0.1806658 0.2174674 0.2944914 0.33195269
# [3,] 0.200 0.2010961 0.2173159 0.2848409 0.32069011
# [4,] 0.225 0.2209668 0.2164004 0.2765359 0.31097486
# [5,] 0.210 0.2239086 0.2564211 0.1187024 0.04408851
# [6,] 0.310 0.3088922 0.2701337 0.1127951 0.04219938
# [7,] 0.410 0.3917367 0.2791581 0.1076727 0.04051487
# [8,] 0.510 0.4728096 0.2846919 0.1032074 0.03900954

#--------------------------------------------------
# C: geosse-like
#--------------------------------------------------

# set nstates <- 3 above

parsCp <- rep(0, length=length(argnames.punctsse(1, k=nstates)))
names(parsCp) <- argnames.punctsse(parsCp, k=nstates)

parsCg <- rep(NA, 7)
names(parsCg) <- argnames.geosse(parsCg)

parsCg['sA'] <- parsCp['lambda222'] <- parsCp['lambda112'] <- 2.22
parsCg['sB'] <- parsCp['lambda333'] <- parsCp['lambda113'] <- 3.33
parsCg['sAB'] <- parsCp['lambda123'] <- 1.23
parsCg['xA'] <- parsCp['mu2'] <- parsCp['q13'] <- 0.8
parsCg['xB'] <- parsCp['mu3'] <- parsCp['q12'] <- 0.9
parsCg['dA'] <- parsCp['q21'] <- 3.0
parsCg['dB'] <- parsCp['q31'] <- 3.5

parsCp2 <- get.ode.pars.punctsse(parsCp, nstates)

y <- c((seq(nstates)+5)/40, seq(nstates)/10 + 0.11)

ans.p <- punctsse.ode(y, c(t0, t0+len), parsCp2, rtol=RTOL, atol=ATOL)[-1,-1]
ans.g <- geosse.ode(y, c(t0, t0+len), parsCg, rtol=RTOL, atol=ATOL)[-1,-1]
ans.p - ans.g < ATOL
