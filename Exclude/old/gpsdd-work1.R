# http://math.acadiau.ca/ACMMaC/howtos/C_R.html
#
# $ gcc -c -I/usr/share/R/include -fPIC gp-bisse-eqs.c
# $ gcc -shared -o gp-bisse-eqs.so gp-bisse-eqs.o
# creates gp-bisse-eqs.so

library(deSolve)
source("../R/util.R")   # now have make.ode()
dyn.load("../src/gp-bisse-eqs.so")
#dyn.unload("../src/gp-bisse-eqs.so")
is.loaded("initmod_gp")
assign("gp.bisse.ode", make.ode("derivs_gp", "gp-bisse-eqs", "initmod_gp", 4, FALSE))

RTOL <- ATOL <- 1e-8
t0 <- 1
len <- c(0.1, 0.5, 1)
y <- c(0.1, 0.2, 0.3, 0.4)

### work below is reverse-chronological ###

# attempt more states
nstates <- 3
assign("gp.bisse.ode", make.ode("derivs_gp", "gp-bisse-eqs", "initmod_gp", nstates*2, FALSE))
y <- seq(from=0.11, by=0.1, length.out=nstates*2)

Lam <- array(seq(from=0.01, by=0.01, length.out=nstates^3), dim=rep(nstates, 3))
Mu <- c(0.2, 0.3, 0.4)
Q <- array(seq(from=0.15, by=0.1, length.out=nstates^2), dim=rep(nstates, 2))
dimnames(Q) <- list(paste("from", seq(nstates), sep="."), paste("to", seq(nstates), sep="."))
pars <- list(lambda=Lam, mu=Mu, q=Q, nstates=as.integer(nstates))
gp.bisse.ode(y, c(t0, t0+len), pars, rtol=RTOL, atol=ATOL)[-1,-1]
# runs fine, but I don't have an answer to compare with


# testing error checking
pars <- list(lambda=Lam, mu=c(0.2, 0.3), q=Q, nstates=as.integer(2))
pars <- list(lambda=Lam, mu=c(0.2, 0.3), q=Q)
pars <- list(blah=Lam, mu=c(0.2, 0.3), q=Q, nstates=as.integer(2))
pars <- list(lambda=Lam, blah=c(0.2, 0.3), q=Q, nstates=as.integer(2))
pars <- list(lambda=Lam, mu=c(0.2, 0.3), blah=Q, nstates=as.integer(2))
pars <- list(lambda=Lam, mu=c(0.2, 0.3), q=Q, blah=as.integer(2))
gp.bisse.ode(y, c(t0, t0+len), pars, rtol=RTOL, atol=ATOL)[-1,-1]

# with non-zero lambda_ijk (i != j, k)
Lam <- array(seq(from=0.01, by=0.01, length.out=8), dim=rep(nstates, 3))
dimnames(Lam) <- list(paste("p", seq(nstates), sep="."), paste("d1", seq(nstates), sep="."), paste("d2", seq(nstates), sep="."))
Lam[1,1,1] <- 1.1
Lam[2,2,2] <- 0.9
Lam[2,,]
pars <- list(lambda=Lam, mu=c(0.2, 0.3), q=Q, nstates=as.integer(2))
gp.bisse.ode(y, c(t0, t0+len), pars, rtol=RTOL, atol=ATOL)[-1,-1]
# should require/check that elements with d1 > d2 are zero
# now get a slightly different answer, as expected

# parse lambda matrix in C code
# with Lam as below, get identical answer!

# include nstates, parse q matrix in C code
pars <- list(lambda=Lam, mu=c(0.2, 0.3), q=Q, nstates=as.integer(2))
gp.bisse.ode(y, c(t0, t0+len), pars, rtol=RTOL, atol=ATOL)[-1,-1]
# works!

# also make lambda an array
#ls <- c(1.1, 0.01, 0.02, 0.03, 0.04, 0.9) # 111, 112, 122, 211, 212, 222
Lam <- array(0, dim=rep(nstates, 3))
dimnames(Lam) <- list(paste("p", seq(nstates), sep="."), paste("d1", seq(nstates), sep="."), paste("d2", seq(nstates), sep="."))
Lam[1,1,1] <- 1.1
Lam[2,2,2] <- 0.9
Lam[2,,]
pars <- list(lambda=Lam, mu=c(0.2, 0.3), q=Q)
gp.bisse.ode(y, c(t0, t0+len), pars, rtol=RTOL, atol=ATOL)[-1,-1]
# works!

# make Q an array
nstates <- 2
qs <- c(0.4, 0.5)   # (q12, q21)
# using weird diags just to test ordering; should eventually use -rowSum()
Q <- array(c(0.01, qs[2], qs[1], 0.02), dim=rep(nstates, 2))
dimnames(Q) <- list(paste("from", seq(nstates), sep="."), paste("to", seq(nstates), sep="."))
pars <- list(lambda=c(1.1, 0.9), mu=c(0.2, 0.3), q=Q)
gp.bisse.ode(y, c(t0, t0+len), pars, rtol=RTOL, atol=ATOL)[-1,-1]
# works!

# test with pars as a list
pars <- list(lambda=c(1.1, 0.9), mu=c(0.2, 0.3), q=c(0.4, 0.5))
gp.bisse.ode(y, c(t0, t0+len), pars, rtol=RTOL, atol=ATOL)[-1,-1]
# works!

pars <- c(1.1, 0.9, 0.2, 0.3, 0.4, 0.5)
gp.bisse.ode(y, c(t0, t0+len), pars, rtol=RTOL, atol=ATOL)[-1,-1]
          [,1]      [,2]      [,3]
[1,] 0.1113393 0.1448683 0.1700954
[2,] 0.2045827 0.2217752 0.2391748
[3,] 0.2732338 0.1915487 0.1260480
[4,] 0.3633701 0.2498325 0.1596702


pars <- c(1, 1, 0.1, 0.1, 0.5, 0.5)
y <- c(0, 0, 1, 0)
     V1          V2         V3         V4
time  1 1.100000000 1.50000000 2.00000000
1     0 0.009472612 0.03870513 0.06185799
2     0 0.009472612 0.03870513 0.06185799
3     1 0.854032509 0.47330355 0.24473155
4     0 0.042666062 0.11592086 0.11309465
