### Testing the node join, initial.conditions.punctsse() ###

#--------------------------------------------------
# Set up a parameter vector
#--------------------------------------------------

argnames.punctsse(1, k=nstates)

# A: bisse-like
nstates <- 2
parsA <- c(1.11, 0, 0, 0, 0, 2.22, 0.8, 0.9, 3, 3.5)
names(parsA) <- argnames.punctsse(parsA, k=nstates)

# B: punct, 2 states
nstates <- 2
parsB <- c(1.11, 1.12, 1.22, 2.11, 2.12, 2.22, 0.8, 0.9, 3, 3.5)
names(parsB) <- argnames.punctsse(parsB, k=nstates)

# C: geosse-like
nstates <- 3
parsC <- rep(0, length=length(argnames.punctsse(1, k=nstates)))
names(parsC) <- argnames.punctsse(parsC, k=nstates)

parsCg <- rep(NA, 7)
names(parsCg) <- argnames.geosse(parsCg)

parsCg['sA'] <- parsC['lambda222'] <- parsC['lambda112'] <- 2.22
parsCg['sB'] <- parsC['lambda333'] <- parsC['lambda113'] <- 3.33
parsCg['sAB'] <- parsC['lambda123'] <- 1.23
parsCg['xA'] <- parsC['mu2'] <- parsC['q13'] <- 0.8
parsCg['xB'] <- parsC['mu3'] <- parsC['q12'] <- 0.9
parsCg['dA'] <- parsC['q21'] <- 3.0
parsCg['dB'] <- parsC['q31'] <- 3.5


#--------------------------------------------------
# Test the node join
#--------------------------------------------------

# init is a list, one element for clade M and one for clade N
#          E_i for M            D_i for M
initM <- c((seq(nstates)+5)/40, seq(nstates)/10 + 0.11)
#          E_i for N            D_i for N
initN <- c(initM[1:nstates],    seq(nstates)/15 + 0.22)
init <- list(initM, initN)

#--------------------
### A: bisse-like ###

ans.musse <- initial.conditions.musse(init, parsA[c(1,6)], t=NA, is.root=FALSE)
# 0.150000  0.175000  0.066822  0.243164
ans.musse == c(init[[1]][1], init[[1]][2], init[[1]][nstates+1] * init[[2]][nstates+1] * parsA[1], init[[1]][nstates+2] * init[[2]][nstates+2] * parsA[6])

ans.punctsse <- initial.conditions.punctsse(init, parsA, t=NA, is.root=FALSE)
# agrees with musse

#-------------------------
### B: punct, 2 states ###

ans.punctsse <- initial.conditions.punctsse(init, parsB, t=NA, is.root=FALSE)
# 0.1500000 0.1750000 0.2917700 0.5430367

DM <- init[[1]][(nstates+1):(2*nstates)]
DN <- init[[2]][(nstates+1):(2*nstates)]
ans.punctsse == c( init[[1]][1], init[[1]][2],
    sum( parsB[1] * (DM[1] * DN[1]),
         parsB[2] * (DM[1] * DN[2] + DM[2] * DN[1]) / 2,
         parsB[3] * (DM[2] * DN[2]) ),
    sum( parsB[4] * (DM[1] * DN[1]),
         parsB[5] * (DM[1] * DN[2] + DM[2] * DN[1]) / 2,
         parsB[6] * (DM[2] * DN[2]) ) )
# agrees

# compare with old, modified to init as list and lambda as matrix
initial.conditions.gpsdd <- function(init, pars, t, is.root=FALSE)
{
  Lam <- array(0, dim=rep(nstates, 3))
  Lam[1,1,1] <- pars['lambda111']
  Lam[1,2,1] <- pars['lambda112']
  Lam[1,2,2] <- pars['lambda122']
  Lam[2,1,1] <- pars['lambda211']
  Lam[2,2,1] <- pars['lambda212']
  Lam[2,2,2] <- pars['lambda222']
  n <- nstates
  # E.1, E.2
  e <- init[[1]][seq(n)]
  # D.1, D.2
  d <- rep(0, n)
  for (i in seq(n)) {
    lambda <- Lam[i,,]
    for (j in seq(n))
      for (k in seq(j))
        d[i] <- d[i] + 0.5 * lambda[j,k] * 
                       sum(init[[1]][c(n+j,n+k)] * init[[2]][c(n+k,n+j)])
  }
  c(e, d)
}
ans.gpsdd <- initial.conditions.gpsdd(init, parsB, t=NA, is.root=FALSE)
# agrees


#---------------------
### C: geosse-like ###

ans.punctsse <- initial.conditions.punctsse(init, parsC, t=NA, is.root=FALSE)
# 0.150000 0.175000 0.200000 0.692716 0.243164 0.573426

ans.geosse <- initial.conditions.geosse(init, parsCg, t=NA, is.root=FALSE)

ans.geosse == ans.punctsse # agrees
