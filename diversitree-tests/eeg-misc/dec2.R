library(diversitree)

fixInNamespace("initial.conditions.mkn", "diversitree")
# function (init, pars, t, idx)
# {
#     browser()
#     init[, 1] * init[, 2]
# }

phy <- read.tree(text="((sp1:1, sp2:1):2, (sp3:1, sp4:1):2);")
states <- c(1, 1, 2, 2)
names(states) <- phy$tip.label
pars <- c(0.1, 0.1)
lik <- make.mkn(phy, states, k=2, control=list(method="ode"))
lik(pars)


phy <- read.tree(text="((sp1:1, sp2:1):2, (sp3:1, sp4:1):2);")
states <- c(1, 1, 2, 3)
names(states) <- phy$tip.label
pars <- rep(0.1, 6)
lik <- make.mkn(phy, states, k=3, control=list(method="ode"))
lik(pars)


### First pass at model.dec is done; ready to try...

library(diversitree)
# from an intermediate ver. 0.9-4
# downgraded deSolve from 1.10-4 to 1.10-3

#--------------------------------------------------
# Test 0: no error messages
#--------------------------------------------------

phy <- read.tree(text="((sp1:1, sp2:1):2, (sp3:1, sp4:1):2);")
states <- c(1, 1, 2, 0)
names(states) <- phy$tip.label
lik <- make.dec2(phy, states)
pars <- rep(0.1, 4)
names(pars) <- argnames(lik)
lik(pars)

#--------------------------------------------------
# Test 1: a simulated tree
#--------------------------------------------------
# A simulated tree with:
# sA = sB = 1.0
# sAB = 0.3
# xA = xB = 0.1
# dA = 0.8, dB = 0.2

source("/home/emma/src/SimTreeSDD/packaged/SimTreeSDD-20101213/misc/ttn.R")
ttn <- read.ttn("dec0.ttn")
phy <- ttn$tree
states <- ttn$tip.states

lik.g0 <- make.geosse(phy, states)
lik.g <- constrain(lik.g0, sB ~ sA)
pars.g <- c(1, 0.3, 0.1, 0.1, 0.5, 0.5)
names(pars.g) <- argnames(lik.g)
ans.g <- find.mle(lik.g, pars.g)
coef(ans.g)
#        sA       sAB        xA        xB        dA        dB 
# 1.0130473 0.3301317 0.6984612 0.2858643 0.4296312 0.8988868 

lik.g <- constrain(lik.g0, sB ~ sA, xB ~ xA)
pars.g <- c(1, 0.3, 0.1, 0.5, 0.5)
#        sA       sAB        xA        dA        dB 
# 0.9972561 0.3569604 0.4460631 0.6306559 0.6864985 

lik.d <- make.dec2(phy, states)
pars.d <- c(0.1, 0.1, 0.5, 0.5)
names(pars.d) <- argnames(lik.d)
ans.d <- find.mle(lik.d, pars.d)
coef(ans.d)
#         xA         xB         dA         dB 
# 23.8389878  0.3094923  0.5705806 12.1966462 

lik.d <- constrain(lik.d, xA ~ xB)
pars.d <- c(0.1, 0.5, 0.5)
#       xB       dA       dB 
# 2.013585 1.130560 1.163226 

### Compare with lagrange ML rates
# max lnL is close, but rates are way off

#   -lnL = 45.42
#   dispersal = 0.4966
#   extinction = 5.054e-07
lik.d <- constrain(lik.d, xA ~ xB, dA ~ dB)
pars.d <- c(0.1, 0.5)
#       xB       dB 
# 2.026340 1.156641 
logLik(ans.d)
# 'log Lik.' -41.73536 (df=2)

# I think lagrange uses a flat root
ans.d <- find.mle(lik.d, pars.d, root=ROOT.FLAT)
lik.d(c(5.054e-07, 0.4966))

#--------------------------------------------------
# Test 2: a larger simulated tree
#--------------------------------------------------
# same params

source("/home/emma/src/SimTreeSDD/packaged/SimTreeSDD-20101213/misc/ttn.R")
ttn <- read.ttn("dec0.ttn")
phy <- ttn$tree
states <- ttn$tip.states

lik.g0 <- make.geosse(phy, states)
lik.g <- constrain(lik.g0, sB ~ sA)
pars.g <- c(1, 0.3, 0.1, 0.1, 0.5, 0.5)
names(pars.g) <- argnames(lik.g)
ans.g <- find.mle(lik.g, pars.g)
coef(ans.g)
#           sA          sAB           xA           xB           dA           dB 
# 1.081867e+00 1.011569e-01 4.936034e-10 2.500478e-02 9.929429e-01 1.900192e-01 

lik.g <- constrain(lik.g0, sB ~ sA, xB ~ xA)
pars.g <- c(1, 0.3, 0.1, 0.5, 0.5)
#           sA          sAB           xA           dA           dB 
# 1.074412e+00 1.096301e-01 2.222895e-08 9.723575e-01 1.937237e-01 

lik.d <- make.dec2(phy, states)
pars.d <- c(0.1, 0.1, 0.5, 0.5)
names(pars.d) <- argnames(lik.d)
ans.d <- find.mle(lik.d, pars.d)
coef(ans.d)
#        xA        xB        dA        dB 
# 1.1826472 3.7303945 3.8204625 0.3066245 

lik.d <- constrain(lik.d, xA ~ xB)
pars.d <- c(0.1, 0.5, 0.5)
#        xB        dA        dB 
# 1.5044522 1.8570890 0.3801415 

