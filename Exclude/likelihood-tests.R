library(diversitreeGSE)
library(help=diversitreeGSE)
packageDescription("diversitreeGSE")

tree <- read.tree(text="((((0:0.461876,(1:0.307717,(2:0.231825,3:0.231825):0.075892):0.154159):0.425922,((4:0.244819,5:0.244819):0.004749,6:0.249568):0.638231):0.142051,7:1.029850):0.038423,(((8:0.510933,(9:0.427929,(10:0.119778,11:0.119778):0.308151):0.083004):0.007428,(12:0.488316,13:0.488316):0.030044):0.100160,14:0.618521):0.449752);")

states <- c(0, 1, 0, 2, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0)
names(states) <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14")

starting.point.geosse(tree)
#        sA        sB       sAB        xA        xB        dA        dB 
# 1.5788510 1.5788510 1.5788510 0.0000000 0.0000000 0.3157702 0.3157702 
# (with possible warning message; gone in 0.4-5)

# 0.4-3 results in parentheses; default is condition.surv=T
# 0.4-5 shown now; condition.surv=F, lambda at root join
# 0.5-2 next
# 0.5-3b now

pars <- c(0.9, 0.8, 0.1, 0.2, 0.3, 0.5, 0.6)
lnL.7par <- make.geosse(tree, states)
names(pars) <- argnames(lnL.7par)
lnL.7par(pars)                    # -23.71574 (-23.39109)
lnL.7par(pars, condition.surv=T)  # -24.10259

lnL <- constrain(lnL.7par, sA ~ sB)
lnL(pars[2:7])                    # -24.34399 (-23.93448)
lnL(pars[2:7], condition.surv=T)  # -24.6769

lnL <- constrain(lnL.7par, sAB ~ 0)
lnL(pars[-3])                     # -23.56325 (-23.22191)
lnL(pars[-3], condition.surv=T)   # -23.90416

lnL <- constrain(lnL.7par, dB ~ dA)
lnL(pars[-6])                     # -23.2386 (-22.91692)

lnL <- constrain(lnL.7par, dB ~ dA, sAB ~ 0)
lnL(pars[-c(3,7)])                # -23.763 (-23.44108)

find.mle(lnL.7par, pars, method="subplex")
# (differs only for the tiny values)
# $par
#           sA           sB          sAB           xA           xB           dA 
# 1.592062e+00 4.066316e-01 7.031295e-10 8.239029e-07 3.233252e-09 1.262807e+00 
#           dB 
# 1.252328e+00 
# $lnLik
# [1] -18.77671
### 0.4.5
# $par
#           sA           sB          sAB           xA           xB           dA 
# 1.592060e+00 4.066320e-01 5.602806e-10 3.026119e-07 7.266436e-07 1.262808e+00 
#           dB 
# 1.252294e+00 
# $lnLik
# [1] -18.77671
### old
# $par
#           sA           sB          sAB           xA           xB           dA 
# 1.446606e+00 4.688561e-01 1.324631e-06 2.720090e-09 3.455288e-08 1.130573e+00 
#           dB 
# 2.564782e+00 
# $lnLik
# [1] -18.70781

find.mle(constrain(lnL.7par, dA ~ dB), pars[-6]) # method="subplex" now default
# $par
#           sA           sB          sAB           xA           xB           dB 
# 1.591789e+00 4.073864e-01 6.526636e-09 2.682779e-06 1.888036e-09 1.262040e+00 
# $lnLik
# [1] -18.77672
### old
# $par
#           sA           sB          sAB           xA           xB           dB 
# 1.480193e+00 3.901173e-01 1.767147e-08 1.594058e-08 9.047256e-07 1.275068e+00 
# $lnLik
# [1] -18.84392

find.mle(constrain(lnL.7par, dA ~ dB, sAB ~ 0), pars[-c(3,7)])
# $par
#           sA           sB           xA           xB           dA 
# 1.591711e+00 4.073796e-01 1.532276e-08 1.519441e-08 1.262004e+00 
# $lnLik
# [1] -18.77672
### old
# $par
#           sA           sB           xA           xB           dA 
# 1.480180e+00 3.901124e-01 2.422951e-06 8.505164e-06 1.275044e+00 
# $lnLik
# [1] -18.84393



lnL.7par(pars, root.p=c(0.9, 0.1, 0))   # warning
lnL.7par(pars, root.p=c(0.9, 0.1), root=diversitreeGSE:::ROOT.GIVEN)  # error
lnL.7par(pars, root.p=c(0.9, 0.1, 0), root=diversitreeGSE:::ROOT.GIVEN)
# -23.56753 (-23.26696)

lnL.7par(pars, root=diversitreeGSE:::ROOT.EQUI)    # -24.32152 (-24.00041)
lnL.7par(pars, root=diversitreeGSE:::ROOT.FLAT)    # -24.25945 (-23.94198)

lnL.7par(pars, condition.surv=FALSE)               # -23.71574 (-23.45394)

set.seed(1)
mcmc(lnL.7par, pars, nsteps=3, lower=0, upper=3, w=rep(3/5, 7), prior=make.prior.exponential(1))
# 1: {1.7079, 0.5958, 0.6397, 0.3059, 0.1992, 0.9057, 1.1472} -> -26.30203
# 2: {1.3404, 1.2697, 0.8503, 0.5133, 0.0089, 1.5693, 0.9488} -> -29.95031
# 3: {1.4314, 1.0871, 0.1138, 0.3022, 0.7059, 0.7149, 2.1334} -> -27.29315
### old:
# 1: {1.7079, 0.5958, 0.6397, 0.3059, 0.1992, 0.9057, 1.1472} -> -26.45726
# 2: {1.3404, 1.2697, 0.8503, 0.6557, 0.5654, 1.1507, 1.5016} -> -29.70205
# 3: {1.1591, 0.5990, 0.8141, 0.7422, 0.8274, 0.9435, 2.3320} -> -29.25331

set.seed(1)
mcmc(lnL.7par, pars, nsteps=3, upper=3, w=rep(3/5, 7), prior=make.prior.exponential(1))
1: {1.7079, 0.8226, 0.0854, 0.3320, 0.3547, 0.8980, 1.3405} -> -25.91188
2: {1.8432, 0.4183, 0.2167, 0.1520, 0.5122, 0.9473, 0.4756} -> -25.02481
3: {1.8248, 0.7531, 0.1291, 0.4015, 0.4574, 1.0977, 1.4606} -> -26.40207

set.seed(1)
mcmc(constrain(lnL.7par, dA ~ dB, sAB ~ 0), pars[-c(3,7)], nsteps=3, lower=0, upper=3, w=rep(3/5, 5), prior=make.prior.exponential(1))
# 1: {1.7079, 0.5958, 0.9004, 0.1992, 1.0660} -> -25.93986
# 2: {1.4176, 0.4126, 1.0722, 0.5315, 2.1264} -> -27.70291
# 3: {2.0392, 0.4766, 0.6001, 0.1083, 1.9912} -> -25.94858
