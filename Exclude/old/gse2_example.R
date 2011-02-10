# July 27, 2010

library(diversitreeGSE)
packageDescription("diversitreeGSE")
# examples here are for version 0.4-5

#--------------------------------------------------
# data
#-------------------------------------------------- 

### an example tree
tree <- read.tree(text="((((0:0.461876,(1:0.307717,(2:0.231825,3:0.231825):0.075892):0.154159):0.425922,((4:0.244819,5:0.244819):0.004749,6:0.249568):0.638231):0.142051,7:1.029850):0.038423,(((8:0.510933,(9:0.427929,(10:0.119778,11:0.119778):0.308151):0.083004):0.007428,(12:0.488316,13:0.488316):0.030044):0.100160,14:0.618521):0.449752);")

### example character states
# note: 0 = AB, 1 = A, 2 = B
states <- c(0, 1, 0, 2, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0)
names(states) <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14")

#--------------------------------------------------
# likelihood functions
#-------------------------------------------------- 

### the likelihood function for the full (unconstrained) model
lnL.7par <- make.gse2(tree, states)

### an initial guess for the parameter values
#         sA   sB   sAB  xA   xB   dA   dB 
pars <- c(0.9, 0.8, 0.1, 0.2, 0.3, 0.5, 0.6)

# can use the class on lnL.7par to name the parameters
names(pars) <- argnames(lnL.7par)

# verify that the likelihood function evaluates with the specified parameters
lnL.7par(pars)
# (I get -23.71574)

### a constrained likelihood function
# remove sAB, make xA and xB the same
lnL.constr1 <- constrain(lnL.7par, sAB ~ 0, xB ~ xA)

# adjust the parameter vector accordingly
pars.constr1 <- pars[-c(3,5)]

# now the constrained likelihood should evaluate okay
lnL.constr1(pars.constr1)
# (I get -23.41976)

### another constrained likelihood function
# make sA and sB the same
lnL.constr2 <- constrain(lnL.7par, sB ~ sA)
pars.constr2 <- pars[-2]

#--------------------------------------------------
# parameter estimation
#-------------------------------------------------- 

### maximum likelihood

# subplex does pretty well, but of course best to try many initial values

find.mle(lnL.7par, pars, method="subplex")
# $par
#           sA           sB          sAB           xA           xB           dA 
# 1.592062e+00 4.066316e-01 7.031295e-10 8.239029e-07 3.233252e-09 1.262807e+00 
#           dB 
# 1.252328e+00 
# $lnLik
# [1] -18.77671

find.mle(lnL.constr1, pars.constr1, method="subplex")
# $par
#           sA           sB           xA           dA           dB 
# 1.592075e+00 4.066262e-01 1.546723e-10 1.262813e+00 1.252262e+00 
# $lnLik
# [1] -18.77671
# (not any worse than the full model)

find.mle(lnL.constr2, pars.constr2, method="subplex")
# $par
#           sA          sAB           xA           xB           dA           dB 
# 1.221990e+00 6.865968e-09 7.784889e-10 7.144535e-01 1.412542e+00 3.674809e+00 
# $lnLik
# [1] -19.36156
# (a bit worse than the full model, though not much power with such a small tree)

### markov chain monte carlo
# this is just an illustration of the syntax; many steps and iterative fiddling with lower/upper/w/prior are obviously required

mcmc(lnL.7par, pars, nsteps=3, lower=0, upper=3, w=rep(3/5, 7), prior=make.prior.exponential(1))
# return value is a data.frame with one line per step

mcmc(lnL.constr1, pars.constr1, nsteps=3, lower=0, upper=3, w=rep(3/5, 5), prior=make.prior.exponential(1))

mcmc(lnL.constr2, pars.constr2, nsteps=3, lower=0, upper=3, w=rep(3/5, 6), prior=make.prior.exponential(1))
