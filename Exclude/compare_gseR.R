dyn.load('gse2-eqs.so')
is.loaded("gse2_initmod")

# from the end of gse2-eqs.R
solve.gse2.C <- function(y, t, pars, rtol=RTOL, atol=ATOL)
  gse2.ode(y, c(0, t), pars, rtol=rtol, atol=atol)[,-1]
solve.gse2 <- solve.gse2.C

# to get make.ode()
source("../R/bisse-eqs.R")

# modified from zzz.R
assign("gse2.ode", make.ode("gse2_derivs", "gse2-eqs", "gse2_initmod", 6))

# to get make.gse2()
source("../R/gse2.R")

# to get get.ordering() and check for naming conflicts
source("../R/bisse.R")

# to get constraining/fixing functions
source("../R/util.R")

# to get find.mle.bisse()
source("../R/mle.R")
source("../R/mle2.R")

library(deSolve)
library(ape)

# for read.ttn()
source("/home/emma/src/miscR/ttn.R")

#--------------------------------------------------
# test 1 (see worksheet.R in gseR)
#-------------------------------------------------- 

tree <- read.tree(text="((((0:0.461876,(1:0.307717,(2:0.231825,3:0.231825):0.075892):0.154159):0.425922,((4:0.244819,5:0.244819):0.004749,6:0.249568):0.638231):0.142051,7:1.029850):0.038423,(((8:0.510933,(9:0.427929,(10:0.119778,11:0.119778):0.308151):0.083004):0.007428,(12:0.488316,13:0.488316):0.030044):0.100160,14:0.618521):0.449752);")

states <- c(0, 1, 0, 2, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0)
names(states) <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14")

lnL.7par = make.gse2(tree, states)
pars = c(0.9, 0.8, 0, 0.2, 0.3, 0.5, 0.6)
lnL.7par(pars)
# here: -23.22191
# in worksheet.R: -23.97088

lnL.6par = fix.par(lnL.7par, c(NA, NA, 0, NA, NA, NA, NA))   # sAB = 0
pars = c(0.9, 0.8, 0.2, 0.3, 0.5, 0.6)
lnL.6par(pars)

#--------------------------------------------------
# test 2
#-------------------------------------------------- 

treefile = "/home/emma/region_phylogeny/simtests2/set9/set9-214.ttn"
treeinfo = read.ttn(treefile)

lnL.7par = make.gse2(treeinfo$tree, treeinfo$states)
lnL.6par = fix.par(lnL.7par, c(NA, NA, 0, NA, NA, NA, NA))   # sAB = 0
class(lnL.6par)

params = c(1.308310, 5.660880e-09, 0.4897199, 0.2917469, 1.328771, 1.883107e-07)
lnL.6par(params)
# here and in worksheet.R: -177.6847

find.mle(lnL.6par, params, lower=0)

# combined constraints
lnL.eqX = constrain.fix.par(lnL.7par, c(NA, NA, NA, NA, 4, NA, NA), c(NA, NA, 0, NA, NA, NA, NA))    # xB = xA, sAB = 0
params = c(1.2, 0.01, 0.3, 1.3, 0.01)   # sA, sB, xA=xB, dA, dB
names(params) = c("sA", "sB", "xAxB", "dA", "dB")
find.mle(lnL.eqX, params, lower=0)

params = c(1.329988, 0.000000, 0, 0.521361, 0.521361, 1.746738, 0.000000)
lnL.7par(params)
# also -177.8142
