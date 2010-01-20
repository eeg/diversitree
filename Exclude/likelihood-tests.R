library(diversitreeGSE)
packageDescription("diversitreeGSE")

tree <- read.tree(text="((((0:0.461876,(1:0.307717,(2:0.231825,3:0.231825):0.075892):0.154159):0.425922,((4:0.244819,5:0.244819):0.004749,6:0.249568):0.638231):0.142051,7:1.029850):0.038423,(((8:0.510933,(9:0.427929,(10:0.119778,11:0.119778):0.308151):0.083004):0.007428,(12:0.488316,13:0.488316):0.030044):0.100160,14:0.618521):0.449752);")

states <- c(0, 1, 0, 2, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0)
names(states) <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14")

starting.point.gse2(tree)
#        sA        sB       sAB        xA        xB        dA        dB 
# 1.5788510 1.5788510 1.5788510 0.0000000 0.0000000 0.3157702 0.3157702 
# with possible warning message


pars <- c(0.9, 0.8, 0.1, 0.2, 0.3, 0.5, 0.6)
lnL.7par <- make.gse2(tree, states)
names(pars) <- argnames(lnL.7par)
lnL.7par(pars)             # -23.39109

lnL <- constrain(lnL.7par, sA ~ sB)
lnL(pars[2:7])             # -23.93448

lnL <- constrain(lnL.7par, sAB ~ 0)
lnL(pars[-3])              # -23.22191

lnL <- constrain(lnL.7par, dB ~ dA)
lnL(pars[-6])              # -22.91692

lnL <- constrain(lnL.7par, dB ~ dA, sAB ~ 0)
lnL(pars[-c(3,7)])         # -23.44108

find.mle(lnL.7par, pars, method="subplex")
# $par
#           sA           sB          sAB           xA           xB           dA 
# 1.446606e+00 4.688561e-01 1.324631e-06 2.720090e-09 3.455288e-08 1.130573e+00 
#           dB 
# 2.564782e+00 
# $lnLik
# [1] -18.70781

find.mle(constrain(lnL.7par, dA ~ dB), pars[-6], method="subplex")
# $par
#           sA           sB          sAB           xA           xB           dB 
# 1.480193e+00 3.901173e-01 1.767147e-08 1.594058e-08 9.047256e-07 1.275068e+00 
# $lnLik
# [1] -18.84392

find.mle(constrain(lnL.7par, dA ~ dB, sAB ~ 0), pars[-c(3,7)], method="subplex")
# $par
#           sA           sB           xA           xB           dA 
# 1.480180e+00 3.901124e-01 2.422951e-06 8.505164e-06 1.275044e+00 
# $lnLik
# [1] -18.84393



lnL.7par(pars, root.p=c(0.9, 0.1, 0))   # warning
lnL.7par(pars, root.p=c(0.9, 0.1), root=diversitreeGSE:::ROOT.GIVEN)  # error
lnL.7par(pars, root.p=c(0.9, 0.1, 0), root=diversitreeGSE:::ROOT.GIVEN)
# -23.26696

lnL.7par(pars, root=diversitreeGSE:::ROOT.EQUI)    # -24.00041
lnL.7par(pars, root=diversitreeGSE:::ROOT.FLAT)    # -23.94198

lnL.7par(pars, condition.surv=FALSE)               # -23.45394
