library(diversitreeGP)

# options(error=recover)

#--------------------------------------------------
# read in test data
#--------------------------------------------------

read.ttn <- function(treefile)
{
    treestr <- readLines(treefile, n=1)
    tree <- read.tree(text=treestr)
    ntips <- Ntip(tree)

    all.states <- read.table(treefile, skip=1)
    tip.states <- all.states[1:ntips, 2]
    names(tip.states) <- all.states[1:ntips, 1]

    return(list(tree = tree, states = tip.states))
}

# bisse-tree.ttn = ancestral_reconstruction/simulation/sim_output/table2/setE7/E7-998.ttn
# sim params = (1.4, 0.6, 0.2, 0.2, 0.5, 0), root state = 0
ttn1 <- read.ttn("bisse-tree.ttn")

#--------------------------------------------------
# prepare likelihood functions and parameter vectors
#--------------------------------------------------

lnL1.C.musse  <- make.musse( ttn1$tree, ttn1$states+1, 2, control=list(backend="CVODES"))
lnL1.C.classe <- make.classe(ttn1$tree, ttn1$states+1, 2, control=list(backend="CVODES"))

pars1.bisse <- c(1.4, 0.6, 0.2, 0.1, 0.5, 0.05)
pars1.musse <- pars1.bisse
pars1.classe <- c(rep(0, 6), pars1.bisse[-seq(2)])
names(pars1.classe) <- argnames(lnL1.C.classe)
pars1.classe[c('lambda111', 'lambda222')] <- pars1.bisse[1:2]

#--------------------------------------------------
# testing
#--------------------------------------------------

lnL1.C.classe(pars1.classe) # should be -284.9406

lnL1.C.musse(pars1.musse)
lnL1.C.classe(pars1.classe)


# insert browser()
fixInNamespace("make.all.branches.C", "diversitreeGP")
fixInNamespace("make.all.branches.C.classe2", "diversitreeGP")

lnL1.classe <- make.classe(ttn1$tree, ttn1$states+1, 2)
lnL1.classe(pars1.classe) # -284.9406 (even when re-run)

lnL1.c.classe <- make.classe(ttn1$tree, ttn1$states+1, 2, control=list(backend="cvodes"))
lnL1.c.classe(pars1.classe) # -284.9406 (even when re-run)


#--------------------------------------------------
# timing
#--------------------------------------------------

system.time(fit2.classe <- find.mle(lnL2.classe, pars2.classe, method="subplex")) # 531.38
system.time(fit2.c.classe <- find.mle(lnL2.c.classe, pars2.classe, method="subplex")) # 544.773
system.time(fit2.C.classe <- find.mle(lnL2.C.classe, pars2.classe, method="subplex")) # 129.936
# but all gave up...maxit


system.time(fit2 <- find.mle(lnL2.classe2, pars2.geosse2, method="subplex"))        # 75.078
system.time(fit2.c <- find.mle(lnL2.c.classe2, pars2.geosse2, method="subplex"))    # 84.885
system.time(fit2.C <- find.mle(lnL2.C.classe2, pars2.geosse2, method="subplex"))    # 13.557
#       sAB        sA        sB        xA        xB        dA        dB 
# 1.0950959 1.4796289 0.4138699 0.3508881 0.8032363 1.7006286 0.4207858 
# [1] -385.89

lnL2.geosse <- make.geosse(ttn2$tree, ttn2$states)
lnL2.c.geosse <- make.geosse(ttn2$tree, ttn2$states, control=list(backend="cvodes"))
lnL2.C.geosse <- make.geosse(ttn2$tree, ttn2$states, control=list(backend="CVODES"))

system.time(fit2 <- find.mle(lnL2.geosse, pars2.geosse, method="subplex"))          # 189.573
#        sA        sB       sAB        xA        xB        dA        dB 
# 1.4540359 0.4188298 1.1107135 0.2841764 0.7553935 1.6754081 0.3790237 
# [1] -385.3282

system.time(fit2.c <- find.mle(lnL2.c.geosse, pars2.geosse, method="subplex"))      # 99.303
system.time(fit2.C <- find.mle(lnL2.C.geosse, pars2.geosse, method="subplex"))      # 23.845
