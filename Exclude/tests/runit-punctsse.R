tol <- 1e-6

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

# geosse-tree.ttn = region_phylogeny/simtests2/set19/set19-97.ttn
# sim params = (1.5, 0.5, 1.0, 0.7, 0.7, 1.5, 1.5), root state = 0
ttn2 <- read.ttn("geosse-tree.ttn")

# punct-tree.ttn = nsf_proposal/prelim_results/punct/set4/tree-100.ttn
# sim params (lam = 6.2, p = 0.2, q = 0; symmetric):
# lam111 lam112 lam122 lam211 lam212 lam222    mu1    mu2    q12    q21 
#   4.96   1.24      0      0   1.24   4.96      3      0      0      0
ttn3 <- read.ttn("punct-tree.ttn")

#--------------------------------------------------
# prepare likelihood functions
#--------------------------------------------------

lnL1.bisse <- make.bisse(ttn1$tree, ttn1$states)
lnL1.punctsse <- make.punctsse(ttn1$tree, ttn1$states+1, 2)
lnL1.punctsse2 <- constrain(lnL1.punctsse, lambda112~0, lambda122~0, lambda211~0, lambda212~0)

lnL2.geosse <- make.geosse(ttn2$tree, ttn2$states)
lnL2.punctsse <- make.punctsse(ttn2$tree, ttn2$states+1, 3)
lnL2.punctsse2 <- constrain(lnL2.punctsse, lambda111~0, lambda122~0,
                            lambda133~0, lambda211~0, lambda212~0, lambda213~0,
                            lambda223~0, lambda233~0, lambda311~0, lambda312~0,
                            lambda313~0, lambda322~0, lambda323~0,
                            lambda112~lambda222, lambda113~lambda333, 
                            mu1~0, q23~0, q32~0, q13~mu2, q12~mu3)

lnL3.punctsse <- make.punctsse(ttn3$tree, ttn3$states+1, 2)
lnL3.punctsse2 <- constrain(lnL3.punctsse, lambda122~0, lambda211~0)

#--------------------------------------------------
# prepare parameter vectors
#--------------------------------------------------

pars1.bisse <- c(1.4, 0.6, 0.2, 0.1, 0.5, 0.05)
pars1.punctsse <- c(rep(0, 6), pars1.bisse[-seq(2)])
names(pars1.punctsse) <- argnames(lnL1.punctsse)
pars1.punctsse[c('lambda111', 'lambda222')] <- pars1.bisse[1:2]

pars2.geosse <- c(1.5, 0.5, 1.0, 0.7, 0.7, 1.4, 1.3)
names(pars2.geosse) <- argnames(lnL2.geosse)
pars2.punctsse <- 0 * starting.point.punctsse(ttn2$tree, 3)
pars2.punctsse['lambda222'] <- pars2.punctsse['lambda112'] <- pars2.geosse['sA']
pars2.punctsse['lambda333'] <- pars2.punctsse['lambda113'] <- pars2.geosse['sB']
pars2.punctsse['lambda123'] <-  pars2.geosse['sAB']
pars2.punctsse['mu2'] <- pars2.punctsse['q13'] <- pars2.geosse['xA']
pars2.punctsse['mu3'] <- pars2.punctsse['q12'] <- pars2.geosse['xB']
pars2.punctsse['q21'] <- pars2.geosse['dA']
pars2.punctsse['q31'] <- pars2.geosse['dB']
pars2.geosse2 <- pars2.geosse[c(3,1,2,4:7)]

pars3.punct <- c(5, 1, 0.1, 0.2, 2, 4, 3, 2.5, 2.1, 2.2)
names(pars3.punct) <- argnames(lnL3.punctsse)

#--------------------------------------------------
# test functions
#--------------------------------------------------

test.lnL1 <- function()
{
    argvals <- list(condition.surv=T)
    ans <- -284.9825
    checkEquals(do.call(lnL1.bisse, c(list(pars1.bisse), argvals)), ans,
                tolerance=tol)
    checkEquals(do.call(lnL1.punctsse, c(list(pars1.punctsse), argvals)), ans,
                tolerance=tol)
    checkEquals(do.call(lnL1.punctsse2, c(list(pars1.bisse), argvals)), ans,
                tolerance=tol)

    argvals <- list(condition.surv=F)
    ans <- -284.9744
    checkEquals(do.call(lnL1.bisse, c(list(pars1.bisse), argvals)), ans,
                tolerance=tol)
    checkEquals(do.call(lnL1.punctsse, c(list(pars1.punctsse), argvals)), ans,
                tolerance=tol)
    checkEquals(do.call(lnL1.punctsse2, c(list(pars1.bisse), argvals)), ans,
                tolerance=tol)

    argvals <- list(condition.surv=T, root=ROOT.GIVEN, root.p=c(0.6, 0.4))
    ans <- -285.2641
    checkEquals(do.call(lnL1.bisse, c(list(pars1.bisse), argvals)), ans,
                tolerance=tol)
    checkEquals(do.call(lnL1.punctsse, c(list(pars1.punctsse), argvals)), ans,
                tolerance=tol)
    checkEquals(do.call(lnL1.punctsse2, c(list(pars1.bisse), argvals)), ans,
                tolerance=tol)
}

test.lnL2 <- function()
{
    argvals <- list(condition.surv=T)
    ans <- -387.1200
    checkEquals(do.call(lnL2.geosse, c(list(pars2.geosse), argvals)), ans,
                tolerance=tol)
    checkEquals(do.call(lnL2.punctsse, c(list(pars2.punctsse), argvals)), ans,
                tolerance=tol)
    checkEquals(do.call(lnL2.punctsse2, c(list(pars2.geosse2), argvals)), ans,
                tolerance=tol)

    argvals <- list(condition.surv=F, root=ROOT.GIVEN, root.p=c(0.5, 0.3, 0.2))
    ans <- -386.9637
    checkEquals(do.call(lnL2.geosse, c(list(pars2.geosse), argvals)), ans,
                tolerance=tol)
    checkEquals(do.call(lnL2.punctsse, c(list(pars2.punctsse), argvals)), ans,
                tolerance=tol)
    checkEquals(do.call(lnL2.punctsse2, c(list(pars2.geosse2), argvals)), ans,
                tolerance=tol)
}

test.lnL3 <- function()
{
    checkEquals(lnL3.punctsse(pars3.punct), 36.75776, tolerance=tol)
}
