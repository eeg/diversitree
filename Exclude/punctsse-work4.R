### Testing PunctSSE MLE and MCMC ###

library(diversitreeGP)

source("/home/emma/src/SimTreeSDD/packaged/SimTreeSDD-20101213/misc/ttn.R")

#--------------------------------------------------
# example trees
#--------------------------------------------------

# bisse-tree.ttn = ancestral_reconstruction/simulation/sim_output/table2/setE7/E7-998.ttn
# sim params = (1.4, 0.6, 0.2, 0.2, 0.5, 0), root state = 0
ttn <- read.ttn("bisse-tree.ttn", nodes=T)

# or

# geosse-tree.ttn = region_phylogeny/simtests2/set19/set19-97.ttn
# sim params = (1.5, 0.5, 1.0, 0.7, 0.7, 1.5, 1.5), root state = 0
ttn <- read.ttn("geosse-tree.ttn", nodes=T)

# or

# punct-tree.tt = nsf_proposal/prelim_results/punct/set4/tree-100.ttn
# sim params (lam = 6.2, p = 0.2, q = 0; symmetric):
# lam111 lam112 lam122 lam211 lam212 lam222    mu1    mu2    q12    q21 
#   4.96   1.24      0      0   1.24   4.96      3      0      0      0
ttn <- read.ttn("punct-tree.ttn", nodes=T)
# ttn$tree$node.label <- NULL

#--------------------------------------------------
# BiSSE MLE
#--------------------------------------------------

lnL.bisse <- make.bisse(ttn$tree, ttn$tip.states)
pars.bisse <- starting.point.bisse(ttn$tree)
#     lambda0     lambda1         mu0         mu1         q01         q10 
# 0.801018555 0.801018555 0.009104681 0.009104681 0.158382775 0.158382775 
lnL.bisse(pars.bisse)   # -313.8406
find.mle(lnL.bisse, pars.bisse)
#      lambda0      lambda1          mu0          mu1          q01          q10 
# 1.612922e+00 4.728302e-01 5.250026e-01 1.064626e-07 5.293502e-01 5.926068e-02 
# -283.0025

lnL.punctsse <- make.punctsse(ttn$tree, ttn$tip.states+1, 2)
pars.punctsse <- c(rep(0, 6), pars.bisse[-seq(2)])
names(pars.punctsse) <- argnames(lnL.punctsse)
pars.punctsse[c('lambda111', 'lambda222')] <- pars.bisse[1:2]
lnL.punctsse(pars.punctsse)

lnL.punct.bi <- constrain(lnL.punctsse, lambda112~0, lambda122~0, lambda211~0, lambda212~0)
pars.punct.bi <- pars.punctsse[-seq(2,5)]
lnL.punct.bi(pars.punct.bi)
find.mle(lnL.punct.bi, pars.punct.bi)
# agrees! (slower, though)

#--------------------------------------------------
# GeoSSE MLE
#--------------------------------------------------

lnL.geosse <- make.geosse(ttn$tree, ttn$tip.states)
pars.geosse <- starting.point.geosse(ttn$tree)
lnL.geosse(pars.geosse, condition.surv=T)   # -444.926
find.mle(lnL.geosse, pars.geosse, condition.surv=T)
#        sA        sB       sAB        xA        xB        dA        dB 
# 1.4898669 0.4136132 1.1075633 0.3826722 0.7675435 1.6632065 0.4956484 
# -385.7713

pars.punctsse <- 0 * starting.point.punctsse(ttn$tree, 3)
pars.punctsse['lambda222'] <- pars.punctsse['lambda112'] <- pars.geosse['sA']
pars.punctsse['lambda333'] <- pars.punctsse['lambda113'] <- pars.geosse['sB']
pars.punctsse['lambda123'] <-  pars.geosse['sAB']
pars.punctsse['mu2'] <- pars.punctsse['q13'] <- pars.geosse['xA']
pars.punctsse['mu3'] <- pars.punctsse['q12'] <- pars.geosse['xB']
pars.punctsse['q21'] <- pars.geosse['dA']
pars.punctsse['q31'] <- pars.geosse['dB']

lnL.punctsse <- make.punctsse(ttn$tree, ttn$tip.states+1, 3)
lnL.punct.geo <- constrain(lnL.punctsse, lambda111~0, lambda122~0, lambda133~0, lambda211~0, lambda212~0, lambda213~0, lambda223~0, lambda233~0, lambda311~0, lambda312~0, lambda313~0, lambda322~0, lambda323~0, mu1~0, q23~0, q32~0, lambda112~lambda222, lambda113~lambda333, q13~mu2, q12~mu3)
pars.punct.geo <- pars.punctsse[c(5,10,18, 20,21, 22,23)]
lnL.punct.geo(pars.punct.geo)
find.mle(lnL.punct.geo, pars.punct.geo)
# lambda123 lambda222 lambda333       mu2       mu3       q12       q13 
# 1.1073530 1.4895624 0.4136894 0.3818640 0.7677768 1.6634969 0.4944137 
# -385.7713

lnL1 <- constrain(lnL.punctsse, lambda111~0, lambda122~0, lambda133~0, lambda211~0, lambda212~0, lambda213~0, lambda223~0, lambda233~0, lambda311~0, lambda312~0, lambda313~0, lambda322~0, lambda323~0, mu1~0, q23~0, q32~0)
lnL2 <- constrain(lnL1, lambda112~lambda222, lambda113~lambda333, q13~mu2, q12~mu3)
lnL2(pars.punct.geo)

# note: constrain.from.pars(likefunc, params): a nice utility function would be to constrain to 0 all params with initial NA values (would need to be able to constrain further afterwards, but seems like that's possible now)

#--------------------------------------------------
# actually punctuated
#--------------------------------------------------

nstates <- 2

lnL.punctsse <- make.punctsse(ttn$tree, ttn$tip.states+1, 2)
lnL.punct <- constrain(lnL.punctsse, lambda122~0, lambda211~0)
pars.punct <- rep(1.0, 8)
names(pars.punct) <- argnames(lnL.punct)
find.mle(lnL.punct, pars.punct)
# lambda111 lambda112 lambda212 lambda222       mu1       mu2       q12       q21 
# 5.7552942 1.0757147 0.3100088 6.3911727 1.5046604 2.5542311 0.4264777 1.0239802 
# 57.29392

# correct that speciation involves character change only ~10-20% of the time
# incorrect that there is substantial gradual evolution
# incorrect that relative extinction rates
# definitely better than the old gpsddd code, at least
