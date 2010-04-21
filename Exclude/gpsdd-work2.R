source("/home/emma/src/miscR/ttn.R")
mycol <- c(rgb(0.832, 0.367, 0), rgb(0, 0.445, 0.695))
ttn <- read.ttn("example-clean.ttn", nodes=F)

plot(ttn$tree, show.tip.label=F, no.margin=T)
tiplabels(text=as.character(seq(ttn$tree$Nnode+1)), bg=mycol[ttn$states+1], adj=0)
dev.off()

#--------------------------------------------------
# straight bisse
#-------------------------------------------------- 

library(diversitree)

lnL.1 <- make.bisse(ttn$tree, ttn$states)
params <- c(1.4, 0.4, 0.2, 0.1, 0.6, 0.35)
names(params) <- argnames(lnL.1)
lnL.1(params)       # -121.3362

#--------------------------------------------------
# gp-bisse
#-------------------------------------------------- 

library(diversitreeGP)

lnL.2 <- make.gpbisse(ttn$tree, ttn$states)

nstates <- 2
Lam <- array(0, dim=rep(nstates, 3))
dimnames(Lam) <- list(paste("p", seq(nstates), sep="."), paste("d1", seq(nstates), sep="."), paste("d2", seq(nstates), sep="."))
Lam[1,1,1] <- 1.4
Lam[2,2,2] <- 0.4
Lam[1,,]

Mu <- c(0.2, 0.1)

Q <- array(0, dim=rep(nstates, 2))
Q[1,2] <- 0.6
Q[2,1] <- 0.35
dimnames(Q) <- list(paste("from", seq(nstates), sep="."), paste("to", seq(nstates), sep="."))

params <- list(lambda=Lam, mu=Mu, q=Q, nstates=as.integer(nstates))

lnL.2(params)       # -121.3362
# agrees with make.bisse!
