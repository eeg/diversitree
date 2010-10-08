### Check my node.fixing-based BiSSE reconstructions with those Rich
### implemented.

### v0.5-3, Oct 8, 2010
### Note: defaults in diversitree's asr() [cf. asr.marginal.bisse()] changed after 0.4-4.

# could instead use /home/emma/src/miscR/ttn.R and plot.ttn()
source("/home/emma/src/miscR/ttn.R")
mycol <- c(rgb(0.832, 0.367, 0), rgb(0, 0.445, 0.695))  # orange, blue
ttn <- read.ttn("example-clean.ttn", nodes=F)

plot(ttn$tree, show.tip.label=F, no.margin=T)
tiplabels(text=as.character(seq(ttn$tree$Nnode+1)), bg=mycol[ttn$states+1], adj=0)
nodelabels(ttn$tree$node.label)

params <- c(1.4, 0.4, 0.2, 0.1, 0.6, 0.35)
names(params) <- c("lambda0", "lambda1", "mu0", "mu1", "q01", "q10")

#--------------------------------------------------
# Part I: asr.marginal (rgf)
#-------------------------------------------------- 

library(diversitree)

lnL <- make.bisse(ttn$tree, ttn$states)

ans <- asr.marginal(lnL, params)
ans <- data.frame(label=ttn$tree$node.label, p1=ans[2,])
i <- order(ans$label)
printme <- data.frame(ans$label[i], round(ans$p1[i], 6))
write.table(printme, file="rgf.anc", quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)

# the assumptions ROOT.EQUI and condition.surv=T are built into asr.marginal.bisse()

#--------------------------------------------------
# Part II: node.fixing (eeg)
#-------------------------------------------------- 

# This code, and the example tree, were taken from
# diversitreeEEG/exclude/reconstruct.R (Oct 2009).

# ROOT.OBS   = 3 (conditional likelihoods)
# ROOT.EQUI  = 2 (stationary)
# ROOT.FLAT  = 1 (equal weights)
# ROOT.GIVEN = 4 (give weights with root.p)

library(diversitreeNF)

# uniform root and conditioning on survival
ans <- reconstruct.bisse(ttn$tree, ttn$states, params, condition.surv=T, root=1)
ans <- data.frame(label=ttn$tree$node.label, ans)
i <- order(ans$label)
printme <- data.frame(ans$label[i], round(ans$p1[i], 6))
write.table(printme, file="uniformT-eeg.anc", quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)
# this agrees perfectly with rgf.anc (but note that it runs a lot slower)

# condlike root and not conditioning on survival
ans <- reconstruct.bisse(ttn$tree, ttn$states, params, condition.surv=F, root=3)
ans <- data.frame(label=ttn$tree$node.label, ans)
i <- order(ans$label)
printme <- data.frame(ans$label[i], round(ans$p1[i], 6))
write.table(printme, file="condlikeF-eeg.anc", quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)

# Why doesn't this agree with my old results?
# This is why:

lnL <- make.bisse(ttn$tree, ttn$states)

# 0.5-3 and 0.4-5
lnL(params, condition.surv=F, root=1)   # -121.9048
lnL(params, condition.surv=F, root=3)   # -121.3429
lnL(params, condition.surv=T, root=1)   # -121.7453
lnL(params, condition.surv=T, root=3)   # -121.3362

# 0.4-4
lnL(params, condition.surv=F, root=1)   # -122.0888
lnL(params, condition.surv=F, root=3)   # -121.7777
lnL(params, condition.surv=T, root=1)   # -121.7453
lnL(params, condition.surv=T, root=3)   # -121.4439

# The real work I've done (geosse, si) is all with 0.4-5 or later, anyway.
