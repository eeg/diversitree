### Check my node.fixing-based BiSSE reconstructions with those Rich
### implemented.
### v0.4-4, Feb 19, 2010

source("/home/emma/src/miscR/ttn.R")
mycol <- c(rgb(0.832, 0.367, 0), rgb(0, 0.445, 0.695))  # orange, blue

ttn <- read.ttn("example-clean.ttn", nodes=F)
plot(ttn$tree, show.tip.label=F, no.margin=T)
tiplabels(text=as.character(seq(ttn$tree$Nnode+1)), bg=mycol[ttn$states+1],
          adj=0)
nodelabels(ttn$tree$node.label)

params <- c(1.4, 0.4, 0.2, 0.1, 0.6, 0.35)
names(params) <- c("lambda0", "lambda1", "mu0", "mu1", "q01", "q10")

# ROOT.OBS   = 3 (conditional likelihoods, the default)
# ROOT.EQUI  = 2 (stationary)
# ROOT.FLAT  = 1 (equal weights)
# ROOT.GIVEN = 4 (give weights with root.p)

#--------------------------------------------------
# Part I: node.fixing (eeg)
#-------------------------------------------------- 

# This code, and the example tree, were taken from
# diversitreeEEG/exclude/reconstruct.R (Oct 2009).

library(diversitreeNF)

# condlike root
ans <- reconstruct.bisse(ttn$tree, ttn$states, params, condition.surv=F)
ans <- data.frame(label=ttn$tree$node.label, ans)
i <- order(ans$label)
printme <- data.frame(ans$label[i], round(ans$p1[i], 6))
write.table(printme, file="condlike-eeg.anc", quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)

# stationary root
ans <- reconstruct.bisse(ttn$tree, ttn$states, params, condition.surv=F, root=2)
ans <- data.frame(label=ttn$tree$node.label, ans)
i <- order(ans$label)
printme <- data.frame(ans$label[i], round(ans$p1[i], 6))
write.table(printme, file="stationary-eeg.anc", quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)

# uniform root
ans <- reconstruct.bisse(ttn$tree, ttn$states, params, condition.surv=F, root=1)
ans <- data.frame(label=ttn$tree$node.label, ans)
i <- order(ans$label)
printme <- data.frame(ans$label[i], round(ans$p1[i], 6))
write.table(printme, file="uniform-eeg.anc", quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)

# root fixed to state 0
ans <- reconstruct.bisse(ttn$tree, ttn$states, params, condition.surv=F,
                         root=4, root.p=c(1,0))
ans <- data.frame(label=ttn$tree$node.label, ans)
i <- order(ans$label)
printme <- data.frame(ans$label[i], round(ans$p1[i], 6))
write.table(printme, file="fix0-eeg.anc", quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)

#--------------------------------------------------
# Part II: Rich (rgf)
#-------------------------------------------------- 

library(diversitree)
?asr.marginal

# Looks like only condlike root is allowed for asr.marginal().
# The root.state argument only gets passed to asr.joint().

lnL <- make.bisse(ttn$tree, ttn$states)

# condlike root
ans <- asr.marginal(lnL, params)
ans <- data.frame(label=ttn$tree$node.label, p1=ans[2,])
i <- order(ans$label)
printme <- data.frame(ans$label[i], round(ans$p1[i], 6))
write.table(printme, file="condlike-rgf.anc", quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)
# agrees perfectly with my eeg

# stationary root
ans <- asr.marginal(lnL, params, root.state=2)
ans <- data.frame(label=ttn$tree$node.label, p1=ans[2,])
i <- order(ans$label)
printme <- data.frame(ans$label[i], round(ans$p1[i], 6))
write.table(printme, file="stationary-rgf.anc", quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)
# root.state appears to have no effect -- results identical to condlike
