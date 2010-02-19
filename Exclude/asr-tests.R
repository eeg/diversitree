### Check my node.fixing-based BiSSE reconstructions with those Rich
### implemented.
### v0.4-4, Feb 19, 2010

source("/home/emma/src/miscR/ttn.R")

ttn <- read.ttn("example-clean.ttn", nodes=F)
params <- c(1.4, 0.4, 0.2, 0.1, 0.6, 0.35)
names(params) <- c("lambda0", "lambda1", "mu0", "mu1", "q01", "q10")

# ROOT.OBS  = 3 (conditional likelihoods, the default)
# ROOT.EQUI = 2 (stationary)
# ROOT.FLAT = 1 (equal weightings)

#--------------------------------------------------
# Part I: node.fixing (eeg)
#-------------------------------------------------- 

# This code, and the example tree, were taken from
# diversitreeEEG/exclude/reconstruct.R (Oct 2009).

library(diversitreeNF)

# condlike root
ans <- reconstruct.bisse(ttn$tree, ttn$states, params)
ans <- data.frame(label=ttn$tree$node.label, ans)
i <- order(ans$label)
printme <- data.frame(ans$label[i], round(ans$p1[i], 6))
write.table(printme, file="condlike-eeg.anc", quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)

# stationary root
ans <- reconstruct.bisse(ttn$tree, ttn$states, params, root=2)
ans <- data.frame(label=ttn$tree$node.label, ans)
i <- order(ans$label)
printme <- data.frame(ans$label[i], round(ans$p1[i], 6))
write.table(printme, file="stationary-eeg.anc", quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)

# uniform root
ans <- reconstruct.bisse(ttn$tree, ttn$states, params, root=1)
ans <- data.frame(label=ttn$tree$node.label, ans)
i <- order(ans$label)
printme <- data.frame(ans$label[i], round(ans$p1[i], 6))
write.table(printme, file="uniform-nf.eeg", quote=FALSE, sep="\t",
            row.names=FALSE, col.names=FALSE)


#--------------------------------------------------
# Part II: Rich (rgf)
#-------------------------------------------------- 

library(diversitree)

# only allows condlike root
