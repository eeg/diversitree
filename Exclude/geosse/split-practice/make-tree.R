library(diversitree)
library(geiger)
set.seed(2)

#           lam1  lam2  lam3  mu1  mu2  mu3  q12  q13  q21, q23, q31, q32
pars.bg <- c(0.1, 0.6, 0.9, 0.01, 0.02, 0.04, 0.3, 0.2, 0.5, 0.1, 0.6, 0.1)
pars.fg <- c(0.1, 0.6, 0.9, 0.01, 0.02, 0.04, 0.5, 0.1, 0.3, 0.2, 0.6, 0.1)  # and will rescale

tree.bg <- tree.musse(pars.bg, max.t=10, x0=2)
max(branching.times(tree.bg))  # 8.471023
tree.fg <- tree.musse(pars.bg*2, max.t=5, x0=2)
max(branching.times(tree.fg))  # 3.558046

tree.fg <- rescaleTree(tree.fg, 2)
tree.fg$tip.label <- paste("fg", tree.fg$tip.label, sep="")
tree.fg$node.label <- paste("fg", tree.fg$node.label, sep="")
names(tree.fg$tip.state) <- tree.fg$tip.label

write.tree(tree.bg, file="back.tre")
write.tree(tree.fg, file="fore.tre")
write.table(tree.bg$tip.state, file="backfore.dat", col.names=F, quote=F, append=F)
write.table(tree.fg$tip.state, file="backfore.dat", col.names=F, quote=F, append=T)

phy <- read.tree(file="bg-fg.tre")
plot(phy)
is.ultrametric(phy)





# AB  A  B
#  0  1  2

# lamAB lamA  lamB    -   xA   xB   xB   xA   dA    -    dB   -
#     0  0.6   0.9    -  0.2  0.3  0.3  0.2   0.5   -   0.6   -

diversitree:::argnames.musse(NULL, 3)

# lam1  lam2  lam3  mu1  mu2  mu3  q12  q13  q21, q23, q31, q32
#  0.1   0.6   0.9  0.1  0.2  0.4  0.3  0.2  0.5  0.1  0.6  0.1

#           lam1  lam2  lam3  mu1  mu2  mu3  q12  q13  q21, q23, q31, q32
pars.bg <- c(0.1, 0.6, 0.9, 0.01, 0.02, 0.04, 0.3, 0.2, 0.5, 0.1, 0.6, 0.1)
tree.bg <- tree.musse(pars.bg, max.t=10, x0=2)
Ntip(tree.bg)   # 215
plot(tree.bg)
max(branching.times(tree.bg))  # 8.471023


# replace a "slow" clade with a rapid radiation

tree.fg <- tree.musse(pars.bg*2, max.t=5, x0=2)
max(branching.times(tree.fg))  # 3.558046
tree.fg2 <- rescaleTree(tree.fg, 2)

tree.fg$tip.label <- paste("fg", tree.fg$tip.label, sep="")
tree.fg$node.label <- paste("fg", tree.fg$node.label, sep="")

write.tree(tree.bg, file="bg.tree")
write.tree(tree.fg2, file="fg.tree")

phy <- read.tree(file="bg-fg.tre")
plot(phy)
