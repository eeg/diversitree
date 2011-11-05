library(ape)
source("/home/emma/src/SimTreeSDD/packaged/SimTreeSDD-20101213/misc/ttn.R")

ttn <- read.ttn("tree.ttn")
plot.ttn(ttn)

max(branching.times(ttn$tree))
# tree0: 3.739662
# tree1: 2.940327
# tree2: 2.796372

ttn <- read.ttn("tree0.ttn")
ttn$tree$tip.label <- paste("0t", ttn$tree$tip.label, sep="")
names(ttn$tip.states) <- paste("0t", names(ttn$tip.states), sep="")
ttn$tree$node.label <- paste("0n", ttn$tree$node.label, sep="")
names(ttn$node.states) <- paste("0n", names(ttn$node.states), sep="")
ttn$states <- ttn$tip.states
write.ttn(ttn, "tree0a.ttn")

ttn <- read.ttn("tree1.ttn")
ttn$tree$tip.label <- paste("1t", ttn$tree$tip.label, sep="")
names(ttn$tip.states) <- paste("1t", names(ttn$tip.states), sep="")
ttn$tree$node.label <- paste("1n", ttn$tree$node.label, sep="")
names(ttn$node.states) <- paste("1n", names(ttn$node.states), sep="")
ttn$states <- ttn$tip.states
write.ttn(ttn, "tree1a.ttn")

ttn <- read.ttn("tree2.ttn")
ttn$tree$tip.label <- paste("2t", ttn$tree$tip.label, sep="")
names(ttn$tip.states) <- paste("2t", names(ttn$tip.states), sep="")
ttn$tree$node.label <- paste("2n", ttn$tree$node.label, sep="")
names(ttn$node.states) <- paste("2n", names(ttn$node.states), sep="")
ttn$states <- ttn$tip.states
write.ttn(ttn, "tree2a.ttn")


phy <- read.tree(file="tree012.tre")
is.ultrametric(phy)
