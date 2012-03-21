################################
# Test the ClaSSE tree simulator
# 19 Mar 2012
################################

# attributes to compare, as distributions across many simulated trees:
#     total number of tips
#     proportion of tips in each state
#     mean of branching times?
#     tree balance statistics
#         ape's balance()
#         apTreeshape -- see http://nunn.rc.fas.harvard.edu/groups/pica/wiki/0c10a/84_Tree_balance_statistics.html

library(diversitree)
# fixInNamespace("make.tree.classe", "diversitree")
library(apTreeshape)

# Use this rather than trees() to get trees recorded in files (set.seed would
# be okay instead, but faster to read in than to simulate), and to avoid memory
# problems from too many unusually large trees.
# Adjust directory names afterwards.
# trees.to.files <- function(pars, type, max.t, ntrees=1000)
trees.to.files <- function(pars, type, ntrees=1000, ...)
{
    sim.func <- switch(type, bisse=tree.bisse, classe=tree.classe, 
                       bisseness=tree.bisseness, musse=tree.musse)
    i <- j <- 0
    while(i < ntrees)
    {
        # phy <- sim.func(pars, max.t=max.t)
        phy <- sim.func(pars, ...)
        if (class(phy) == "phylo")
        {
            if(Ntip(phy) > 4)
            {
                i <- i + 1
                filename <- paste(type, "/tree-", i, ".ttn", sep="")
                write.tree(phy, filename)
                write.table(phy$tip.state, file=filename, append=TRUE, quote=FALSE, 
                            sep="\t", row.names=TRUE, col.names=FALSE)
            }
        } else j <- j + 1
    }
    print(paste(i, "trees written,", j, "trees discarded"))
}

read.treefile <- function(filename)
{
    treestr <- readLines(filename, n=1)
    phy <- read.tree(text=treestr)

    states <- read.table(filename, skip=1)
    phy$tip.state <- states[,2]
    names(phy$tip.state) <- states[,1]

    return(phy)
}

get.tree.stats <- function(filelist, state1)
{
    ntrees <- length(filelist)
    empty <- rep(NA, ntrees)
    stats <- data.frame(ntips=empty, freq1=empty, meanbrt=empty, varbrt=empty, colless=empty, sackin=empty)
    for (i in seq(ntrees))
    {
        phy <- read.treefile(filelist[i])
        phy$tip.state <- phy$tip.state[1:Ntip(phy)] # to fix SimTreeSDD trees
        stats$ntips[i] <- Ntip(phy)
        stats$freq1[i] <- sum(phy$tip.state == state1) / Ntip(phy)
        stats$meanbrt[i] <- mean(branching.times(phy))
        stats$varbrt[i] <- var(branching.times(phy))

        aphy <- as.treeshape(phy)
        # "yule" is probably better than not
        stats$colless[i] <- colless(aphy, norm="yule")
        stats$sackin[i] <- sackin(aphy, norm="yule")
    }
    return(stats)
}

plot.tree.stats <- function(filename, dat1, dat2)
{
    mycolors <- c(rgb(0, 0.445, 0.695), rgb(0.832, 0.367, 0)) # blue, orange
    stats <- names(dat1)
    pdf(filename)
    for (i in seq(length(stats)))
    {
        x1 <- dat1[,i]
        x2 <- dat2[,i]
        plot(density(x1), xlim=range(c(x1, x2)), main=stats[i], xlab="", col=mycolors[1], lwd=2)
        lines(density(x2), col=mycolors[2], lwd=2)
    }
    dev.off()
}

#--------------------------------------------------
# 1. Compare with BiSSE
#--------------------------------------------------

### simulate

pars1.bisse <- c(1.4, 0.6, 0.2, 0.1, 0.5, 0.05)
pars1.classe <- c(rep(0, 6), pars1.bisse[-seq(2)])
names(pars1.classe) <- diversitree:::argnames.classe(NULL, 2)
pars1.classe[c('lambda111', 'lambda222')] <- pars1.bisse[1:2]

set.seed(1)
trees.to.files(pars1.bisse, "bisse", max.t=8)
set.seed(1)
trees.to.files(pars1.classe, "classe", max.t=8)

### summarize

bisse.stats <- get.tree.stats(sort(dir("1-bisse", full.names=T)), 0)
classe.stats <- get.tree.stats(sort(dir("1-classe", full.names=T)), 1)
plot.tree.stats("1-bisse-classe.pdf", bisse.stats, classe.stats)
# agreement is good on all fronts

#--------------------------------------------------
# 2-5. Compare with BiSSE-ness
#--------------------------------------------------

### simulate

bn.to.cl <- function(pars)
{
    pars.cl <- pars * NA
    names(pars.cl) <- diversitree:::argnames.classe(NULL, 2)
    pars.cl["lambda111"] <- pars["lambda0"] * (1 - pars["p0c"])
    pars.cl["lambda112"] <- pars["lambda0"] * pars["p0c"] * pars["p0a"]
    pars.cl["lambda122"] <- pars["lambda0"] * pars["p0c"] * (1 - pars["p0a"])
    pars.cl["lambda211"] <- pars["lambda1"] * pars["p1c"] * (1 - pars["p1a"])
    pars.cl["lambda212"] <- pars["lambda1"] * pars["p1c"] * pars["p1a"]
    pars.cl["lambda222"] <- pars["lambda1"] * (1 - pars["p1c"])
    pars.cl[c("mu1", "mu2", "q12", "q21")] <- pars[c("mu0", "mu1", "q01", "q10")]
    return(pars.cl)
}

pars2.bness <- c(1.4, 0.6, 0.2, 0.1, 0.5, 0.05, 0.2, 0.15, 0.22, 0.08)
names(pars2.bness) <- diversitree:::argnames.bisseness(NULL)
pars2.classe <- bn.to.cl(pars2.bness)

set.seed(1)
trees.to.files(pars2.bness, "bisseness", max.t=7)
set.seed(1)
trees.to.files(pars2.classe, "classe", max.t=7)

### summarize

bness.stats <- get.tree.stats(sort(dir("2-bisseness", full.names=T)), 0)
classe.stats <- get.tree.stats(sort(dir("2-classe", full.names=T)), 1)
plot.tree.stats("2-bness-classe.pdf", bness.stats, classe.stats)
# agreement on balance stats doesn't seem quite good enough

### To try to isolate a possible problem, remove some processes.

# 3. Ensure that all change is "asymmetric".
pars3.bness <- c(1.4, 0.6, 0.2, 0.1, 0.5, 0.05, 0.2, 1, 0.22, 1)
names(pars3.bness) <- diversitree:::argnames.bisseness(NULL)
pars3.classe <- bn.to.cl(pars3.bness)

set.seed(1)
trees.to.files(pars3.bness, "bisseness", max.t=7)
set.seed(1)
trees.to.files(pars3.classe, "classe", max.t=7)

bness.stats <- get.tree.stats(sort(dir("3-bisseness", full.names=T)), 0)
classe.stats <- get.tree.stats(sort(dir("3-classe", full.names=T)), 1)
plot.tree.stats("3-bness-classe.pdf", bness.stats, classe.stats)

# 4. Ensure that all change is "symmetric". (also increased p1c)
pars4.bness <- c(1.4, 0.6, 0.2, 0.1, 0.5, 0.05, 0.2, 0, 0.7, 0)
names(pars4.bness) <- diversitree:::argnames.bisseness(NULL)
pars4.classe <- bn.to.cl(pars4.bness)

set.seed(1)
trees.to.files(pars4.bness, "bisseness", max.t=7)
set.seed(1)
trees.to.files(pars4.classe, "classe", max.t=7)

bness.stats <- get.tree.stats(sort(dir("4-bisseness", full.names=T)), 0)
classe.stats <- get.tree.stats(sort(dir("4-classe", full.names=T)), 1)
plot.tree.stats("4-bness-classe.pdf", bness.stats, classe.stats)

# 5. Repeat original, but with increased p1c

pars5.bness <- c(1.4, 0.6, 0.2, 0.1, 0.5, 0.05, 0.2, 0.15, 0.7, 0.08)
names(pars5.bness) <- diversitree:::argnames.bisseness(NULL)
pars5.classe <- bn.to.cl(pars5.bness)

set.seed(1)
trees.to.files(pars5.bness, "bisseness", max.t=7)
set.seed(1)
trees.to.files(pars5.classe, "classe", max.t=7)

bness.stats <- get.tree.stats(sort(dir("5-bisseness", full.names=T)), 0)
classe.stats <- get.tree.stats(sort(dir("5-classe", full.names=T)), 1)
plot.tree.stats("5-bness-classe.pdf", bness.stats, classe.stats)

# Doesn't look so bad overall.
# Also simulated a second set of trees under pars2, but with set.seed(2).
#   Results are appended to 2-bness-classe.pdf

#--------------------------------------------------
# 6. Compare with MuSSE
#--------------------------------------------------

# 4 states
mu.to.cl <- function(pars)
{
    lam <- array(0, dim=c(4,4,4))
    lam[1,1,1] <- pars[1]
    lam[2,2,2] <- pars[2]
    lam[3,3,3] <- pars[3]
    lam[4,4,4] <- pars[4]
    parlist <- list(lambda=lam, mu=pars[5:8], q=matrix(0,4,4), nstates=4)

    pars.cl <- diversitree:::flatten.pars.classe(parlist)
    pars.cl[45:56] <- pars[9:20]
    return(pars.cl)
}

pars6.musse <- c( 0.2,0.4,0.6,0.8, 0.08,0.06,0.04,0.02, 0.1,0.15,0.2, 0.25,0.3,0.35, 0.4,0.45,0.5, 0.11,0.22,0.33 )
names(pars6.musse) <- diversitree:::argnames.musse(NULL, 4)
pars6.classe <- mu.to.cl(pars6.musse)

set.seed(1)
trees.to.files(pars6.musse, "musse", max.t=10, x0=1)
set.seed(1)
trees.to.files(pars6.classe, "classe", max.t=10, x0=1)

musse.stats <- get.tree.stats(sort(dir("6-musse", full.names=T)), 1)
classe.stats <- get.tree.stats(sort(dir("6-classe", full.names=T)), 1)
plot.tree.stats("6-musse-classe.pdf", musse.stats, classe.stats)

#--------------------------------------------------
# 7. Compare with GeoSSE
#--------------------------------------------------

ge.to.cl <- function(pars)
{
    pars.cl <- rep(0, 27)
    names(pars.cl) <- diversitree:::argnames.classe(NULL, 3)
    pars.cl['lambda222'] <- pars.cl['lambda112'] <- pars['sA']
    pars.cl['lambda333'] <- pars.cl['lambda113'] <- pars['sB']
    pars.cl['lambda123'] <-  pars['sAB']
    pars.cl['mu2'] <- pars.cl['q13'] <- pars['xA']
    pars.cl['mu3'] <- pars.cl['q12'] <- pars['xB']
    pars.cl['q21'] <- pars['dA']
    pars.cl['q31'] <- pars['dB']
    return(pars.cl)
}

pars7.geosse <- c(1.5, 0.5, 1.0, 0.7, 0.6, 1.4, 0.9)
names(pars7.geosse) <- diversitree:::argnames.geosse(NULL)
pars7.classe <- ge.to.cl(pars7.geosse)
diversitree:::stationary.freq.classe(pars7.classe, 3)   # adjust params so freq1 is more different from freq3

set.seed(1)
trees.to.files(pars7.classe, "classe", max.t=4.5, x0=1)
# geosse trees simulated with SimTreeSDD

geosse.stats <- get.tree.stats(sort(dir("7-geosse", full.names=T)), 0)
classe.stats <- get.tree.stats(sort(dir("7-classe", full.names=T)), 1)
plot.tree.stats("7-geosse-classe.pdf", geosse.stats, classe.stats)

# plotting
phy <- tree.classe(pars7.classe, max.t=4, x0=1, include.extinct=T)
plot(history.from.sim.discrete(phy, 1:3), phy)

