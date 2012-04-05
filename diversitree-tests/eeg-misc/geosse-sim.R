library(diversitreeEEG)

fixInNamespace("tree.geosse", "diversitreeEEG")
fixInNamespace("tree.classe", "diversitreeEEG")

pars.geosse <- c(1.5, 0.5, 1.0, 0.7, 0.6, 1.4, 0.9)
names(pars.geosse) <- diversitreeEEG:::default.argnames.geosse()
pars.classe <- diversitreeEEG:::pars.ge.to.cl(pars.geosse)

diversitreeEEG:::stationary.freq.geosse(pars.geosse)
diversitreeEEG:::stationary.freq.classe(pars.classe, 3)

set.seed(1)
phy.g <- trees(pars.geosse, type="geosse", n=2, max.t=3, x0=0)
phy.g <- trees(pars.geosse, type="geosse", n=2, max.t=3)
set.seed(1)
phy.c <- trees(pars.classe, type="classe", n=2, max.t=3, x0=1)
phy.c <- trees(pars.classe, type="classe", n=2, max.t=3)

set.seed(3)
phy <- tree.geosse(pars.geosse, max.t=3, x0=1, include.extinct=T)
plot(history.from.sim.discrete(phy, 0:2), phy)
set.seed(3)
phy <- tree.classe(pars.classe, max.t=3, x0=2, include.extinct=T)
plot(history.from.sim.discrete(phy, 1:3), phy)
