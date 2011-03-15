### Profiling PunctSSE MLE ###

library(diversitreeGP)

source("/home/emma/src/SimTreeSDD/packaged/SimTreeSDD-20101213/misc/ttn.R")
ttn <- read.ttn("tests/geosse-tree.ttn", nodes=T)

lnL.geosse <- make.geosse(ttn$tree, ttn$tip.states)
pars.geosse <- starting.point.geosse(ttn$tree)
lnL.geosse(pars.geosse, condition.surv=T)   # -444.926

Rprof("geosse.out")
find.mle(lnL.geosse, pars.geosse, condition.surv=T)
Rprof(NULL)

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
pars.punct.geo <- pars.punctsse[c(5,10,18, 20,21, 24,26)]  # oops -- was 22,23!
lnL.punct.geo(pars.punct.geo)

Rprof("punctsse-geosse.out")
find.mle(lnL.punct.geo, pars.punct.geo)
Rprof(NULL)


prof.geosse <- summaryRprof("geosse.out")

print(str(prof.geosse))
# List of 4
#  $ by.self        :'data.frame':        47 obs. of  4 variables:
#   ..$ self.time : num [1:47] 39.6 37.6 30.5 24.6 18.5 ...
#   ..$ self.pct  : num [1:47] 16.9 16 13 10.5 7.9 ...
#   ..$ total.time: num [1:47] 52.3 234.3 75.1 154.7 233.6 ...
#   ..$ total.pct : num [1:47] 22.3 100 32.1 66 99.7 ...
#  $ by.total       :'data.frame':        58 obs. of  4 variables:
#   ..$ total.time: num [1:58] 234 234 234 234 234 ...
#   ..$ total.pct : num [1:58] 100 100 100 100 100 ...
#   ..$ self.time : num [1:58] 0 0 0 0 0 ...
#   ..$ self.pct  : num [1:58] 0 0 0 0 0 ...
#  $ sample.interval: num 0.02
#  $ sampling.time  : num 234

library(profr)

prof.geosse <- parse_rprof("geosse.out")
plot(prof.geosse)
sort(table(prof.geosse$f))
# tail end:
#            inherits     matrix.to.list      is.data.frame              .Call 
#                 224                233                261                298 
#          geosse.ode                sum                  t            rowSums 
#                 318                412                553               1146 
#  initial.conditions           branches 
#                1969               3761 

prof.punctsse <- parse_rprof("punctsse-geosse.out")
sort(table(prof.punctsse$f))
# tail end:
#        diag<-     unique.default          is.vector             unlist 
#           920               1295               1628               2536 
#       rowSums              .Call             unique               pmax 
#          2704               3004               3119               3359 
#  punctsse.ode                FUN                  t        seq.default 
#          4385               5289               5344               6752 
#        lapply initial.conditions             sapply                seq 
#          7950               8401               9782               9831 
#      branches 
#         15453 



# TODO
# time-trial just lnL computation, but for a large tree
# try reducing seq and l/sapply, especially in initial.conditions.punctsse -- see if that reduces runtime

pars.geosse[1:7] <- c(2, 3, 1.5, 0.7, 0.8, 1.1, 1.4)
# use other definitions from above

# -476.1513
lnL.geosse(pars.geosse, condition.surv=T)
lnL.punct.geo(pars.punct.geo, condition.surv=T)

system.time(lnL.geosse(pars.geosse, condition.surv=T))
  #  user  system elapsed 
  # 0.060   0.000   0.061
system.time(lnL.punct.geo(pars.punct.geo, condition.surv=T))
  #  user  system elapsed 
  # 0.170   0.000   0.177

fixInNamespace("initial.conditions.punctsse", "diversitreeGP")
#
#   old: d <- sapply(seq(n), get.di)
#   new: d <- sapply(nseq, get.di)
#   0.170   0.000   0.168
#
#   old: nseq <- seq(n)
#   new: nseq <- seq_len(n)
#   0.160   0.010   0.165
#
#   sapply isn't really vectorized
#   old: d <- sapply(nseq, get.di)
#   new: d <- unlist(lapply(nseq, get.di))
#   0.150   0.000   0.154 

fixInNamespace("make.branches.punctsse", "diversitreeGP")
#   old: idx.lm <- seq(x)
#   new: idx.lm <- seq_len(x)
#   no real difference

### This next looks like a big improvement, but messier to implement.
# store icp.* in cache? then pass to initial.conditions during call from all.branches
# consider also for make.branches.punctsse

# initial.conditions.punctsse
# pull out all the index computation that depends only on n
#   0.110   0.000   0.111 

icp.n <- 3
icp.nseq <- seq_len(icp.n)
icp.nlam <- icp.n*(icp.n+1)/2
icp.Didx <- (icp.n+1):(2*icp.n)
icp.j <- rep(icp.nseq, times=seq(icp.n,1,-1))
icp.k <- unlist(lapply(1:icp.n, function(i) icp.nseq[i:icp.n]))

initial.conditions.punctsse <- function(init, pars, t, is.root=FALSE) {
  e <- init[[1]][icp.nseq]
  DM <- init[[1]][icp.Didx]
  DN <- init[[2]][icp.Didx]
  get.di <- function(i)
  {
    x <- seq((i-1)*icp.nlam+1, i*icp.nlam)
    sum(pars[x] * 0.5 * (DM[icp.j] * DN[icp.k] + DM[icp.k] * DN[icp.j]))
  }
  d <- unlist(lapply(icp.nseq, get.di))
  c(e, d)
}
