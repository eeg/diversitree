### Profiling PunctSSE MLE ###

library(diversitreeGP)

source("/home/emma/src/SimTreeSDD/packaged/SimTreeSDD-20101213/misc/ttn.R")
ttn <- read.ttn("tests/geosse-tree.ttn", nodes=T)

lnL.geosse <- make.geosse(ttn$tree, ttn$tip.states)
pars.geosse <- starting.point.geosse(ttn$tree)
pars.geosse[1:7] <- c(1.4, 0.4, 1.1, 0.35, 0.75, 1.6, 0.5)  # to speed things up
lnL.geosse(pars.geosse, condition.surv=T)   # -386.0813

system.time(find.mle(lnL.geosse, pars.geosse, condition.surv=T))
#   user  system elapsed 
# 34.220   0.000  34.241 

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
pars.punct.geo <- pars.punctsse[c(5,10,18, 20,21, 24,26)]
lnL.punct.geo(pars.punct.geo)   # -386.0813

system.time(find.mle(lnL.punct.geo, pars.punct.geo, condition.surv=T))
# 49.280   0.010  49.347
#  user  system elapsed 
# 64.55    0.06   64.71

Rprof("punctsse-geosse.out")
find.mle(lnL.punct.geo, pars.punct.geo, condition.surv=T)
Rprof(NULL)

library(profr)

prof.geosse <- parse_rprof("geosse.out")
sort(table(prof.geosse$f))
#     matrix.to.list           inherits      is.data.frame              .Call 
#                 28                 38                 42                 52 
#         geosse.ode                sum                  t            rowSums 
#                 55                 57                 88                167 
# initial.conditions           branches 
#                281                527 

prof.punctsse <- parse_rprof("punctsse-geosse.out")
sort(table(prof.punctsse$f))
#     sum             diag<-              .Call                seq 
#      71                124                267                300 
# rowSums       punctsse.ode                FUN                  t 
#     399                422                462                535 
#  lapply             unlist initial.conditions           branches 
#     582                646                705               1404 



fixInNamespace("make.initial.conditions.punctsse", "diversitreeGP")

# try reshaping pars[lambdas] to allow real vectorization
lambdas <- matrix(pars[1:(n*n*(n+1)/2)], nrow=n, byrow=T)
d <- rowSums(lambdas * DM.DN)

lambda111 lambda112 lambda113 lambda122 lambda123 lambda133 lambda211 lambda212 
      0.0       1.4       0.4       0.0       1.1       0.0       0.0       0.0 
lambda213 lambda222 lambda223 lambda233 lambda311 lambda312 lambda313 lambda322 
      0.0       1.4       0.0       0.0       0.0       0.0       0.0       0.0 
lambda323 lambda333 
      0.0       0.4 
     [,1] [,2] [,3] [,4] [,5] [,6]
[1,]    0  1.4  0.4  0.0  1.1  0.0
[2,]    0  0.0  0.0  1.4  0.0  0.0
[3,]    0  0.0  0.0  0.0  0.0  0.4

[1] 0.11634555 0.12461659 0.10015749 0.13347296 0.10730019 0.08603294
     [,1]      [,2]       [,3]      [,4]      [,5]       [,6]
[1,]    0 0.1868621 0.04653822 0.0000000 0.1279801 0.00000000
[2,]    0 0.0000000 0.00000000 0.1502203 0.0000000 0.00000000
[3,]    0 0.0000000 0.00000000 0.0000000 0.0000000 0.03441318

i = 1
lams = pars[1:6]
DM.DN = c(0.11634555, 0.12461659, 0.10015749, 0.13347296, 0.10730019, 0.08603294)
sum(lams * DM.DN)

sum(lambdas[1,] * DM.DN)

rowSums(lambdas * DM.DN)

rowSums(lambdas * rbind(DM.DN, DM.DN, DM.DN))

matrix(DM.DN, byrow=T, nrow=n, ncol=6)
rowSums(lambdas * matrix(DM.DN, byrow=T, nrow=n, ncol=6))
d <- rowSums(lambdas * matrix(DM.DN, byrow=T, nrow=n, ncol=nlambdas/n)) # 52

d <- apply(lambdas, 1, function(i) sum(i * DM.DN))       # 62
for (i in nseq) d[i] <- sum(pars[lam.idx[i,]] * DM.DN)   # 48
d <- apply(lam.idx, 1, function(i) sum(pars[i] * DM.DN)) # 62
