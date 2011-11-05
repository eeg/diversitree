library(diversitree)

phy <- read.tree(file="backfore.tre")
tmp <- read.table(file="backfore.dat")
states <- tmp[,2]
names(states) <- tmp[,1]

#--------------------------------------------------
# musse
#--------------------------------------------------

lik.m <- make.musse(phy, states, 3)
lik.s <- make.musse.split(phy, states, 3, c("fgnd1"), split.t=c(3.5))

pars.bg <- c(0.1, 0.6, 0.9, 0.01, 0.02, 0.04, 0.3, 0.2, 0.5, 0.1, 0.6, 0.1)
pars.fg <- c(0.1, 0.6, 0.9, 0.01, 0.02, 0.04, 0.5, 0.1, 0.3, 0.2, 0.6, 0.1) * 1.5

lik.m(pars.bg)              # -648.0743
lik.s(c(pars.bg, pars.bg))  # -648.0743
lik.s(c(pars.bg, pars.fg))  # -635.2749

ans.m <- find.mle(lik.m, pars.bg, method="subplex")
# 'log Lik.' -640.7726 (df=12)
#      lambda1      lambda2      lambda3          mu1          mu2          mu3 
# 1.162621e-01 5.358947e-01 1.312242e+00 4.773151e-01 2.534687e-08 1.662715e-01 
#          q12          q13          q21          q23          q31          q32 
# 1.215177e-01 3.446764e-01 4.781790e-01 2.603423e-02 1.000956e+00 1.630480e-01 

ans.s <- find.mle(lik.s, c(pars.bg, pars.fg), method="subplex")
# 'log Lik.' -616.7279 (df=24)
#    lambda1.1    lambda2.1    lambda3.1        mu1.1        mu2.1        mu3.1 
# 7.841177e-02 5.236306e-01 1.221682e+00 8.231421e-09 4.998364e-08 4.896992e-01 
#        q12.1        q13.1        q21.1        q23.1        q31.1        q32.1 
# 1.939652e-01 3.275805e-01 3.359429e-01 4.440101e-02 6.282391e-01 2.833790e-03 
#    lambda1.2    lambda2.2    lambda3.2        mu1.2        mu2.2        mu3.2 
# 4.240050e-09 1.611740e+00 3.047389e+00 1.237229e-10 8.419144e-07 2.868245e-08 
#        q12.2        q13.2        q21.2        q23.2        q31.2        q32.2 
# 1.005687e+00 1.048794e+01 2.267067e+00 3.459829e-08 1.558793e+01 2.847778e-08 

anova(ans.m, split=ans.s)
#       Df   lnLik    AIC  ChiSq Pr(>|Chi|)    
# full  12 -640.77 1305.5                      
# split 24 -616.73 1281.5 48.089  3.016e-06 ***

#--------------------------------------------------
# geosse
#--------------------------------------------------

lik.g <- make.geosse(phy, states-1)
lik.s <- make.geosse.split(phy, states-1, c("fgnd1"), split.t=c(3.5))

pars1 <- c(0.1, 0.6, 0.6, 0.02, 0.04, 0.5, 0.6)
pars2 <- c(0.1, 0.6, 0.6, 0.02, 0.04, 0.6, 0.5) * 1.5

lik.g(pars1)    # -730.9623
ans.g <- find.mle(lik.g, pars1, method="subplex")
# 'log Lik.' -679.2592 (df=7)
#           sA           sB          sAB           xA           xB           dA           dB 
# 2.944270e-01 6.582528e-01 2.423275e-09 3.337100e-01 3.082441e-01 2.856723e-01 9.543262e-01 

lik.s(c(pars1, pars1))    # -730.9623
lik.s(c(pars1, pars2))    # -719.3296

ans.s <- find.mle(lik.s, c(pars1, pars2), method="subplex")
ans.s <- find.mle(lik.s, coef(ans.s), method="subplex")
# 'log Lik.' -656.3848 (df=14)
#         sA.1         sB.1        sAB.1         xA.1         xB.1         dA.1          dB.1
# 2.555688e-01 5.064676e-01 1.074906e-02 1.433254e-01 1.483200e-01 1.964381e-01  6.035578e-01
#         sA.2         sB.2        sAB.2         xA.2         xB.2         dA.2         dB.2 
# 8.695977e-01 1.338836e+00 2.694903e-08 1.271626e+00 4.583919e-01 1.252704e+00 3.115052e+00 

#--------------------------------------------------
# more geosse
#--------------------------------------------------

source("/home/emma/src/SimTreeSDD/packaged/SimTreeSDD-20101213/misc/ttn.R")
ttn <- read.ttn("tree012.ttn")

pars0 <- c(1, 1, 0.5, 0.1, 0.1, 1, 1)
pars1 <- c(1.5, 0.5, 0.5, 0.1, 0.1, 1, 1)
pars2 <- c(1, 1, 0.5, 0.1, 0.1, 0.5, 1.5)

lik.0 <- make.geosse(ttn$tree, ttn$tip.states)
lik.01 <- make.geosse.split(ttn$tree, ttn$tip.states, nodes=c("1n130"), split.t=c(0))
lik.02 <- make.geosse.split(ttn$tree, ttn$tip.states, nodes=c("2n234"), split.t=c(0))
lik.012 <- make.geosse.split(ttn$tree, ttn$tip.states, nodes=c("1n130", "2n234"), split.t=c(0, 0))

lik.0(pars0)                    # -734.7009
lik.01(c(pars0, pars1))         # -727.531
lik.02(c(pars0, pars2))         # -733.0873
lik.012(c(pars0, pars1, pars2)) # -725.9247

lik.01c <- constrain(lik.01, sAB.2~sAB.1, xA.2~xA.1, xB.2~xB.1)
lik.02c <- constrain(lik.02, sAB.2~sAB.1, xA.2~xA.1, xB.2~xB.1)
lik.012c <- constrain(lik.012, sAB.2~sAB.1, xA.2~xA.1, xB.2~xB.1, sAB.3~sAB.1, xA.3~xA.1, xB.3~xB.1)

pars.01c <- c(pars0, pars1[c(1,2,6,7)])
names(pars.01c) <- argnames(lik.01c)

pars.02c <- c(pars0, pars2[c(1,2,6,7)])
names(pars.02c) <- argnames(lik.02c)

pars.012c <- c(pars0, pars1[c(1,2,6,7)], pars2[c(1,2,6,7)])
names(pars.012c) <- argnames(lik.012c)

ans.0 <- find.mle(lik.0, pars0, method="subplex")
# 'log Lik.' -731.3538 (df=7)
#           sA           sB          sAB           xA           xB           dA           dB 
# 1.031820e+00 8.628284e-01 4.137763e-01 3.258811e-09 1.266916e-07 9.664977e-01 9.880026e-01 

ans.01c <- find.mle(lik.01c, pars.01c, method="subplex")
# 'log Lik.' -724.5755 (df=11)
#         sA.1         sB.1        sAB.1         xA.1         xB.1         dA.1         dB.1 
# 9.211448e-01 9.400531e-01 3.953276e-01 5.623124e-02 1.648901e-08 8.729154e-01 1.116992e+00 
#         sA.2         sB.2         dA.2         dB.2 
# 1.623283e+00 3.296529e-01 1.034546e+00 8.502274e-01 

ans.02c <- find.mle(lik.02c, pars.02c, method="subplex")
# 'log Lik.' -725.1294 (df=11)
#         sA.1         sB.1        sAB.1         xA.1         xB.1         dA.1         dB.1 
# 1.079823e+00 7.797680e-01 3.999258e-01 9.085485e-09 1.556507e-10 1.232690e+00 9.373133e-01 
#         sA.2         sB.2         dA.2         dB.2 
# 8.765820e-01 1.077263e+00 2.798255e-01 1.103185e+00 

ans.012c <- find.mle(lik.012c, pars.012c, method="subplex")
# 'log Lik.' -718.4374 (df=15)
#         sA.1         sB.1        sAB.1         xA.1         xB.1         dA.1         dB.1 
# 9.264594e-01 8.472284e-01 4.122328e-01 1.498948e-07 4.467635e-10 1.296034e+00 9.997090e-01 
#         sA.2         sB.2         dA.2         dB.2 
# 1.587718e+00 3.505630e-01 1.032757e+00 8.139111e-01 
#         sA.3          sB.3         dA.3         dB.3 
# 8.760642e-01  1.076841e+00 2.807598e-01 1.103829e+00 

anova(ans.0, split01=ans.01c, split02=ans.02c, split012=ans.012c)
#          Df   lnLik    AIC  ChiSq Pr(>|Chi|)   
# full      7 -731.35 1476.7                     
# split01  11 -724.58 1471.2 13.556   0.008854 **
# split02  11 -725.13 1472.3 12.449   0.014308 * 
# split012 15 -718.44 1466.9 25.833   0.001122 **
