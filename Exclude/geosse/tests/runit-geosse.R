tol <- 1e-6

#--------------------------------------------------
# test data for everything below
#--------------------------------------------------

tree <- read.tree(text="((((0:0.461876,(1:0.307717,(2:0.231825,3:0.231825):0.075892):0.154159):0.425922,((4:0.244819,5:0.244819):0.004749,6:0.249568):0.638231):0.142051,7:1.029850):0.038423,(((8:0.510933,(9:0.427929,(10:0.119778,11:0.119778):0.308151):0.083004):0.007428,(12:0.488316,13:0.488316):0.030044):0.100160,14:0.618521):0.449752);")
states <- c(0, 1, 0, 2, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0)
names(states) <- as.character(seq(Ntip(tree))-1)

lnL.full <- make.geosse(tree, states)

rate.names <- c("sA", "sB", "sAB", "xA", "xB", "dA", "dB")
pars0 <- c(0.9, 0.8, 0.1, 0.2, 0.3, 0.5, 0.6)

#--------------------------------------------------
# test functions
#--------------------------------------------------

test.starting <- function()
{
    ans <- starting.point.geosse(tree)
    ans0 <- c( rep(3.339799, 3), rep(1.669899, 4) )
    names(ans0) <- rate.names
    checkEquals(ans, ans0, tolerance=tol)

    ans <- starting.point.geosse(tree, yule=TRUE)
    ans0 <- c( rep(1.8861312, 3), rep(0, 2), rep(0.1886131, 2) )
    names(ans0) <- rate.names
    checkEquals(ans, ans0, tolerance=tol)
}

test.lnL <- function()
{
    pars <- pars0
    checkEquals(lnL.full(pars), -23.71574, tolerance=tol)
    checkEquals(lnL.full(pars, condition.surv=T), -24.10259, tolerance=tol)
    checkEquals(lnL.full(pars, root=ROOT.EQUI), -24.32152, tolerance=tol)
    checkEquals(lnL.full(pars, root=ROOT.GIVEN, root.p=c(0.6,0.4,0.2)), 
                -23.78353, tolerance=tol)
    checkEquals(lnL.full(pars, root=ROOT.GIVEN, root.p=rep(1,3)/3), 
                lnL.full(pars, root=ROOT.FLAT), tolerance=tol)

    lnL <- constrain(lnL.full, sAB ~ 0)
    checkEquals(lnL(pars[-3]), -23.56325, tolerance=tol)

    lnL <- constrain(lnL.full, dB ~ dA, sAB ~ 0)
    checkEquals(lnL(pars[-c(3,7)]), -23.763, tolerance=tol)
}

test.mle <- function()
{
    pars <- pars0
    names(pars) <- argnames(lnL.full)
    ans <- find.mle(lnL.full, pars)
    ans0.par <- c(1.592065e+00, 4.066657e-01, 7.598541e-06, 6.827111e-07, 
                  1.675306e-04, 1.262941e+00, 1.252454e+00)
    names(ans0.par) <- rate.names
    checkEquals(ans$par, ans0.par, tolerance=tol)
    checkEquals(ans$lnLik, -18.77676, tolerance=tol)

    lnL <- constrain(lnL.full, dB ~ dA, sAB ~ 0)
    pars <- pars0[-c(3,7)]
    names(pars) <- argnames(lnL)
    ans <- find.mle(lnL, pars)
    ans0.par <- c(1.591802e+00, 4.074054e-01, 2.284247e-06, 2.727531e-05, 
                  1.262041e+00)
    names(ans0.par) <- rate.names[-c(3,7)]
    checkEquals(ans$par, ans0.par, tolerance=tol)
    checkEquals(ans$lnLik, -18.77673, tolerance=tol)
}

test.mcmc <- function()
{
    set.seed(1)
    pars <- pars0
    ans <- mcmc(lnL.full, pars, nsteps=3, lower=0, upper=3, w=3/5, 
                prior=make.prior.exponential(seq(0.8, by=0.1, length.out=7)),
                print.every=0) 
    print(ans)
    checkEquals(ans$p[3], -29.78957, tolerance=tol)

    set.seed(1)
    lnL <- constrain(lnL.full, dB ~ dA, sAB ~ 0)
    pars <- pars0[-c(3,7)]
    ans <- mcmc(lnL, pars, nsteps=3, lower=0, upper=3, w=seq(0.8, by=0.1,
                length.out=5), prior=make.prior.exponential(1), print.every=0)
    print(ans)
    checkEquals(ans$p[3], -24.58257, tolerance=tol)
}
