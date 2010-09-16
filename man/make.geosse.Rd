\name{make.geosse}
\alias{make.geosse}
\alias{starting.point.geosse}

\title{Geographic State Speciation and Extinction Model}

\description{Prepare to run GeoSSE (Geographic State Speciation and
  Extinction) on a phylogenetic tree and character distribution.  This
  function creates a likelihood function that can be used in
  \link[=find.mle]{maximum likelihood} or \link[=mcmc]{Bayesian}
  inference.}

\usage{
make.geosse(tree, states, sampling.f=NULL, safe=FALSE, strict=TRUE)
starting.point.geosse(tree, q.div=5, yule=FALSE)
}

\arguments{
  \item{tree}{A phylogenetic tree, in \code{ape} \dQuote{phylo} format.}
  \item{states}{A vector of character states, each of which must be 0
    (in both regions/widespread), 1 or 2 (endemic to one region), or
    \code{NA} if the state is unknown.  This vector must have names that
    correspond to the tip labels in the phylogenetic tree
    (\code{tree$tip.label}).}
%  \item{unresolved}{Unresolved clades are not yet implemented for geosse.}
  \item{sampling.f}{Vector of length 3 with the estimated proportion of
    extant species in states 0, 1 and 2 that are included in the
    phylogeny.  A value of \code{c(0.5, 0.75, 1)} means that half of
    species in state 0, three quarters of species in state 1, and all the
    species in state 2 are included in the phylogeny.  By default all
    species are assumed to be known.}
  \item{safe}{Use a safer, slower, integrator?  Leave this alone unless
    crashes are happening.}
  \item{strict}{The \code{states} vector is always checked to make sure
    that the values are 0, 1 and 2 only.  If \code{strict} is \code{TRUE}
    (the default), then the additional check is made that \emph{every}
    state is present.  The likelihood models tend to be poorly behaved
    where states are missing.}
  \item{q.div}{Ratio of diversification rate to state change rate.}
  \item{yule}{Logical: should starting parameters be Yule estimates
    rather than birth-death estimates?}
}

\details{
  \code{make.geosse} returns a function of class \code{geosse}.  This
  function has argument list (and default values)

  \preformatted{
    f(pars, prior=NULL, condition.surv=FALSE, root=ROOT.OBS,
    root.p=NULL, intermediates=FALSE)
  }
  
  The arguments are interpreted as
  \itemize{
    \item \code{pars} A vector of seven parameters, in the order
      \code{sA}, \code{sB}, \code{sAB}, \code{xA}, \code{xB},
      \code{dA}, \code{dB}.
    \item \code{condition.surv} (logical): should the likelihood
      calculation condition on survival of two lineages and the
      speciation event subtending them?  This is currently not done.
      %This is done by default, following Nee et al. 1994.
    \item \code{root}: Behaviour at the root (see FitzJohn et al. 2009).
      The possible options are
      \itemize{
        \item \code{ROOT.FLAT}: A flat prior, weighting
        \eqn{D_0}{D0}, \eqn{D_1}{D1} and \eqn{D_2}{D2} equally.
        \item \code{ROOT.EQUI}: Use the equilibrium distribution
        of the model.
        \item \code{ROOT.OBS}: Weight \eqn{D_0}{D0}, \eqn{D_1}{D1} and
        \eqn{D_2}{D2} by their relative probability of observing the
        data, following FitzJohn et al. 2009.
        \item \code{ROOT.GIVEN}: Root will be in state i
        with probability \code{root.p[i+1]}.
        \item \code{ROOT.BOTH}: Don't do anything at the root,
        and return both values.  (Note that this will not give you a
        likelihood!)
    }
    \item \code{root.p}: Root weightings for use when
      \code{root=ROOT.GIVEN}.  Must be a vector of length 3 with sum 1.
    \item \code{intermediates}: Add intermediates to the returned value as
      attributes:
      \itemize{
        \item \code{cache}: Cached tree traversal information.
        \item \code{intermediates}: Mostly branch end information.
        \item \code{vals}: Root \eqn{D} values.
    }
    At this point, you will have to poke about in the source for more
    information on these.
  }
  
  \code{starting.point.geosse} produces a not-that-great starting point
  based on \code{starting.point.bisse}.
}

\section{Unresolved clade information}{
  This is not (yet?) available for GeoSSE.
}


\seealso{
  \code{\link{constrain}} for making submodels, \code{\link{find.mle}}
  for ML parameter estimation, \code{\link{mcmc}} for MCMC integration.

  The help pages for \code{\link{find.mle}} has further examples of ML
  searches on full and constrained BiSSE models.  Things work similarly
  for GeoSSE, just with different parameters.
}

% FIXME
\examples{
pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
set.seed(4)
phy <- tree.bisse(pars, max.t=30, x0=0)

## Here is the 52 species tree with the true character history coded.
## Red is state '1', which has twice the speciation rate of black (state
## '0').
h <- history.from.sim.discrete(phy, 0:1)
plot(h, phy)

lik <- make.bisse(phy, phy$tip.state)
lik(pars) # -159.71

## Heuristic guess at a starting point, based on the constant-rate
## birth-death model (uses \link{make.bd}).
p <- starting.point.bisse(phy)

%% This is code for better getting at the q's
%% system.time(find.mle(constrain(make.mk2(phy, phy$tip.state), q10 ~
%% q01), .1)$par)

## Start an ML search from this point.  This takes some time (~7s)
fit <- find.mle(lik, p, method="subplex")
logLik(fit) # -158.69

## The estimated parameters aren't too far away from the real ones, even
## with such a small tree
rbind(real=pars,
      estimated=round(coef(fit), 2))

## Test a constrained model where the speciation rates are set equal
## (takes ~4s).
lik.l <- constrain(lik, lambda1 ~ lambda0)
fit.l <- find.mle(lik.l, p[-1], method="subplex")
logLik(fit.l) # -158.74

## Despite the difference in the estimated parameters, there is no
## statistical support for this difference:
anova(fit, equal.lambda=fit.l)

## Run an MCMC.  Because we are fitting six parameters to a tree with
## only 50 species, priors will be needed.  I will use an exponential
## prior with rate 1/(2r), where r is the character independent
## diversificiation rate:
prior <- make.prior.exponential(1 / (2 * (p[1] - p[3])))

## This takes quite a while to run, so is not run by default
\dontrun{
tmp <- mcmc(lik, fit$par, nsteps=100, prior=prior,
            lower=rep(0, 6), upper=rep(Inf, 6), w=rep(.1, 6))
w <- diff(sapply(tmp[2:7], range))
samples <- mcmc(lik, fit$par, nsteps=1000, prior=prior,
                lower=rep(0, 6), upper=rep(Inf, 6), w=rep(.1, 6))

## See \link{profiles.plot} for more information on plotting these
## profiles.
col <- c("blue", "red")
profiles.plot(samples[c("lambda0", "lambda1")], col.line=col, las=1,
              xlab="Speciation rate",
              ylab="Posterior probability density")
legend("topright", c("lambda0", "lambda1"), col=col, lty=1)
}

## BiSSE reduces to the birth-death model and Mk2 when diversification
## is state independent (i.e., lambda0 ~ lambda1 and mu0 ~ mu1).
lik.mk2 <- make.mk2(phy, phy$tip.state)
lik.bd <- make.bd(phy)

## 1. BiSSE / Birth-Death
## Set the q01 and q10 parameters to arbitrary numbers (need not be
## symmetric), and constrain the lambdas and mus to be the same for each
## state.  The likelihood function now has just two parameters and
## will be proprtional to Nee's birth-death based likelihood:
lik.bisse.bd <- constrain(lik,
                          lambda1 ~ lambda0, mu1 ~ mu0,
                          q01 ~ .01, q10 ~ .02)
pars <- c(.1, .03)
## These differ by -167.3861 for both parameter sets:
lik.bisse.bd(pars)   - lik.bd(pars)
lik.bisse.bd(2*pars) - lik.bd(2*pars)

## 2. BiSSE / Mk2
## Same idea as above: set all diversification parameters to arbitrary
## values (but symmetric this time):
lik.bisse.mk2 <- constrain(lik,
                           lambda0 ~ .1, lambda1 ~ .1,
                           mu0 ~ .03, mu1 ~ .03)
## Differ by -150.4740 for both parameter sets.
lik.bisse.mk2(pars)   - lik.mk2(pars)
lik.bisse.mk2(2*pars) - lik.mk2(2*pars)

## 3. Sampled BiSSE / Birth-Death
## Pretend that the tree is only .6 sampled:
lik.bd2 <- make.bd(phy, sampling.f=.6)
lik.bisse2 <- make.bisse(phy, phy$tip.state, sampling.f=c(.6, .6))
lik.bisse2.bd <- constrain(lik.bisse2,
                           lambda1 ~ lambda0, mu1 ~ mu0,
                           q01 ~ .01, q10 ~ .01)

## Difference of -167.2876
lik.bisse2.bd(pars)   - lik.bd2(pars)
lik.bisse2.bd(2*pars) - lik.bd2(2*pars)
}

\references{
  FitzJohn R.G., Maddison W.P., and Otto S.P. 2009. Estimating
  trait-dependent speciation and extinction rates from incompletely
  resolved phylogenies. Syst. Biol. 58:595-611.
  
  Goldberg E.E., Lancaster L.T., and Ree R.H. Phylogenetic inference of
  reciprocal effects between geographic range evolution and
  diversification. in review at Syst Biol.

  Maddison W.P., Midford P.E., and Otto S.P. 2007. Estimating
  a binary character's effect on speciation and
  extinction. Syst. Biol. 56:701-710.

  Nee S., May R.M., and Harvey P.H. 1994. The reconstructed
  evolutionary process. Philos.  Trans. R. Soc. Lond. B
  Biol. Sci. 344:305-311.
}

\author{Emma E. Goldberg}
\keyword{models}

