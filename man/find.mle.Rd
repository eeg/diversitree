\name{find.mle}
\alias{find.mle}
\alias{find.mle.bisse}
\alias{find.mle.default}
\alias{starting.point}

\title{Maximimum Likelihood Inference}

\description{Find the maximum likelihood point of a model by nonlinear
  optimisation.}

\usage{
find.mle(func, x.init, ...)
\method{find.mle}{default}(func, x.init, ...)
\method{find.mle}{bisse}(func, x.init, control=list(),
                         lower=NULL, fail.value=NULL, ...)
starting.point(tree, q.div=5)
}

\arguments{
  \item{func}{A likelihood function.  This is assumed to return the log
    likelihood (see Details).  The function must take a vector of
    parameters as the first argument.}
  \item{x.init}{Initial starting point for the optimisation.}
  \item{control}{List of options to \code{\link{optim}}}
  \item{lower}{Lower limits for the parameters.  If omitted, the
    \code{bisse} method adds sensible lower limits in the absence of
    additional information (see Details).}
  \item{fail.value}{Value to use if the likelihood calculation fails.
    \code{optim} requires finite values for all evaluations, so if you
    have a problem with failed evaluations set this
    to be an unlikely value (but not too much smaller than observed
    numbers).}
  \item{...}{Arguments passed to the function \code{f}}
  \item{tree}{A phylogeny}
  \item{q.div}{Ratio of diversification rate to character change rate}
}

\details{
  In general, it will not matter if this not return the natural log
  likelihood (raw likelihood or a different log base would work fine),
  so long as the likelihood is scaled so that larger values are more
  likely.  This is in contrast with the \code{mle} function in the
  package \code{stats4}.

  Arguments for constrained BiSSE models are processed a little
  differently.  Constrained BiSSE models have fewer than six parameters;
  if a constrained model is detected and \code{x.init} has six
  arguments, then the arguments \code{lower} and \code{x.init} are
  subset so that only the appropriate elements are used.  This does not
  work for other vectorised arguments to \code{\link{optim}} (e.g.,
  \code{upper}).

  Lower limits are required for optimisation to work well.  Being rates,
  all parameters have a natural lower limit of zero, but the likelihood
  will be undefined with polymorphic state data if both \code{q} values
  are zero, and undefined if both \code{lambda} values are zero.  The
  lower limit is set to be \code{1e-4} for both \code{lambda} and
  \code{q} values.  A warning will be given if the final point lies on
  any given constraint; it would be worth rerunning the analysis with
  contraints relaxed if this warning is given.

  \code{starting.point} produces a poor starting point to start from,
  based on the character-independent birth-death model.
}

\keyword{models}

\examples{
pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
set.seed(2)
phy <- tree.bisse(pars, max.t=60, x0=0)

plot(phy)

f <- make.bisse(phy, phy$tip.state)
f(pars)

## This will take about 20s.
\dontrun{
fit.mle <- find.mle(f, pars)
fit.mle
}
}