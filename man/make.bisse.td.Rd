\name{make.bisse.td}
\alias{make.bisse.td}

\title{Binary State Speciation and Extinction Model: Time Chunks}

\description{Create a likelihood function for a BiSSE model where
  different chunks of time have different parameters.  This code is
  experimental!}

\usage{
make.bisse.td(tree, states, n.epoch, unresolved=NULL,
             sampling.f=NULL, nt.extra=10, safe=FALSE)
}

\arguments{
  \item{tree}{A phylogenetic tree, in \code{ape} \dQuote{phylo} format.}
  \item{states}{A vector of character states, each of which must be 0 or
    1, or \code{NA} if the state is unknown.  This vector must have
    names that correspond to the tip labels in the phylogenetic tree
    (\code{tree$tip.label}).  For tips
    corresponding to unresolved clades, the state should be \code{NA}.}
  \item{n.epoch}{Number of epochs.  1 corresponds to plain BiSSE, so
    this will generally be an integer at least 2.}
  \item{unresolved}{Unresolved clade information: see
    \code{\link{make.bisse}}. (Currently this is not supported.)}
  \item{sampling.f}{Vector of length 2 with the estimated proportion of
    extant species in state 0 and 1 that are included in the phylogeny.
    See \code{\link{make.bisse}}.}
  \item{nt.extra}{The number of species modelled in unresolved clades
    (this is in addition to the largest observed clade).}
  \item{safe}{Use a safer, slower, integrator?  Leave this alone unless
    crashes are happening.}
}

\section{Warning}{
  This will change soon!  In particular, the likelihood function
  requires that the "time" of the final parameter set is \code{Inf};
  all optimisations and parameter searches will have to be on
  constrained models that remove this so that the ML/MCMC functions do
  not change this parameter!
}


\examples{
set.seed(4)
pars <- c(0.1, 0.2, 0.03, 0.03, 0.01, 0.01)
phy <- tree.bisse(pars, max.t=30, x0=0)

## Show that the calculations are correct for 
lik <- make.bisse(phy, phy$tip.state)
lik.td <- make.bisse.td(phy, phy$tip.state, 2)

argnames(lik.td)
pars.td <- c(2,   pars,
             Inf, pars)
names(pars.td) <- argnames(lik.td)

## Calculations are identical to a reasonable tolerance.
lik(pars) - lik.td(pars.td)

## Constrain the times
lik.td.use <- constrain(lik.td, t.t1 ~ 2, t.t2 ~ Inf)

## This is marked "don't run" because it takes a minute or too.
\dontrun{
fit <- find.mle(lik, pars)
fit.td <- find.mle(lik.td.use, pars.td[argnames(lik.td.use)],
                   control=list(maxit=20000, reltol=0.0001))
}
}