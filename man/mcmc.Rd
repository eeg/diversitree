\name{mcmc}
\alias{mcmc}
\title{Simple Markov Chain Monte Carlo with Slice Sampling}

\description{Run a simple-minded MCMC using slice samples for
  independent updating of each variable.}

\usage{
mcmc(f, x0, y0, nsteps, w, lower, upper, type=TYPE.LOG, ...)
}

\arguments{
  \item{f}{Function to run MCMC on.}
  \item{x0}{Initial x location (vector).}
  \item{y0}{Function value at \code{x0}.}
  \item{nsteps}{Number of MCMC steps to take.}
  \item{w}{(Approximate) width of the high probability region; this is
    the step size.  Must be same length as \code{x0}.}
  \item{lower}{Lower bounds on parameter space (vector).}
  \item{upper}{Upper bounds on parameter space (vector).}
  \item{type}{Type of value that the function \code{f} returns.  The
    default is to assume it returns a log-value, \code{bisse:::TYPE.RAW}
    assumes a raw probability and \code{bisse:::TYPE.NLOG} assumes a
    negative log-probability.}
  \item{...}{Arguments passed to the function \code{f}}
}

\keyword{models}
