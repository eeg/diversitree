\name{constrain.par}
\alias{constrain.par}
\alias{fix.par}
\title{Constrain Parameters of a Model}

\description{Constrain a function so that some arguments are always
  set.  If \code{f} is a function that takes a vector \code{x} as its
  first argument, this function returns a new function that takes a
  shorter vector \code{x} with some elements constrained to the value of
  other elements of \code{x} (\code{constrain.par}) or fixed to
  particular values (\code{fix.par}).
}

\usage{
constrain.par(f, rel)
fix.par(f, rel)
}
\arguments{
  \item{f}{A function to constrain.}
  \item{rel}{An integer vector to constrain parameters (see Details).}
}

\details{
  The argument \code{rel} should be a vector of integers.  Arguments
  that are not constrained should be \code{NA}, and arguments that are
  should be an integer referring to the unconstrained argument \emph{in
  the original list} that they are constrained to.
}

\section{Warning}{
  No checks are done that the resulting function makes any sense; some
  sanity checking needs to be added to the \code{rel} argument to make
  sure that \code{pars[rel]} returns something useful.  In particular,
  constraining an already constrained function is fraught with danger.
}

\keyword{programming}
