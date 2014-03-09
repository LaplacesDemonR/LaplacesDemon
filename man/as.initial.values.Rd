\name{as.initial.values}
\alias{as.initial.values}
\title{Initial Values}
\description{
  This function returns the most recent posterior samples from an object
  of class \code{demonoid} or \code{demonoid.hpc}, the posterior means
  of an object of class \code{iterquad}, the posterior modes of an
  object of class \code{laplace} or \code{vb}, the posterior means of an
  object of class \code{pmc} with one mixture component, or the latest
  means of the importance sampling distribution of an object of class
  \code{pmc} with multiple mixture components. The returned values are
  intended to be the initial values for future updates.
}
\usage{
as.initial.values(x)
}
\arguments{
  \item{x}{This is an object of class \code{demonoid},
    \code{demonoid.hpc}, \code{iterquad}, \code{laplace}, \code{pmc}, or
    \code{vb}.}
}
\details{
  Unless it is known beforehand how many iterations are required for
  \code{\link{IterativeQuadrature}}, \code{\link{LaplaceApproximation}},
  or \code{\link{VariationalBayes}} to converge, MCMC in
  \code{\link{LaplacesDemon}} to appear converged, or the normalized
  perplexity to stabilize in \code{\link{PMC}}, multiple updates are
  necessary. An additional update, however, should not begin with the
  same initial values as the original update, because it will have to
  repeat the work already accomplished. For this reason, the
  \code{as.initial.values} function may be used at the end of an update
  to change the previous initial values to the latest values.

  When using \code{\link{LaplacesDemon.hpc}}, \code{as.initial.values}
  should be used when the output is of class \code{demonoid.hpc}, before
  the \code{\link{Combine}} function is used to combine the multiple
  chains for use with \code{\link{Consort}} and other functions, because
  the \code{\link{Combine}} function returns an object of class
  \code{demonoid}, and the number of chains will become unknown. The
  \code{\link{Consort}} function may suggest using
  \code{as.initial.values}, but when applied to an object of class
  \code{demonoid}, it will return the latest values as if there were
  only one chain.
}
\value{
  The returned value is a vector (or matrix in the case of an object of
  class \code{demonoid.hpc}, or \code{pmc} with multiple mixture
  components) of the latest values, which may now be used as initial
  values for a future update.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{Combine}},
  \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{LaplacesDemon.hpc}},
  \code{\link{PMC}}, and
  \code{\link{VariationalBayes}}.
}
\keyword{Initial Values}