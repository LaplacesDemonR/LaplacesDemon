\name{as.covar}
\alias{as.covar}
\title{Proposal Covariance}
\description{
  This function returns the most recent covariance matrix or a list of
  blocking covariance matrices from an object of class \code{demonoid},
  the most recent covariance matrix from \code{iterquad},
  \code{laplace}, or \code{vb}, the most recent covariance matrix from
  the chain with the lowest deviance in an object of class
  \code{demonoid.hpc}, and a number of covariance matrices of an object
  of class \code{pmc} equal to the number of mixture components. The
  returned covariance matrix or matrices are intended to be the initial
  proposal covariance matrix or matrices for future updates. A variance
  vector from an object of class \code{demonoid} or \code{demonoid.hpc}
  is converted to a covariance matrix.
}
\usage{
as.covar(x)
}
\arguments{
  \item{x}{This is an object of class \code{demonoid},
    \code{demonoid.hpc}, \code{iterquad}, \code{laplace}, \code{pmc}, or
    \code{vb}.}
}
\details{
  Unless it is known beforehand how many iterations are required for
  iterative quadrature, Laplace Approximation, or Variational Bayes to
  converge, MCMC to appear converged, or the normalized perplexity to
  stabilize in PMC, multiple updates are necessary. An additional
  update, however, should not begin with the same proposal covariance
  matrix or matrices as the original update, because it will have to
  repeat the work already accomplished. For this reason, the
  \code{as.covar} function may be used at the end of an update to change
  the previous initial values to the latest values.

  The \code{as.covar} function is most helpful with objects of class
  \code{pmc} that have multiple mixture components. For more
  information, see \code{\link{PMC}}.
}
\value{
  The returned value is a matrix (or array in the case of PMC with
  multiple mixture components) of the latest observed or proposal
  covariance, which may now be used as an initial proposal covariance
  matrix or matrices for a future update.
}
\author{Statisticat, LLC \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{LaplacesDemon.hpc}},
  \code{\link{PMC}}, and
  \code{\link{VariationalBayes}}.
}
\keyword{Utility}