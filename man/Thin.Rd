\name{Thin}
\alias{Thin}
\title{Thin}
\description{
  This function reduces the number of posterior samples by retaining
  every \eqn{k}th sample.
}
\usage{
Thin(x, By=1)
}
\arguments{
  \item{x}{This is a vector or matrix of posterior samples to be
    thinned.}
  \item{By}{This argument specifies that every \eqn{k}th posterior
    sample will be retained, and \code{By} defaults to 1, meaning that
    thinning will not occur, because every sample will be retained.}
}
\details{
  A thinned matrix of posterior samples is a matrix in which only every
  \eqn{k}th posterior sample (or row) in the original matrix is
  retained. The act of thinning posterior samples has been criticized as
  throwing away information, which is correct. However, it is common
  practice to thin posterior samples, usually associated with MCMC such
  as \code{\link{LaplacesDemon}}, for two reasons. First, Each chain
  (column vector) in a matrix of posterior samples probably has higher
  autocorrelation than desired, which reduces the effective sample size
  (see \code{\link{ESS}} for more information). Therefore, a thinned
  matrix usually contains posterior samples that are closer to
  independent than an un-thinned matrix. The other reason for the
  popularity of thinning is that it a user may not have the
  random-access memory (RAM) to store large, un-thinned matrices of
  posterior samples.

  \code{\link{LaplacesDemon}} and \code{\link{PMC}} automatically thin
  posterior samples, deviance samples, and samples of monitored
  variables, according to its own user-specified argument. The
  \code{Thin} function is made available here, should it be necessary to
  thin posterior samples outside of objects of class \code{demonoid} or
  \code{pmc}.
}
\value{
  The \code{Thin} argument returns a thinned matrix. When \code{x} is a
  vector, the returned object is a matrix with 1 column.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{ESS}},
  \code{\link{LaplacesDemon}}, and
  \code{\link{PMC}}.
}
\examples{
library(LaplacesDemon)
x <- matrix(runif(100), 10, 10)
Thin(x, By=2)
}
\keyword{Diagnostic}
\keyword{MCMC}
\keyword{Monte Carlo}
\keyword{Utility}
