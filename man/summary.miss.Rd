\name{summary.miss}
\alias{summary.miss}
\title{MISS Summary}
\description{
  This function summarizes posterior predictive distributions from
  an object of class \code{miss}.
}
\usage{\method{summary}{miss}(object, \dots)}
\arguments{
  \item{object}{An object of class \code{miss} is required.}
  \item{\dots}{Additional arguments are unused.}
}
\details{
  This function summarizes the posterior predictive distributions from
  an object of class \code{miss}.
}
\value{
  This function returns a \eqn{M \times 7}{M x 7} matrix, in which each
  row is the posterior predictive distribution of one of \eqn{M} missing
  values. Columns are Mean, SD, MCSE, ESS, LB, Median, and UB.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{MISS}}.
}
\examples{### See the MISS function for an example.}
\keyword{Imputation}
\keyword{summary}
