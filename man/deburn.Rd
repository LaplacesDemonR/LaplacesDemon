\name{deburn}
\alias{deburn}
\title{De-Burn}
\description{
  The \code{deburn} function discards or removes a user-specified number
  of burn-in iterations from an object of class \code{demonoid}.
}
\usage{
deburn(x, BurnIn=0)
}
\arguments{
  \item{x}{This is an object of class \code{demonoid}.}
  \item{BurnIn}{This argument defaults to \code{BurnIn=0}, and accepts
    an integer that indicates the number of iterations to discard as
    burn-in.}
}
\details{
  Documentation for the \code{\link{burnin}} function provides an
  introduction to the concept of burn-in as it relates to Markov chains.

  The \code{deburn} function discards a number of the first posterior
  samples, as specified by the \code{BurnIn} argument. Stationarity is
  not checked, because it is assumed the user has a reason for using the
  \code{deburn} function, rather than using the results from the object
  of class \code{demonoid}. Therefore, the posterior samples in
  \code{Posterior1} and \code{Posterior2} are identical, as are
  \code{Summary1} and \code{Summary2}.
}
\value{
  The \code{deburn} function returns an object of class \code{demonoid}.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{burnin}} and
  \code{\link{LaplacesDemon}}.
}
\examples{
### Assuming the user has Fit which is an object of class demonoid:
#library(LaplacesDemon)
#Fit2 <- deburn(Fit, BurnIn=100)
}
\keyword{Stationarity}
\keyword{Utility}
