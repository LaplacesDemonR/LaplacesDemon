\name{is.stationary}
\alias{is.stationary}
\title{Logical Check of Stationarity}
\description{
  This function returns \code{TRUE} if the object is stationary
  according to the \code{\link{Geweke.Diagnostic}} function, and
  \code{FALSE} otherwise.
}
\usage{
is.stationary(x)
}
\arguments{
    \item{x}{This is a vector, matrix, or object of class
      \code{demonoid}.}
}
\details{
  Stationarity, here, refers to the limiting distribution in a Markov
  chain. A series of samples from a Markov chain, in which each sample
  is the result of an iteration of a Markov chain Monte Carlo (MCMC)
  algorithm, is analyzed for stationarity, meaning whether or not the
  samples trend or its moments change across iterations. A stationary
  posterior distribution is an equilibrium distribution, and assessing
  stationarity is an important diagnostic toward inferring Markov chain
  convergence.

  In the cases of a matrix or an object of class \code{demonoid}, all
  Markov chains (as column vectors) must be stationary for
  \code{is.stationary} to return \code{TRUE}.

  Alternative ways to assess stationarity of chains are to use the
  \code{\link{BMK.Diagnostic}} or \code{\link{Heidelberger.Diagnostic}}
  functions.
}
\value{
  \code{is.stationary} returns a logical value indicating whether or not
  the supplied object is stationary according to the
  \code{\link{Geweke.Diagnostic}} function.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{BMK.Diagnostic}},
  \code{\link{Geweke.Diagnostic}},
  \code{\link{Heidelberger.Diagnostic}}, and
  \code{\link{LaplacesDemon}}.
}
\examples{
library(LaplacesDemon)
is.stationary(rnorm(100))
is.stationary(matrix(rnorm(100),10,10))
}
\keyword{Diagnostic}
\keyword{Stationarity}
\keyword{Utility}