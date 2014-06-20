\name{IAT}
\alias{IAT}
\title{Integrated Autocorrelation Time}
\description{
  The \code{IAT} function estimates integrated autocorrelation time,
  which is the computational inefficiency of a continuous chain or MCMC
  sampler. IAT is also called the IACT, ACT, autocorrelation time,
  autocovariance time, correlation time, or inefficiency factor. A lower
  value of \code{IAT} is better. \code{IAT} is a MCMC diagnostic that is
  an estimate of the number of iterations, on average, for an
  independent sample to be drawn, given a continuous chain or Markov
  chain. Put another way, \code{IAT} is the number of correlated samples
  with the same variance-reducing power as one independent sample.

  IAT is a univariate function. A multivariate form is not included.
}
\usage{
IAT(x)
}
\arguments{
  \item{x}{This requried argument is a vector of samples from a chain.}
}
\details{
  \code{IAT} is a MCMC diagnostic that is often used to compare
  continuous chains of MCMC samplers for computational inefficiency,
  where the sampler with the lowest \code{IAT}s is the most efficient
  sampler. Otherwise, chains may be compared within a model, such as
  with the output of \code{\link{LaplacesDemon}} to learn about the
  inefficiency of the continuous chain. For more information on
  comparing MCMC algorithmic inefficiency, see the
  \code{\link{Juxtapose}} function.
  
  \code{IAT} is also estimated in the \code{\link{PosteriorChecks}}
  function. \code{IAT} is usually applied to a stationary, continuous
  chain after discarding burn-in iterations (see \code{\link{burnin}}
  for more information). The \code{IAT} of a continuous chain correlates
  with the variability of the mean of the chain, and relates to
  Effective Sample Size (\code{\link{ESS}}) and Monte Carlo Standard
  Error (\code{\link{MCSE}}).

  \code{IAT} and \code{\link{ESS}} are inversely related, though not
  perfectly, because each is estimated a little differently. Given
  \eqn{N}{N} samples and taking autocorrelation into account,
  \code{\link{ESS}} estimates a reduced number of \eqn{M}{M} samples.
  Conversely, \code{IAT} estimates the number of autocorrelated samples,
  on average, required to produce one independently drawn sample.
  
  The \code{IAT} function is similar to the \code{IAT} function in the
  \code{Rtwalk} package of Christen and Fox (2010), which is currently
  unavailabe on CRAN.
}
\value{
  The \code{IAT} function returns the integrated autocorrelation time of
  a chain.
}
\references{
  Christen, J.A. and Fox, C. (2010). "A General Purpose Sampling
  Algorithm for Continuous Distributions (the t-walk)". \emph{Bayesian
  Analysis}, 5(2), p. 263--282.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{burnin}},
  \code{\link{Compare}},
  \code{\link{ESS}},
  \code{\link{LaplacesDemon}},
  \code{\link{MCSE}}, and
  \code{\link{PosteriorChecks}}.
}
\examples{
library(LaplacesDemon)
theta <- rnorm(100)
IAT(theta)
}
\keyword{Diagnostic}
\keyword{MCMC}
