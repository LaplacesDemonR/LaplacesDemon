\name{Hangartner.Diagnostic}
\alias{Hangartner.Diagnostic}
\title{Hangartner's Convergence Diagnostic}
\description{
  Hangartner et al. (2011) proposed a convergence diagnostic for
  discrete Markov chains. A simple Pearson's Chi-squared test for
  two or more non-overlapping periods of a discrete Markov chain
  is a reliable diagnostic of convergence. It does not rely upon the
  estimation of spectral density, on suspect normality assumptions, or
  determining overdispersion within a small number of outcomes, all of
  which can be problematic with discrete measures. A discrete Markov
  chain is split into two or more non-overlapping windows. Two windows
  are recommended, and results may be sensitive to the number of
  selected windows, as well as sample size. As such, a user may try
  several window configurations before concluding there is no evidence
  of non-convergence.

  As the number of discrete events in the sample space increases, this
  diagnostic becomes less appropriate and standard diagnostics become
  more appropriate.
}
\usage{
Hangartner.Diagnostic(x, J=2)
}
\arguments{
  \item{x}{This required argument is a vector of marginal posterior
    samples of a discrete Markov chain, such as selected from the output
    of \code{\link{LaplacesDemon}}.}
  \item{J}{This argument specifies the number \eqn{J} of windows to be
    used, and defaults to \eqn{J=2}.}
}
\value{
  The \code{Hangartner.Diagnostic} returns an object of class
  \code{hangartner}, including the output from a Pearson's Chi-squared
  test. A frequentist p-value less than or equal to 0.05 is usually
  considered to be indicative of non-convergence.
}
\references{
  Hangartner, D., Gill, J., and Cranmer, S., (2011). "An MCMC Diagnostic
  for Purely Discrete Parameters". Paper presented at the annual meeting
  of the Southern Political Science Association, Hotel InterContinental,
  New Orleans, Louisiana Online.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{LaplacesDemon}} and
  \code{\link{TransitionMatrix}}.
}
\examples{
library(LaplacesDemon)
N <- 1000
K <- 3
x <- rcat(N, rep(1/K,K))
hd <- Hangartner.Diagnostic(x, J=2)
hd
}
\keyword{Diagnostic}
\keyword{MCMC}
