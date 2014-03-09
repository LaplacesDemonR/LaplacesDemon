\name{plotSamples}
\alias{plotSamples}
\title{Plot Samples}
\description{
  This function provides basic plots that are extended to include
  samples.
}
\usage{
plotSamples(X, Style="KDE", LB=0.025, UB=0.975, Title=NULL)
}
\arguments{
  \item{X}{This required argument is a \eqn{N \times S}{N x S} numerical
    matrix of \eqn{N} records and \eqn{S} samples.}
  \item{Style}{This argument accepts the following quoted strings:
    "barplot", "dotchart", "hist", "KDE", or "Time-Series". It defaults
    to \code{Style="KDE"}.}
  \item{LB}{This argument accepts the lower bound of a probability
    interval, which must be in the interval [0,0.5).}
  \item{UB}{This argument accepts the upper bound of a probability
    interval, which must be in the interval (0.5,1].}
  \item{Title}{This argument defaults to \code{NULL}, and otherwise
    accepts a quoted string that will be the title of the plot.}
}
\details{
  The \code{plotSamples} function extends several basic plots from
  points to samples. For example, it is common to use the \code{hist}
  function to plot a histogram from a column vector. However, the user
  may desire to plot a histogram of a column vector that was sampled
  numerous times, rather than a simple column vector, in which a
  (usually 95\%) probability interval is also plotted to show the
  uncertainty around the sampled median of each bin in the histogram.

  The \code{plotSamples} function extends the \code{barplot},
  \code{dotchart}, and \code{hist} functions to include uncertainty due
  to samples. The \code{KDE} style of plot is added so that a
  probability interval is shown around a sampled kernel density estimate
  of a distribution, and the \code{Time-Series} style of plot is added
  so that a probability interval is shown around a sampled univariate
  time-series.

  For each style of plot, three quantiles are plotted: the lower bound
  (LB), median, and upper bound (UB).

  One of many potential Bayesian applications is to examine the
  uncertainty in a predictive distribution.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\examples{
#library(LaplacesDemon)
#N <- 100
#S <- 100
#X <- matrix(rnorm(N*S),N,S)
#rownames(X) <- 1:100
#plotSamples(X, Style="barplot", LB=0.025, UB=0.975)
#plotSamples(X[1:10,], Style="dotchart", LB=0.025, UB=0.975)
#plotSamples(X, Style="hist", LB=0.025, UB=0.975)
#plotSamples(X, Style="KDE", LB=0.025, UB=0.975)
#plotSamples(X, Style="Time-Series", LB=0.025, UB=0.975)
}
\keyword{Plot}