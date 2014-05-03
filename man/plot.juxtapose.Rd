\name{plot.juxtapose}
\alias{plot.juxtapose}
\title{Plot MCMC Juxtaposition}
\description{
  This may be used to plot a juxtaposition of MCMC algorithms according
  either to \code{\link{IAT}} or ISM (Independent Samples per Minute).
}
\usage{\method{plot}{juxtapose}(x, Style="ISM", \dots)}
\arguments{
  \item{x}{This required argument is an object of class
    \code{juxtapose}.}
  \item{Style}{This argument accepts either \code{IAT} or \code{ISM},
    and defaults to \code{ISM}.}
  \item{\dots}{Additional arguments are unused.}
}
\details{
  When \code{Style="IAT"}, the medians and 95\% probability intervals of
  the integrated autocorrelation times (IATs) of MCMC algorithms are
  displayed in a caterpillar plot. The best, or least inefficient, MCMC
  algorithm is the algorithm with the lowest IAT.

  When \code{Style="ISM"}, the medians and 95\% probability intervals of
  the numbers of independent samples per minute (ISM) of MCMC algorithms
  are displayed in a caterpillar plot. The best, or least inefficient,
  MCMC algorithm is the algorithm with the highest ISM.

  For more information, see the \code{\link{Juxtapose}} function.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{Juxtapose}}}
\keyword{Plot}