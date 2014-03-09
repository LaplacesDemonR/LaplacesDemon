\name{plot.miss}
\alias{plot.miss}
\title{Plot samples from the output of MISS}
\description{
  This may be used to plot, or save plots of, samples in an object of
  class \code{miss}. Plots include a trace plot, density plot, and
  autocorrelation or ACF plot.
}
\usage{
\method{plot}{miss}(x, PDF=FALSE, \dots)
}
\arguments{
  \item{x}{
    This required argument is an object of class \code{miss}.}
  \item{PDF}{
    This logical argument indicates whether or not the user wants
    Laplace's Demon to save the plots as a .pdf file.}
  \item{\dots}{Additional arguments are unused.}
}
\details{
  The plots are arranged in a \eqn{3 \times 3}{3 x 3} matrix. Each row
  represents the predictive distribution of a missing value. The
  left column displays trace plots, the middle column displays kernel
  density plots, and the right column displays autocorrelation (ACF)
  plots.

  Trace plots show the thinned history of the predictive distribution,
  with its value in the y-axis moving by iteration across the x-axis.
  Simulations of a predictive distribution with good properties do not
  suggest a trend upward or downward as it progresses across the x-axis
  (it should appear stationary), and it should mix well, meaning it
  should appear as though random samples are being taken each time from
  the same target distribution. Visual inspection of a trace plot cannot
  verify convergence, but apparent non-stationarity or poor mixing can
  certainly suggest non-convergence. A red, smoothed line also appears
  to aid visual inspection.

  Kernel density plots depict the marginal posterior distribution.
  There is no distributional assumption about this density.

  Autocorrelation plots show the autocorrelation or serial correlation
  between sampled values at nearby iterations. Samples with
  autocorrelation do not violate any assumption, but are inefficient
  because they reduce the effective sample size (\code{\link{ESS}}), and
  indicate that the chain is not mixing well, since each value is
  influenced by values that are previous and nearby. The x-axis
  indicates lags with respect to samples by iteration, and the y-axis
  represents autocorrelation. The ideal autocorrelation plot shows
  perfect correlation at zero lag, and quickly falls to zero
  autocorrelation for all other lags.
}
\author{Statisticat, LLC \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{MISS}}.}
\examples{### See the MISS function for an example.}
\keyword{Plot}