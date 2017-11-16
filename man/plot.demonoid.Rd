\name{plot.demonoid}
\alias{plot.demonoid}
\alias{plot.demonoid.hpc}
\title{Plot samples from the output of Laplace's Demon}
\description{
  This may be used to plot, or save plots of, samples in an object of
  class \code{demonoid} or \code{demonoid.hpc}. Plots include a trace
  plot, density plot, autocorrelation or ACF plot, and if an adaptive
  algorithm was used, the absolute difference in the proposal variance,
  or the value of epsilon, across adaptations.
}
\usage{
\method{plot}{demonoid}(x, BurnIn=0, Data, PDF=FALSE, Parms, FileName, \dots)
\method{plot}{demonoid.hpc}(x, BurnIn=0, Data, PDF=FALSE, Parms, FileName, \dots)
}
\arguments{
  \item{x}{
    This required argument is an object of class \code{demonoid} or
    \code{demonoid.hpc}.}
  \item{BurnIn}{
    This argument requires zero or a positive integer that indicates the
    number of thinned samples to discard as burn-in for the purposes of
    plotting. For more information on burn-in, see \code{\link{burnin}}.}
  \item{Data}{
    This required argument must receive the list of data that was
    supplied to \code{\link{LaplacesDemon}} to create the object of
    class \code{demonoid}.}
  \item{PDF}{
    This logical argument indicates whether or not the user wants
    Laplace's Demon to save the plots as a .pdf file.}
  \item{Parms}{
    This argument accepts a vector of quoted strings to be matched for
    selecting parameters for plotting. This argument defaults to
    \code{NULL} and selects every parameter for plotting. Each quoted
    string is matched to one or more parameter names with the
    \code{grep} function. For example, if the user specifies
    \code{Parms=c("eta", "tau")}, and if the parameter names
    are beta[1], beta[2], eta[1], eta[2], and tau, then all parameters
    will be selected, because the string \code{eta} is within
    \code{beta}. Since \code{grep} is used, string matching uses
    regular expressions, so beware of meta-characters, though these are
    acceptable: ".", "[", and "]".}
  \item{FileName}{
    This argument accepts a string and save the plot under the specified name. If \code{PDF=FALSE} this argument in unused. By default, \code{FileName = paste0("laplacesDemon-plot_", format(Sys.time(), "yyyy-mm-dd_h:m:s"), ".pdf")}
  }
  \item{\dots}{Additional arguments are unused.}
}
\details{
  The plots are arranged in a \eqn{3 \times 3}{3 x 3} matrix. Each row
  represents a parameter, the deviance, or a monitored variable. The
  left column displays trace plots, the middle column displays kernel
  density plots, and the right column displays autocorrelation (ACF)
  plots.

  Trace plots show the thinned history of the chain or Markov chain,
  with its value in the y-axis moving by thinned sample across the x-axis.
  A chain or Markov chain with good properties does not suggest a trend
  upward or downward as it progresses across the x-axis (it should
  appear stationary), and it should mix well, meaning it should appear
  as though random samples are being taken each time from the same
  target distribution. Visual inspection of a trace plot cannot
  verify convergence, but apparent non-stationarity or poor mixing can
  certainly suggest non-convergence. A red, smoothed line also appears
  to aid visual inspection.

  Kernel density plots depict the marginal posterior distribution.
  Although there is no distributional assumption about this density,
  kernel density estimation uses Gaussian basis functions.

  Autocorrelation plots show the autocorrelation or serial correlation
  between values of thinned samples at nearby thinned samples. Samples with
  autocorrelation do not violate any assumption, but are inefficient
  because they reduce the effective sample size (\code{\link{ESS}}), and
  indicate that the chain is not mixing well, since each value is
  influenced by values that are previous and nearby. The x-axis
  indicates lags with respect to thinned samples, and the y-axis
  represents autocorrelation. The ideal autocorrelation plot shows
  perfect correlation at zero lag, and quickly falls to zero
  autocorrelation for all other lags.

  If an adaptive algorithm was used, then the distribution of absolute
  differences in the proposal variances, or the value of epsilon, is
  plotted across adaptations. The proposal variance, or epsilon, should
  change less as the adaptive algorithm approaches the target
  distributions. The absolute differences in the proposal variance plot
  should approach zero. This is called the condition of diminishing
  adaptation. If it is not approaching zero, then consider using a
  different adaptive MCMC algorithm. The following quantiles are plotted
  for absolute changes proposal variance: 0.025, 0.500, and 0.975.
}
\author{Statisticat, LLC \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{burnin}},
  \code{\link{ESS}},
  \code{\link{LaplacesDemon}}, and
  \code{\link{LaplacesDemon.hpc}}.}
\examples{### See the LaplacesDemon function for an example.}
\keyword{Plot}