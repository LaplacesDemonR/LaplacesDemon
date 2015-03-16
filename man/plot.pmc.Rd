\name{plot.pmc}
\alias{plot.pmc}
\title{Plot samples from the output of PMC}
\description{
  This may be used to plot, or save plots of, samples in an object of
  class \code{pmc}. Plots include a trace plot and density plot for
  parameters, a density plot for deviance and monitored variables, and
  convergence plots.
}
\usage{
\method{plot}{pmc}(x, BurnIn=0, Data, PDF=FALSE, Parms, \dots)
}
\arguments{
  \item{x}{This required argument is an object of class \code{pmc}.}
  \item{BurnIn}{This argument requires zero or a positive integer that
    indicates the number of iterations to discard as burn-in for the
    purposes of plotting.}
  \item{Data}{This required argument must receive the list of data that
    was supplied to \code{\link{PMC}} to create the object of class
    \code{pmc}.}
  \item{PDF}{This logical argument indicates whether or not the user
    wants Laplace's Demon to save the plots as a .pdf file.}
  \item{Parms}{This argument accepts a vector of quoted strings to be
    matched for selecting parameters for plotting. This argument
    defaults to \code{NULL} and selects every parameter for
    plotting. Each quoted string is matched to one or more parameter
    names with the \code{grep} function. For example, if the user
    specifies \code{Parms=c("eta", "tau")}, and if the parameter names
    are beta[1], beta[2], eta[1], eta[2], and tau, then all parameters
    will be selected, because the string \code{eta} is within
    \code{beta}. Since \code{grep} is used, string matching uses
    regular expressions, so beware of meta-characters, though these are
    acceptable: ".", "[", and "]".}
  \item{\dots}{Additional arguments are unused.}
}
\details{
  The plots are arranged in a \eqn{2 \times 2}{2 x 2} matrix. Each row
  represents a parameter, the deviance, or a monitored variable. For
  parameters, the left column displays trace plots and the right column
  displays kernel density plots.

  Trace plots show the history of the distribution of independent
  importance samples. When multiple mixture components are used, each
  mixture component has a different color. These plots are unavailable
  for the deviance and monitored variables.
  
  Kernel density plots depict the marginal posterior distribution.
  Although there is no distributional assumption about this density,
  kernel density estimation uses Gaussian basis functions.

  Following these plots are three plots for convergence. First, ESSN
  (red) and perplexity (black) are plotted by iteration. Convergence
  occurs when both of these seem to stabilize, and higher is
  better. The second plot shows the distribution of the normalized
  importance weights by iteration. The third plot appears only when
  multiple mixture components are used. The third plot displays the
  probabilities of each mixture component by iteration. Although the
  last two plots are not formally convergence plots, they are provided
  so the user can verify the distribution of importance weights and the
  mixture probabilities have become stable.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{ESS}} and
  \code{\link{PMC}}.}
\examples{### See the PMC function for an example.}
\keyword{Plot}