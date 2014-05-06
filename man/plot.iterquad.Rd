\name{plot.iterquad}
\alias{plot.iterquad}
\title{Plot the output of \code{\link{IterativeQuadrature}}}
\description{
  This may be used to plot, or save plots of, the iterated history of
  the parameters and, if posterior samples were taken, density plots of
  parameters and monitors in an object of class \code{iterquad}.
}
\usage{\method{plot}{iterquad}(x, Data, PDF=FALSE, Parms, \dots)}
\arguments{
  \item{x}{
    This required argument is an object of class \code{iterquad}.}
  \item{Data}{
    This required argument must receive the list of data that was
    supplied to \code{\link{IterativeQuadrature}} to create the object
    of class \code{iterquad}.}
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
  \item{\dots}{Additional arguments are unused.}
}
\details{
  The plots are arranged in a \eqn{2 \times 2}{2 x 2} matrix. The
  purpose of the iterated history plots is to show how the value of each
  parameter and the deviance changed by iteration as the
  \code{\link{IterativeQuadrature}} attempted to fit a normal
  distribution to the marginal posterior distributions.

  The plots on the right show several densities, described below.

  \itemize{
    \item The transparent black density is the normalized quadrature
    weights for non-standard normal distributions, \eqn{M}. For
    multivariate quadrature, there are often multiple weights at a given
    node, and the average \eqn{M} is shown. Vertical black lines
    indicate the nodes.
    \item The transparent red density is the normalized LP weights. For
    multivariate quadrature, there are often multiple weights at a given
    node, and the average normalized and weighted LP is shown. Vertical
    red lines indicate the nodes.
    \item The transparent green density is the normal density implied
    given the conditional mean and conditional variance.
    \item The transparent blue density is the kernel density estimate
    of posterior samples generated with Sampling Importance Resampling.
    This is plotted only if the algorithm converged, and if
    \code{sir=TRUE}.
  }
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{\code{\link{IterativeQuadrature}}}
\examples{### See the IterativeQuadrature function for an example.}
\keyword{Plot}