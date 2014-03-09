\name{plot.vb}
\alias{plot.vb}
\title{Plot the output of \code{\link{VariationalBayes}}}
\description{
  This may be used to plot, or save plots of, the iterated history of
  the parameters and variances, and if posterior samples were taken,
  density plots of parameters and monitors in an object of class
  \code{vb}.
}
\usage{\method{plot}{vb}(x, Data, PDF=FALSE, Parms, \dots)}
\arguments{
  \item{x}{
    This required argument is an object of class \code{vb}.}
  \item{Data}{
    This required argument must receive the list of data that was
    supplied to \code{\link{VariationalBayes}} to create the object
    of class \code{vb}.}
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
  The plots are arranged in a \eqn{3 \times 3}{3 x 3} matrix. The
  purpose of the iterated history plots is to show how the value of each
  parameter, variance, and the deviance changed by iteration as the
  \code{\link{VariationalBayes}} attempted to maximize the logarithm
  of the unnormalized joint posterior density. If the algorithm
  converged, and if \code{sir=TRUE} in
  \code{\link{VariationalBayes}}, then plots are produced of
  selected parameters and all monitored variables.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{\code{\link{VariationalBayes}}}
\examples{### See the VariationalBayes function for an example.}
\keyword{Plot}