\name{Gelfand.Diagnostic}
\alias{Gelfand.Diagnostic}
\title{Gelfand's Convergence Diagnostic}
\description{
  Gelfand et al. (1990) proposed a convergence diagnostic for Markov
  chains. The \code{Gelfand.Diagnostic} function is an interpretation of
  Gelfand's ``thick felt-tip pen'' MCMC convergence diagnostic. This
  diagnostic plots a series of kernel density plots at \eqn{k}
  intervals of cumulative samples. Given a vector of \eqn{S} samples
  from a marginal posterior distribution, \eqn{\theta}{theta}, multiple
  kernel density lines are plotted together, where each includes samples
  from a different interval. It is assumed that \code{\link{burnin}}
  iterations have been discarded.

  Gelfand et al. (1990) assert that convergence is violated when the
  plotted lines are farther apart than the width of a thick, felt-tip
  pen. This depends on the size of the plot, and, of course, the
  pen. The estimated width of a ``thick felt-tip pen'' is included as a
  black, vertical line. The pen in \code{Gelfand.Diagnostic} is included
  for historical reasons. This diagnostic requires numerous samples.
}
\usage{
Gelfand.Diagnostic(x, k=3, pen=FALSE)
}
\arguments{
  \item{x}{This required argument is a vector of marginal posterior
    samples, such as selected from the output of
    \code{\link{LaplacesDemon}}.}
  \item{k}{This argument specifies the number \eqn{k} of kernel
    density plots given cumulative intervals of samples. This argument
    defaults to \eqn{k=3}{k=3}.}
  \item{pen}{Logical. This argument defaults to \code{pen=FALSE}. When
    \code{pen=TRUE}, the thick felt-tip pen is included as a black,
    vertical line.}
}
\value{
  The \code{Gelfand.Diagnostic} returns a plot.
}
\references{
  Gelfand, A.E., Hills, S., Racine-Poon, A., and Smith,
  A.F.M. (1990). "Illustration of Bayesian Inference in Normal Data
  Models Using Gibbs Sampling". \emph{Journal of the American
  Statistical Association}, 85, p. 972--985.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{burnin}} and
  \code{\link{LaplacesDemon}}.
}
\examples{
library(LaplacesDemon)
x <- rnorm(1000)
Gelfand.Diagnostic(x)
}
\keyword{Diagnostic}
\keyword{MCMC}
