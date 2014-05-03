\name{joint.pr.plot}
\alias{joint.pr.plot}
\title{Joint Probability Region Plot}
\description{
  Given two vectors, the \code{joint.pr.plot} function creates a
  scatterplot with ellipses of probability regions.
}
\usage{
joint.pr.plot(x, y, quantiles=c(0.25,0.50,0.75,0.95))
}
\arguments{
  \item{x}{This required argument is a vector.}
  \item{y}{This required argument is a vector.}
  \item{quantiles}{These are the quantiles for which probability regions
    are estimated with ellipses. The center of the ellipse is plotted by
    default. The 0.95 quantile creates a probability region that
    contains approximately 95\% of the data or samples of \code{x} and
    \code{y}. By default, four quantiles are included.}
}
\details{
  A probability region is also commonly called a credible region. For
  more information on probability regions, see \code{\link{p.interval}}.

  Joint probability regions are plotted only for two variables, and the
  regions are estimated with functions modified from the \code{car}
  package. The internal ellipse functions assume bivariate normality.

  This function is often used to plot posterior distributions of
  samples, such as from the \code{\link{LaplacesDemon}} function.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{LaplacesDemon}} and
  \code{\link{p.interval}}
}
\examples{
library(LaplacesDemon)
x <- rnorm(100)
y <- rnorm(100)
joint.pr.plot(x, y)
}
\keyword{Plot}
