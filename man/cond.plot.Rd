\name{cond.plot}
\alias{cond.plot}
\title{Conditional Plots}
\description{
  This function provides several styles of conditional plots with base
  graphics.
}
\usage{
cond.plot(x, y, z, Style="smoothscatter")
}
\arguments{
  \item{x}{This required argument accepts a numeric vector.}
  \item{y}{This argument accepts a numeric vector, and is only used with
    some styles.}
  \item{z}{This required argument accepts a discrete vector.}
  \item{Style}{This argument specifies the style of plot, and accepts
    "boxplot", "densover" (density overlay), "hist", "scatter", or
    "smoothscatter".}
}
\details{
  The \code{cond.plot} function provides simple conditional plots with
  base graphics.  All plot styles are conditional upon \code{z}. Up to
  nine conditional plots are produced in a panel.

  Plots include:

  boxplot: y ~ x | z
  densover: f(x | z)
  hist: x | z
  scatter: x, y | z
  smoothscatter: x, y | z

  The \code{cond.plot} function is not intended to try to compete with
  some of the better graphics packages, but merely to provide simple
  functionality.
}
\value{
  Conditional plots are returned.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{joint.density.plot}} and
  \code{\link{joint.pr.plot}}.
}
\examples{
library(LaplacesDemon)
x <- rnorm(1000)
y <- runif(1000)
z <- rcat(1000, rep(1/4,4))
cond.plot(x, y, z, Style="smoothscatter")
}
\keyword{Plot}