\name{joint.density.plot}
\alias{joint.density.plot}
\title{Joint Density Plot}
\description{
  This function plots the joint kernel density from samples of two
  marginal posterior distributions.
}
\usage{
joint.density.plot(x, y, Title=NULL, contour=TRUE, color=FALSE, Trace=NULL)
}
\arguments{
  \item{x,y}{These are vectors consisting of samples from two marginal
    posterior distributions, such as those output by
    \code{\link{LaplacesDemon}} in components \code{Posterior1} (all
    samples) or \code{Posterior2} (stationary samples).}
  \item{Title}{This is the title of the joint posterior density plot.}
  \item{contour}{This logical argument indicates whether or not contour
    lines will be added to the plot. \code{contour} defaults to
    \code{TRUE}.}
  \item{color}{This logical argument indicates whether or not color will
    be added to the plot. \code{color} defaults to \code{FALSE}.}
  \item{Trace}{This argument defaults to \code{NULL}, in which case it
    does not trace the exploration of the joint density. To trace the
    exploration of the joint density, specify \code{Trace} with the
    beginning and ending iteration or sample. For example, to view the
    trace of the first ten iterations or samples, specify
    \code{Trace=c(1,10)}.}
}
\details{
  This function produces either a bivariate scatterplot that may have
  kernel density contour lines added, or a bivariate plot with kernel 
  density-influenced colors, which may also have kernel density contour
  lines added. A joint density plot may be more informative than two
  univariate density plots.

  The \code{Trace} argument allows the user to view the exploration of
  the joint density, such as from MCMC chain output. An efficient
  algorithm jumps to random points of the joint density, and an
  inefficient algorithm explores more slowly. The initial point of the
  trace (which is the first element passed to \code{Trace}) is plotted
  with a green dot. The user should consider plotting the joint density of
  the two marginal posterior distributions with the highest
  \code{\link{IAT}}, as identified with the
  \code{\link{PosteriorChecks}} function, since these are the two least
  efficient MCMC chains. Different sequences of iterations may be
  plotted. This `joint trace plot' may show behavior of the MCMC
  algorithm to the user.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{\code{\link{IAT}},
  \code{\link{LaplacesDemon}}, and
  \code{\link{PosteriorChecks}}
}
\examples{
library(LaplacesDemon)
X <- rmvn(1000, runif(2), diag(2))
joint.density.plot(X[,1], X[,2], Title="Joint Density Plot",
     contour=TRUE, color=FALSE)
joint.density.plot(X[,1], X[,2], Title="Joint Density Plot",
     contour=FALSE, color=TRUE)
joint.density.plot(X[,1], X[,2], Title="Joint Density Plot",
     contour=TRUE, color=TRUE)
joint.density.plot(X[,1], X[,2], Title="Joint Trace Plot",
     contour=FALSE, color=TRUE, Trace=c(1,10))
}
\keyword{Plot}