\name{print.demonoid}
\alias{print.demonoid}
\title{Print an object of class \code{demonoid} to the screen.}
\description{
  This may be used to print the contents of an object of class
  \code{demonoid} to the screen.
}
\usage{\method{print}{demonoid}(x, \dots)}
\arguments{
  \item{x}{An object of class \code{demonoid} is required.}
  \item{\dots}{Additional arguments are unused.}
}
\details{
  If the user has an object of class \code{demonoid.hpc}, then the
  \code{print} function may still be used by specifying the chain as a
  component in a list, such as printing the second chain with
  \code{print(Fit[[2]])} when the \code{demonoid.hpc} object is named
  \code{Fit}, for example.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{Consort}},
  \code{\link{LaplacesDemon}}, and
  \code{\link{LaplacesDemon.hpc}}.
}
\examples{### See the LaplacesDemon function for an example.}
\keyword{print}
