\name{is.appeased}
\alias{is.appeased}
\title{Appeased}
\description{
  This function returns \code{TRUE} if Laplace's Demon is appeased by the
  object of class \code{demonoid}, and \code{FALSE} otherwise. If
  appeased, then the object passes several tests that indicate potential
  convergence of the Markov chains.
}
\usage{
is.appeased(x)
}
\arguments{
  \item{x}{This is an object of class \code{demonoid}.}
}
\details{
  After updating a model with the \code{\link{LaplacesDemon}} function,
  an output object is created. The output object is of class
  \code{demonoid}. The object may be passed to the \code{\link{Consort}}
  function, which will apply several criteria regarding the potential
  convergence of its Markov chains. If all criteria are met, then
  Laplace's Demon is appeased. Otherwise, Laplace's Demon suggests R
  code to be copy/pasted and executed. The \code{\link{Consort}}
  function prints a large amount of information to the screen. The
  \code{is.appeased} function may be applied as an alternative, though
  it only informs the user as to whether or not Laplace's Demon was
  appeased, as \code{TRUE} or \code{FALSE}.
}
\value{
  The \code{is.appeased} function returns a logical value indicating
  whether or not the supplied object passes several potential Markov
  chain convergence criteria. If the object passes all criteria, then
  Laplace's Demon is appeased, and the logical value returned is
  \code{TRUE}.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{Consort}} and
  \code{\link{LaplacesDemon}}.
}
\keyword{Diagnostic}
\keyword{Utility}