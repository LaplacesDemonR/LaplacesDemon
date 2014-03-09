\name{is.bayesian}
\alias{is.bayesian}
\title{Logical Check of a Bayesian Model}
\description{
  This function provides a logical test of whether or not a \code{Model}
  specification function is Bayesian.
}
\usage{
is.bayesian(Model, Initial.Values, Data)
}
\arguments{
  \item{Model}{This is a model specification function. For more
    information, see the \code{\link{LaplacesDemon}} function.}
  \item{Initial.Values}{This is a vector of initial values, or current
    parameter values. For more information, see the
    \code{\link{LaplacesDemon}} function.}
  \item{Data}{This is a list of data. For more information, see the
    \code{\link{LaplacesDemon}} function.}
}
\details{
  This function tests whether or not a model is Bayesian by comparing
  the first two returned arguments: the logarithm of the unnormalized
  joint posterior density (\code{LP}) and deviance (\code{Dev}). The
  deviance (D) is
  
  \deqn{\mathrm{D} = -2 \mathrm{LL}}{D = -2 LL},

  where LL is the log-likelihood. Consequently,

  \deqn{\mathrm{LL} = \mathrm{D} / -2}{LL = D / -2},

  and LP is the sum of LL and prior probability densities. If LP = LL,
  then the model is not Bayesian, because prior densities are absent.
}
\value{
  The \code{is.bayesian} function returns a logical value of \code{TRUE}
  when the model is Bayesian, and \code{FALSE} otherwise.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{LaplacesDemon}}.
}
\keyword{Utility}