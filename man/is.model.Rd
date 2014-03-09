\name{is.model}
\alias{is.model}
\title{Logical Check of a Model}
\description{
  This function provides a logical test of whether or not a \code{Model}
  specification function meets mininum requirements to be considered as
  such.
}
\usage{
is.model(Model, Initial.Values, Data)
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
  This function tests for minimum criteria for \code{Model} to be
  considered a model specification function. Specifically, it tests:

  \itemize{
    \item \code{Model} must be a function
    \item \code{Model} must execute without errors
    \item \code{Model} must return a list
    \item \code{Model} must have five components in the list
    \item The first component must be named LP and have length 1
    \item The second component must be named Dev and have length 1
    \item The third component must be named Monitor
    \item The lengths of Monitor and mon.names must be equal
    \item The fourth component must be named yhat
    \item The fifth component must be named parm
    \item The lengths of parm and parm.names must be equal
  }

  This function is not extensive, and checks only for these minimum
  criteria. Additional checks are conducted in
  \code{\link{IterativeQuadrature}}, \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}}, \code{\link{PMC}}, and
  \code{\link{VariationalBayes}}.
}
\value{
  The \code{is.model} function returns a logical value of \code{TRUE}
  when \code{Model} meets minimum criteria of a model specification
  function, and \code{FALSE} otherwise.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{PMC}}, and
  \code{\link{VariationalBayes}}.
}
\keyword{Utility}