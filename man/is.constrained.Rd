\name{is.constrained}
\alias{is.constrained}
\title{Logical Check of Constraints}
\description{
  This function provides a logical test of constraints for each initial
  value or parameter for a model specification, given data.
}
\usage{
is.constrained(Model, Initial.Values, Data)
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
  This function is useful for testing whether or not initial values
  changed due to constraints when being passed through a \code{Model}
  specification function. If any initial value changes, then the
  constrained values that are ouput in the fifth component of the
  \code{Model} specification are suitable as initial values, not the
  tested initial values.

  A parameter may be constrained and this function may not discover the
  constraint, since the discovery depends on the initial values and
  whether or not they change as they are passed through the model.
}
\value{
  The \code{is.constrained} function returns a logical vector, equal in
  length to the number of initial values. Each element receives
  \code{TRUE} if the corresponding initial value changed due to a
  constraint, or \code{FALSE} if it did not.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{LaplacesDemon}}.
}
\keyword{Utility}