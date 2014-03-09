\name{is.constant}
\alias{is.constant}
\title{Logical Check of a Constant}
\description{
  This function provides a logical test of whether or not a vector is a
  constant.
}
\usage{
is.constant(x)
}
\arguments{
  \item{x}{This is a vector.}
}
\details{
  As opposed to a variable, a constant is a vector in which the elements
  contain less than or equal to one unique value.
}
\value{
  The \code{is.constant} function returns a logical result, reporting
  \code{TRUE} when a vector is a constant, or \code{FALSE} otherwise.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{unique}}
}
\examples{
library(LaplacesDemon)
is.constant(rep(1,10)) #TRUE
is.constant(1:10) #FALSE
}
\keyword{Utility}