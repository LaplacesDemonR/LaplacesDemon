\name{as.ppc}
\alias{as.ppc}
\title{As Posterior Predictive Check}
\description{
  This function converts an object of class \code{demonoid.val} to an
  object of class \code{demonoid.ppc}.
}
\usage{
as.ppc(x, set=3)
}
\arguments{
  \item{x}{This is an object of class \code{demonoid.val}.}
  \item{set}{This is an integer that indicates which list component is
    to be used. When \code{set=1}, the modeled data set is used. When
    \code{set=2}, the validation data set is used. When \code{set=3},
    both data sets are used.}
}
\details{
  After using the \code{\link{Validate}} function for holdout
  validation, it is often suggested to perform posterior predictive
  checks. The \code{as.ppc} function converts the output object of
  \code{\link{Validate}}, which is an object of class
  \code{demonoid.val}, to an object of class \code{demonoid.ppc}. The
  returned object is the same as if it were created with the
  \code{\link{predict.demonoid}} function, rather than the
  \code{\link{Validate}} function.

  After this conversion, the user may use posterior predictive checks,
  as usual, with the \code{\link{summary.demonoid.ppc}} function.
}
\value{
  The returned object is an object of class \code{demonoid.ppc}.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{predict.demonoid}},
  \code{\link{summary.demonoid.ppc}}, and
  \code{\link{Validate}}.
}
\keyword{Initial Values}