\name{interval}
\alias{interval}
\title{Constrain to Interval}
\description{
  This function constrains the value(s) of a scalar, vector, matrix, or
  array to a specified interval, \eqn{[a,b]}. In Bayesian inference, it
  is often used both to truncate a parameter to an interval, such as
  \eqn{p(\theta) \in [a,b]}{p(theta) in [a,b]}. The \code{interval}
  function is often used in conjunction with the \code{\link{dtrunc}} 
  function to truncate the prior probability distribution associated
  with the constrained parameter. While \code{\link{dtrunc}} prevents
  assigning density outside of its interval and re-estimates density
  within the interval, the \code{interval} function is used to prevent
  the parameter from moving outside of the interval in the first place.

  After the parameter is constrained to an interval in
  \code{\link{IterativeQuadrature}}, \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}}, \code{\link{PMC}}, or
  \code{\link{VariationalBayes}}, the constrained parameter should be
  updated back into the \code{parm} vector, so the algorithm knows it
  has been constrained.

  This is unrelated to the probability interval (see
  \code{\link{p.interval}} and \code{\link{LPL.interval}}).
}
\usage{
interval(x, a=-Inf, b=Inf, reflect=TRUE)
}
\arguments{
  \item{x}{This required argument is a scalar, vector, matrix or array,
    and its elements will be constrained to the interval
    [\code{a},\code{b}].}
  \item{a}{This optional argument allows the specification of the lower
    bound of the interval, and defaults to \code{-Inf}.}
  \item{b}{This optional argument allows the specification of the upper
    bound of the interval, and defaults to \code{Inf}.}
  \item{reflect}{Logical. When \code{TRUE}, a value outside of the
  constrained interval is reflected or bounced back into the
  interval. When \code{FALSE}, a value outside of the interval is
  assigned the nearest boundary of the interval. This argument defaults
  to \code{TRUE}.}
}
\details{
  It is common for a parameter to be constrained to an interval. The
  \code{interval} function provides two methods of constraining
  proposals. The default is to reflect an out-of-bounds proposal off of
  the boundaries until the proposal is within the specified
  interval. This is rare in the literature but works very well in
  practice. The other method does not reflect off of boundaries, but
  sets the value equal to the violated boundary. This is also rare in
  the literature and is not generally recommended.

  If the \code{interval} function is unacceptable, then there are
  several alternatives.

  It is common to re-parameterize by transforming the constrained
  parameter to the real line. For example, a positive-only scale
  parameter may be log-transformed. A parameter that is re-parameterized
  to the real line often mixes better in MCMC, exhibiting a higher
  effective sample size (\code{\link{ESS}}), and each evaluation of the
  model specification function is faster as well. However, without a
  hard constraint, it remains possible for the transformed parameter
  still become problematic, such as a log-transformed scale parameter
  that reaches negative infinity. This is much more common in the
  literature.

  Another method is to allow the parameters to move outside of the
  desired, constrained interval in MCMC during the model update, and
  when the model update is finished, to discard any samples outside of
  the constraint boundaries. This is a method of rejecting unacceptable
  proposals in regions of zero probability. However, it is possible for
  parameters to remain outside of acceptable bounds long enough to be
  problematic.

  In \code{\link{LaplacesDemon}}, the Gibbs sampler allows more control
  in the FC function, where a user can customize how constraints are
  handled.
}
\value{
  The \code{interval} function returns a scalar, vector, matrix, or
  array in accord with its argument, \code{x}. Each element is
  constrained to the interval [\code{a},\code{b}].
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{dtrunc}},
  \code{\link{ESS}},
  \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{LPL.interval}},
  \code{\link{PMC}},
  \code{\link{p.interval}},
  \code{\link{VariationalBayes}}.
}
\examples{
#See the Examples vignette for numerous examples.
library(LaplacesDemon)
x <- 2
interval(x,0,1)
X <- matrix(runif(25,-2,2),5,5)
interval(X,-1,1)
}
\keyword{Utility}