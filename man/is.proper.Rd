\name{is.proper}
\alias{is.proper}
\title{Logical Check of Propriety}
\description{
  This function provides a logical check of the propriety of a
  univariate prior probability distribution or the joint posterior
  distribution.
}
\usage{
is.proper(f, a, b, tol=1e-5)
}
\arguments{
  \item{f}{This is either a probability density function or an object of
    class \code{demonoid}, \code{laplace}, \code{pmc}, or \code{vb}.}
  \item{a}{This is the lower limit of integration, and may be negative
    infinity.}
  \item{b}{This is the upper limit of integration, and may be positive
    infinity.}
  \item{tol}{This is the tolerance, and indicates the allowable
    difference from one.}
}
\details{
  A proper probability distribution is a probability distribution that
  integrates to one, and an improper probability distribution does not
  integrate to one. If a probability distribution integrates to any
  positive and finite value other than one, then it is an improper
  distribution, but is merely unnormalized. An unnormalized distribution
  may be multiplied by a constant so that it integrates to one.

  In Bayesian inference, the posterior probability distribution should
  be proper. An improper prior distribution can cause an improper
  posterior distribution. When the posterior distribution is improper,
  inferences are invalid, it is non-integrable, and Bayes factors cannot
  be used (though there are exceptions).

  To avoid these problems, it is suggested that the prior probability
  distribution should be proper, though it is possible to use an
  improper prior distribution and have it result in a proper posterior
  distribution.

  To check the propriety of a univariate prior probability distribution,
  create a function \code{f}. For example, to check the propriety of a
  vague normal distribution, such as

  \deqn{\theta \sim \mathcal{N}(0,1000)}{theta ~ N(0,1000)}

  the function is \code{function(x){dnormv(x,0,1000)}}. Next, set the
  lower and upper limits of integration, \code{a} and
  \code{b}. Internally, this function calls \code{integrate} from base
  R, which uses adaptive quadrature. By using \eqn{f(x)} as shorthand
  for the specified function, \code{is.proper} will check to see if the
  area of the following integral is one:

  \deqn{\int^b_a f(x)dx}{integral from a to b of f(x)dx}

  Multivariate prior probability distributions currently cannot be
  checked for approximate propriety. This is currently unavailable in
  this package.

  To check the propriety of the joint posterior distribution, the only
  argument to be supplied is an object of class \code{demonoid},
  \code{iterquad}, \code{laplace}, \code{pmc}, or \code{vb}. The
  \code{is.proper} function checks the logarithm of the marginal
  likelihood (see \code{\link{LML}}) for a finite value, and returns
  \code{TRUE} when the LML is finite. This indicates that the marginal
  likelihood is finite for all observed \eqn{\textbf{y}}{y} in the model
  data set. This implies:

  \deqn{\int p(\theta|\textbf{y})p(\theta)d\theta < \infty}{integral of
  p(theta|y)p(theta)dtheta < Inf}

  If the object is of class \code{demonoid} and the algorithm was
  adaptive, or if the object is of class \code{iterquad},
  \code{laplace}, or \code{vb} and the algorithm did not converge, then
  \code{is.proper} will return \code{FALSE} because LML was not
  estimated. In this case, it is possible for the joint posterior to be
  proper, but \code{is.proper} will be unable to determine propriety
  without the estimate of LML. If desired, the \code{\link{LML}} may be
  estimated by the user, and if it is finite, then the joint posterior
  distribution is proper.
}
\value{
  The \code{is.proper} function returns a logical value indicating
  whether or not the univariate prior or joint posterior probability
  distribution integrates to one within its specified limits.
  \code{TRUE} is returned for a proper univariate probability
  distribution.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{dnormv}},
  \code{\link{integrate}},
  \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{LML}},
  \code{\link{PMC}}, and
  \code{\link{VariationalBayes}}.
}
\examples{
library(LaplacesDemon)
### Prior Probability Distribution
is.proper(function(x) {dnormv(x,0,1000)}, -Inf, Inf) #x ~ N(0,1000)
is.proper(function(x) {dhalfcauchy(x,25)}, 0, Inf) #x ~ HC(25)
is.proper(function(x) {dunif(x,0,1)}, 0, 1) #x ~ U(0,1)
is.proper(function(x) {dunif(x,-Inf,Inf)}, -Inf, Inf) #x ~ U(-Inf,Inf)
### Joint Posterior Distribution
##This assumes that Fit is an object of class demonoid, iterquad,
##  laplace, or pmc
#is.proper(Fit)
}
\keyword{Diagnostic}
\keyword{Distribution}
\keyword{Utility}