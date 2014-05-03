\name{dist.Categorical}
\alias{dcat}
\alias{qcat}
\alias{rcat}
\title{Categorical Distribution}
\description{
  This is the density and random deviates function for the categorical
  distribution with probabilities parameter \eqn{p}.
}
\usage{
dcat(x, p, log=FALSE)
qcat(pr, p, lower.tail=TRUE, log.pr=FALSE)
rcat(n, p)
}
\arguments{
  \item{x}{This is a vector of discrete data with \eqn{k} discrete
    categories, and is of length \eqn{n}. This function also accepts
    \eqn{x} after it has been converted to an \eqn{n \times k}{n x k}
    indicator matrix, such as with the \code{as.indicator.matrix} function.}
  \item{n}{This is the number of observations, which must be a positive
    integer that has length 1. When \code{p} is supplied to \code{rcat}
    as a matrix, \code{n} must equal the number of rows in \code{p}.}
  \item{p}{This is a vector of length \eqn{k} or \eqn{n \times k}{n x k}
    matrix of probabilities. The \code{qcat} function requires a
    vector.}
  \item{pr}{This is a vector of probabilities, or log-probabilities.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
  \item{log.pr}{Logical. if \code{TRUE}, probabilities \eqn{pr} are
    given as \eqn{\log(pr)}{log(pr)}.}
  \item{lower.tail}{Logical. if \code{TRUE} (default), probabilities
    are \eqn{Pr[X \le x]}{Pr[X <= x]}, otherwise,
    \eqn{Pr[X > x]}{Pr[X > x]}.}
}
\details{
    \itemize{
  \item Application: Discrete Univariate
  \item Density: \eqn{p(\theta) = \sum \theta p}{p(theta) = Sum (theta * p)}
  \item Inventor: Unknown (to me, anyway)
  \item Notation 1: \eqn{\theta \sim \mathcal{CAT}(p)}{theta ~ CAT(p)}
  \item Notation 2: \eqn{p(\theta) = \mathcal{CAT}(\theta | p)}{p(theta) = CAT(theta | p)}
  \item Parameter 1: probabilities \eqn{p}
  \item Mean: \eqn{E(\theta)}{E(theta)} = Unknown
  \item Variance: \eqn{var(\theta)}{var(theta)} = Unknown
  \item Mode: \eqn{mode(\theta)}{mode(theta)} = Unknown
}
Also called the discrete distribution, the categorical distribution
describes the result of a random event that can take on one of \eqn{k}
possible outcomes, with the probability \eqn{p} of each outcome
separately specified. The vector \eqn{p} of probabilities for each
event must sum to 1. The categorical distribution is often used, for
example, in the multinomial logit model. The conjugate prior is the
Dirichlet distribution.
}
\value{
  \code{dcat} gives the density and
  \code{rcat} generates random deviates.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{as.indicator.matrix}},
  \code{\link{ddirichlet}}, and
  \code{\link{dmultinom}}.
}
\examples{
library(LaplacesDemon)
dcat(x=1, p=c(0.3,0.3,0.4))
rcat(n=10, p=c(0.1,0.3,0.6))
}
\keyword{Distribution}
