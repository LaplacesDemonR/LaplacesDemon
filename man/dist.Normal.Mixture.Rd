\name{dist.Normal.Mixture}
\alias{dnormm}
\alias{pnormm}
\alias{rnormm}
\title{Mixture of Normal Distributions}
\description{
  These functions provide the density, cumulative, and random generation
  for the mixture of univariate normal distributions with probability
  \eqn{p}, mean \eqn{\mu}{mu} and standard deviation \eqn{\sigma}{sigma}.
}
\usage{
dnormm(x, p, mu, sigma, log=FALSE)
pnormm(q, p, mu, sigma, lower.tail=TRUE, log.p=FALSE)
rnormm(n, p, mu, sigma)
}
\arguments{
  \item{x,q}{This is vector of values at which the density will be evaluated.}
  \item{p}{This is a vector of length \eqn{M} of probabilities for \eqn{M}
    components. The sum of the vector must be one.}
  \item{n}{This is the number of observations, which must be a positive
    integer that has length 1.}
  \item{mu}{This is a vector of length \eqn{M} that is the mean
    parameter \eqn{\mu}{mu}.}
  \item{sigma}{This is a vector of length \eqn{M} that is the standard
    deviation parameter \eqn{\sigma}{sigma}, which must be positive.}
  \item{lower.tail}{Logical. This defaults to \code{TRUE}.}
  \item{log,log.p}{Logical. If \code{TRUE}, then probabilities
    \eqn{p} are given as \eqn{\log(p)}{log(p)}.}
}
\details{
  \itemize{
  \item Application: Continuous Univariate
  \item Density: \eqn{p(\theta) = \sum p_i \mathcal{N}(\mu_i,
    \sigma^2_i)}{p(theta) = sum p[i] N(mu[i], sigma[i]^2)}
  \item Inventor: Unknown
  \item Notation 1: \eqn{\theta \sim \mathcal{N}(\mu, \sigma^2)}{theta ~ N(mu, sigma^2)}
  \item Notation 2: \eqn{p(\theta) = \mathcal{N}(\theta | \mu,
    \sigma^2)}{p(theta) = N(theta | mu, sigma^2)}
  \item Parameter 1: mean parameters \eqn{\mu}{mu}
  \item Parameter 2: standard deviation parameters \eqn{\sigma > 0}{sigma > 0}
  \item Mean: \eqn{E(\theta) = \sum p_i \mu_i}{E(theta) = sum p[i] mu[i]}
  \item Variance: \eqn{var(\theta) = \sum p_i \sigma^{0.5}_i}{var(theta)
    = sum p[i] sigma[i]^(0.5)}
  \item Mode:
}

A mixture distribution is a probability distribution that is a
combination of other probability distributions, and each distribution is
called a mixture component, or component. A probability (or weight)
exists for each component, and these probabilities sum to one. A
mixture distribution (though not these functions here in particular) may
contain mixture components in which each component is a different
probability distribution. Mixture distributions are very flexible, and
are often used to represent a complex distribution with an unknown
form. When the number of mixture components is unknown, Bayesian
inference is the only sensible approach to estimation.

A normal mixture, or Gaussian mixture, distribution is a combination of
normal probability distributions.
}
\value{
  \code{dnormm} gives the density,
  \code{pnormm} returns the CDF, and
  \code{rnormm} generates random deviates.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{ddirichlet}} and
  \code{\link{dnorm}}.
}
\examples{
library(LaplacesDemon)
p <- c(0.3,0.3,0.4)
mu <- c(-5, 1, 5)
sigma <- c(1,2,1)
x <- seq(from=-10, to=10, by=0.1)
plot(x, dnormm(x, p, mu, sigma, log=FALSE), type="l") #Density
plot(x, pnormm(x, p, mu, sigma), type="l") #CDF
plot(density(rnormm(10000, p, mu, sigma))) #Random Deviates
}
\keyword{Distribution}