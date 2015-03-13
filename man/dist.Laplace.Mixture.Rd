\name{dist.Laplace.Mixture}
\alias{dlaplacem}
\alias{plaplacem}
\alias{rlaplacem}
\title{Mixture of Laplace Distributions}
\description{
  These functions provide the density, cumulative, and random generation
  for the mixture of univariate Laplace distributions with probability
  \eqn{p}, location \eqn{\mu}{mu} and scale \eqn{\sigma}{sigma}.
}
\usage{
dlaplacem(x, p, location, scale, log=FALSE)
plaplacem(q, p, location, scale)
rlaplacem(n, p, location, scale)
}
\arguments{
  \item{x,q}{This is vector of values at which the density will be evaluated.}
  \item{p}{This is a vector of length \eqn{M} of probabilities for \eqn{M}
    components. The sum of the vector must be one.}
  \item{n}{This is the number of observations, which must be a positive
    integer that has length 1.}
  \item{location}{This is a vector of length \eqn{M} that is the location
    parameter \eqn{\mu}{mu}.}
  \item{scale}{This is a vector of length \eqn{M} that is the scale
    parameter \eqn{\sigma}{sigma}, which must be positive.}
  \item{log}{Logical. If \code{TRUE}, then the logarithm of the density is
    returned.}
}
\details{
  \itemize{
  \item Application: Continuous Univariate
  \item Density: \eqn{p(\theta) = \sum p_i \mathcal{L}(\mu_i,
    \sigma_i)}{p(theta) = sum p[i] L(mu[i], sigma[i])}
  \item Inventor: Unknown
  \item Notation 1: \eqn{\theta \sim \mathcal{L}(\mu, \sigma)}{theta ~ L(mu, sigma)}
  \item Notation 2: \eqn{p(\theta) = \mathcal{L}(\theta | \mu,
    \sigma)}{p(theta) = L(theta | mu, sigma)}
  \item Parameter 1: location parameters \eqn{\mu}{mu}
  \item Parameter 2: scale parameters \eqn{\sigma > 0}{sigma > 0}
  \item Mean: \eqn{E(\theta) = \sum p_i \mu_i}{E(theta) = sum p[i] mu[i]}
  \item Variance: 
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

A Laplace mixture distribution is a combination of Laplace probability
distributions.

One of many applications of Laplace mixture distributions is the Laplace
Mixture Model (LMM).
}
\value{
  \code{dlaplacem} gives the density,
  \code{plaplacem} returns the CDF, and
  \code{rlaplacem} generates random deviates.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{ddirichlet}} and
  \code{\link{dlaplace}}.
}
\examples{
library(LaplacesDemon)
p <- c(0.3,0.3,0.4)
mu <- c(-5, 1, 5)
sigma <- c(1,2,1)
x <- seq(from=-10, to=10, by=0.1)
plot(x, dlaplacem(x, p, mu, sigma, log=FALSE), type="l") #Density
plot(x, plaplacem(x, p, mu, sigma), type="l") #CDF
plot(density(rlaplacem(10000, p, mu, sigma))) #Random Deviates
}
\keyword{Distribution}