\name{dist.Asymmetric.Log.Laplace}
\alias{dallaplace}
\alias{pallaplace}
\alias{qallaplace}
\alias{rallaplace}
\title{Asymmetric Log-Laplace Distribution}
\description{
  These functions provide the density, distribution function, quantile
  function, and random generation for the univariate, asymmetric,
  log-Laplace distribution with location parameter \eqn{\mu}{mu},
  scale parameter \eqn{\lambda}{\lambda}, and asymmetry or skewness
  parameter \eqn{\kappa}{kappa}.
}
\usage{
dallaplace(x, location=0, scale=1, kappa=1, log=FALSE)
pallaplace(q, location=0, scale=1, kappa=1)
qallaplace(p, location=0, scale=1, kappa=1)
rallaplace(n, location=0, scale=1, kappa=1)
}
\arguments{
  \item{x, q}{These are each a vector of quantiles.}
  \item{p}{This is a vector of probabilities.}
  \item{n}{This is the number of observations, which must be a positive
    integer that has length 1.}
  \item{location}{This is the location parameter \eqn{\mu}{mu}.}
  \item{scale}{This is the scale parameter \eqn{\lambda}{lambda}, which
    must be positive.}
  \item{kappa}{This is the asymmetry or skewness parameter
    \eqn{\kappa}{kappa}, which must be positive.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
  \item Application: Continuous Univariate
  \item Density 1: \eqn{p(\theta) = \exp(-\mu)\frac{(\sqrt(2)\kappa / \lambda)(\sqrt(2) / \lambda\kappa)}{(\sqrt(2)\kappa / \lambda)+(\sqrt(2) / (\lambda\kappa))} \exp(-(\frac{\sqrt(2)\kappa}{\lambda})+1), \quad \theta \ge \exp(\mu)}{p(theta) = exp(-mu) * (sqrt(2)*kappa/lambda) * (sqrt(2)/(lambda*kappa)) / ((sqrt(2)*kappa/lambda)+(sqrt(2)/(lambda*kappa))) * exp(-((sqrt(2)*kappa/lambda)+1)), theta >= exp(mu)}
  \item Density 2: \eqn{p(\theta) = \exp(-\mu) \frac{(\sqrt(2)\kappa / \lambda) (\sqrt(2) / (\lambda\kappa))}{(\sqrt(2)\kappa / \lambda) + (\sqrt(2) / (\lambda\kappa))} \exp(\frac{\sqrt(2)(\log(\theta)-\mu)}{\lambda\kappa} - (\log(\theta)-\mu)), \quad \theta < \exp(\mu)}{p(theta) = exp(-mu) * (sqrt(2)*kappa/lambda) * (sqrt(2)/(lambda*kappa)) / ((sqrt(2)*kappa/lambda)+(sqrt(2)/(lambda*kappa))) * exp(((sqrt(2)*(log(theta)-mu)) / (lambda*kappa)) - (log(theta)-mu)), theta < exp(mu)}
  \item Inventor: Pierre-Simon Laplace
  \item Notation 1: \eqn{\theta \sim \mathcal{ALL}(\mu, \lambda,
    \kappa)}{theta ~ ALL(mu, lambda, kappa)}
  \item Notation 2: \eqn{p(\theta) = \mathcal{ALL}(\theta | \mu,
    \lambda, \kappa)}{p(theta) = ALL(theta | mu, lambda, kappa)}
  \item Parameter 1: location parameter \eqn{\mu}{mu}
  \item Parameter 2: scale parameter \eqn{\lambda > 0}{lambda > 0}
  \item Mean: \eqn{E(\theta) = }{E(theta) = }
  \item Variance: \eqn{var(\theta) = }{var(theta) = }
  \item Mode: \eqn{mode(\theta) = }{mode(theta) = }
}

The univariate, asymmetric log-Laplace distribution is derived from the
Laplace distribution. Multivariate and symmetric versions also exist.

These functions are similar to those in the \code{VGAM} package.
}
\value{
  \code{dallaplace} gives the density,
  \code{pallaplace} gives the distribution function,
  \code{qallaplace} gives the quantile function, and
  \code{rallaplace} generates random deviates.
}
\references{
  Kozubowski, T. J. and Podgorski, K. (2003). "Log-Laplace Distributions".
  \emph{International Mathematical Journal}, 3, p. 467--495.
}
\seealso{
  \code{\link{dalaplace}},
  \code{\link{dexp}},
  \code{\link{dlaplace}},
  \code{\link{dlaplacep}},
  \code{\link{dllaplace}},
  \code{\link{dmvl}},
  \code{\link{dnorm}},
  \code{\link{dnormp}},
  \code{\link{dnormv}}.
}
\examples{
library(LaplacesDemon)
x <- dallaplace(1,0,1,1)
x <- pallaplace(1,0,1,1)
x <- qallaplace(0.5,0,1,1)
x <- rallaplace(100,0,1,1)

#Plot Probability Functions
x <- seq(from=0.1, to=10, by=0.1)
plot(x, dallaplace(x,0,1,0.5), ylim=c(0,1), type="l", main="Probability Function",
     ylab="density", col="red")
lines(x, dallaplace(x,0,1,1), type="l", col="green")
lines(x, dallaplace(x,0,1,5), type="l", col="blue")
legend(5, 0.9, expression(paste(mu==0, ", ", lambda==1, ", ", kappa==0.5),
     paste(mu==0, ", ", lambda==1, ", ", kappa==1),
     paste(mu==0, ", ", lambda==1, ", ", kappa==5)),
     lty=c(1,1,1), col=c("red","green","blue"))
}
\keyword{Distribution}