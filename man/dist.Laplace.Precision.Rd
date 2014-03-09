\name{dist.Laplace.Precision}
\alias{dlaplacep}
\alias{plaplacep}
\alias{qlaplacep}
\alias{rlaplacep}
\title{Laplace Distribution: Precision Parameterization}
\description{
  These functions provide the density, distribution function, quantile
  function, and random generation for the univariate, symmetric, Laplace
  distribution with location parameter \eqn{\mu}{mu} and precision
  parameter \eqn{\tau}{tau}, which is the inverse of the usual scale
  parameter, \eqn{\lambda}{lambda}.
}
\usage{
dlaplacep(x, mu=0, tau=1, log=FALSE)
plaplacep(q, mu=0, tau=1)
qlaplacep(p, mu=0, tau=1)
rlaplacep(n, mu=0, tau=1)
}
\arguments{
  \item{x, q}{These are each a vector of quantiles.}
  \item{p}{This is a vector of probabilities.}
  \item{n}{This is the number of observations, which must be a positive
    integer that has length 1.}
  \item{mu}{This is the location parameter \eqn{\mu}{mu}.}
  \item{tau}{This is the precision parameter \eqn{\tau}{tau}, which
    must be positive.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
  \item Application: Continuous Univariate
  \item Density: \eqn{p(\theta) = \frac{\tau}{2} \exp(-\tau |\theta-\mu|)}{p(theta) = (tau/2) * exp(-tau*abs(theta-mu))}
  \item Inventor: Pierre-Simon Laplace (1774)
  \item Notation 1: \eqn{\theta \sim \mathrm{Laplace}(\mu,\tau^{-1})}{theta ~
  Laplace(mu, tau^(-1))}
  \item Notation 2: \eqn{\theta \sim \mathcal{L}(\mu, \tau^{-1})}{theta ~
    L(mu, tau^(-1))}
  \item Notation 3: \eqn{p(\theta) =
    \mathrm{Laplace}(\mu,\tau^{-1})}{p(theta) = Laplace(mu, tau^(-1))}
  \item Notation 4: \eqn{p(\theta) = \mathcal{L}(\theta | \mu,
    \tau^{-1})}{p(theta) = L(theta | mu, tau(-1))}
  \item Parameter 1: location parameter \eqn{\mu}{mu}
  \item Parameter 2: precision parameter \eqn{\tau > 0}{tau > 0}
  \item Mean: \eqn{E(\theta) = \mu}{E(theta) = mu}
  \item Variance: \eqn{var(\theta) = 2 \tau^{-2}}{var(theta) = 2*tau^{-2}}
  \item Mode: \eqn{mode(\theta) = \mu}{mode(theta) = mu}
}

The Laplace distribution is also called the double exponential 
distribution, because it looks like two exponential distributions back to
back with respect to location \eqn{\mu}{mu}. It is also called the
``First Law of Laplace'', just as the normal distribution is referred to
as the ``Second Law of Laplace''. The Laplace distribution is
symmetric with respect to \eqn{\mu}{mu}, though there are asymmetric
versions of the Laplace distribution. The PDF of the Laplace
distribution is reminiscent of the normal distribution; however,
whereas the normal distribution is expressed in terms of the squared
difference from the mean \eqn{\mu}{mu}, the Laplace density is
expressed in terms of the absolute difference from the mean,
\eqn{\mu}{mu}. Consequently, the Laplace distribution has fatter
tails than the normal distribution. It has been argued that the Laplace
distribution fits most things in nature better than the normal
distribution. Elsewhere, there are a large number of extensions to the
Laplace distribution, including asymmetric versions and
multivariate versions, among many more. These functions provide the
precision parameterization for convenience and familiarity in Bayesian
inference.
}
\value{
  \code{dlaplacep} gives the density,
  \code{plaplacep} gives the distribution function,
  \code{qlaplacep} gives the quantile function, and
  \code{rlaplacep} generates random deviates.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{dalaplace}},
  \code{\link{dexp}},
  \code{\link{dlaplace}},
  \code{\link{dmvl}},
  \code{\link{dnorm}},
  \code{\link{dnormp}}, and
  \code{\link{dnormv}}.
}
\examples{
library(LaplacesDemon)
x <- dlaplacep(1,0,1)
x <- plaplacep(1,0,1)
x <- qlaplacep(0.5,0,1)
x <- rlaplacep(100,0,1)

#Plot Probability Functions
x <- seq(from=-5, to=5, by=0.1)
plot(x, dlaplacep(x,0,0.5), ylim=c(0,1), type="l", main="Probability Function",
     ylab="density", col="red")
lines(x, dlaplacep(x,0,1), type="l", col="green")
lines(x, dlaplacep(x,0,2), type="l", col="blue")
legend(2, 0.9, expression(paste(mu==0, ", ", tau==0.5),
     paste(mu==0, ", ", tau==1), paste(mu==0, ", ", tau==2)),
     lty=c(1,1,1), col=c("red","green","blue"))
}
\keyword{Distribution}