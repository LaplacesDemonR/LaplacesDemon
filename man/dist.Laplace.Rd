\name{dist.Laplace}
\alias{dlaplace}
\alias{plaplace}
\alias{qlaplace}
\alias{rlaplace}
\title{Laplace Distribution: Univariate Symmetric}
\description{
  These functions provide the density, distribution function, quantile
  function, and random generation for the univariate, symmetric, Laplace
  distribution with location parameter \eqn{\mu}{mu} and scale
  parameter \eqn{\lambda}{\lambda}.
}
\usage{
dlaplace(x, location=0, scale=1, log=FALSE)
plaplace(q, location=0, scale=1)
qlaplace(p, location=0, scale=1)
rlaplace(n, location=0, scale=1)
}
\arguments{
  \item{x, q}{These are each a vector of quantiles.}
  \item{p}{This is a vector of probabilities.}
  \item{n}{This is the number of observations, which must be a positive
    integer that has length 1.}
  \item{location}{This is the location parameter \eqn{\mu}{mu}.}
  \item{scale}{This is the scale parameter \eqn{\lambda}{lambda}, which
    must be positive.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
  \item Application: Continuous Univariate
  \item Density: \eqn{p(\theta) = \frac{1}{2 \lambda} \exp(-\frac{|\theta - \mu|}{\lambda})}{p(theta) = (1 / (2*lambda)) * exp(-(abs(theta - mu) / lambda))}
  \item Inventor: Pierre-Simon Laplace (1774)
  \item Notation 1: \eqn{\theta \sim
    \mathrm{Laplace}(\mu,\lambda)}{theta ~ Laplace(mu, lambda)}
  \item Notation 2: \eqn{\theta \sim \mathcal{L}(\mu, \lambda)}{theta ~
    L(mu, lambda)}
  \item Notation 3: \eqn{p(\theta) = \mathrm{Laplace}(\theta | \mu,
    \lambda)}{p(theta) = Laplace(theta | mu, lambda)}
  \item Notation 4: \eqn{p(\theta) = \mathcal{L}(\theta | \mu,
    \lambda)}{p(theta) = L(theta | mu, lambda)}
  \item Parameter 1: location parameter \eqn{\mu}{mu}
  \item Parameter 2: scale parameter \eqn{\lambda > 0}{lambda > 0}
  \item Mean: \eqn{E(\theta) = \mu}{E(theta) = mu}
  \item Variance: \eqn{var(\theta) = 2 \lambda^2}{var(theta) = 2*lambda^2}
  \item Mode: \eqn{mode(\theta) = \mu}{mode(theta) = mu}
}

The Laplace distribution (Laplace, 1774) is also called the double
exponential distribution, because it looks like two exponential
distributions back to back with respect to location \eqn{\mu}{mu}. It is
also called the ``First Law of Laplace'', just as the normal
distribution is referred to as the ``Second Law of Laplace''. The
Laplace distribution is symmetric with respect to \eqn{\mu}{mu}, though
there are asymmetric versions of the Laplace distribution. The PDF of
the Laplace distribution is reminiscent of the normal distribution;
however, whereas the normal distribution is expressed in terms of the
squared difference from the mean \eqn{\mu}{mu}, the Laplace density is
expressed in terms of the absolute difference from the mean,
\eqn{\mu}{mu}. Consequently, the Laplace distribution has fatter
tails than the normal distribution. It has been argued that the Laplace
distribution fits most things in nature better than the normal
distribution.

There are many extensions to the Laplace distribution, such as the
asymmetric Laplace, asymmetric log-Laplace, Laplace (re-parameterized
for precision), log-Laplace, multivariate Laplace, and skew-Laplace,
among many more.

These functions are similar to those in the \code{VGAM} package.
}
\value{
  \code{dlaplace} gives the density,
  \code{plaplace} gives the distribution function,
  \code{qlaplace} gives the quantile function, and
  \code{rlaplace} generates random deviates.
}
\references{
  Laplace, P. (1774). "Memoire sur la Probabilite des Causes par les
  Evenements." l'Academie Royale des Sciences, 6, 621--656. English
  translation by S.M. Stigler in 1986 as "Memoir on the Probability
  of the Causes of Events" in \emph{Statistical Science}, 1(3),
  p. 359--378.
}
\seealso{
  \code{\link{dalaplace}},
  \code{\link{dallaplace}},
  \code{\link{dexp}},
  \code{\link{dlaplacep}},
  \code{\link{dllaplace}},
  \code{\link{dmvl}},
  \code{\link{dnorm}},
  \code{\link{dnormp}},
  \code{\link{dnormv}},
  \code{\link{dsdlaplace}}, and
  \code{\link{dslaplace}}.
}
\examples{
library(LaplacesDemon)
x <- dlaplace(1,0,1)
x <- plaplace(1,0,1)
x <- qlaplace(0.5,0,1)
x <- rlaplace(100,0,1)

#Plot Probability Functions
x <- seq(from=-5, to=5, by=0.1)
plot(x, dlaplace(x,0,0.5), ylim=c(0,1), type="l", main="Probability Function",
     ylab="density", col="red")
lines(x, dlaplace(x,0,1), type="l", col="green")
lines(x, dlaplace(x,0,2), type="l", col="blue")
legend(2, 0.9, expression(paste(mu==0, ", ", lambda==0.5),
     paste(mu==0, ", ", lambda==1), paste(mu==0, ", ", lambda==2)),
     lty=c(1,1,1), col=c("red","green","blue"))
}
\keyword{Distribution}