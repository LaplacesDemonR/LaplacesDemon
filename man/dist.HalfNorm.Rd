\name{dist.HalfNormal}
\alias{dhalfnorm}
\alias{phalfnorm}
\alias{qhalfnorm}
\alias{rhalfnorm}
\title{Half-Normal Distribution}
\description{
  These functions provide the density, distribution function, quantile
  function, and random generation for the half-normal distribution.
}
\usage{
dhalfnorm(x, scale=sqrt(pi/2), log=FALSE)
phalfnorm(q, scale=sqrt(pi/2), lower.tail=TRUE, log.p=FALSE)
qhalfnorm(p, scale=sqrt(pi/2), lower.tail=TRUE, log.p=FALSE)
rhalfnorm(n, scale=sqrt(pi/2))
}
\arguments{
  \item{x,q}{These are each a vector of quantiles.}
  \item{p}{This is a vector of probabilities.}
  \item{n}{This is the number of observations, which must be a positive
    integer that has length 1.}
  \item{scale}{This is the scale parameter \eqn{\sigma}{sigma}, which
    must be positive.}
  \item{log,log.p}{Logical. If \code{log=TRUE}, then the logarithm of the
    density or result is returned.}
  \item{lower.tail}{Logical. If \code{lower.tail=TRUE} (default),
    probabilities are \eqn{Pr[X \le x]}{Pr[X <= x]}, otherwise,
    \eqn{Pr[X > x]}{Pr[X > x]}.}
}
\details{
  \itemize{
  \item Application: Continuous Univariate
  \item Density: \eqn{p(\theta) = \frac{2 \sigma}{\pi}
    \exp(-\frac{\theta^2 \sigma^2}{\pi}), \quad \theta \ge 0}{p(theta) = 2*sigma/pi e^-(theta^2*sigma^2/pi), theta >= 0}
  \item Inventor: Derived from the normal or Gaussian
  \item Notation 1: \eqn{\theta \sim \mathcal{HN}(\sigma)}{theta ~ HALF-N(sigma)}
  \item Notation 2: \eqn{p(\theta) = \mathcal{HN}(\theta | \sigma)}{p(theta) = HN(theta | sigma)}
  \item Parameter 1: scale parameter \eqn{\sigma > 0}{sigma > 0}
  \item Mean: \eqn{E(\theta) = \frac{1}{\sigma}}{E(theta) = 1 / sigma}
  \item Variance: \eqn{var(\theta) = \frac{\pi-2}{2 \sigma^2}}{var(theta) = (pi-2)/(2*sigma^2)}
  \item Mode: \eqn{mode(\theta) = 0}{mode(theta) = 0}
}

The half-normal distribution is recommended as a weakly informative prior
distribution for a scale parameter that may be useful as an alternative
to the half-Cauchy, half-t, or vague gamma.
}
\value{
  \code{dhalfnorm} gives the density,
  \code{phalfnorm} gives the distribution function,
  \code{qhalfnorm} gives the quantile function, and
  \code{rhalfnorm} generates random deviates.
}
\seealso{
  \code{\link{dnorm}},
  \code{\link{dnormp}}, and
  \code{\link{dnormv}}.
}
\examples{
library(LaplacesDemon)
x <- dhalfnorm(1)
x <- phalfnorm(1)
x <- qhalfnorm(0.5)
x <- rhalfnorm(10)

#Plot Probability Functions
x <- seq(from=0.1, to=20, by=0.1)
plot(x, dhalfnorm(x,0.1), ylim=c(0,1), type="l", main="Probability Function",
     ylab="density", col="red")
lines(x, dhalfnorm(x,0.5), type="l", col="green")
lines(x, dhalfnorm(x,1), type="l", col="blue")
legend(2, 0.9, expression(sigma==0.1, sigma==0.5, sigma==1),
     lty=c(1,1,1), col=c("red","green","blue"))
}
\keyword{Distribution}