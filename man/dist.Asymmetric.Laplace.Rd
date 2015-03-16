\name{dist.Asymmetric.Laplace}
\alias{dalaplace}
\alias{palaplace}
\alias{qalaplace}
\alias{ralaplace}
\title{Asymmetric Laplace Distribution: Univariate}
\description{
  These functions provide the density, distribution function, quantile
  function, and random generation for the univariate, asymmetric Laplace
  distribution with location parameter \code{location}, scale parameter
  \code{scale}, and asymmetry or skewness parameter \code{kappa}.
}
\usage{
dalaplace(x, location=0, scale=1, kappa=1, log=FALSE)
palaplace(q, location=0, scale=1, kappa=1)
qalaplace(p, location=0, scale=1, kappa=1)
ralaplace(n, location=0, scale=1, kappa=1)
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
  \item Density: \eqn{p(\theta) = \frac{\kappa \sqrt{2}}{\lambda
      (1+\kappa^2)} \exp(-|\theta-\mu| \frac{\sqrt{2}}{\lambda}
    \kappa^{|\theta-\mu|} |\theta-\mu|)}{((kappa*sqrt(2)) / (lambda*(1+kappa^2))) * exp(-|x-mu| * (sqrt(2)/lambda) * kappa^(|x-mu|) * |x-mu|)}
  \item Inventor: Kotz, Kozubowski, and Podgorski (2001)
  \item Notation 1: \eqn{\theta \sim \mathcal{AL}(\mu, \lambda,
    \kappa)}{theta ~ AL(mu, lambda, kappa)}
  \item Notation 2: \eqn{p(\theta) = \mathcal{AL}(\theta | \mu, \lambda,
  \kappa)}{p(theta) = AL(theta | mu, lambda, kappa)}
  \item Parameter 1: location parameter \eqn{\mu}{mu}
  \item Parameter 2: scale parameter \eqn{\lambda > 0}{lambda > 0}
  \item Parameter 3: skewness parameter \eqn{\kappa > 0}{kappa > 0}
  \item Mean: \eqn{E(\theta) = \mu + \lambda \frac{1/\kappa - \kappa}{\sqrt{2}}}{E(theta) = mu + lambda*(1/kappa - kappa) / sqrt(2)}
  \item Variance: \eqn{var(\theta) = \lambda^2 \frac{1 + \kappa^4}{2 \kappa^2}}{var(theta) = lambda^2 * (1 + kappa^4) / (2*kappa^2)}
  \item Mode: \eqn{mode(\theta) = \mu}{mode(theta) = mu}
}

The asymmetric Laplace of Kotz, Kozubowski, and Podgorski (2001), also
referred to as AL, is an extension of the univariate, symmetric Laplace
distribution to allow for skewness. It is parameterized according to
three parameters: location parameter \eqn{\mu}{mu}, scale parameter
\eqn{\lambda}{lambda}, and asymmetry or skewness parameter
\eqn{\kappa}{kappa}. The special case of \eqn{\kappa=1}{kappa=1} is the
symmetric Laplace distribution. Values of \eqn{\kappa}{kappa} in the
intervals \eqn{(0, 1)}{(0,1)} and \eqn{(1, \infty)}{(1, Inf)},
correspond to positive (right) and negative (left) skewness,
respectively. The AL distribution is leptokurtic, and its kurtosis
ranges from 3 to 6 as \eqn{\kappa}{kappa} ranges from 1 to infinity. The
skewness of the AL has been useful in engineering and finance. As an
example, the AL distribution has been used as a replacement for
Gaussian-distributed GARCH residuals. There is also an extension to the
asymmetric multivariate Laplace distribution.

The asymmetric Laplace distribution is demonstrated in Kozubowski and
Podgorski (2001) to be well-suited for financial modeling, specifically
with currency exchange rates.

These functions are similar to those in the \code{VGAM} package.
}
\value{
  \code{dalaplace} gives the density,
  \code{palaplace} gives the distribution function,
  \code{qalaplace} gives the quantile function, and
  \code{ralaplace} generates random deviates.
}
\references{
  Kotz, S., Kozubowski, T.J., and Podgorski, K. (2001). "The Laplace
  Distribution and Generalizations: a Revisit with Applications to
  Communications, Economics, Engineering, and Finance". Boston:
  Birkhauser.

  Kozubowski, T.J. and Podgorski, K. (2001). "Asymmetric Laplace
  Laws and Modeling Financial Data". \emph{Mathematical and Computer
  Modelling}, 34, p. 1003-1021.
}
\seealso{\code{\link{dlaplace}} and
  \code{\link{dallaplace}}
}
\examples{
library(LaplacesDemon)
x <- dalaplace(1,0,1,1)
x <- palaplace(1,0,1,1)
x <- qalaplace(0.5,0,1,1)
x <- ralaplace(100,0,1,1)

#Plot Probability Functions
x <- seq(from=-5, to=5, by=0.1)
plot(x, dalaplace(x,0,1,0.5), ylim=c(0,1), type="l", main="Probability Function",
     ylab="density", col="red")
lines(x, dalaplace(x,0,1,1), type="l", col="green")
lines(x, dalaplace(x,0,1,5), type="l", col="blue")
legend(1, 0.9, expression(paste(mu==0, ", ", lambda==1, ", ", kappa==0.5),
     paste(mu==0, ", ", lambda==1, ", ", kappa==1),
     paste(mu==0, ", ", lambda==1, ", ", kappa==5)),
     lty=c(1,1,1), col=c("red","green","blue"))
}
\keyword{Distribution}