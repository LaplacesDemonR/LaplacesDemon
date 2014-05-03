\name{dist.HalfCauchy}
\alias{dhalfcauchy}
\alias{phalfcauchy}
\alias{qhalfcauchy}
\alias{rhalfcauchy}
\title{Half-Cauchy Distribution}
\description{
  These functions provide the density, distribution function, quantile
  function, and random generation for the half-Cauchy distribution.
}
\usage{
dhalfcauchy(x, scale=25, log=FALSE)
phalfcauchy(q, scale=25)
qhalfcauchy(p, scale=25)
rhalfcauchy(n, scale=25)
}
\arguments{
  \item{x,q}{These are each a vector of quantiles.}
  \item{p}{This is a vector of probabilities.}
  \item{n}{This is the number of observations, which must be a positive
    integer that has length 1.}
  \item{scale}{This is the scale parameter \eqn{\alpha}{alpha}, which
    must be positive.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
  \item Application: Continuous Univariate
  \item Density: \eqn{p(\theta) = \frac{2 \alpha}{\pi(\theta^2 +
      \alpha^2)}, \quad \theta > 0}{p(theta) = 2alpha / pi(theta^2 + alpha^2), theta >= 0}
  \item Inventor: Derived from Cauchy
  \item Notation 1: \eqn{\theta \sim \mathcal{HC}(\alpha)}{theta ~ HC(alpha)}
  \item Notation 2: \eqn{p(\theta) = \mathcal{HC}(\theta | \alpha)}{p(theta) = HC(theta | alpha)}
  \item Parameter 1: scale parameter \eqn{\alpha > 0}{alpha > 0}
  \item Mean: \eqn{E(\theta)}{E(theta)} = does not exist
  \item Variance: \eqn{var(\theta)}{var(theta)} = does not exist
  \item Mode: \eqn{mode(\theta) = 0}{mode(theta) = 0}
}

The half-Cauchy distribution with scale \eqn{\alpha=25}{alpha=25} is a
recommended, default, weakly informative prior distribution for a scale
parameter. Otherwise, the scale, \eqn{\alpha}{alpha}, is recommended to
be set to be just a little larger than the expected standard deviation,
as a weakly informative prior distribution on a standard deviation
parameter.

The Cauchy distribution is known as a pathological distribution because
its mean and variance are undefined, and it does not satisfy the central
limit theorem.
}
\value{
  \code{dhalfcauchy} gives the density,
  \code{phalfcauchy} gives the distribution function,
  \code{qhalfcauchy} gives the quantile function, and
  \code{rhalfcauchy} generates random deviates.
}
\seealso{
  \code{\link{dcauchy}}
}
\examples{
library(LaplacesDemon)
x <- dhalfcauchy(1,25)
x <- phalfcauchy(1,25)
x <- qhalfcauchy(0.5,25)
x <- rhalfcauchy(1,25)

#Plot Probability Functions
x <- seq(from=0, to=20, by=0.1)
plot(x, dhalfcauchy(x,1), ylim=c(0,1), type="l", main="Probability Function",
     ylab="density", col="red")
lines(x, dhalfcauchy(x,5), type="l", col="green")
lines(x, dhalfcauchy(x,10), type="l", col="blue")
legend(2, 0.9, expression(alpha==1, alpha==5, alpha==10),
     lty=c(1,1,1), col=c("red","green","blue"))
}
\keyword{Distribution}