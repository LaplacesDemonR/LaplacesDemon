\name{dist.Inverse.Gaussian}
\alias{dinvgaussian}
\alias{rinvgaussian}
\title{Inverse Gaussian Distribution}
\description{
  This is the density function and random generation from the inverse
  gaussian distribution.
}
\usage{
dinvgaussian(x, mu, lambda, log=FALSE)
rinvgaussian(n, mu, lambda)
}
\arguments{
  \item{n}{This is the number of draws from the distribution.}
  \item{x}{This is the scalar location to evaluate density.}
  \item{mu}{This is the mean parameter, \eqn{\mu}{mu}.}
  \item{lambda}{This is the inverse-variance parameter,
    \eqn{\lambda}{lambda}.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
    \item Application: Continuous Univariate
    \item Density: \eqn{p(\theta) = \frac{\lambda}{(2 \pi
	\theta^3)^{1/2}} \exp(-\frac{\lambda (\theta - \mu)^2}{2 \mu^2
	\theta}), \theta > 0}{p(theta) = (lambda / (2*pi*theta^3))^(1/2) *
        exp(-((lambda*(theta-mu)^2) / (2*mu^2*theta))), theta > 0}
    \item Inventor: Schrodinger (1915)
    \item Notation 1: \eqn{\theta \sim \mathcal{N}^{-1}(\mu,
      \lambda)}{theta ~ N^-1(mu, lambda)}
    \item Notation 2: \eqn{p(\theta) = \mathcal{N}^{-1}(\theta | \mu,
      \lambda)}{p(theta) = N^-1(theta | mu, lambda)}
    \item Parameter 1: shape \eqn{\mu > 0}{mu > 0}
    \item Parameter 2: scale \eqn{\lambda > 0}{lambda > 0}
    \item Mean: \eqn{E(\theta) = \mu}{E(theta) = mu}
    \item Variance: \eqn{var(\theta) = \frac{\mu^3}{\lambda}}{var(theta)
      = mu^3/lambda}
    \item Mode: \eqn{mode(\theta) = \mu((1 + \frac{9 \mu^2}{4
      \lambda^2})^{1/2} - \frac{3 \mu}{2 \lambda})}{mode(theta) = mu*((1
    + ((9*mu^2)/(4*lambda^2)))^(1/2) - \frac{3*mu}{2*lambda})}
}

The inverse-Gaussian distribution, also called the Wald distribution, is
used when modeling dependent variables that are positive and
continuous. When
\eqn{\lambda \rightarrow \infty}{lambda tends to infinity} (or variance
to zero), the inverse-Gaussian distribution becomes similar to a normal
(Gaussian) distribution. The name, inverse-Gaussian, is misleading,
because it is not the inverse of a Gaussian distribution, which is
obvious from the fact that \eqn{\theta}{theta} must be positive.
}
\value{
  \code{dinvgaussian} gives the density and
  \code{rinvgaussian} generates random deviates.
}
\references{
  Schrodinger E. (1915). "Zur Theorie der Fall-und Steigversuche an
  Teilchenn mit Brownscher Bewegung". \emph{Physikalische Zeitschrift},
  16, p. 289--295.
}
\seealso{
  \code{\link{dnorm}},
  \code{\link{dnormp}}, and
  \code{\link{dnormv}}.
}
\examples{
library(LaplacesDemon)
x <- dinvgaussian(2, 1, 1)
x <- rinvgaussian(10, 1, 1)

#Plot Probability Functions
x <- seq(from=1, to=20, by=0.1)
plot(x, dinvgaussian(x,1,0.5), ylim=c(0,1), type="l", main="Probability Function",
     ylab="density", col="red")
lines(x, dinvgaussian(x,1,1), type="l", col="green")
lines(x, dinvgaussian(x,1,5), type="l", col="blue")
legend(2, 0.9, expression(paste(mu==1, ", ", sigma==0.5),
     paste(mu==1, ", ", sigma==1), paste(mu==1, ", ", sigma==5)),
     lty=c(1,1,1), col=c("red","green","blue"))
}
\keyword{Distribution}

