\name{dist.Inverse.ChiSquare}
\alias{dinvchisq}
\alias{rinvchisq}
\title{(Scaled) Inverse Chi-Squared Distribution}
\description{
  This is the density function and random generation for the (scaled)
  inverse chi-squared distribution.
}
\usage{
dinvchisq(x, df, scale, log=FALSE)
rinvchisq(n, df, scale=1/df)
}
\arguments{
  \item{x}{This is a vector of quantiles.}
  \item{n}{This is the number of observations. If \code{length(n) > 1},
    then the length is taken to be the number required.}
  \item{df}{This is the degrees of freedom parameter, usually
    represented as \eqn{\nu}{nu}.}
  \item{scale}{This is the scale parameter, usually represented as
    \eqn{\lambda}{lambda}.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
  \item Application: Continuous Univariate
  \item Density: \deqn{p(\theta) = \frac{{\nu/2}^{\nu/2}}{\Gamma(\nu/2)}
    \lambda^\nu \frac{1}{\theta}^{\nu/2+1} \exp(-\frac{\nu
      \lambda^2}{2\theta}), \theta \ge 0}{p(theta) =
    ((nu/2)^(nu/2))/(\Gamma(nu/2)) lambda^nu (1/theta)^((nu/2)+1)
    exp(-(nu lambda^2)/(2*theta)), theta >= 0}
  \item Inventor: Derived from the chi-squared distribution
  \item Notation 1: \eqn{\theta \sim \chi^{-2}(\nu, \lambda)}{theta ~
    chi^(-2)(nu, lambda)}
  \item Notation 2: \eqn{p(\theta) = \chi^{-2}(\theta | \nu,
    \lambda)}{p(theta) = chi^(-2)(theta | nu, lambda)}
  \item Parameter 1: degrees of freedom parameter \eqn{\nu > 0}{nu > 0}
  \item Parameter 2: scale parameter \eqn{\lambda}{lambda}
  \item Mean: \eqn{E(\theta)}{E(theta)} = unknown
  \item Variance: \eqn{var(\theta)}{var(theta)} = unknown
  \item Mode: \eqn{mode(\theta) = }{mode(theta) = }
}

The inverse chi-squared distribution, also called the
inverted chi-square distribution, is the multiplicate inverse of the
chi-squared distribution. If \eqn{x} has the chi-squared distribution
with \eqn{\nu}{nu} degrees of freedom, then \eqn{1 / x} has the
inverse chi-squared distribution with \eqn{\nu}{nu} degrees of freedom,
and \eqn{\nu / x}{nu / x} has the inverse chi-squared distribution with
\eqn{\nu}{nu} degrees of freedom.

These functions are similar to those in the GeoR package.
}
\value{
  \code{dinvchisq} gives the density and
  \code{rinvchisq} generates random deviates.
}
\seealso{
  \code{\link{dchisq}}
}
\examples{
library(LaplacesDemon)
x <- dinvchisq(1,1,1)
x <- rinvchisq(10,1)

#Plot Probability Functions
x <- seq(from=0.1, to=5, by=0.01)
plot(x, dinvchisq(x,0.5,1), ylim=c(0,1), type="l", main="Probability Function",
     ylab="density", col="red")
lines(x, dinvchisq(x,1,1), type="l", col="green")
lines(x, dinvchisq(x,5,1), type="l", col="blue")
legend(3, 0.9, expression(paste(nu==0.5, ", ", lambda==1),
     paste(nu==1, ", ", lambda==1), paste(nu==5, ", ", lambda==1)),
     lty=c(1,1,1), col=c("red","green","blue"))
}
\keyword{Distribution}
