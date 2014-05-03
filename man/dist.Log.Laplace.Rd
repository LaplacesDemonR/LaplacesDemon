\name{dist.Log.Laplace}
\alias{dllaplace}
\alias{pllaplace}
\alias{qllaplace}
\alias{rllaplace}
\title{Log-Laplace Distribution: Univariate Symmetric}
\description{
  These functions provide the density, distribution function, quantile
  function, and random generation for the univariate, symmetric,
  log-Laplace distribution with location parameter \code{location} and
  scale parameter \code{scale}.
}
\usage{
dllaplace(x, location=0, scale=1, log=FALSE)
pllaplace(q, location=0, scale=1)
qllaplace(p, location=0, scale=1)
rllaplace(n, location=0, scale=1)
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
  \item Density 1: \eqn{p(\theta) =
    \frac{(\sqrt{2}/\lambda)^2}{2(\sqrt{2}/\lambda)}
    \exp(-(\sqrt{2}/\lambda)(\theta - \mu)), \theta \ge
    \exp(\mu)}{p(theta) = ((sqrt(2)/lambda)^2 / (2*(sqrt(2)/lambda))) *
    exp(-(sqrt(2)/lambda)*(theta - mu)), theta >= exp(mu)}
  \item Density 2: \eqn{p(\theta) =
    \frac{(\sqrt{2}/\lambda)^2}{2(\sqrt{2}/\lambda)}
    \exp((\sqrt{2}/\lambda)(\theta - \mu)), \theta < \exp(\mu)}{p(theta)
    = ((sqrt(2)/lambda)^2 / (2*(sqrt(2)/lambda))) *
    exp((sqrt(2)/lambda)*(theta - mu)),theta < exp(mu)}
  \item Inventor: Pierre-Simon Laplace
  \item Notation 1: \eqn{\theta \sim \mathcal{LL}(\mu, \lambda)}{theta ~
    LL(mu, lambda)}
  \item Notation 2: \eqn{p(\theta) = \mathcal{LL}(\theta | \mu,
    \lambda)}{p(theta) = LL(theta | mu, lambda)}
  \item Parameter 1: location parameter \eqn{\mu}{mu}
  \item Parameter 2: scale parameter \eqn{\lambda > 0}{lambda > 0}
  \item Mean: \eqn{E(\theta) = }{E(theta) = }
  \item Variance: \eqn{var(\theta) = }{var(theta) = }
  \item Mode: \eqn{mode(\theta) = }{mode(theta) = }
}

The univariate, symmetric log-Laplace distribution is derived from the
Laplace distribution. Multivariate and asymmetric versions also exist.

These functions are similar to those in the \code{VGAM} package.
}
\value{
  \code{dllaplace} gives the density,
  \code{pllaplace} gives the distribution function,
  \code{qllaplace} gives the quantile function, and
  \code{rllaplace} generates random deviates.
}
\references{
  Kozubowski, T. J. and Podgorski, K. (2003). "Log-Laplace Distributions".
  \emph{International Mathematical Journal}, 3, p. 467--495.
}
\seealso{
  \code{\link{dalaplace}},
  \code{\link{dallaplace}},
  \code{\link{dexp}},
  \code{\link{dlaplace}},
  \code{\link{dlaplacep}},
  \code{\link{dmvl}},
  \code{\link{dnorm}},
  \code{\link{dnormp}}, and
  \code{\link{dnormv}}.
}
\examples{
library(LaplacesDemon)
x <- dllaplace(1,0,1)
x <- pllaplace(1,0,1)
x <- qllaplace(0.5,0,1)
x <- rllaplace(100,0,1)

#Plot Probability Functions
x <- seq(from=0.1, to=20, by=0.1)
plot(x, dllaplace(x,0,0.1), ylim=c(0,1), type="l", main="Probability Function",
     ylab="density", col="red")
lines(x, dllaplace(x,0,0.5), type="l", col="green")
lines(x, dllaplace(x,0,1.5), type="l", col="blue")
legend(2, 0.9, expression(paste(mu==0, ", ", lambda==0.1),
     paste(mu==0, ", ", lambda==0.5), paste(mu==0, ", ", lambda==1.5)),
     lty=c(1,1,1), col=c("red","green","blue"))
}
\keyword{Distribution}