\name{dist.Generalized.Pareto}
\alias{dgpd}
\alias{rgpd}
\title{Generalized Pareto Distribution}
\description{
  These are the density and random generation functions for the
  generalized Pareto distribution.
}
\usage{
dgpd(x, mu, sigma, xi, log=FALSE)
rgpd(n, mu, sigma, xi)
}
\arguments{
  \item{x}{This is a vector of data.}
  \item{n}{This is a positive scalar integer, and is the number of
    observations to generate randomly.}
  \item{mu}{This is a scalar or vector location parameter
    \eqn{\mu}{mu}. When \eqn{\xi}{xi} is non-negative,
    \eqn{\mu}{mu} must not be greater than \eqn{\textbf{x}}{x}. When
    \eqn{\xi}{xi} is negative, \eqn{\mu}{mu} must be less than
    \eqn{\textbf{x} + \sigma / \xi}{x + sigma/xi}.}
  \item{sigma}{This is a positive-only scalar or vector of scale
    parameters \eqn{\sigma}{sigma}.}
  \item{xi}{This is a scalar or vector of shape parameters
    \eqn{\xi}{xi}.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
    \item Application: Continuous Univariate
    \item Density: \eqn{p(\theta) = \frac{1}{\sigma}(1 +
      \xi\textbf{z})^(-1/\xi + 1)}{1/sigma (1 + xi z)^(-1/xi +
      1)} where \eqn{\textbf{z} = \frac{\theta - \mu}{\sigma}}{z =
      (theta - mu)/sigma}
    \item Inventor: Pickands (1975)
    \item Notation 1: \eqn{\theta \sim \mathcal{GPD}(\mu, \sigma,
      \xi)}{theta ~ GPD(mu, sigma, xi)}
    \item Notation 2: \eqn{p(\theta) \sim \mathcal{GPD}(\theta |
      \mu, \sigma, \xi)}{p(theta) ~ GPD(theta | mu, sigma, xi)}
    \item Parameter 1: location \eqn{\mu}{mu}, where \eqn{\mu \le
      \theta}{mu <= theta} when \eqn{\xi \ge 0}{xi >= 0}, and
      \eqn{\mu \ge \theta + \sigma / \xi}{mu >= theta + sigma/xi}
      when \eqn{\xi < 0}{xi < 0}
    \item Parameter 2: scale \eqn{\sigma > 0}{sigma > 0}
    \item Parameter 3: shape \eqn{\xi}{xi}
    \item Mean: \eqn{\mu + \frac{\sigma}{1 - \xi}}{mu + sigma / (1 -
      xi)} when \eqn{\xi < 1}{xi < 1}
    \item Variance: \eqn{\frac{\sigma^2}{(1 - \xi)^2 (1 -
      2\xi)}}{sigma^2 / [(1 - xi)^2 (1 - 2 xi)]} when \eqn{\xi
      < 0.5}{xi < 0.5}
    \item Mode: 
}

The generalized Pareto distribution (GPD) is a more flexible extension
of the Pareto (\code{\link{dpareto}}) distribution. It is equivalent to
the exponential distribution when both \eqn{\mu = 0}{mu = 0} and
\eqn{\xi = 0}{xi = 0}, and it is equivalent to the Pareto
distribution when \eqn{\mu = \sigma / \xi}{mu = sigma / xi} and
\eqn{\xi > 0}{xi > 0}.

The GPD is often used to model the tails of another distribution, and
the shape parameter \eqn{\xi}{xi} relates to
tail-behavior. Distributions with tails that decrease exponentially are
modeled with shape \eqn{\xi = 0}{xi = 0}. Distributions with tails that
decrease as a polynomial are modeled with a positive shape
parameter. Distributions with finite tails are modeled with a negative
shape parameter.
}
\value{
  \code{dgpd} gives the density, and
  \code{rgpd} generates random deviates.
}
\references{
  Pickands J. (1975). "Statistical Inference Using Extreme Order
  Statistics". \emph{The Annals of Statistics}, 3, p. 119--131.
}
\seealso{
  \code{\link{dpareto}}
}
\examples{
library(LaplacesDemon)
x <- dgpd(0,0,1,0,log=TRUE)
x <- rgpd(10,0,1,0)
}
\keyword{Distribution}
