\name{dist.Power.Exponential}
\alias{dpe}
\alias{ppe}
\alias{qpe}
\alias{rpe}
\title{Power Exponential Distribution: Univariate Symmetric}
\description{
  These functions provide the density, distribution function, quantile
  function, and random generation for the univariate, symmetric, power
  exponential distribution with location parameter \eqn{\mu}{mu}, scale
  parameter \eqn{\sigma}{sigma}, and kurtosis parameter
  \eqn{\kappa}{kappa}.
}
\usage{
dpe(x, mu=0, sigma=1, kappa=2, log=FALSE)
ppe(q, mu=0, sigma=1, kappa=2, lower.tail=TRUE, log.p=FALSE)
qpe(p, mu=0, sigma=1, kappa=2, lower.tail=TRUE, log.p=FALSE)
rpe(n, mu=0, sigma=1, kappa=2)
}
\arguments{
  \item{x, q}{These are each a vector of quantiles.}
  \item{p}{This is a vector of probabilities.}
  \item{n}{This is the number of observations, which must be a positive
    integer that has length 1.}
  \item{mu}{This is the location parameter \eqn{\mu}{mu}.}
  \item{sigma}{This is the scale parameter \eqn{\sigma}{sigma}, which
    must be positive.}
  \item{kappa}{This is the kurtosis parameter \eqn{\kappa}{kappa}, which
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
  \item Density: \eqn{p(\theta) = \frac{1}{2 \kappa^{1/\kappa}
    \Gamma(1+\frac{1}{\kappa}) \sigma}
  \exp(-\frac{|\theta-\mu|^{\kappa}}{\kappa \sigma^\kappa})}{p(theta) = 1/(2
  kappa^(1/kappa) Gamma(1+1/kappa) sigma) * exp(-(|theta -
  mu|^kappa)/(kappa sigma^kappa))}
  \item Inventor: Subbotin, M.T. (1923)
  \item Notation 1: \eqn{\theta \sim \mathcal{PE}(\mu, \sigma, \kappa)}{theta ~
    PE(mu, sigma, kappa)}
  \item Notation 2: \eqn{p(\theta) = \mathcal{PE}(\theta | \mu,
    \sigma, \kappa)}{p(theta) = PE(theta | mu, sigma, kappa)}
  \item Parameter 1: location parameter \eqn{\mu}{mu}
  \item Parameter 2: scale parameter \eqn{\sigma > 0}{sigma > 0}
  \item Parameter 3: kurtosis parameter \eqn{\kappa > 0}{kappa > 0}
  \item Mean: \eqn{E(\theta) = \mu}{E(theta) = mu}
  \item Variance: \eqn{var(\theta) = }{var(theta) = }
  \item Mode: \eqn{mode(\theta) = \mu}{mode(theta) = mu}
}

The power exponential distribution is also called the exponential power
distribution, generalized error distribution, generalized Gaussian
distribution, and generalized normal distribution. The original form was
introduced by Subbotin (1923) and re-parameterized by Lunetta
(1963). These functions use the more recent parameterization by Lunetta
(1963). A shape parameter, \eqn{\kappa > 0}{kappa > 0}, is added to the
normal distribution. When \eqn{\kappa=1}{kappa=1}, the power exponential
distribution is the same as the Laplace distribution. When
\eqn{\kappa=2}{kappa=2}, the power exponential distribution is the same
as the normal distribution. As \eqn{\kappa \rightarrow \infty}{kappa ->
  infinity}, this becomes a uniform distribution \eqn{\in (\mu-\sigma,
  \mu+\sigma)}{in (mu-sigma, mu+sigma)}. Tails that are heavier than
normal occur when \eqn{\kappa < 2}{kappa < 2}, or lighter than normal
when \eqn{\kappa > 2}{kappa > 2}. This distribution is univariate and
symmetric, and there exist multivariate and asymmetric versions.

These functions are similar to those in the \code{normalp} package.
}
\value{
  \code{dpe} gives the density,
  \code{ppe} gives the distribution function,
  \code{qpe} gives the quantile function, and
  \code{rpe} generates random deviates.
}
\references{
  Lunetta, G. (1963). "Di una Generalizzazione dello Schema della Curva
  Normale". \emph{Annali della Facolt`a di Economia e Commercio di
  Palermo}, 17, p. 237--244.

  Subbotin, M.T. (1923). "On the Law of Frequency of Errors".
  \emph{Matematicheskii Sbornik}, 31, p. 296--301.
}
\seealso{
  \code{\link{dlaplace}},
  \code{\link{dlaplacep}},
  \code{\link{dmvpe}},
  \code{\link{dnorm}},
  \code{\link{dnormp}},
  \code{\link{dnormv}}, and
  \code{\link{dunif}}.
}
\examples{
library(LaplacesDemon)
x <- dpe(1,0,1,2)
x <- ppe(1,0,1,2)
x <- qpe(0.5,0,1,2)
x <- rpe(100,0,1,2)

#Plot Probability Functions
x <- seq(from=0.1, to=3, by=0.01)
plot(x, dpe(x,0,1,0.1), ylim=c(0,1), type="l", main="Probability Function",
     ylab="density", col="red")
lines(x, dpe(x,0,1,2), type="l", col="green")
lines(x, dpe(x,0,1,5), type="l", col="blue")
legend(1.5, 0.9, expression(paste(mu==0, ", ", sigma==1, ", ", kappa==0.1),
     paste(mu==0, ", ", sigma==1, ", ", kappa==2),
     paste(mu==0, ", ", sigma==1, ", ", kappa==5)),
     lty=c(1,1,1), col=c("red","green","blue"))
}
\keyword{Distribution}