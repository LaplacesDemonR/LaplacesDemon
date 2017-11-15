\name{dist.Log.Normal.Precision}
\alias{dlnormp}
\alias{plnormp}
\alias{qlnormp}
\alias{rlnormp}
\title{Log-Normal Distribution: Precision Parameterization}
\description{
  These functions provide the density, distribution function, quantile
  function, and random generation for the univariate log-normal
  distribution with mean \eqn{\mu}{mu} and precision \eqn{\tau}{tau}.
}
\usage{
dlnormp(x, mu, tau=NULL, var=NULL, log=FALSE)
plnormp(q, mu, tau, lower.tail=TRUE, log.p=FALSE)
qlnormp(p, mu, tau, lower.tail=TRUE, log.p=FALSE)
rlnormp(n, mu, tau=NULL, var=NULL)
}
\arguments{
  \item{x, q}{These are each a vector of quantiles.}
  \item{p}{This is a vector of probabilities.}
  \item{n}{This is the number of observations, which must be a positive
    integer that has length 1.}
  \item{mu}{This is the mean parameter \eqn{\mu}{mu}.}
  \item{tau}{This is the precision parameter \eqn{\tau}{tau}, which
    must be positive. Tau and var cannot be used together}
  \item{var}{This is the variance parameter, which must be positive. Tau and var cannot be used together}
  \item{log, log.p}{Logical. If \code{TRUE}, then probabilities
    \eqn{p} are given as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{Logical. If \code{TRUE} (default), then probabilities
    are \eqn{Pr[X \le x]}{Pr[X <= x]}, otherwise,
    \eqn{Pr[X > x]}{Pr[X > x]}.}
}
\details{
  \itemize{
  \item Application: Continuous Univariate
  \item Density: \eqn{p(\theta) = \sqrt{\frac{\tau}{2\pi}}
    \frac{1}{\theta} \exp(-\frac{\tau}{2} (\log(\theta -
    \mu))^2)}{p(theta) = sqrt(tau/(2*pi)) * (1/theta) *
    exp(-(tau/2)*(log(theta-mu))^2)}
  \item Inventor: Carl Friedrich Gauss or Abraham De Moivre
  \item Notation 1: \eqn{\theta \sim \mathrm{Log-}\mathcal{N}(\mu,
    \tau^{-1})}{theta ~ Log-N(mu, tau^(-1))}
  \item Notation 2: \eqn{p(\theta) = \mathrm{Log-}\mathcal{N}(\theta | \mu,
    \tau^{-1})}{p(theta) = Log-N(theta | mu, tau^(-1))}
  \item Parameter 1: mean parameter \eqn{\mu}{mu}
  \item Parameter 2: precision parameter \eqn{\tau > 0}{tau > 0}
  \item Mean: \eqn{E(\theta) = \exp(\mu + \tau^{-1} / 2)}{E(theta) =
    exp(mu + tau^(-1) / 2)}
  \item Variance: \eqn{var(\theta) = (\exp(\tau^{-1}) - 1)\exp(2\mu +
    \tau^{-1})}{var(theta) = exp(tau^(-1) - 1) * exp(2*mu + tau^(-1))}
  \item Mode: \eqn{mode(\theta) = \exp(\mu - \tau^{-1})}{mode(theta) =
    exp(mu - tau^(-1))}
}

The log-normal distribution, also called the Galton distribution, is
applied to a variable whose logarithm is normally-distributed. The
distribution is usually parameterized with mean and variance, or in
Bayesian inference, with mean and precision, where precision is the
inverse of the variance. In contrast, \code{Base R} parameterizes the
log-normal distribution with the mean and standard deviation. These
functions provide the precision parameterization for convenience and
familiarity.

A flat distribution is obtained in the limit as
\eqn{\tau \rightarrow 0}{tau -> 0}.

These functions are similar to those in \code{base R}.
}
\value{
  \code{dlnormp} gives the density,
  \code{plnormp} gives the distribution function,
  \code{qlnormp} gives the quantile function, and
  \code{rlnormp} generates random deviates.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{dnorm}},
  \code{\link{dnormp}},
  \code{\link{dnormv}}, and
  \code{\link{prec2var}}.
}
\examples{
library(LaplacesDemon)
x <- dlnormp(1,0,1)
x <- plnormp(1,0,1)
x <- qlnormp(0.5,0,1)
x <- rlnormp(100,0,1)

#Plot Probability Functions
x <- seq(from=0.1, to=3, by=0.01)
plot(x, dlnormp(x,0,0.1), ylim=c(0,1), type="l", main="Probability Function",
     ylab="density", col="red")
lines(x, dlnormp(x,0,1), type="l", col="green")
lines(x, dlnormp(x,0,5), type="l", col="blue")
legend(2, 0.9, expression(paste(mu==0, ", ", tau==0.1),
     paste(mu==0, ", ", tau==1), paste(mu==0, ", ", tau==5)),
     lty=c(1,1,1), col=c("red","green","blue"))
}
\keyword{Distribution}