\name{dist.Student.t.Precision}
\alias{dstp}
\alias{pstp}
\alias{qstp}
\alias{rstp}
\title{Student t Distribution: Precision Parameterization}
\description{
  These functions provide the density, distribution function, quantile
  function, and random generation for the univariate Student t
  distribution with location parameter \eqn{\mu}{mu}, precision parameter
  \eqn{\tau}{tau}, and degrees of freedom parameter \eqn{\nu}{nu}.
}
\usage{
dstp(x, mu=0, tau=1, nu=10, log=FALSE)
pstp(q, mu=0, tau=1, nu=10, lower.tail=TRUE, log.p=FALSE)
qstp(p, mu=0, tau=1, nu=10, lower.tail=TRUE, log.p=FALSE)
rstp(n, mu=0, tau=1, nu=10)
}
\arguments{
  \item{x, q}{These are each a vector of quantiles.}
  \item{p}{This is a vector of probabilities.}
  \item{n}{This is the number of observations, which must be a positive
    integer that has length 1.}
  \item{mu}{This is the location parameter \eqn{\mu}{mu}.}
  \item{tau}{This is the precision parameter \eqn{\tau}{tau}, which
    must be positive.}
  \item{nu}{This is the degrees of freedom parameter \eqn{\nu}{nu},
    which must be positive.}
  \item{lower.tail}{Logical. If \code{lower.tail=TRUE}, then
    probabilities are \eqn{Pr[X \le x]}{Pr[X <= x]}, otherwise,
    \eqn{Pr[X > x]}{Pr[X > x]}.}
  \item{log, log.p}{Logical. If \code{log=TRUE}, then the logarithm of the
    density or probability is returned.}
}
\details{
  \itemize{
  \item Application: Continuous Univariate
  \item Density: \eqn{p(\theta) =
    \frac{\Gamma((\nu+1)/2)}{\Gamma(\nu/2)} \sqrt{\frac{\tau}{\nu\pi}} (1
    + \frac{\tau}{\nu} (\theta-\mu)^2)^{-(\nu+1)/2}}{p(theta) =
    (Gamma((nu+1)/2) / Gamma(nu/2)) * sqrt(tau/(nu*pi)) *
     (1 + (tau/nu)*(theta-mu)^2)^(-(nu+1)/2)}
  \item Inventor: William Sealy Gosset (1908)
  \item Notation 1: \eqn{\theta \sim \mathrm{t}(\mu, \sqrt{\tau^{-1}},
    \nu)}{theta ~ t(mu, sqrt(tau^(-1)),nu)}
  \item Notation 2: \eqn{p(\theta) = \mathrm{t}(\theta | \mu,
    \sqrt{\tau^{-1}}, \nu)}{p(theta) = t(theta | mu, sqrt(tau^(-1)),
    nu)}
  \item Parameter 1: location parameter \eqn{\mu}{mu}
  \item Parameter 2: precision parameter \eqn{\tau > 0}{tau > 0}
  \item Parameter 3: degrees of freedom \eqn{\nu > 0}{nu > 0}
  \item Mean: \eqn{E(\theta) = \mu}{E(theta) = mu}, for \eqn{\nu > 1}{nu
    > 1}, otherwise undefined
  \item Variance: \eqn{var(\theta) = \frac{1}{\tau}\frac{\nu}{\nu -
      2}}{var(theta) = (1/tau)[nu / (nu - 2)]}, for \eqn{\nu > 2}{nu >
    2}
  \item Mode: \eqn{mode(\theta) = \mu}{mode(theta) = mu}
}

The Student t-distribution is often used as an alternative to the normal
distribution as a model for data. It is frequently the case that real
data have heavier tails than the normal distribution allows for. The
classical approach was to identify outliers and exclude or downweight
them in some way. However, it is not always easy to identify outliers
(especially in high dimensions), and the Student t-distribution is a
natural choice of model-form for such data. It provides a parametric
approach to robust statistics.

The degrees of freedom parameter, \eqn{\nu}{nu}, controls the kurtosis
of the distribution, and is correlated with the precision parameter
\eqn{\tau}{tau}. The likelihood can have multiple local maxima
and, as such, it is often necessary to fix \eqn{\nu}{nu} at a fairly
low value and estimate the other parameters taking this as given.
Some authors report that values between 3 and 9 are often good
choices, and some authors suggest 5 is often a good choice.

In the limit \eqn{\nu \rightarrow \infty}{nu -> infinity}, the
Student t-distribution approaches
\eqn{\mathcal{N}(\mu, \sigma^2)}{N(mu, sigma^2)}. The
case of \eqn{\nu = 1}{nu = 1} is the Cauchy distribution.
}
\value{
  \code{dstp} gives the density,
  \code{pstp} gives the distribution function,
  \code{qstp} gives the quantile function, and
  \code{rstp} generates random deviates.
}
\seealso{
  \code{\link{dcauchy}},
  \code{\link{dmvt}},
  \code{\link{dmvtp}},
  \code{\link{dnorm}},
  \code{\link{dnormp}},
  \code{\link{dnormv}},
  \code{\link{dst}},
  \code{\link{dt}}.
}
\examples{
library(LaplacesDemon)
x <- dstp(1,0,1,10)
x <- pstp(1,0,1,10)
x <- qstp(0.5,0,1,10)
x <- rstp(100,0,1,10)

#Plot Probability Functions
x <- seq(from=-5, to=5, by=0.1)
plot(x, dstp(x,0,1,0.1), ylim=c(0,1), type="l", main="Probability Function",
     ylab="density", col="red")
lines(x, dstp(x,0,1,1), type="l", col="green")
lines(x, dstp(x,0,1,10), type="l", col="blue")
legend(1, 0.9, expression(paste(mu==0, ", ", tau==1, ", ", nu==0.5),
     paste(mu==0, ", ", tau==1, ", ", nu==1),
     paste(mu==0, ", ", tau==1, ", ", nu==10)),
     lty=c(1,1,1), col=c("red","green","blue"))
}
\keyword{Distribution}