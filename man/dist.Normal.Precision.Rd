\name{dist.Normal.Precision}
\alias{dnormp}
\alias{pnormp}
\alias{qnormp}
\alias{rnormp}
\title{Normal Distribution: Precision Parameterization}
\description{
  These functions provide the density, distribution function, quantile
  function, and random generation for the univariate normal distribution
  with mean \eqn{\mu}{mu} and precision \eqn{\tau}{tau}.
}
\usage{
dnormp(x, mean=0, prec=1, log=FALSE)
pnormp(q, mean=0, prec=1, lower.tail=TRUE, log.p=FALSE)
qnormp(p, mean=0, prec=1, lower.tail=TRUE, log.p=FALSE)
rnormp(n, mean=0, prec=1)
}
\arguments{
  \item{x, q}{These are each a vector of quantiles.}
  \item{p}{This is a vector of probabilities.}
  \item{n}{This is the number of observations, which must be a positive
    integer that has length 1.}
  \item{mean}{This is the mean parameter \eqn{\mu}{mu}.}
  \item{prec}{This is the precision parameter \eqn{\tau}{tau}, which
    must be positive.}
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
    \exp(-\frac{\tau}{2} (\theta-\mu)^2)}{p(theta) = sqrt(tau/(2*pi)) * exp(-(tau/2)*(x-mu)^2)}
  \item Inventor: Carl Friedrich Gauss or Abraham De Moivre
  \item Notation 1: \eqn{\theta \sim \mathcal{N}(\mu, \tau^{-1})}{theta ~ N(mu, tau^(-1))}
  \item Notation 2: \eqn{p(\theta) = \mathcal{N}(\theta | \mu,
    \tau^{-1})}{p(theta) = N(theta | mu, tau^(-1))}
  \item Parameter 1: mean parameter \eqn{\mu}{mu}
  \item Parameter 2: precision parameter \eqn{\tau > 0}{tau > 0}
  \item Mean: \eqn{E(\theta) = \mu}{E(theta) = mu}
  \item Variance: \eqn{var(\theta) = \tau^{-1}}{var(theta) = tau^(-1)}
  \item Mode: \eqn{mode(\theta) = \mu}{mode(theta) = mu}
}

The normal distribution, also called the Gaussian distribution and the
Second Law of Laplace, is usually parameterized with mean and variance,
or in Bayesian inference, with mean and precision, where precision is
the inverse of the variance. In contrast, \code{Base R} parameterizes
the normal distribution with the mean and standard deviation. These
functions provide the precision parameterization for convenience and
familiarity.

Some authors attribute credit for the normal distribution to Abraham
de Moivre in 1738. In 1809, Carl Friedrich Gauss published his
monograph ``Theoria motus corporum coelestium in sectionibus conicis
solem ambientium'', in which he introduced the method of least squares,
method of maximum likelihood, and normal distribution, among many other
innovations.

Gauss, himself, characterized this distribution according to mean and
precision, though his definition of precision differed from the modern
one. The modern Bayesian use of precision \eqn{\tau}{tau} developed
because it was more straightforward to estimate \eqn{\tau}{tau} with a
gamma distribution as a conjugate prior, than to estimate
\eqn{\sigma^2}{sigma^2} with an inverse-gamma distribution as a
conjugate prior.

Although the normal distribution is very common, it often does not fit
data as well as more robust alternatives with fatter tails, such as the
Laplace or Student t distribution.

A flat distribution is obtained in the limit as
\eqn{\tau \rightarrow 0}{tau -> 0}.

For models where the dependent variable, y, is specified to be
normally distributed given the model, the Jarque-Bera test (see
\code{\link{plot.demonoid.ppc}} or \code{\link{plot.laplace.ppc}}) may
be used to test the residuals.

These functions are similar to those in \code{base R}.
}
\value{
  \code{dnormp} gives the density,
  \code{pnormp} gives the distribution function,
  \code{qnormp} gives the quantile function, and
  \code{rnormp} generates random deviates.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{dlaplace}},
  \code{\link{dnorm}},
  \code{\link{dnormv}},
  \code{\link{prec2var}},
  \code{\link{dst}},
  \code{\link{dt}},
  \code{\link{plot.demonoid.ppc}}, and
  \code{\link{plot.laplace.ppc}}.
}
\examples{
library(LaplacesDemon)
x <- dnormp(1,0,1)
x <- pnormp(1,0,1)
x <- qnormp(0.5,0,1)
x <- rnormp(100,0,1)

#Plot Probability Functions
x <- seq(from=-5, to=5, by=0.1)
plot(x, dnormp(x,0,0.5), ylim=c(0,1), type="l", main="Probability Function",
     ylab="density", col="red")
lines(x, dnormp(x,0,1), type="l", col="green")
lines(x, dnormp(x,0,5), type="l", col="blue")
legend(2, 0.9, expression(paste(mu==0, ", ", tau==0.5),
     paste(mu==0, ", ", tau==1), paste(mu==0, ", ", tau==5)),
     lty=c(1,1,1), col=c("red","green","blue"))
}
\keyword{Distribution}