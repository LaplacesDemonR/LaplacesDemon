\name{dist.Normal.Variance}
\alias{dnormv}
\alias{pnormv}
\alias{qnormv}
\alias{rnormv}
\title{Normal Distribution: Variance Parameterization}
\description{
  These functions provide the density, distribution function, quantile
  function, and random generation for the univariate normal distribution
  with mean \eqn{\mu}{mu} and variance \eqn{\sigma^2}{sigma^2}.
}
\usage{
dnormv(x, mean=0, var=1, log=FALSE)
pnormv(q, mean=0, var=1, lower.tail=TRUE, log.p=FALSE)
qnormv(p, mean=0, var=1, lower.tail=TRUE, log.p=FALSE)
rnormv(n, mean=0, var=1)
}
\arguments{
  \item{x, q}{These are each a vector of quantiles.}
  \item{p}{This is a vector of probabilities.}
  \item{n}{This is the number of observations, which must be a positive
    integer that has length 1.}
  \item{mean}{This is the mean parameter \eqn{\mu}{mu}.}
  \item{var}{This is the variance parameter \eqn{\sigma^2}{sigma^2}, which
    must be positive.}
  \item{log, log.p}{Logical. If \code{TRUE}, then probabilities \eqn{p}
    are given as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{Logical. If \code{TRUE} (default), then probabilities
    are \eqn{Pr[X \le x]}{Pr[X <= x]}, otherwise,
    \eqn{Pr[X > x]}{Pr[X > x]}.}
}
\details{
  \itemize{
  \item Application: Continuous Univariate
  \item Density: \eqn{p(\theta) = \frac{1}{\sqrt{2\pi\sigma^2}}
    \exp(-\frac{(\theta-\mu)^2}{2\sigma^2})}{p(theta) = (1/(sqrt(2*pi*sigma^2))) * exp(-((theta-mu)^2/(2*sigma^2)))}
  \item Inventor: Carl Friedrich Gauss or Abraham De Moivre
  \item Notation 1: \eqn{\theta \sim \mathcal{N}(\mu, \sigma^2)}{theta ~ N(mu, sigma^2)}
  \item Notation 2: \eqn{p(\theta) = \mathcal{N}(\theta | \mu,
    \sigma^2)}{p(theta) = N(theta | mu, sigma^2)}
  \item Parameter 1: mean parameter \eqn{\mu}{mu}
  \item Parameter 2: variance parameter \eqn{\sigma^2 > 0}{sigma^2 > 0}
  \item Mean: \eqn{E(\theta) = \mu}{E(theta) = mu}
  \item Variance: \eqn{var(\theta) = \sigma^2}{var(theta) = sigma^2}
  \item Mode: \eqn{mode(\theta) = \mu}{mode(theta) = mu}
}

The normal distribution, also called the Gaussian distribution and the
Second Law of Laplace, is usually parameterized with mean and variance.
\code{Base R} uses the mean and standard deviation. These functions
provide the variance parameterization for convenience and familiarity.
For example, it is easier to code \code{dnormv(1,1,1000)} than
\code{dnorm(1,1,sqrt(1000))}.

Some authors attribute credit for the normal distribution to Abraham
de Moivre in 1738. In 1809, Carl Friedrich Gauss published his
monograph ``Theoria motus corporum coelestium in sectionibus conicis
solem ambientium'', in which he introduced the method of least squares,
method of maximum likelihood, and normal distribution, among many other
innovations.

Gauss, himself, characterized this distribution according to mean and
precision, though his definition of precision differed from the modern
one.

Although the normal distribution is very common, it often does not fit
data as well as more robust alternatives with fatter tails, such as the
Laplace or Student t distribution.

A flat distribution is obtained in the limit as
\eqn{\sigma^2 \rightarrow \infty}{sigma^2 -> infinity}.

For models where the dependent variable, y, is specified to be
normally distributed given the model, the Jarque-Bera test (see
\code{\link{plot.demonoid.ppc}} or \code{\link{plot.laplace.ppc}}) may
be used to test the residuals.

These functions are similar to those in \code{base R}.
}
\value{
  \code{dnormv} gives the density,
  \code{pnormv} gives the distribution function,
  \code{qnormv} gives the quantile function, and
  \code{rnormv} generates random deviates.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{dlaplace}},
  \code{\link{dnorm}},
  \code{\link{dnormp}},
  \code{\link{dst}},
  \code{\link{dt}},
  \code{\link{plot.demonoid.ppc}}, and
  \code{\link{plot.laplace.ppc}}.
}
\examples{
library(LaplacesDemon)
x <- dnormv(1,0,1)
x <- pnormv(1,0,1)
x <- qnormv(0.5,0,1)
x <- rnormv(100,0,1)

#Plot Probability Functions
x <- seq(from=-5, to=5, by=0.1)
plot(x, dnormv(x,0,0.5), ylim=c(0,1), type="l", main="Probability Function",
     ylab="density", col="red")
lines(x, dnormv(x,0,1), type="l", col="green")
lines(x, dnormv(x,0,5), type="l", col="blue")
legend(2, 0.9, expression(paste(mu==0, ", ", sigma^2==0.5),
     paste(mu==0, ", ", sigma^2==1), paste(mu==0, ", ", sigma^2==5)),
     lty=c(1,1,1), col=c("red","green","blue"))
}
\keyword{Distribution}