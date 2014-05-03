\name{dist.Bernoulli}
\alias{dbern}
\alias{pbern}
\alias{qbern}
\alias{rbern}
\title{Bernoulli Distribution}
\description{
  These functions provide the density, distribution function, quantile
  function, and random generation for the Bernoulli distribution.
}
\usage{
dbern(x, prob, log=FALSE)
pbern(q, prob, lower.tail=TRUE, log.p=FALSE)
qbern(p, prob, lower.tail=TRUE, log.p=FALSE)
rbern(n, prob)
}
\arguments{
  \item{x, q}{These are each a vector of quantiles.}
  \item{p}{This is a vector of probabilities.}
  \item{n}{This is the number of observations. If \code{length(n) > 1},
    then the length is taken to be the number required.}
  \item{prob}{This is the probability of success on each trial.}
  \item{log, log.p}{Logical. if \code{TRUE}, probabilities \eqn{p} are
    given as \eqn{\log(p)}{log(p)}.}
  \item{lower.tail}{Logical. if \code{TRUE} (default), probabilities
    are \eqn{Pr[X \le x]}{Pr[X <= x]}, otherwise,
    \eqn{Pr[X > x]}{Pr[X > x]}.}
}
\details{
  \itemize{
  \item Application: Continuous Univariate
  \item Density: \eqn{p(\theta) = {p}^{\theta}
    {(1-p)}^{1-\theta}}{p(theta) = p^theta (1-p)^(1-theta)}, \eqn{\theta = 0,1}{theta = 0,1}
  \item Inventor: Jacob Bernoulli
  \item Notation 1: \eqn{\theta \sim \mathcal{BERN}(p)}{theta ~ Bern(p)}
  \item Notation 2: \eqn{p(\theta) = \mathcal{BERN}(\theta | p)}{p(theta) = Bern(theta | p)}
  \item Parameter 1: probability parameter \eqn{0 \le p \le 1}{0 <= p <=
  1}
  \item Mean: \eqn{E(\theta) = p}{E(theta) = p}
  \item Variance: \eqn{var(\theta) = \frac{p}{1-p}}{var(theta) = p / (1-p)}
  \item Mode: \eqn{mode(\theta) =}{mode(theta) =}
}

The Bernoulli distribution is a binomial distribution with
\eqn{n=1}{n=1}, and one instance of a Bernoulli distribution is called a
Bernoulli trial. One coin flip is a Bernoulli trial, for example. The
categorical distribution is the generalization of the Bernoulli
distribution for variables with more than two discrete values. The
beta distribution is the conjugate prior distribution of the Bernoulli
distribution. The geometric distribution is the number of Bernoulli
trials needed to get one success.
}
\value{
  \code{dbern} gives the density,
  \code{pbern} gives the distribution function,
  \code{qbern} gives the quantile function, and
  \code{rbern} generates random deviates.
}
\seealso{
  \code{\link{dbinom}}
}
\examples{
library(LaplacesDemon)
dbern(1, 0.7)
rbern(10, 0.5)
}
\keyword{Distribution}
