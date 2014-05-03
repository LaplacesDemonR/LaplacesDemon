\name{dist.Skew.Laplace}
\alias{dslaplace}
\alias{pslaplace}
\alias{qslaplace}
\alias{rslaplace}
\title{Skew-Laplace Distribution: Univariate}
\description{
  These functions provide the density, distribution function, quantile
  function, and random generation for the univariate, skew-Laplace
  distribution with location parameter \eqn{\mu}{mu}, and two mixture
  parameters: \eqn{\alpha}{alpha} and \eqn{\beta}{beta}.
}
\usage{
dslaplace(x, mu, alpha, beta, log=FALSE)
pslaplace(q, mu, alpha, beta)
qslaplace(p, mu, alpha, beta)
rslaplace(n, mu, alpha, beta)
}
\arguments{
  \item{x, q}{These are each a vector of quantiles.}
  \item{p}{This is a vector of probabilities.}
  \item{n}{This is the number of observations, which must be a positive
    integer that has length 1.}
  \item{mu}{This is the location parameter \eqn{\mu}{mu}.}
  \item{alpha}{This is a mixture parameter \eqn{\alpha}{alpha}, which
    must be positive.}
  \item{beta}{This is a mixture parameter \eqn{\beta}{beta}, which
    must be positive.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
  \item Application: Continuous Univariate
  \item Density 1: \eqn{p(\theta) = \frac{1}{\alpha + \beta}
    \exp(\frac{\theta - \mu}{\alpha}), \theta \le \mu}{p(theta) =
    (1/(alpha+beta)) exp((theta-mu)/alpha), theta <= mu}
  \item Density 2: \eqn{p(\theta) = \frac{1}{\alpha + \beta}
    \exp(\frac{\mu - \theta}{\beta}), \theta > \mu}{p(theta) =
    (1/(alpha+beta)) exp((mu-theta)/beta), theta > mu}
  \item Inventor: Fieller, et al. (1992)
  \item Notation 1: \eqn{\theta \sim \mathcal{SL}(\mu, \alpha, \beta)}{theta ~
    SL(mu, alpha, beta)}
  \item Notation 2: \eqn{p(\theta) = \mathcal{SL}(\theta | \mu,
    \alpha, \beta)}{p(theta) = SL(theta | mu, alpha, beta)}
  \item Parameter 1: location parameter \eqn{\mu}{mu}
  \item Parameter 2: mixture parameter \eqn{\alpha > 0}{alpha > 0}
  \item Parameter 3: mixture parameter \eqn{\beta > 0}{beta > 0}
  \item Mean: \eqn{E(\theta) = \mu + \beta - \alpha}{E(theta) = mu +
    beta - alpha}
  \item Variance: \eqn{var(\theta) = \alpha^2 + \beta^2}{var(theta) =
    alpha^2 + beta^2}
  \item Mode: \eqn{mode(\theta) = \mu}{mode(theta) = mu}
}

This is the three-parameter general skew-Laplace distribution, which is
an extension of the two-parameter central skew-Laplace distribution. The
general form allows the mode to be shifted along the real line with
parameter \eqn{\mu}{mu}. In contrast, the central skew-Laplace has mode
zero, and may be reproduced here by setting \eqn{\mu=0}{mu=0}.

The general skew-Laplace distribution is a mixture of a negative
exponential distribution with mean \eqn{\beta}{beta}, and the negative
of an exponential distribution with mean \eqn{\alpha}{alpha}. The
weights of the positive and negative components are proportional to
their means. The distribution is symmetric when
\eqn{\alpha=\beta}{alpha=beta}, in which case the mean is \eqn{\mu}{mu}.

These functions are similar to those in the \code{HyperbolicDist} package.
}
\value{
  \code{dslaplace} gives the density,
  \code{pslaplace} gives the distribution function,
  \code{qslaplace} gives the quantile function, and
  \code{rslaplace} generates random deviates.
}
\references{
  Fieller, N.J., Flenley, E.C., and Olbricht, W. (1992). "Statistics of
  Particle Size Data". \emph{Applied Statistics}, 41, p. 127--146.
}
\seealso{
  \code{\link{dalaplace}},
  \code{\link{dexp}},
  \code{\link{dlaplace}}, 
  \code{\link{dlaplacep}}, and
  \code{\link{dsdlaplace}}.
}
\examples{
library(LaplacesDemon)
x <- dslaplace(1,0,1,1)
x <- pslaplace(1,0,1,1)
x <- qslaplace(0.5,0,1,1)
x <- rslaplace(100,0,1,1)

#Plot Probability Functions
x <- seq(from=0.1, to=3, by=0.01)
plot(x, dslaplace(x,0,1,1), ylim=c(0,1), type="l", main="Probability Function",
     ylab="density", col="red")
lines(x, dslaplace(x,0,0.5,2), type="l", col="green")
lines(x, dslaplace(x,0,2,0.5), type="l", col="blue")
legend(1.5, 0.9, expression(paste(mu==0, ", ", alpha==1, ", ", beta==1),
     paste(mu==0, ", ", alpha==0.5, ", ", beta==2),
     paste(mu==0, ", ", alpha==2, ", ", beta==0.5)),
     lty=c(1,1,1), col=c("red","green","blue"))
}
\keyword{Distribution}