\name{dist.Pareto}
\alias{dpareto}
\alias{ppareto}
\alias{qpareto}
\alias{rpareto}
\title{Pareto Distribution}
\description{
  These functions provide the density, distribution function, quantile
  function, and random generation for the pareto distribution.
}
\usage{
dpareto(x, alpha, log=FALSE)
ppareto(q, alpha)
qpareto(p, alpha)
rpareto(n, alpha)
}
\arguments{
  \item{x,q}{These are each a vector of quantiles.}
  \item{p}{This is a vector of probabilities.}
  \item{n}{This is the number of observations, which must be a positive
    integer that has length 1.}
  \item{alpha}{This is the shape parameter \eqn{\alpha}{alpha}, which
    must be positive.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density or result is returned.}
}
\details{
  \itemize{
  \item Application: Continuous Univariate
  \item Density: \eqn{p(\theta) =
    \frac{\alpha}{\theta^{\alpha+1}}, \theta \ge 1}{p(theta) = alpha /
    theta^(alpha + 1), theta >= 1}
  \item Inventor: Vilfredo Pareto (1848-1923)
  \item Notation 1: \eqn{\theta \sim \mathcal{PA}(\alpha)}{theta ~ PA(alpha)}
  \item Notation 2: \eqn{p(\theta) = \mathcal{PA}(\theta |
    \alpha)}{p(theta) = PA(theta | alpha)}
  \item Parameter 1: shape parameter \eqn{\alpha > 0}{alpha > 0}
  \item Mean: \eqn{E(\theta) = \frac{\alpha}{\alpha - 1}}{E(theta) =
    alpha / (alpha - 1)}
  \item Variance: \eqn{var(\theta) =
    \frac{\alpha}{(\alpha-1)^2(\alpha-2)}, \alpha > 2}{var(theta) =
    alpha / (alpha - 1)^2 (alpha - 2), alpha > 2}
  \item Mode: \eqn{mode(\theta) = 1}{mode(theta) = 1}
}

The Pareto distribution, sometimes called the Bradford distribution, is
related to the exponential distribution. The gamma distribution is the
conjugate prior distribution for the shape parameter \eqn{\alpha}{alpha}
in the Pareto distribution. The Pareto distribution is the conjugate
prior distribution for the range parameters of a uniform
distribution. An extension, elsewhere, is the symmetric Pareto
distribution.
}
\value{
  \code{dpareto} gives the density,
  \code{ppareto} gives the distribution function,
  \code{qpareto} gives the quantile function, and
  \code{rpareto} generates random deviates.
}
\seealso{
  \code{\link{dexp}},
  \code{\link{dlnorm}},
  \code{\link{dlnormp}},
  \code{\link{dnorm}},
  \code{\link{dnormp}},
  \code{\link{dnormv}}.
}
\examples{
library(LaplacesDemon)
x <- dpareto(1,1)
x <- ppareto(0.5,1)
x <- qpareto(0.5,1)
x <- rpareto(10,1)

#Plot Probability Functions
x <- seq(from=1, to=5, by=0.01)
plot(x, dpareto(x,0.1), ylim=c(0,1), type="l", main="Probability Function",
     ylab="density", col="red")
lines(x, dpareto(x,0.5), type="l", col="green")
lines(x, dpareto(x,1), type="l", col="blue")
legend(2, 0.9, expression(alpha==0.1, alpha==0.5, alpha==1),
     lty=c(1,1,1), col=c("red","green","blue"))
}
\keyword{Distribution}