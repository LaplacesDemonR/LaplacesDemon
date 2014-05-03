\name{dist.Inverse.Beta}
\alias{dinvbeta}
\alias{rinvbeta}
\title{Inverse Beta Distribution}
\description{
  This is the density function and random generation from the inverse
  beta distribution.
}
\usage{
dinvbeta(x, a, b, log=FALSE)
rinvbeta(n, a, b)
}
\arguments{
  \item{n}{This is the number of draws from the distribution.}
  \item{x}{This is a location vector at which to evaluate density.}
  \item{a}{This is the scalar shape parameter \eqn{\alpha}{alpha}.}
  \item{b}{This is the scalar shape parameter \eqn{\beta}{beta}}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
    \item Application: Continuous Univariate
    \item Density: \eqn{p(\theta) = \frac{\theta^{\alpha - 1} (1 +
	\theta)^{-\alpha - \beta}}{\beta(\alpha, \beta)}}{(theta^(alpha -
      1) * (1 + theta)^(-alpha - beta)) / beta(alpha, beta)}
    \item Inventor: Dubey (1970)
    \item Notation 1: \eqn{\theta \sim \mathcal{B}^{-1}(\alpha, \beta)}{theta ~ B^-1(alpha, beta)}
    \item Notation 2: \eqn{p(\theta) = \mathcal{B}^{-1}(\theta | \alpha,
    \beta)}{p(theta) = B^-1(theta | alpha, beta)}
    \item Parameter 1: shape \eqn{\alpha > 0}{alpha > 0}
    \item Parameter 2: shape \eqn{\beta > 0}{beta > 0}
    \item Mean: \eqn{E(\theta) = \frac{\alpha}{\beta - 1}}{E(theta) =
      alpha / (beta - 1)}, for \eqn{\beta > 1}{beta > 1}
    \item Variance: \eqn{var(\theta) = \frac{\alpha(\alpha + \beta -
	1)}{(\beta - 1)^2 (\beta - 2)}}{var(theta) = (alpha * (alpha + beta
      - 1)) / ((beta - 1)^2 * (beta - 2))}
    \item Mode: \eqn{mode(\theta) = \frac{\alpha - 1}{\beta +
	1}}{mode(theta) = (alpha - 1) / (beta + 1)}
}

The inverse-beta, also called the beta prime distribution, applies to
variables that are continuous and positive. The inverse beta is the
conjugate prior distribution of a parameter of a Bernoulli distribution
expressed in odds.

The inverse-beta distribution has also been extended to the generalized
beta prime distribution, though it is not (yet) included here.

}
\value{
  \code{dinvbeta} gives the density and
  \code{rinvbeta} generates random deviates.
}
\references{
  Dubey, S.D. (1970). "Compound Gamma, Beta and F Distributions".
  \emph{Metrika}, 16, p. 27--31.
}
\seealso{
  \code{\link{dbeta}}
}
\examples{
library(LaplacesDemon)
x <- dinvbeta(5:10, 2, 3)
x <- rinvbeta(10, 2, 3)

#Plot Probability Functions
x <- seq(from=0.1, to=20, by=0.1)
plot(x, dinvbeta(x,2,2), ylim=c(0,1), type="l", main="Probability Function",
     ylab="density", col="red")
lines(x, dinvbeta(x,2,3), type="l", col="green")
lines(x, dinvbeta(x,3,2), type="l", col="blue")
legend(2, 0.9, expression(paste(alpha==2, ", ", beta==2),
     paste(alpha==2, ", ", beta==3), paste(alpha==3, ", ", beta==2)),
     lty=c(1,1,1), col=c("red","green","blue"))
}
\keyword{Distribution}

