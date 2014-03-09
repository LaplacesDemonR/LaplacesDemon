\name{dist.Inverse.Gamma}
\alias{dinvgamma}
\alias{rinvgamma}
\title{Inverse Gamma Distribution}
\description{
  This is the density function and random generation from the inverse
  gamma distribution.
}
\usage{
dinvgamma(x, shape=1, scale=1, log=FALSE)
rinvgamma(n, shape=1, scale=1)
}
\arguments{
  \item{n}{This is the number of draws from the distribution.}
  \item{x}{This is the scalar location to evaluate density.}
  \item{shape}{This is the scalar shape parameter \eqn{\alpha}{alpha},
    which defaults to one.}
  \item{scale}{This is the scalar scale parameter \eqn{\beta}{beta},
    which defaults to one.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
    \item Application: Continuous Univariate
    \item Density: \eqn{p(\theta) = \frac{\beta^\alpha}{\Gamma(\alpha)}
      \theta^{-(\alpha + 1)} \exp(-\frac{\beta}{\theta}), \quad \theta > 0}{p(theta) = (beta^alpha / Gamma(alpha)) *
      theta^(-(alpha + 1)) * exp(-beta / theta), theta > 0}
    \item Inventor: Unknown (to me, anyway)
    \item Notation 1: \eqn{\theta \sim \mathcal{G}^{-1}(\alpha, \beta)}{theta ~ G^-1(alpha, beta)}
    \item Notation 2: \eqn{p(\theta) = \mathcal{G}^{-1}(\theta | \alpha,
    \beta)}{p(theta) = G^-1(theta | alpha, beta)}
    \item Parameter 1: shape \eqn{\alpha > 0}{alpha > 0}
    \item Parameter 2: scale \eqn{\beta > 0}{beta > 0}
    \item Mean: \eqn{E(\theta) = \frac{\beta}{\alpha - 1}}{E(theta) =
      beta / (alpha - 1)}, for \eqn{\alpha > 1}{alpha > 1}
    \item Variance: \eqn{var(\theta) = \frac{\beta^2}{(\alpha - 1)^2
	(\alpha - 2)}, \alpha > 2}{var(theta) = beta^2 / ((alpha - 1)^2 * (alpha - 2)), alpha > 2}
    \item Mode: \eqn{mode(\theta) = \frac{\beta}{\alpha + 1}}{mode(theta) = beta / (alpha + 1)}
}

The inverse-gamma is the conjugate prior distribution for the normal
or Gaussian variance, and has been traditionally specified as a vague
prior in that application. The density is always finite; its integral is
finite if \eqn{\alpha > 0}{alpha > 0}. Prior information decreases as
\eqn{\alpha, \beta \rightarrow 0}{alpha, beta -> 0}.

These functions are similar to those in the \code{MCMCpack} package.
}
\value{
  \code{dinvgamma} gives the density and
  \code{rinvgamma} generates random deviates. The parameterization
  is consistent with the Gamma Distribution in the stats package.
}
\seealso{
  \code{\link{dgamma}},
  \code{\link{dnorm}},
  \code{\link{dnormp}}, and
  \code{\link{dnormv}}.
}
\examples{
library(LaplacesDemon)
x <- dinvgamma(4.3, 1.1)
x <- rinvgamma(10, 3.3)

#Plot Probability Functions
x <- seq(from=0.1, to=20, by=0.1)
plot(x, dinvgamma(x,1,1), ylim=c(0,1), type="l", main="Probability Function",
     ylab="density", col="red")
lines(x, dinvgamma(x,1,0.6), type="l", col="green")
lines(x, dinvgamma(x,0.6,1), type="l", col="blue")
legend(2, 0.9, expression(paste(alpha==1, ", ", beta==1),
     paste(alpha==1, ", ", beta==0.6), paste(alpha==0.6, ", ", beta==1)),
     lty=c(1,1,1), col=c("red","green","blue"))
}
\keyword{Distribution}

