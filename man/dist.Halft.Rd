\name{dist.Halft}
\alias{dhalft}
\alias{phalft}
\alias{qhalft}
\alias{rhalft}
\title{Half-t Distribution}
\description{
  These functions provide the density, distribution function, quantile function, and random generation for the half-t distribution.
}
\usage{
dhalft(x, scale=25, nu=1, log=FALSE)
phalft(q, scale=25, nu=1)
qhalft(p, scale=25, nu=1)
rhalft(n, scale=25, nu=1)
}
\arguments{
  \item{x,q}{These are each a vector of quantiles.}
  \item{p}{This is a vector of probabilities.}
  \item{n}{This is the number of observations, which must be a positive
    integer that has length 1.}
  \item{scale}{This is the scale parameter \eqn{\alpha}{alpha}, which
    must be positive.}
  \item{nu}{This is the scalar degrees of freedom parameter, which is
    usually represented as \eqn{\nu}{nu}.}
  \item{log}{Logical. If \code{log=TRUE} then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
  \item Application: Continuous Univariate
  \item Density: \eqn{p(\theta) = (1 + \frac{1}{\nu} (\theta /
    \alpha)^2)^{-(\nu+1)/2}, \quad \theta \ge 0}{p(theta) = (1 +
    (1/nu)*(theta/alpha)^2)^(-(nu+1)/2), theta >= 0}
  \item Inventor: Derived from the Student t
  \item Notation 1: \eqn{\theta \sim \mathcal{HT}(\alpha, \nu)}{theta ~
    HT(alpha,nu)}
  \item Notation 2: \eqn{p(\theta) = \mathcal{HT}(\theta | \alpha,
    \nu)}{p(theta) = HT(theta | alpha,nu)}
  \item Parameter 1: scale parameter \eqn{\alpha > 0}{alpha > 0}
  \item Parameter 2: degrees of freedom parameter \eqn{\nu}{nu}
  \item Mean: \eqn{E(\theta)}{E(theta)} = unknown
  \item Variance: \eqn{var(\theta)}{var(theta)} = unknown
  \item Mode: \eqn{mode(\theta) = 0}{mode(theta) = 0}
}

The half-t distribution is derived from the Student t distribution, and
is useful as a weakly informative prior distribution for a scale
parameter. It is more adaptable than the default recommended
half-Cauchy, though it may also be more difficult to estimate due to its
additional degrees of freedom parameter, \eqn{\nu}{nu}. When
\eqn{\nu=1}{nu=1}, the density is proportional to a proper half-Cauchy
distribution. When \eqn{\nu=-1}{nu=-1}, the density becomes an improper,
uniform prior distribution. For more information on propriety, see
\code{is.proper}.

Wand et al. (2011) demonstrated that the half-t distribution may be
represented as a scale mixture of inverse-gamma distributions. This
representation is useful for conjugacy.
}
\value{
  \code{dhalft} gives the density,
  \code{phalft} gives the distribution function,
  \code{qhalft} gives the quantile function, and
  \code{rhalft} generates random deviates.
}
\references{
  Wand, M.P., Ormerod, J.T., Padoan, S.A., and Fruhwirth, R. (2011).
  "Mean Field Variational Bayes for Elaborate Distributions".
  \emph{Bayesian Analysis}, 6: p. 847--900.
}
\seealso{
  \code{\link{dhalfcauchy}},
  \code{\link{dst}},
  \code{\link{dt}},
  \code{\link{dunif}}, and
  \code{\link{is.proper}}.
}
\examples{
library(LaplacesDemon)
x <- dhalft(1,25,1)
x <- phalft(1,25,1)
x <- qhalft(0.5,25,1)
x <- rhalft(10,25,1)

#Plot Probability Functions
x <- seq(from=0.1, to=20, by=0.1)
plot(x, dhalft(x,1,-1), ylim=c(0,1), type="l", main="Probability Function",
     ylab="density", col="red")
lines(x, dhalft(x,1,0.5), type="l", col="green")
lines(x, dhalft(x,1,500), type="l", col="blue")
legend(2, 0.9, expression(paste(alpha==1, ", ", nu==-1),
     paste(alpha==1, ", ", nu==0.5), paste(alpha==1, ", ", nu==500)),
     lty=c(1,1,1), col=c("red","green","blue"))
}
\keyword{Distribution}