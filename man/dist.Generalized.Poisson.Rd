\name{dist.Generalized.Poisson}
\alias{dgpois}
\title{Generalized Poisson Distribution}
\description{
  The density function is provided for the univariate, discrete,
  generalized Poisson distribution with location parameter
  \eqn{\lambda}{lambda} and scale parameter \eqn{\omega}{omega}.
}
\usage{
dgpois(x, lambda=0, omega=0, log=FALSE)
}
\arguments{
  \item{x}{This is a vector of quantiles.}
  \item{lambda}{This is the parameter \eqn{\lambda}{lambda}.}
  \item{omega}{This is the parameter \eqn{\omega}{omega}, which
    should be in the interval [0,1) for positive counts.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
  \item Application: Discrete Univariate
  \item Density: \eqn{p(\theta) = (1 - \omega) \lambda \frac{[(1 -
      \omega) \lambda + \omega \theta]^{\theta - 1}}{\theta!} \exp{-[(1
      - \omega) \lambda + \omega \theta]}}{(1 - omega) lambda (((1 -
    omega) lambda + omega theta)^(y-1) / theta!) e(-((1 - omega) lambda
    + omega theta))}
  \item Inventor: Consul (1989) and Ntzoufras et al. (2005)
  \item Notation 1: \eqn{\theta \sim
    \mathrm{GP}(\lambda,\omega)}{theta ~ GP(lambda, omega)}
  \item Notation 2: \eqn{p(\theta) = \mathrm{GP}(\theta | \lambda,
    \omega)}{p(theta) = GP(theta | lambda, omega)}
  \item Parameter 1: location parameter \eqn{\lambda}{lambda}
  \item Parameter 2: scale parameter \eqn{\omega \in [0,1)}{omega in [0,1)}
  \item Mean: \eqn{E(\theta) = \lambda}{E(theta) = lambda}
  \item Variance: \eqn{var(\theta) = \lambda(1 -
    \omega)^{-2}}{var(theta) = lambda(1 - omega)^(-2)}
}

The generalized Poisson distribution (Consul, 1989) is also called the
Lagrangian Poisson distribution. The simple Poisson distribution is a
special case of the generalized Poisson distribution. The generalized
Poisson distribution is used in generalized Poisson regression as an
extension of Poisson regression that accounts for overdispersion.

The \code{dgpois} function is parameterized according to Ntzoufras et
al. (2005), which is easier to interpret and estimates better with MCMC.

Valid values for omega are in the interval [0,1) for positive counts.
For \eqn{\omega = 0}{omega = 0}, the generalized Poisson reduces to a
simple Poisson with mean \eqn{\lambda}{lambda}. Note that it is possible
for \eqn{\omega < 0}{omega < 0}, but this implies underdispersion in
count data, which is uncommon. The \code{dgpois} function returns
warnings or errors, so \eqn{\omega}{omega} should be non-negative here.

The dispersion index (DI) is a variance-to-mean ratio, and is \eqn{DI =
  (1 - \omega)^{-2}}{DI = (1 - omega)^(-2)}. A simple Poisson has DI=1.
When DI is far from one, the assumption that the variance equals the
mean of a simple Poisson is violated.
}
\value{
  \code{dgpois} gives the density.
}
\references{
  Consul, P. (1989). `"Generalized Poisson Distribution: Properties and
  Applications". Marcel Decker: New York, NY.

  Ntzoufras, I., Katsis, A., and Karlis, D. (2005). "Bayesian Assessment
  of the Distribution of Insurance Claim Counts using Reversible Jump
  MCMC", \emph{North American Actuarial Journal}, 9, p. 90--108.
}
\seealso{
  \code{\link{dnbinom}} and
  \code{\link{dpois}}.
}
\examples{
library(LaplacesDemon)
y <- rpois(100, 5)
lambda <- rpois(100, 5)
x <- dgpois(y, lambda, 0.5)

#Plot Probability Functions
x <- seq(from=0, to=20, by=1)
plot(x, dgpois(x,1,0.5), ylim=c(0,1), type="l", main="Probability Function",
     ylab="density", col="red")
lines(x, dlaplace(x,1,0.6), type="l", col="green")
lines(x, dlaplace(x,1,0.7), type="l", col="blue")
legend(2, 0.9, expression(paste(lambda==1, ", ", omega==0.5),
     paste(lambda==1, ", ", omega==0.6), paste(lambda==1, ", ", omega==0.7)),
     lty=c(1,1,1), col=c("red","green","blue"))
}
\keyword{Distribution}