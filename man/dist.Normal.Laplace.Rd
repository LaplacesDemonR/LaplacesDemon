\name{dist.Normal.Laplace}
\alias{dnormlaplace}
\alias{rnormlaplace}
\title{Normal-Laplace Distribution: Univariate Asymmetric}
\description{
  These functions provide the density and random generation for the
  univariate, asymmetric, normal-Laplace distribution with location
  parameter \eqn{\mu}{mu}, scale parameter \eqn{\sigma}{\sigma}, and
  tail-behavior parameters \eqn{\alpha}{alpha} and \eqn{\beta}{beta}.
}
\usage{
dnormlaplace(x, mu=0, sigma=1, alpha=1, beta=1, log=FALSE)
rnormlaplace(n, mu=0, sigma=1, alpha=1, beta=1)
}
\arguments{
  \item{x}{This is a vector of data.}
  \item{n}{This is the number of observations, which must be a positive
    integer that has length 1.}
  \item{mu}{This is the location parameter \eqn{\mu}{mu}.}
  \item{sigma}{This is the scale parameter \eqn{\sigma}{sigma}, which
    must be positive.}
  \item{alpha}{This is shape parameter \eqn{\alpha}{alpha} for
    left-tail behavior.}
  \item{beta}{This is shape parameter \eqn{\beta}{beta} for
    right-tail behavior.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
    \item Application: Continuous Univariate
    \item Density: \eqn{p(\theta) = \frac{\alpha\beta}{\alpha +
	\beta}\phi\frac{\theta - \mu}{\sigma} [R(\alpha\sigma -
      \frac{\theta - \mu}{\sigma}) + R(\beta\sigma + \frac{\theta -
	\mu}{\sigma})]}{p(theta) = [(alpha*beta)/(alpha + beta)] phi
      [(theta - mu) / sigma] [R(alpha*sigma - [(theta - mu) / sigma]) +
      R(beta*sigma + [(theta - mu) / sigma])]}
  \item Inventor: Reed (2006)
  \item Notation 1: \eqn{\theta \sim
    \mathrm{NL}(\mu,\sigma,\alpha,\beta)}{theta ~ NL(mu, sigma, alpha,
    beta)}
  \item Notation 2: \eqn{p(\theta) = \mathrm{NL}(\theta | \mu,
    \sigma, \alpha, \beta)}{p(theta) = NL(theta | mu, sigma, alpha,
    beta)}
  \item Parameter 1: location parameter \eqn{\mu}{mu}
  \item Parameter 2: scale parameter \eqn{\sigma > 0}{sigma > 0}
  \item Parameter 3: shape parameter \eqn{\alpha > 0}{alpha > 0}
  \item Parameter 4: shape parameter \eqn{\beta > 0}{beta > 0}
  \item Mean:
  \item Variance:
  \item Mode:
}

The normal-Laplace (NL) distribution is the convolution of a normal
distribution and a skew-Laplace distribution. When the NL distribution
is symmetric (when \eqn{\alpha = \beta}{alpha = beta}), it behaves
somewhat like the normal distribution in the middle of its range,
somewhat like the Laplace distribution in its tails, and functions
generally between the normal and Laplace distributions. Skewness is
parameterized by including a skew-Laplace component. It may be applied,
for example, to the logarithmic price of a financial instrument.

Parameters \eqn{\alpha}{alpha} and \eqn{\beta}{beta} determine the
behavior in the left and right tails, respectively. A small value
corresponds to heaviness in the corresponding tail. As
\eqn{\sigma}{sigma} approaches zero, the NL distribution approaches a
skew-Laplace distribution. As \eqn{\beta}{beta} approaches infinity,
the NL distribution approaches a normal distribution, though it never
quite reaches it.
}
\value{
  \code{dnormlaplace} gives the density, and
  \code{rnormlaplace} generates random deviates.
}
\references{
  Reed, W.J. (2006). "The Normal-Laplace Distribution and Its Relatives".
  In \emph{Advances in Distribution Theory, Order Statistics and
  Inference}, p. 61--74, Birkhauser, Boston.
}
\seealso{
  \code{\link{dalaplace}},
  \code{\link{dallaplace}},
  \code{\link{daml}},
  \code{\link{dlaplace}}, and
  \code{\link{dnorm}}
}
\examples{
library(LaplacesDemon)
x <- dnormlaplace(1,0,1,0.5,2)
x <- rnormlaplace(100,0,1,0.5,2)

#Plot Probability Functions
x <- seq(from=-5, to=5, by=0.1)
plot(x, dlaplace(x,0,0.5), ylim=c(0,1), type="l", main="Probability Function",
     ylab="density", col="red")
lines(x, dlaplace(x,0,1), type="l", col="green")
lines(x, dlaplace(x,0,2), type="l", col="blue")
legend(2, 0.9, expression(paste(mu==0, ", ", lambda==0.5),
     paste(mu==0, ", ", lambda==1), paste(mu==0, ", ", lambda==2)),
     lty=c(1,1,1), col=c("red","green","blue"))
}
\keyword{Distribution}