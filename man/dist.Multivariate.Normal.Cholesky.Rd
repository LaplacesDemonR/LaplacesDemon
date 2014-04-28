\name{dist.Multivariate.Normal.Cholesky}
\alias{dmvnc}
\alias{rmvnc}
\title{Multivariate Normal Distribution: Cholesky Parameterization}
\description{
  These functions provide the density and random number generation for
  the multivariate normal distribution, given the Cholesky
  parameterization.
}
\usage{
dmvnc(x, mu, U, log=FALSE) 
rmvnc(n=1, mu, U)
}
\arguments{
  \item{x}{This is data or parameters in the form of a vector of length
       \eqn{k} or a matrix with \eqn{k} columns.}
  \item{n}{This is the number of random draws.}
  \item{mu}{This is mean vector \eqn{\mu}{mu} with length \eqn{k} or
       matrix with \eqn{k} columns.}
  \item{U}{This is the \eqn{k \times k}{k x k} upper-triangular matrix
    that is Cholesky factor \eqn{\textbf{U}}{U} of covariance matrix
    \eqn{\Sigma}{Sigma}.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
       density is returned.}
}
\details{
  \itemize{
  \item Application: Continuous Multivariate
  \item Density: \eqn{p(\theta) = \frac{1}{(2\pi)^{k/2}|\Sigma|^{1/2}}
    \exp(-\frac{1}{2}(\theta - \mu)'\Sigma^{-1}(\theta - \mu))}{p(theta)
    = (1/((2*pi)^(k/2)*|Sigma|^(1/2))) *
    exp(-(1/2)*(theta-mu)'*Sigma^(-1)*(theta-mu))}
  \item Inventor: Unknown (to me, anyway)
  \item Notation 1: \eqn{\theta \sim \mathcal{MVN}(\mu, \Sigma)}{theta ~ MVN(mu, Sigma)}
  \item Notation 2: \eqn{\theta \sim \mathcal{N}_k(\mu, \Sigma)}{theta ~ N[k](mu, Sigma)}
  \item Notation 3: \eqn{p(\theta) = \mathcal{MVN}(\theta | \mu, \Sigma)}{p(theta) = MVN(theta | mu, Sigma)}
  \item Notation 4: \eqn{p(\theta) = \mathcal{N}_k(\theta | \mu, \Sigma)}{p(theta) = N[k](theta | mu, Sigma)}
  \item Parameter 1: location vector \eqn{\mu}{mu}
  \item Parameter 2: \eqn{k \times k}{k x k} positive-definite matrix \eqn{\Sigma}{Sigma}
  \item Mean: \eqn{E(\theta) = \mu}{E(theta) = mu}
  \item Variance: \eqn{var(\theta) = \Sigma}{var(theta) = Sigma}
  \item Mode: \eqn{mode(\theta) = \mu}{mode(theta) = mu}
}

The multivariate normal distribution, or multivariate Gaussian
distribution, is a multidimensional extension of the one-dimensional
or univariate normal (or Gaussian) distribution. A random vector is
considered to be multivariate normally distributed if every linear
combination of its components has a univariate normal distribution.
This distribution has a mean parameter vector \eqn{\mu}{mu} of length
\eqn{k} and an upper-triangular \eqn{k \times k}{k x k} matrix that is
Cholesky factor \eqn{\textbf{U}}{U}, as per the \code{\link{chol}}
function for Cholesky decomposition.

In practice, \eqn{\textbf{U}}{U} is fully unconstrained for proposals
when its diagonal is log-transformed. The diagonal is exponentiated
after a proposal and before other calculations. Overall, the Cholesky
parameterization is faster than the traditional parameterization.
Compared with \code{dmvn}, \code{dmvnc} must additionally
matrix-multiply the Cholesky back to the covariance matrix, but it
does not have to check for or correct the covariance matrix to
positive-definiteness, which overall is slower. Compared with
\code{rmvn}, \code{rmvnc} is faster because the Cholesky decomposition
has already been performed.

For models where the dependent variable, Y, is specified to be
distributed multivariate normal given the model, the Mardia test (see
\code{\link{plot.demonoid.ppc}}, \code{\link{plot.laplace.ppc}}, or
\code{\link{plot.pmc.ppc}}) may be used to test the residuals.
}
\value{
  \code{dmvnc} gives the density and 
  \code{rmvnc} generates random deviates.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{chol}},
  \code{\link{dinvwishartc}},
  \code{\link{dmvn}},
  \code{\link{dmvnp}},
  \code{\link{dmvnpc}},
  \code{\link{dnorm}},
  \code{\link{dnormp}},
  \code{\link{dnormv}},
  \code{\link{plot.demonoid.ppc}},
  \code{\link{plot.laplace.ppc}}, and
  \code{\link{plot.pmc.ppc}}.
}
\examples{
library(LaplacesDemon)
Sigma <- diag(3)
U <- chol(Sigma)
x <- dmvnc(c(1,2,3), c(0,1,2), U)
X <- rmvnc(1000, c(0,1,2), U)
joint.density.plot(X[,1], X[,2], color=TRUE)
}
\keyword{Distribution}