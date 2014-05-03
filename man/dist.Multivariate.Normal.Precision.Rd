\name{dist.Multivariate.Normal.Precision}
\alias{dmvnp}
\alias{rmvnp}
\title{Multivariate Normal Distribution: Precision Parameterization}
\description{
  These functions provide the density and random number generation for
  the multivariate normal distribution, given the precision
  parameterization.
}
\usage{
dmvnp(x, mu, Omega, log=FALSE) 
rmvnp(n=1, mu, Omega)
}
\arguments{
  \item{x}{This is data or parameters in the form of a vector of length
       \eqn{k} or a matrix with \eqn{k} columns.}
  \item{n}{This is the number of random draws.}
  \item{mu}{This is mean vector \eqn{\mu}{mu} with length \eqn{k} or
       matrix with \eqn{k} columns.}
  \item{Omega}{This is the \eqn{k \times k}{k x k} precision matrix
       \eqn{\Omega}{Omega}.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
       density is returned.}
}
\details{
  \itemize{
  \item Application: Continuous Multivariate
  \item Density: \eqn{p(\theta) = (2\pi)^{-p/2} |\Omega|^{1/2}
    \exp(-\frac{1}{2} (\theta-\mu)^T \Omega (\theta-\mu))}{p(theta) = (2*pi)^(-p/2) * |Omega|^(1/2) * exp(-(1/2)*(theta-mu)^t Omega (theta-mu))}
  \item Inventor: Unknown (to me, anyway)
  \item Notation 1: \eqn{\theta \sim \mathcal{MVN}(\mu, \Omega^{-1})}{theta ~ MVN(mu, Omega^(-1))}
  \item Notation 2: \eqn{\theta \sim \mathcal{N}_k(\mu, \Omega^{-1})}{theta ~ N[k](mu, Omega^(-1))}
  \item Notation 3: \eqn{p(\theta) = \mathcal{MVN}(\theta | \mu, \Omega^{-1})}{p(theta) = MVN(theta | mu, Omega^(-1))}
  \item Notation 4: \eqn{p(\theta) = \mathcal{N}_k(\theta | \mu, \Omega^{-1})}{p(theta) = N[k](theta | mu, Omega^(-1))}
  \item Parameter 1: location vector \eqn{\mu}{mu}
  \item Parameter 2: positive-definite \eqn{k \times k}{k x k} precision matrix \eqn{\Omega}{Omega}
  \item Mean: \eqn{E(\theta) = \mu}{E(theta) = mu}
  \item Variance: \eqn{var(\theta) = \Omega^{-1}}{var(theta) = Omega^(-1)}
  \item Mode: \eqn{mode(\theta) = \mu}{mode(theta) = mu}
}

The multivariate normal distribution, or multivariate Gaussian
distribution, is a multidimensional extension of the one-dimensional
or univariate normal (or Gaussian) distribution. It is usually
parameterized with mean and a covariance matrix, or in Bayesian
inference, with mean and a precision matrix, where the precision matrix
is the matrix inverse of the covariance matrix. These functions
provide the precision parameterization for convenience and
familiarity. It is easier to calculate a multivariate normal density
with the precision parameterization, because a matrix inversion can be
avoided.

A random vector is considered to be multivariate normally distributed
if every linear combination of its components has a univariate normal
distribution. This distribution has a mean parameter vector
\eqn{\mu}{mu} of length \eqn{k} and a \eqn{k \times k}{k x k}
precision matrix \eqn{\Omega}{Omega}, which must be positive-definite.

The conjugate prior of the mean vector is another multivariate normal
distribution. The conjugate prior of the precision matrix is the
Wishart distribution (see \code{\link{dwishart}}).

When applicable, the alternative Cholesky parameterization should be
preferred. For more information, see \code{\link{dmvnpc}}.

For models where the dependent variable, Y, is specified to be
distributed multivariate normal given the model, the Mardia test (see
\code{\link{plot.demonoid.ppc}}, \code{\link{plot.laplace.ppc}}, or
\code{\link{plot.pmc.ppc}}) may be used to test the residuals.
}
\value{
  \code{dmvnp} gives the density and 
  \code{rmvnp} generates random deviates.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{dmvn}},
  \code{\link{dmvnc}},
  \code{\link{dmvnpc}},
  \code{\link{dnorm}},
  \code{\link{dnormp}},
  \code{\link{dnormv}},
  \code{\link{dwishart}},
  \code{\link{plot.demonoid.ppc}},
  \code{\link{plot.laplace.ppc}}, and
  \code{\link{plot.pmc.ppc}}.
}
\examples{
library(LaplacesDemon)
x <- dmvnp(c(1,2,3), c(0,1,2), diag(3))
X <- rmvnp(1000, c(0,1,2), diag(3))
joint.density.plot(X[,1], X[,2], color=TRUE)
}
\keyword{Distribution}