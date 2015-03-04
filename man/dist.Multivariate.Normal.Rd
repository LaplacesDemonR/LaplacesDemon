\name{dist.Multivariate.Normal}
\alias{dmvn}
\alias{rmvn}
\title{Multivariate Normal Distribution}
\description{
  These functions provide the density and random number generation for
  the multivariate normal distribution.
}
\usage{
dmvn(x, mu, Sigma, log=FALSE) 
rmvn(n=1, mu, Sigma)
}
\arguments{
  \item{x}{This is data or parameters in the form of a vector of length
       \eqn{k} or a matrix with \eqn{k} columns.}
  \item{n}{This is the number of random draws.}
  \item{mu}{This is mean vector \eqn{\mu}{mu} with length \eqn{k} or
       matrix with \eqn{k} columns.}
  \item{Sigma}{This is the \eqn{k \times k}{k x k} covariance matrix
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
  \item Inventors: Robert Adrain (1808), Pierre-Simon Laplace (1812), and Francis Galton (1885)
  \item Notation 1: \eqn{\theta \sim \mathcal{MVN}(\mu, \Sigma)}{theta ~ MVN(mu, Sigma)}
  \item Notation 2: \eqn{\theta \sim \mathcal{N}_k(\mu, \Sigma)}{theta ~ N[k](mu, Sigma)}
  \item Notation 3: \eqn{p(\theta) = \mathcal{MVN}(\theta | \mu, \Sigma)}{p(theta) = MVN(theta | mu, Sigma)}
  \item Notation 4: \eqn{p(\theta) = \mathcal{N}_k(\theta | \mu, \Sigma)}{p(theta) = N[k](theta | mu, Sigma)}
  \item Parameter 1: location vector \eqn{\mu}{mu}
  \item Parameter 2: positive-definite \eqn{k \times k}{k x k} covariance matrix \eqn{\Sigma}{Sigma}
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
\eqn{k} and a \eqn{k \times k}{k x k} covariance matrix
\eqn{\Sigma}{Sigma}, which must be positive-definite.

The conjugate prior of the mean vector is another multivariate normal
distribution. The conjugate prior of the covariance matrix is the
inverse Wishart distribution (see \code{\link{dinvwishart}}).

When applicable, the alternative Cholesky parameterization should be
preferred. For more information, see \code{\link{dmvnc}}.

For models where the dependent variable, Y, is specified to be
distributed multivariate normal given the model, the Mardia test (see
\code{\link{plot.demonoid.ppc}}, \code{\link{plot.laplace.ppc}}, or
\code{\link{plot.pmc.ppc}}) may be used to test the residuals.
}
\value{
  \code{dmvn} gives the density and 
  \code{rmvn} generates random deviates.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{dinvwishart}},
  \code{\link{dmatrixnorm}},
  \code{\link{dmvnc}},
  \code{\link{dmvnp}},
  \code{\link{dnorm}},
  \code{\link{dnormp}},
  \code{\link{dnormv}},
  \code{\link{plot.demonoid.ppc}},
  \code{\link{plot.laplace.ppc}}, and
  \code{\link{plot.pmc.ppc}}.
}
\examples{
library(LaplacesDemon)
x <- dmvn(c(1,2,3), c(0,1,2), diag(3))
X <- rmvn(1000, c(0,1,2), diag(3))
joint.density.plot(X[,1], X[,2], color=TRUE)
}
\keyword{Distribution}