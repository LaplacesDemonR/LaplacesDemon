\name{dist.Multivariate.Normal.Precision.Cholesky}
\alias{dmvnpc}
\alias{rmvnpc}
\title{Multivariate Normal Distribution: Precision-Cholesky Parameterization}
\description{
  These functions provide the density and random number generation for
  the multivariate normal distribution, given the precision-Cholesky
  parameterization.
}
\usage{
dmvnpc(x, mu, U, log=FALSE) 
rmvnpc(n=1, mu, U)
}
\arguments{
  \item{x}{This is data or parameters in the form of a vector of length
       \eqn{k} or a matrix with \eqn{k} columns.}
  \item{n}{This is the number of random draws.}
  \item{mu}{This is mean vector \eqn{\mu}{mu} with length \eqn{k} or
       matrix with \eqn{k} columns.}
  \item{U}{This is the \eqn{k \times k}{k x k} upper-triangular of the
    precision matrix that is Cholesky factor \eqn{\textbf{U}}{U} of
    precision matrix \eqn{\Omega}{Omega}.}
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
provide the precision-Cholesky parameterization for convenience and
familiarity. It is easier to calculate a multivariate normal density
with the precision parameterization, because a matrix inversion can be
avoided. The precision matrix is replaced with an upper-triangular
\eqn{k \times k}{k x k} matrix that is Cholesky factor
\eqn{\textbf{U}}{U}, as per the \code{\link{chol}} function for Cholesky
decomposition.

A random vector is considered to be multivariate normally distributed
if every linear combination of its components has a univariate normal
distribution. This distribution has a mean parameter vector
\eqn{\mu}{mu} of length \eqn{k} and a \eqn{k \times k}{k x k}
precision matrix \eqn{\Omega}{Omega}, which must be positive-definite.

In practice, \eqn{\textbf{U}}{U} is fully unconstrained for proposals
when its diagonal is log-transformed. The diagonal is exponentiated
after a proposal and before other calculations. Overall, Cholesky
parameterization is faster than the traditional parameterization.
Compared with \code{dmvnp}, \code{dmvnpc} must additionally
matrix-multiply the Cholesky back to the covariance matrix, but it
does not have to check for or correct the precision matrix to
positive-definiteness, which overall is slower. Compared with
\code{rmvnp}, \code{rmvnpc} is faster because the Cholesky decomposition
has already been performed.

For models where the dependent variable, Y, is specified to be
distributed multivariate normal given the model, the Mardia test (see
\code{\link{plot.demonoid.ppc}}, \code{\link{plot.laplace.ppc}}, or
\code{\link{plot.pmc.ppc}}) may be used to test the residuals.
}
\value{
  \code{dmvnpc} gives the density and 
  \code{rmvnpc} generates random deviates.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{chol}},
  \code{\link{dmvn}},
  \code{\link{dmvnc}},
  \code{\link{dmvnp}},
  \code{\link{dnorm}},
  \code{\link{dnormp}},
  \code{\link{dnormv}},
  \code{\link{dwishartc}},
  \code{\link{plot.demonoid.ppc}},
  \code{\link{plot.laplace.ppc}}, and
  \code{\link{plot.pmc.ppc}}.
}
\examples{
library(LaplacesDemon)
Omega <- diag(3)
U <- chol(Omega)
x <- dmvnpc(c(1,2,3), c(0,1,2), U)
X <- rmvnpc(1000, c(0,1,2), U)
joint.density.plot(X[,1], X[,2], color=TRUE)
}
\keyword{Distribution}