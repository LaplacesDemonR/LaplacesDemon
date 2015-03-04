\name{dist.Normal.Inverse.Wishart}
\alias{dnorminvwishart}
\alias{rnorminvwishart}
\title{Normal-Inverse-Wishart Distribution}
\description{
  These functions provide the density and random number generation for
  the normal-inverse-Wishart distribution.
}
\usage{
dnorminvwishart(mu, mu0, lambda, Sigma, S, nu, log=FALSE) 
rnorminvwishart(n=1, mu0, lambda, S, nu)
}
\arguments{
  \item{mu}{This is data or parameters in the form of a vector of length
    \eqn{k} or a matrix with \eqn{k} columns.}
  \item{mu0}{This is mean vector \eqn{\mu_0}{mu[0]} with length \eqn{k} or
    matrix with \eqn{k} columns.}
  \item{lambda}{This is a positive-only scalar.}
  \item{n}{This is the number of random draws.}
  \item{nu}{This is the scalar degrees of freedom \eqn{\nu}{nu}.}
  \item{Sigma}{This is a \eqn{k \times k}{k x k} covariance matrix
    \eqn{\Sigma}{Sigma}.}
  \item{S}{This is the symmetric, positive-semidefinite, \eqn{k \times
      k}{k x k} scale matrix \eqn{\textbf{S}}{S}.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
       density is returned.}
}
\details{
  \itemize{
    \item Application: Continuous Multivariate
    \item Density: \eqn{p(\mu, \Sigma) = \mathcal{N}(\mu | \mu_0,
      \frac{1}{\lambda}\Sigma) \mathcal{W}^{-1}(\Sigma | \nu,
      \textbf{S})}{p(mu, Sigma) = N(mu | mu[0], (1/lambda) Sigma)
      W^(-1)(Sigma | nu, S)}
  \item Inventors: Unknown
  \item Notation 1: \eqn{(\mu, \Sigma) \sim \mathcal{NIW}(\mu_0, \lambda,
    \textbf{S}, \nu)}{(mu, Sigmaa) ~ NIW(mu[0], lambda, S, nu)}
  \item Notation 2: \eqn{p(\mu, \Sigma) = \mathcal{NIW}(\mu, \Sigma |
    \mu_0, \lambda, \textbf{S}, \nu)}{p(mu, Sigma) = NIW(mu, Sigma |
    mu[0], lambda, S, nu)}
  \item Parameter 1: location vector \eqn{\mu_0}{mu[0]}
  \item Parameter 2: \eqn{\lambda > 0}{lambda > 0}
  \item Parameter 3: symmetric, positive-semidefinite
    \eqn{k \times k}{k x k} scale matrix \eqn{\textbf{S}}{S}
  \item Parameter 4: degrees of freedom \eqn{\nu \ge k}{nu >= k}
  \item Mean: Unknown
  \item Variance: Unknown
  \item Mode: Unknown
}

The normal-inverse-Wishart distribution, or Gaussian-inverse-Wishart
distribution, is a multivariate four-parameter continuous probability
distribution. It is the conjugate prior of a multivariate normal
distribution with unknown mean and covariance matrix.
}
\value{
  \code{dnorminvwishart} gives the density and 
  \code{rnorminvwishart} generates random deviates and returns a list
    with two components.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{dmvn}} and
  \code{\link{dinvwishart}}.
}
\examples{
library(LaplacesDemon)
K <- 3
mu <- rnorm(K)
mu0 <- rnorm(K)
nu <- K + 1
S <- diag(K)
lambda <- runif(1) #Real scalar
Sigma <- as.positive.definite(matrix(rnorm(K^2),K,K))
x <- dnorminvwishart(mu, mu0, lambda, Sigma, S, nu, log=TRUE)
out <- rnorminvwishart(n=10, mu0, lambda, S, nu)
joint.density.plot(out$mu[,1], out$mu[,2], color=TRUE)
}
\keyword{Distribution}