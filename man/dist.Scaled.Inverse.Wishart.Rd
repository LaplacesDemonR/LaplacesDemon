\name{dist.Scaled.Inverse.Wishart}
\alias{dsiw}
\alias{rsiw}
\title{Scaled Inverse Wishart Distribution}
\description{
  These functions provide the density and random number generation
     for the scaled inverse Wishart distribution.
}
\usage{
   dsiw(Q, nu, S, zeta, mu, delta, log=FALSE)
   rsiw(nu, S, mu, delta)
}
\arguments{
  \item{Q}{This is the symmetric, positive-definite
    \eqn{k \times k}{k x k} matrix \eqn{\textbf{Q}}.}
  \item{nu}{This is the scalar degrees of freedom, \eqn{\nu}{nu}
    regarding \eqn{\textbf{Q}}. The default recommendation is
    \code{nu=k+1}.}
  \item{S}{This is the symmetric, positive-semidefinite
    \eqn{k \times k}{k x k} scale matrix \eqn{\textbf{S}}{S} regarding
    \eqn{\textbf{Q}}. The default recommendation is \code{S=diag(k)}.}
  \item{zeta}{This is a positive-only vector of length \eqn{k} of
    auxiliary scale parameters \eqn{\zeta}{zeta}.}
  \item{mu}{This is a vector of length \eqn{k} of location
    hyperparameters \eqn{\mu}{mu} regarding \eqn{\zeta}{zeta}.}
  \item{delta}{This is a positive-only vector of length \eqn{k} of
    scale hyperparameters \eqn{\delta}{delta} regarding
    \eqn{\zeta}{zeta}.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
    \item Application: Continuous Multivariate
    \item Density: (see below)
    \item Inventor: O'Malley and Zaslavsky (2005)
    \item Notation 1: \eqn{p(\Sigma) \sim \mathcal{SIW}(\textbf{Q},
      \nu, \textbf{S}, \zeta, \mu, \delta)}{p(Sigma) ~ SIW(Q, nu,
      S, zeta, mu, delta)}
    \item Notation 2: \eqn{p(\Sigma) = \mathcal{SIW}(\Sigma |
      \textbf{Q}, \nu, \textbf{S}, \zeta, \mu, \delta}{p(Sigma) =
      SIW(Sigma | Q, nu, S, zeta, mu, delta)}
    \item Parameter 1: symmetric, positive-definite
      \eqn{k \times k}{k x k} matrix \eqn{\textbf{Q}}{Q}
    \item Parameter 2: degrees of freedom \eqn{\nu}{nu}
    \item Parameter 3: symmetric, positive-semidefinite
      \eqn{k \times k}{k x k} scale matrix \eqn{\textbf{S}}{S}
    \item Parameter 4: Auxiliary scale parameter vector
      \eqn{\zeta}{zeta}
    \item Parameter 5: Hyperparameter location vector
      \eqn{\mu}{mu}
    \item Parameter 6: Hyperparameter scale vector
      \eqn{\delta}{delta}
    \item Mean: 
    \item Variance: 
    \item Mode: 
}

The scaled inverse Wishart (SIW) distribution is a prior probability
distribution for a covariance matrix, and is an alternative to the
inverse Wishart distribution.

While the inverse Wishart distribution is applied directly to covariance
matrix \eqn{\Sigma}{Sigma}, the SIW distribution is applied to a
decomposed matrix \eqn{\textbf{Q}} and diagonal scale matrix
\eqn{\zeta}{zeta}. For information on how to apply it to
\eqn{\textbf{Q}}, see the example below.

SIW is more flexible than the inverse Wishart distribution because it
has additional, and some say somewhat redundant, scale parameters. This
makes up for one limitation of the inverse Wishart, namely that all
uncertainty about posterior variances is represented in one parameter.
The SIW prior may somewhat alleviate the dependency in the inverse
Wishart between variances and correlations, though the SIW prior still
retains some of this relationship.

The Huang-Wand (\code{\link{dhuangwand}}) prior is a hierarchical
alternative.
}
\value{
  \code{dsiw} gives the density and
  \code{rsiw} generates random deviates.
}
\references{
  O'Malley, A.J. and Zaslavsky, A.M. (2005), "Domain-Level Covariance
  Analysis for Survey Data with Structured Nonresponse".
}
\seealso{
  \code{\link{dhuangwand}},
  \code{\link{dinvwishartc}},
  \code{\link{dmvn}}, and
  \code{\link{dwishart}}.
}
\examples{
library(LaplacesDemon)
### In the model specification function, input U and zeta, then:
# Q <- t(U) %*% U
# Zeta <- diag(zeta)
# Sigma <- Zeta %*% Q %*% Zeta
# Sigma.prior <- dsiw(Q, nu=Data$K+1, S=diag(Data$K), zeta, mu=0, delta=1)
### Examples
x <- dsiw(diag(3), 4, diag(3), runif(3), rep(0,3), rep(1,3), log=TRUE)
x <- rsiw(4, diag(3), rep(0,3), rep(1,3))
}
\keyword{Distribution}
