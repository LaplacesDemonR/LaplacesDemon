\name{dist.HuangWand}
\alias{dhuangwand}
\alias{dhuangwandc}
\alias{rhuangwand}
\alias{rhuangwandc}
\title{Huang-Wand Distribution}
\description{
  These are the density and random generation functions for the
  Huang-Wand prior distribution for a covariance matrix.
}
\usage{
dhuangwand(x, nu=2, a, A, log=FALSE)
dhuangwandc(x, nu=2, a, A, log=FALSE)
rhuangwand(nu=2, a, A)
rhuangwandc(nu=2, a, A)
}
\arguments{
  \item{x}{This is a \eqn{k \times k}{k x k} positive-definite
    covariance matrix \eqn{\Sigma}{Sigma} for \code{dhuangwand}, or the
    Cholesky factor \eqn{\textbf{U}}{U} of the covariance matrix for
    \code{dhuangwandc}.}
  \item{nu}{This is a scalar degrees of freedom parameter
    \eqn{\nu}{nu}. The default is \code{nu=2}, which is an
    uninformative prior, resulting in marginal uniform distributions
    on the correlation matrix.}
  \item{a}{This is a positive-only vector of scale parameters
    \eqn{a} of length \eqn{k}.}
  \item{A}{This is a positive-only vector of of scale hyperparameters
    \eqn{A} of length \eqn{k}. Larger values result in a more
    uninformative prior. A default, uninformative prior is
    \code{A=rep(1e6,k)}.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
    \item Application: Continuous Multivariate
    \item Density: \eqn{p(\theta) =
      \mathcal{W}^{-1}_{\nu+k-1}(2 \nu diag(1/a))
      \mathcal{G}^{-1}(1/2, 1/A^2)}{p(theta) =
      W^(-1)[nu+k-1](2*nu*diag(1/a)) G^(-1)(1/2, 1/A^2)}
    \item Inventor: Huang and Wand (2013)
    \item Notation 1: \eqn{\theta \sim \mathcal{HW}_\nu(\textbf{a},
      \textbf{A})}{theta ~ HW[nu](a, A)}
    \item Notation 2: \eqn{p(\theta) \sim \mathcal{HW}_\nu(\theta |
      \textbf{a}, \textbf{A})}{p(theta) ~ HW[nu](theta | a, A)}
    \item Parameter 1: degrees of freedom \eqn{\nu}{nu}
    \item Parameter 2: scale \eqn{a > 0}
    \item Parameter 3: scale \eqn{A > 0}
    \item Mean: 
    \item Variance: 
    \item Mode: 
}

Huang and Wand (2013) proposed a prior distribution for a covariance
matrix that uses a hierarchical inverse Wishart. This is a more flexible
alternative to the inverse Wishart distribution, and the Huang-Wand
prior retains conjugacy. The Cholesky parameterization is also provided
here.

The Huang-Wand prior distribution alleviates two main limitations of an
inverse Wishart distribution. First, the uncertainty in the diagonal
variances of a covariance matrix that is inverse Wishart distributed is
represented with only one degrees of freedom parameter, which may be too
restrictive. The Huang-Wand prior overcomes this limitation. Second, the
inverse Wishart distribution imposes a dependency between variance and
correlation. The Huang-Wand prior lessens, but does not fully remove,
this dependency.

The standard deviations of a Huang-Wand distributed covariance matrix
are half-t distributed, as \eqn{\mathcal{HT}(\nu, \textbf{A})}{HT(nu,
  A)}. This is in accord with modern assumptions about distributions of
scale parameters, and is also useful for sparse covariance matrices.

The \code{rhuangwand} function allows either \code{a} or \code{A} to be
missing. When \code{a} is missing, the covariance matrix is generated
from the hyperparameters. When \code{A} is missing, the covariance
matrix is generated from the parameters.
}
\value{
  \code{dhuangwand} and \code{dhuangwandc} give the density, and
  \code{rhuangwand} and \code{rhuangwandc} generate random deviates.
}
\references{
  Huang, A., Wand, M., et al. (2013), "Simple Marginally Noninformative
  Prior Distributions for Covariance Matrices". \emph{Bayesian
  Analysis}, 8, p. 439--452.
}
\seealso{
  \code{\link{dhalft}} and
  \code{\link{dinvwishart}}
}
\examples{
library(LaplacesDemon)
dhuangwand(diag(3), nu=2, a=runif(3), A=rep(1e6,3), log=TRUE)
rhuangwand(nu=2, A=rep(1e6, 3)) #Missing a
rhuangwand(nu=2, a=runif(3)) #Missing A
}
\keyword{Distribution}
