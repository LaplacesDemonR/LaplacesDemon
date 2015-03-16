\name{dist.Matrix.Normal}
\alias{dmatrixnorm}
\alias{rmatrixnorm}
\title{Matrix Normal Distribution}
\description{
  These functions provide the density and random number generation for
  the matrix normal distribution.
}
\usage{
dmatrixnorm(X, M, U, V, log=FALSE) 
rmatrixnorm(M, U, V)
}
\arguments{
  \item{X}{This is data or parameters in the form of a matrix with
    \eqn{n} rows and \eqn{k} columns.}
  \item{M}{This is mean matrix with \eqn{n} rows and \eqn{k} columns.}
  \item{U}{This is a \eqn{n \times n}{n x n} positive-definite scale
    matrix.}
  \item{V}{This is a \eqn{k \times k}{k x k} positive-definite scale
    matrix.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
       density is returned.}
}
\details{
  \itemize{
    \item Application: Continuous Multivariate Matrix
    \item Density: \eqn{p(\theta) =
      \frac{\exp(-0.5tr[V^{-1}(X-M)'U^{-1}(X-M)])}{(2\pi)^{nk/2}|V|^{n/2}|U|^{k/2}}}{p(theta) = exp(-0.5tr[V^(-1)(X-M)'U^(-1)(X-M)])/[(2*pi)^(nk/2)|V|^(n/2)|U|(k/2)]}
  \item Inventors: Unknown
  \item Notation 1: \eqn{\theta \sim \mathcal{MN}_{n \times k}(M, U,
    V)}{theta ~ MN[n x k](M, U, V)}
  \item Notation 2: \eqn{p(\theta) = \mathcal{MN}_{n \times k}(\theta |
    M, U, V)}{p(theta) = MN[n x k](theta | M, U, V)}
  \item Parameter 1: location \eqn{n \times k}{n x k} matrix \eqn{M}
  \item Parameter 2: positive-definite \eqn{n \times n}{n x n} scale matrix \eqn{U}
  \item Parameter 3: positive-definite \eqn{k \times k}{k x k} scale matrix \eqn{V}
  \item Mean: \eqn{E(\theta) = M}{E(theta) = M}
  \item Variance: Unknown
  \item Mode: Unknown
}

The matrix normal distribution is also called the matrix Gaussian,
matrix-variate normal, or matrix-variate Gaussian distribution. It is a
generalization of the multivariate normal distribution to matrix-valued
random variables.

An example of the use of a matrix normal distribution is multivariate
regression, in which there is a \eqn{j \times k}{j x k} matrix of
regression effects of \eqn{j} predictors for \eqn{k} dependent
variables. For univariate regression, having only one dependent
variable, the \eqn{j} regression effects may be multivariate normally
distributed. For multivariate regression, this multivariate normal
distribution may be extended to a matrix normal distribution to account
for relationships of the regression effects across \eqn{k} dependent
variables. In this example, the matrix normal distribution is the
conjugate prior distribution for these regression effects.

The matrix normal distribution has two covariance matrices, one for the
rows and one for the columns. When \eqn{U} is diagonal, the rows are
independent. When \eqn{V} is diagonal, the columns are independent.
}
\value{
  \code{dmatrixnorm} gives the density and 
  \code{rmatrixnorm} generates random deviates.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{dinvmatrixgamma}},
  \code{\link{dmatrixgamma}}, and
  \code{\link{dmvn}}.
}
\examples{
library(LaplacesDemon)
N <- 10
K <- 4
U <- as.positive.definite(matrix(rnorm(N*N),N,N))
V <- as.positive.definite(matrix(rnorm(K*K),K,K))
x <- dmatrixnorm(matrix(0,N,K), matrix(0,N,K), U, V)
X <- rmatrixnorm(matrix(0,N,K), U, V)
joint.density.plot(X[,1], X[,2], color=TRUE)
}
\keyword{Distribution}