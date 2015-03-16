\name{dist.Inverse.Matrix.Gamma}
\alias{dinvmatrixgamma}
\title{Inverse Matrix Gamma Distribution}
\description{
  This function provides the density for the inverse matrix gamma
  distribution.
}
\usage{
dinvmatrixgamma(X, alpha, beta, Psi, log=FALSE)
}
\arguments{
  \item{X}{This is a \eqn{k \times k}{k x k} positive-definite
    covariance matrix.}
  \item{alpha}{This is a scalar shape parameter (the degrees of freedom),
    \eqn{\alpha}{alpha}.}
  \item{beta}{This is a scalar, positive-only scale parameter,
    \eqn{\beta}{beta}.}
  \item{Psi}{This is a \eqn{k \times k}{k x k} positive-definite scale
    matrix.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
       density is returned.}
}
\details{
  \itemize{
    \item Application: Continuous Multivariate Matrix
    \item Density: \eqn{p(\theta) = \frac{|\Psi|^\alpha}{\beta^{k
	  \alpha} \Gamma_k(\alpha)}
      |\theta|^{-\alpha-(k+1)/2}\exp(tr(-\frac{1}{\beta}\Psi\theta^{-1}))}{p(theta) = {|Psi|^alpha / [beta^(k alpha) Gamma[k](alpha)]} |theta|^[-alpha-(k+1)/2] exp(tr(-(1/beta)Psi theta^(-1)))}
  \item Inventors: Unknown
  \item Notation 1: \eqn{\theta \sim \mathcal{IMG}_k(\alpha, \beta,
    \Psi)}{theta ~ IMG[k](alpha, beta, Psi)}
  \item Notation 2: \eqn{p(\theta) = \mathcal{IMG}_k(\theta | \alpha,
    \beta, \Psi)}{p(theta) = IMG[k](theta | alpha, beta, Psi)}
  \item Parameter 1: shape \eqn{\alpha > 2}{alpha > 2}
  \item Parameter 2: scale \eqn{\beta > 0}{beta > 0}
  \item Parameter 3: positive-definite \eqn{k \times k}{k x k} scale
    matrix \eqn{\Psi}{Psi}
  \item Mean: 
  \item Variance: 
  \item Mode: 
}

The inverse matrix gamma (IMG), also called the inverse matrix-variate
gamma, distribution is a generalization of the inverse gamma
distribution to positive-definite matrices. It is a more general and
flexible version of the inverse Wishart distribution
(\code{\link{dinvwishart}}), and is a conjugate prior of the covariance
matrix of a multivariate normal distribution (\code{\link{dmvn}}) and
matrix normal distribution (\code{\link{dmatrixnorm}}).

The compound distribution resulting from compounding a matrix normal
with an inverse matrix gamma prior over the covariance matrix is a
generalized matrix t-distribution.

The inverse matrix gamma distribution is identical to the inverse
Wishart distribution when \eqn{\alpha = \nu / 2}{alpha = nu / 2} and
\eqn{\beta = 2}{beta = 2}.
}
\value{
  \code{dinvmatrixgamma} gives the density.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{dinvgamma}}
  \code{\link{dmatrixnorm}},
  \code{\link{dmvn}}, and
  \code{\link{dinvwishart}}
}
\examples{
library(LaplacesDemon)
k <- 10
dinvmatrixgamma(X=diag(k), alpha=(k+1)/2, beta=2, Psi=diag(k), log=TRUE)
dinvwishart(Sigma=diag(k), nu=k+1, S=diag(k), log=TRUE)
}
\keyword{Distribution}