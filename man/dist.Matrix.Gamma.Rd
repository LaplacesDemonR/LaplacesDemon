\name{dist.Matrix.Gamma}
\alias{dmatrixgamma}
\title{Matrix Gamma Distribution}
\description{
  This function provides the density for the matrix gamma distribution.
}
\usage{
dmatrixgamma(X, alpha, beta, Sigma, log=FALSE)
}
\arguments{
  \item{X}{This is a \eqn{k \times k}{k x k} positive-definite precision
    matrix.}
  \item{alpha}{This is a scalar shape parameter (the degrees of freedom),
    \eqn{\alpha}{alpha}.}
  \item{beta}{This is a scalar, positive-only scale parameter,
    \eqn{\beta}{beta}.}
  \item{Sigma}{This is a \eqn{k \times k}{k x k} positive-definite scale
    matrix.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
       density is returned.}
}
\details{
  \itemize{
    \item Application: Continuous Multivariate Matrix
    \item Density: \eqn{p(\theta) = \frac{|\Sigma|^{-\alpha}}{\beta^{k
	  \alpha} \Gamma_k(\alpha)}
      |\theta|^{\alpha-(k+1)/2}\exp(tr(-\frac{1}{\beta}\Sigma^{-1}\theta))}{p(theta) = {|Sigma|^(-alpha) / [beta^(k alpha) Gamma[k](alpha)]} |theta|^[alpha-(k+1)/2] exp(tr(-(1/beta)Sigma^(-1)theta))}
  \item Inventors: Unknown
  \item Notation 1: \eqn{\theta \sim \mathcal{MG}_k(\alpha, \beta,
    \Sigma)}{theta ~ MG[k](alpha, beta, Sigma)}
  \item Notation 2: \eqn{p(\theta) = \mathcal{MG}_k(\theta | \alpha,
    \beta, \Sigma)}{p(theta) = MG[k](theta | alpha, beta, Sigma)}
  \item Parameter 1: shape \eqn{\alpha > 2}{alpha > 2}
  \item Parameter 2: scale \eqn{\beta > 0}{beta > 0}
  \item Parameter 3: positive-definite \eqn{k \times k}{k x k} scale matrix \eqn{\Sigma}{Sigma}
  \item Mean: 
  \item Variance: 
  \item Mode: 
}

The matrix gamma (MG), also called the matrix-variate gamma,
distribution is a generalization of the gamma distribution to
positive-definite matrices. It is a more general and flexible version of
the Wishart distribution (\code{\link{dwishart}}), and is a conjugate
prior of the precision matrix of a multivariate normal distribution
(\code{\link{dmvnp}}) and matrix normal distribution
(\code{\link{dmatrixnorm}}).

The compound distribution resulting from compounding a matrix normal
with a matrix gamma prior over the precision matrix is a generalized
matrix t-distribution.

The matrix gamma distribution is identical to the Wishart distribution
when \eqn{\alpha = \nu / 2}{alpha = nu / 2} and
\eqn{\beta = 2}{beta = 2}.
}
\value{
  \code{dmatrixgamma} gives the density.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{dgamma}}
  \code{\link{dmatrixnorm}},
  \code{\link{dmvnp}}, and
  \code{\link{dwishart}}
}
\examples{
library(LaplacesDemon)
k <- 10
dmatrixgamma(X=diag(k), alpha=(k+1)/2, beta=2, Sigma=diag(k), log=TRUE)
dwishart(Omega=diag(k), nu=k+1, S=diag(k), log=TRUE)
}
\keyword{Distribution}