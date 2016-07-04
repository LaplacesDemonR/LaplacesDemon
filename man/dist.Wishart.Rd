\name{dist.Wishart}
\alias{dwishart}
\alias{rwishart}
\title{Wishart Distribution}
\description{
  These functions provide the density and random number generation
  for the Wishart distribution.
}
\usage{
   dwishart(Omega, nu, S, log=FALSE)
   rwishart(nu, S)
}
\arguments{
  \item{Omega}{This is the symmetric, positive-definite \eqn{k \times
      k}{k x k} matrix \eqn{\Omega}{Omega}.}
  \item{nu}{This is the scalar degrees of freedom \eqn{\nu}{nu}.}
  \item{S}{This is the symmetric, positive-semidefinite, \eqn{k \times
      k}{k x k} scale matrix \eqn{\textbf{S}}{S}.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
    \item Application: Continuous Multivariate
    \item Density: \eqn{p(\theta) = (2^{\nu k/2} \pi^{k(k-1)/4}
      \prod^k_{i=1} \Gamma(\frac{\nu+1-i}{2}))^{-1} |\textbf{S}|^{-nu/2}
      |\Omega|^{(nu-k-1)/2} \exp(-\frac{1}{2} tr(\textbf{S}^{-1}
      \Omega))}{p(theta) = (2^(nu*k/2) * pi^(k(k-1)/4) *
      [Gamma((nu+1-i)/2) * ... * Gamma((nu+1-k)/2)])^(-1) * |S|^(-nu/2)
      * |Omega|^((nu-k-1)/2) * exp(-(1/2) * tr(S^(-1) Omega))}
    \item Inventor: John Wishart (1928)
    \item Notation 1: \eqn{\Omega \sim
      \mathcal{W}_{\nu}(\textbf{S})}{Omega ~ W[nu](S)}
    \item Notation 2: \eqn{p(\Omega) = \mathcal{W}_{\nu}(\Omega |
      \textbf{S})}{p(Omega) = W[nu](Omega | S)}
    \item Parameter 1: degrees of freedom \eqn{\nu \ge k}{nu >= k}
    \item Parameter 2: symmetric, positive-semidefinite
    \eqn{k \times k}{k x k} scale matrix \eqn{\textbf{S}}{S}
    \item Mean: \eqn{E(\Omega) = \nu \textbf{S}}{E(Omega) = nuS}
    \item Variance: \eqn{var(\Omega) = \nu (\textbf{S}^2_{i,j} +
      \textbf{S}_{i,i} \textbf{S}_{j,j})}{var(Omega) = nu(S[i,j]^2 +
      S[i,i]S[j,j])}
    \item Mode: \eqn{mode(\Omega) = (\nu - k - 1)
      \textbf{S}}{mode(Omega) = (nu-k-1)S}, for \eqn{\nu \ge k + 1}{nu
      >= k + 1}
}

The Wishart distribution is a generalization to multiple dimensions of
the chi-square distribution, or, in the case of non-integer degrees of
freedom, of the gamma distribution. However, the Wishart distribution is
not called the multivariate chi-squared distribution because the
marginal distribution of the off-diagonal elements is not chi-squared.

The Wishart is the conjugate prior distribution for the precision matrix
\eqn{\Omega}{Omega}, the inverse of which (covariance matrix
\eqn{\Sigma}{Sigma}) is used in a multivariate normal distribution.

The integral is finite when \eqn{\nu \ge k}{nu >= k}, where
\eqn{\nu}{nu} is the scalar degrees of freedom parameter, and \eqn{k} is
the dimension of scale matrix \eqn{\textbf{S}}{S}. The density is finite
when \eqn{\nu ge k + 1}{nu >= k + 1}, which is recommended.

The degrees of freedom, \eqn{\nu}{nu}, is equivalent to specifying a
prior sample size, indicating the confidence in \eqn{\textbf{S}}{S},
where \eqn{\textbf{S}}{S} is a prior guess at the order of covariance
matrix \eqn{\Sigma}{Sigma}. A flat prior distribution is
obtained as \eqn{\nu \rightarrow 0}{nu -> 0}.

When applicable, the alternative Cholesky parameterization should be
preferred. For more information, see \code{\link{dwishartc}}.

The Wishart prior lacks flexibility, having only one parameter,
\eqn{\nu}{nu}, to control the variability for all \eqn{k(k + 1)/2}
elements. Popular choices for the scale matrix \eqn{\textbf{S}}{S}
include an identity matrix or sample covariance matrix. When the model
sample size is small, the specification of the scale matrix can be
influential.

Although the related inverse Wishart distribution has a dependency
between variance and correlation, the Wishart distribution does not
have this dependency.

The matrix gamma (\code{\link{dmatrixgamma}}) distribution is a more
general version of the Wishart distribution, and the Yang-Berger
(\code{\link{dyangberger}}) distribution is an alterative that is a
least informative prior (LIP).
}
\value{
  \code{dwishart} gives the density and
  \code{rwishart} generates random deviates.
}
\references{
  Wishart, J. (1928). "The Generalised Product Moment Distribution in
  Samples from a Normal Multivariate Population". \emph{Biometrika},
  20A(1-2), p. 32--52.
}
\seealso{
  \code{\link{dchisq}},
  \code{\link{dgamma}},
  \code{\link{dinvwishart}},
  \code{\link{dmatrixgamma}},
  \code{\link{dmvnp}},
  \code{\link{dwishartc}},
  \code{\link{Prec2Cov}}, and
  \code{\link{dyangberger}}.
}
\examples{
library(LaplacesDemon)
x <- dwishart(matrix(c(2,-.3,-.3,4),2,2), 3, matrix(c(1,.1,.1,1),2,2))
x <- rwishart(3, matrix(c(1,.1,.1,1),2,2))
}
\keyword{Distribution}
