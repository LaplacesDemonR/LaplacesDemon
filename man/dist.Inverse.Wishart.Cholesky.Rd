\name{dist.Inverse.Wishart.Cholesky}
\alias{dinvwishartc}
\alias{rinvwishartc}
\title{Inverse Wishart Distribution: Cholesky Parameterization}
\description{
  These functions provide the density and random number generation
  for the inverse Wishart distribution with the Cholesky
  parameterization.
}
\usage{
   dinvwishartc(U, nu, S, log=FALSE)
   rinvwishartc(nu, S)
}
\arguments{
  \item{U}{This is the upper-triangular \eqn{k \times k}{k x k} matrix
    for the Cholesky factor \eqn{\textbf{U}}{U} of covariance matrix
    \eqn{\Sigma}{Sigma}.}
  \item{nu}{This is the scalar degrees of freedom, \eqn{\nu}{nu}.}
  \item{S}{This is the symmetric, positive-semidefinite
    \eqn{k \times k}{k x k} scale matrix \eqn{\textbf{S}}{S}.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
    \item Application: Continuous Multivariate
    \item Density: \eqn{p(\theta) = (2^{\nu k/2} \pi^{k(k-1)/4}
      \prod^k_{i=1} \Gamma(\frac{\nu+1-i}{2}))^{-1} |\textbf{S}|^{nu/2}
      |\Omega|^{-(nu-k-1)/2} \exp(-\frac{1}{2} tr(\textbf{S}
      \Omega^{-1}))}{p(theta) = (2^(nu*k/2) * pi^(k(k-1)/4) *
      [Gamma((nu+1-i)/2) * ... * Gamma((nu+1-k)/2)])^(-1) * |S|^(nu/2) *
      |Omega|^(-(nu-k-1)/2) * exp(-(1/2) * tr(S Omega^(-1)))}
    \item Inventor: John Wishart (1928)
    \item Notation 1: \eqn{\Sigma \sim
      \mathcal{W}^{-1}_{\nu}(\textbf{S}^{-1})}{Sigma ~
      W^(-1)[nu](S^(-1))}
    \item Notation 2: \eqn{p(\Sigma) = \mathcal{W}^{-1}_{\nu}(\Sigma |
      \textbf{S}^{-1})}{p(Sigma) = W^-1[nu](Sigma | S^(-1))}
    \item Parameter 1: degrees of freedom \eqn{\nu}{nu}
    \item Parameter 2: symmetric, positive-semidefinite
      \eqn{k \times k}{k x k} scale matrix \eqn{\textbf{S}}{S}
    \item Mean: \eqn{E(\Sigma) = \frac{\textbf{S}}{\nu - k -
	1}}{E(Sigma) = S / (nu - k - 1)}
    \item Variance: 
    \item Mode: \eqn{mode(\Sigma) = \frac{\textbf{S}}{\nu + k +
	1}}{mode(Sigma) = S / (nu + k + 1)}
}

The inverse Wishart distribution is a probability distribution defined on
real-valued, symmetric, positive-definite matrices, and is used as the
conjugate prior for the covariance matrix, \eqn{\Sigma}{Sigma}, of a
multivariate normal distribution. In this parameterization,
\eqn{\Sigma}{Sigma} has been decomposed to the upper-triangular Cholesky
factor \eqn{\textbf{U}}{U}, as per \code{\link{chol}}. The
inverse-Wishart density is always finite, and the integral is always
finite. A degenerate form occurs when \eqn{\nu < k}{nu < k}.

In practice, \eqn{\textbf{U}}{U} is fully unconstrained for proposals
when its diagonal is log-transformed. The diagonal is exponentiated
after a proposal and before other calculations. Overall, the
Cholesky parameterization is faster than the traditional
parameterization. Compared with \code{dinvwishart}, \code{dinvwishartc}
must additionally matrix-multiply the Cholesky back to the covariance
matrix, but it does not have to check for or correct the covariance
matrix to positive-semidefiniteness, which overall is slower. Compared
with \code{rinvwishart}, \code{rinvwishartc} must additionally
calculate a Cholesky decomposition, and is therefore slower.

The inverse Wishart prior lacks flexibility, having only one parameter,
\eqn{\nu}{nu}, to control the variability for all \eqn{k(k + 1)/2}
elements. Popular choices for the scale matrix \eqn{\textbf{S}}{S}
include an identity matrix or sample covariance matrix. When the model
sample size is small, the specification of the scale matrix can be
influential.

One of many alternatives is to use hierarchical priors,
in which the main diagonal of the (identity) scale matrix and the
degrees of freedom are treated as unknowns (Bouriga and Feron, 2011;
Daniels and Kass, 1999). A hierarchical inverse Wishart prior provides
shrinkage toward diagonality. Another alternative is to abandon the
inverse Wishart distribution altogether for the more flexible method of
Barnard et al. (2000) or the horseshoe distribution (\code{\link{dhs}})
for sparse covariance matrices.
}
\value{
  \code{dinvwishartc} gives the density and
  \code{rinvwishartc} generates random deviates.
}
\references{
  Barnard, J., McCulloch, R., and Meng, X. (2000). "Modeling Covariance
  Matrices in Terms of Standard Deviations and Correlations, with
  Application to Shrinkage". \emph{Statistica Sinica}, 10, p. 1281--1311.
  
  Bouriga, M. and Feron, O. (2011). "Estimation of Covariance Matrices
  Based on Hierarchical Inverse-Wishart Priors". URL:
  http://www.citebase.org/abstract?id=oai:arXiv.org:1106.3203.
  
  Daniels, M., and Kass, R. (1999). "Nonconjugate Bayesian Estimation of
  Covariance Matrices and its use in Hierarchical Models". \emph{Journal
  of the American Statistical Association}, 94(448), p. 1254--1263.

  Wishart, J. (1928). "The Generalised Product Moment Distribution in
  Samples from a Normal Multivariate Population". \emph{Biometrika},
  20A(1-2), p. 32--52.
}
\seealso{
  \code{\link{chol}},
  \code{\link{Cov2Prec}},
  \code{\link{dhs}},
  \code{\link{dmvn}},
  \code{\link{dmvnc}},
  \code{\link{dmvtc}},
  \code{\link{dwishart}}, and
  \code{\link{dwishartc}}.
}
\examples{
library(LaplacesDemon)
Sigma <- matrix(c(2,-.3,-.3,4),2,2)
U <- chol(Sigma)
x <- dinvwishartc(U, 3, matrix(c(1,.1,.1,1),2,2))
x <- rinvwishartc(3, matrix(c(1,.1,.1,1),2,2))
}
\keyword{Distribution}
