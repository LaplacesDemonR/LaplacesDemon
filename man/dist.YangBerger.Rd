\name{dist.YangBerger}
\alias{dyangberger}
\alias{dyangbergerc}
\title{Yang-Berger Distribution}
\description{
  This is the density function for the Yang-Berger prior distribution
  for a covariance matrix or precision matrix.
}
\usage{
dyangberger(x, log=FALSE)
dyangbergerc(x, log=FALSE)
}
\arguments{
  \item{x}{This is the \eqn{k \times k}{k x k} positive-definite
    covariance matrix or precision matrix for \code{dyangberger} or the
    Cholesky factor \eqn{\textbf{U}}{U} of the covariance matrix or
    precision matrix for \code{dyangbergerc}.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
    \item Application: Continuous Multivariate
    \item Density: \eqn{p(\theta) = \frac{1}{|\theta|^{\prod (d_j -
      d_{j-1})}}}{p(theta) = 1 / |theta|^(prod (d[j] - d[j-1]))},
      where \eqn{d} are increasing eigenvalues. See equation 13 in
      Yang and Berger (1994).
    \item Inventor: Yang and Berger (1994)
    \item Notation 1: \eqn{\theta \sim \mathcal{YB}}{p(theta) ~ YB}
    \item Mean: 
    \item Variance: 
    \item Mode: 
}

Yang and Berger (1994) derived a least informative prior (LIP) for a
covariance matrix or precision matrix. The Yang-Berger (YB) distribution
does not have any parameters. It is a reference prior for objective
Bayesian inference. The Cholesky parameterization is also provided here.

The YB prior distribution results in a proper posterior. It involves an
eigendecomposition of the covariance matrix or precision matrix. It is
difficult to interpret a model that uses the YB prior, due to a lack of
intuition regarding the relationship between eigenvalues and
correlations.

Compared to Jeffreys prior for a covariance matrix, this reference prior
encourages equal eigenvalues, and therefore results in a covariance
matrix or precision matrix with a better shrinkage of its
eigenstructure.
}
\value{
  \code{dyangberger} and \code{dyangbergerc} give the density.
}
\references{
  Yang, R. and Berger, J.O. (1994). "Estimation of a Covariance Matrix
  using the Reference Prior". \emph{Annals of Statistics}, 2,
  p. 1195-1211.
}
\seealso{
  \code{\link{dinvwishart}} and
  \code{\link{dwishart}}
}
\examples{
library(LaplacesDemon)
X <- matrix(c(1,0.8,0.8,1), 2, 2)
dyangberger(X, log=TRUE)
}
\keyword{Distribution}
