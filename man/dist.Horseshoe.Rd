\name{dist.Horseshoe}
\alias{dhs}
\alias{rhs}
\title{Horseshoe Distribution}
\description{
  This is the density function and random generation from the horseshoe
  distribution.
}
\usage{
dhs(x, lambda, tau, log=FALSE)
rhs(n, lambda, tau)
}
\arguments{
  \item{n}{This is the number of draws from the distribution.}
  \item{x}{This is a location vector at which to evaluate density.}
  \item{lambda}{This vector is a positive-only local parameter
    \eqn{\lambda}{lambda}.}
  \item{tau}{This scalar is a positive-only global parameter
    \eqn{\tau}{tau}.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
    \item Application: Multivariate Scale Mixture
    \item Density: (see below)
    \item Inventor: Carvalho et al. (2008)
    \item Notation 1: \eqn{\theta \sim \mathcal{HS}(\lambda,
      \tau)}{theta ~ HS(lambda, tau)}
    \item Notation 2: \eqn{p(\theta) = \mathcal{HS}(\theta | \lambda,
      \tau)}{p(theta) = HS(theta | lambda, tau)}
    \item Parameter 1: local scale \eqn{\lambda > 0}{lambda > 0}
    \item Parameter 2: global scale \eqn{\tau > 0}{tau > 0}
    \item Mean: \eqn{E(\theta)}{E(theta)}
    \item Variance: \eqn{var(\theta)}{var(theta)}
    \item Mode: \eqn{mode(\theta)}{mode(theta)}
}

The horseshoe distribution (Carvalho et al., 2008) is a heavy-tailed
mixture distribution that can be considered a variance mixture,
and it is in the family of multivariate scale mixtures of normals.

The horseshoe distribution was proposed as a prior distribution, and
recommended as a default choice for shrinkage priors in the presence of
sparsity. Horseshoe priors are most appropriate in large-p models where
dimension reduction is necessary to avoid overly complex models that
predict poorly, and also perform well in estimating a sparse covariance
matrix via Cholesky decomposition (Carvalho et al., 2009).

When the number of parameters in variable selection is assumed to be
sparse, meaning that most elements are zero or nearly zero, a horseshoe
prior is a desirable alternative to the Laplace-distributed parameters
in the LASSO, or the parameterization in ridge regression. When the true
value is far from zero, the horseshoe prior leaves the parameter
unshrunk. Yet, the horseshoe prior is accurate in shrinking parameters
that are truly zero or near-zero. Parameters near zero are shrunk more
than parameters far from zero. Therefore, parameters far from zero
experience less shrinkage and are closer to their true values. The
horseshoe prior is valuable in discriminating signal from noise.

By replacing the Laplace-distributed parameters in LASSO with
horseshoe-distributed parameters and including a global scale, the
result is called horseshoe regression.
}
\value{
  \code{dhs} gives the density and
  \code{rhs} generates random deviates.
}
\references{
  Carvalho, C.M., Polson, N.G., and Scott, J.G. (2008). "The Horseshoe
  Estimator for Sparse Signals". \emph{Discussion Paper 2008-31}. Duke
  University Department of Statistical Science.

  Carvalho, C.M., Polson, N.G., and Scott, J.G. (2009). "Handling
  Sparsity via the Horseshoe". \emph{Journal of Machine Learning
  Research}, 5, p. 73--80.
}
\seealso{
  \code{\link{dlaplace}}
}
\examples{
library(LaplacesDemon)
x <- rnorm(100)
lambda <- rhalfcauchy(100, 5)
tau <- 5
x <- dhs(x, lambda, tau, log=TRUE)
x <- rhs(100, lambda=lambda, tau=tau)
plot(density(x))
}
\keyword{Distribution}
