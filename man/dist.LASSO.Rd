\name{dist.LASSO}
\alias{dlasso}
\title{LASSO Distribution}
\description{
  This is the density function for the Bayesian LASSO prior distribution.
}
\usage{
dlasso(x, sigma, tau, lambda, a, b, v, w, log=FALSE)
}
\arguments{
  \item{x}{This is a location vector of length \eqn{J} at which to
    evaluate density.}
  \item{sigma}{This is a positive-only scalar parameter
    \eqn{\sigma}{sigma}.}
  \item{tau}{This is a positive-only vector, \eqn{\tau}{tau}, of length
    \eqn{J} of local sparsity.}
  \item{lambda}{This isa positive-only scalar, \eqn{\lambda}{lambda}, of
    global sparsity.}
  \item{a, b}{These are scalar hyperparameters for gamma distributed
    \eqn{\lambda}{lambda}.}
  \item{v, w}{These are scalar hyperparameters for inverse-gamma
    distributed \eqn{\sigma}{sigma}.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
    \item Application: Multivariate Scale Mixture
    \item Density: (see below)
    \item Inventor: Parks and Casella (2008)
    \item Notation 1: \eqn{\theta \sim \mathcal{LASSO}(\sigma, \tau,
      \lambda, a, b, v, w)}{theta ~ LASSO(sigma, tau, lambda, a, b, v, w)}
    \item Notation 2: \eqn{p(\theta) = \mathcal{LASSO}(\theta | \sigma,
      \tau, \lambda, a, b, v, w)}{p(theta) = LASSO(theta | sigma, tau,
      lambda, a, b, v, w)}
    \item Parameter 1: global scale \eqn{\sigma > 0}{sigma > 0}
    \item Parameter 2: local scale \eqn{\tau > 0}{tau > 0}
    \item Parameter 3: global scale \eqn{\lambda > 0}{lambda > 0}
    \item Parameter 4: hyperparameter scale \eqn{a > 0}
    \item Parameter 5: hyperparameter scale \eqn{b > 0}
    \item Parameter 6: hyperparameter scale \eqn{v > 0}
    \item Parameter 7: hyperparameter scale \eqn{w > 0}
    \item Mean: \eqn{E(\theta)}{E(theta)}
    \item Variance: 
    \item Mode: 
}

The Bayesian LASSO distribution (Parks and Casella, 2008) is a
heavy-tailed mixture distribution that can be considered a variance
mixture, and it is in the family of multivariate scale mixtures of
normals.

The LASSO distribution was proposed as a prior distribution, as a
Bayesian version of the frequentist LASSO. It is applied as a shrinkage
prior in the presence of sparsity for \eqn{J} regression effects. LASSO
priors are most appropriate in large-dimensional models where dimension
reduction is necessary to avoid overly complex models that predict
poorly.

The Bayesian LASSO is an alternative to horseshoe regression and ridge
regression.
}
\value{
  \code{dlasso} gives the density
}
\references{
  Park, T. and Casella, G. (2008). "The Bayesian Lasso", \emph{Journal
  of the American Statistical Association}, 103, p. 672--680.
}
\seealso{
  \code{\link{dhs}}
}
\examples{
library(LaplacesDemon)
x <- rnorm(100)
sigma <- rhalfcauchy(1, 5)
tau <- rhalfcauchy(100, 5)
lambda <- rhalfcauchy(1, 5)
a <- rhalfcauchy(1, 5)
b <- rhalfcauchy(1, 5)
v <- rhalfcauchy(1, 5)
w <- rhalfcauchy(1, 5)
x <- dlasso(x, sigma, tau, lambda, a, b, v, w, log=TRUE)
}
\keyword{Distribution}
