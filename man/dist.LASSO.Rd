\name{dist.LASSO}
\alias{dlasso}
\alias{rlasso}
\title{LASSO Distribution}
\description{
  These functions provide the density and random generation for the
  Bayesian LASSO prior distribution.
}
\usage{
dlasso(x, sigma, tau, lambda, a=1, b=1, log=FALSE)
rlasso(n, sigma, tau, lambda, a=1, b=1)
}
\arguments{
  \item{x}{This is a location vector of length \eqn{J} at which to
    evaluate density.}
  \item{n}{This is the number of observations, which must be a
    positive integer that has length 1.}
  \item{sigma}{This is a positive-only scalar hyperparameter
    \eqn{\sigma}{sigma}, which is also the residual standard deviation.}
  \item{tau}{This is a positive-only vector of hyperparameters,
    \eqn{\tau}{tau}, of length \eqn{J} regarding local sparsity.}
  \item{lambda}{This is a positive-only scalar hyperhyperparameter,
    \eqn{\lambda}{lambda}, of global sparsity.}
  \item{a, b}{These are positive-only scalar hyperhyperhyperparameters
    for gamma distributed \eqn{\lambda}{lambda}.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
    \item Application: Multivariate Scale Mixture
    \item Density: \eqn{p(\theta) \sim \mathcal{N}_k(0, \sigma^2
      diag(\tau^2))(\frac{1}{sigma^2})
      \mathcal{EXP}(\frac{\lambda^2}{2}) \mathcal{G}(a,b)}{p(theta) ~
      N[k](0, sigma^2 diag(tau^2))(1/sigma^2) EXP(lambda^2/2) G(a,b)}    
    \item Inventor: Parks and Casella (2008)
    \item Notation 1: \eqn{\theta \sim \mathcal{LASSO}(\sigma, \tau,
      \lambda, a, b)}{theta ~ LASSO(sigma, tau, lambda, a, b)}
    \item Notation 2: \eqn{p(\theta) = \mathcal{LASSO}(\theta | \sigma,
      \tau, \lambda, a, b)}{p(theta) = LASSO(theta | sigma, tau,
      lambda, a, b)}
    \item Parameter 1: hyperparameter global scale
      \eqn{\sigma > 0}{sigma > 0}
    \item Parameter 2: hyperparameter local scale
      \eqn{\tau > 0}{tau > 0}
    \item Parameter 3: hyperhyperparameter global scale
      \eqn{\lambda > 0}{lambda > 0}
    \item Parameter 4: hyperhyperhyperparameter scale \eqn{a > 0}
    \item Parameter 5: hyperhyperhyperparameter scale \eqn{b > 0}
    \item Mean: \eqn{E(\theta)}{E(theta)}
    \item Variance: 
    \item Mode: 
}

The Bayesian LASSO distribution (Parks and Casella, 2008) is a
heavy-tailed mixture distribution that can be considered a variance
mixture, and it is in the family of multivariate scale mixtures of
normals.

The LASSO distribution was proposed as a prior distribution, as a
Bayesian version of the frequentist LASSO, introduced by Tibshirani
(1996). It is applied as a shrinkage prior in the presence of sparsity
for \eqn{J} regression effects. LASSO priors are most appropriate in
large-dimensional models where dimension reduction is necessary to
avoid overly complex models that predict poorly.

The Bayesian LASSO results in regression effects that are a compromise
between regression effects in the frequentist LASSO and ridge
regression. The Bayesian LASSO applies more shrinkage to weak regression
effects than ridge regression.

The Bayesian LASSO is an alternative to horseshoe regression and ridge
regression.
}
\value{
  \code{dlasso} gives the density and
  \code{rlasso} generates random deviates.
}
\references{
  Park, T. and Casella, G. (2008). "The Bayesian Lasso". \emph{Journal
  of the American Statistical Association}, 103, p. 672--680.

  Tibshirani, R. (1996). "Regression Shrinkage and Selection via the
  Lasso". \emph{Journal of the Royal Statistical Society}, Series B, 58,
  p. 267--288.
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
x <- dlasso(x, sigma, tau, lambda, log=TRUE)
x <- rlasso(length(tau), sigma, tau, lambda)
}
\keyword{Distribution}
