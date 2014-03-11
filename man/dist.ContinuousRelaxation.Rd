\name{dist.ContinuousRelaxation}
\alias{dcrmrf}
\alias{rcrmrf}
\title{Continuous Relaxation of a Markov Random Field Distribution}
\description{
  This is the density function and random generation from the continuous
  relaxation of a Markov random field (MRF) distribution.
}
\usage{
dcrmrf(x, alpha, Omega, log=FALSE)
rcrmrf(n, alpha, Omega)
}
\arguments{
  \item{x}{This is a vector of length \eqn{k}.}
  \item{n}{This is the number of random deviates to generate.}
  \item{alpha}{This is a vector of length \eqn{k} of shape parameters.}
  \item{Omega}{This is the \eqn{k \times k}{k x k} precision matrix
       \eqn{\Omega}{Omega}.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
    \item Application: Continuous Multivariate
    \item Density: \deqn{p(\theta) \propto \exp(-\frac{1}{2} \theta^T
      \Omega^{-1} \theta) \prod_i (1 + \exp(\theta_i +
      alpha_i))}{p(theta) = exp(-0.5 theta^T Omega^(-1) theta) prod i=1
      (1 + exp(theta[i] + alpha[i]))}
    \item Inventor: Zhang et al. (2012)
    \item Notation 1: \eqn{\theta \sim \mathcal{CRMRF}(\alpha,
      \Omega)}{theta ~ CRMRF(alpha, Omega)}
    \item Notation 2: \eqn{p(\theta) = \mathcal{CRMRF}(\theta | \alpha,
      \Omega)}{p(theta) = CRMRF(theta | alpha, Omega)}
    \item Parameter 1: shape vector \eqn{\alpha}{alpha}
    \item Parameter 2: positive-definite \eqn{k \times k}{k x k} matrix
      \eqn{\Omega}{Omega}
    \item Mean: \eqn{E(\theta)}{E(theta)}
    \item Variance: \eqn{var(\theta)}{var(theta)}
    \item Mode: \eqn{mode(\theta)}{mode(theta)}
}

It is often easier to solve or optimize a problem with continuous
variables rather than a problem that involves discrete variables. A
continuous variable may also have a gradient, contour, and curvature
that may be useful for optimization or sampling. Continuous MCMC
samplers are far more common.

Zhang et al. (2012) introduced a generalized form of the Gaussian
integral trick from statistical physics to transform a discrete variable
so that it may be estimated with continuous variables. An auxiliary
Gaussian variable is added to a discrete Markov random field (MRF) so
that discrete dependencies cancel out, allowing the discrete variable to
be summed away, and leaving a continuous problem. The resulting
continuous representation of the problem allows the model to be updated
with a continuous MCMC sampler, and may benefit from a MCMC sampler that
uses derivatives. Another advantage of continuous MCMC is that
stationarity of discrete Markov chains is problematic to assess.

A disadvantage of solving a discrete problem with continuous parameters
is that the continuous solution requires more parameters.
}
\value{
  \code{dcrmrf} gives the density and
  \code{rcrmrf} generates random deviates.
}
\references{
  Zhang, Y., Ghahramani, Z., Storkey, A.J., and Sutton, C.A. (2012).
  "Continuous Relaxations for Discrete Hamiltonian Monte Carlo".
  \emph{Advances in Neural Information Processing Systems}, 25,
  p. 3203--3211.
}
\seealso{
  \code{\link{dmvn}}
}
\examples{
library(LaplacesDemon)
x <- dcrmrf(rnorm(5), rnorm(5), diag(5))
x <- rcrmrf(10, rnorm(5), diag(5))
}
\keyword{Distribution}
