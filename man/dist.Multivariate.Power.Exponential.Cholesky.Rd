\name{dist.Multivariate.Power.Exponential.Cholesky}
\alias{dmvpec}
\alias{rmvpec}
\title{Multivariate Power Exponential Distribution: Cholesky Parameterization}
\description{
  These functions provide the density and random number generation for
  the multivariate power exponential distribution, given the Cholesky
  parameterization.
}
\usage{
dmvpec(x=c(0,0), mu=c(0,0), U, kappa=1, log=FALSE)
rmvpec(n, mu=c(0,0), U, kappa=1)
}
\arguments{
  \item{x}{This is data or parameters in the form of a vector of length
    \eqn{k} or a matrix with \eqn{k} columns.}
  \item{n}{This is the number of random draws.}
  \item{mu}{This is mean vector \eqn{\mu}{mu} with length \eqn{k} or
    matrix with \eqn{k} columns.}
  \item{U}{This is the \eqn{k \times k}{k x k} upper-triangular matrix
    that is Cholesky factor \eqn{\textbf{U}}{U} of covariance matrix
    \eqn{\Sigma}{Sigma}.}
  \item{kappa}{This is the kurtosis parameter, \eqn{\kappa}{kappa}, and
    must be positive.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
  \item Application: Continuous Multivariate
  \item Density: \deqn{p(\theta) = \frac{k\Gamma(k/2)}{\pi^{k/2}
      \sqrt{|\Sigma|} \Gamma(1 + k/(2\kappa)) 2^{1 + k/(2\kappa)}}
    \exp(-\frac{1}{2}(\theta-\mu)^T \Sigma
    (\theta-\mu))^{\kappa}}{p(theta) = ((k*Gamma(k/2)) / (pi^(k/2) * 
    sqrt(|Sigma|) * Gamma(1 + k/(2*kappa)) * 2^(1 + k/(2*kappa)))) *
          exp(-(1/2)*(theta-mu)^T Sigma (theta-mu))^kappa}
  \item Inventor: Gomez, Gomez-Villegas, and Marin (1998)
  \item Notation 1: \eqn{\theta \sim \mathcal{MPE}(\mu, \Sigma,
    \kappa)}{theta ~ MPE(mu, Sigma, kappa)}
  \item Notation 2: \eqn{\theta \sim \mathcal{PE}_k(\mu, \Sigma,
    \kappa)}{theta ~ PE[k](mu, Sigma, kappa)}
  \item Notation 3: \eqn{p(\theta) = \mathcal{MPE}(\theta | \mu, \Sigma,
    \kappa)}{p(theta) = MPE(theta | mu, Sigma, kappa)}
  \item Notation 4: \eqn{p(\theta) = \mathcal{PE}_k(\theta | \mu,
    \Sigma, \kappa)}{p(theta) = PE[k](theta | mu, Sigma, kappa)}
  \item Parameter 1: location vector \eqn{\mu}{mu}
  \item Parameter 2: positive-definite \eqn{k \times k}{k x k}
  covariance matrix \eqn{\Sigma}{Sigma}
  \item Parameter 3: kurtosis parameter \eqn{\kappa}{kappa}
  \item Mean: \eqn{E(\theta) = }{E(theta) = }
  \item Variance: \eqn{var(\theta) =}{var(theta) = }
  \item Mode: \eqn{mode(\theta) = }{mode(theta) = }
}

The multivariate power exponential distribution, or multivariate
exponential power distribution, is a multidimensional extension of
the one-dimensional or univariate power exponential distribution.
Gomez-Villegas (1998) and Sanchez-Manzano et al. (2002) proposed
multivariate and matrix generalizations of the PE family of
distributions and studied their properties in relation to multivariate
Elliptically Contoured (EC) distributions.

The multivariate power exponential distribution includes the
multivariate normal distribution (\eqn{\kappa = 1}{kappa = 1}) and
multivariate Laplace distribution (\eqn{\kappa = 0.5}{kappa = 0.5}) as
special cases, depending on the kurtosis or \eqn{\kappa}{kappa}
parameter. A multivariate uniform occurs as
\eqn{\kappa \rightarrow \infty}{kappa -> infinity}.

If the goal is to use a multivariate Laplace distribution, the
\code{dmvlc} function will perform faster and more accurately.

In practice, \eqn{\textbf{U}}{U} is fully unconstrained for proposals
when its diagonal is log-transformed. The diagonal is exponentiated
after a proposal and before other calculations. Overall, the Cholesky
parameterization is faster than the traditional parameterization.
Compared with \code{dmvpe}, \code{dmvpec} must additionally
matrix-multiply the Cholesky back to the covariance matrix, but it
does not have to check for or correct the covariance matrix to
positive-definiteness, which overall is slower. Compared with
\code{rmvpe}, \code{rmvpec} is faster because the Cholesky decomposition
has already been performed.

The \code{rmvpec} function is a modified form of the rmvpowerexp function
in the MNM package.
}
\value{
  \code{dmvpec} gives the density and
  \code{rmvpec} generates random deviates.
}
\references{
  Gomez, E., Gomez-Villegas, M.A., and Marin, J.M. (1998). "A
  Multivariate Generalization of the Power Exponential Family of
  Distributions". \emph{Communications in Statistics-Theory and
  Methods}, 27(3), p. 589--600.

  Sanchez-Manzano, E.G., Gomez-Villegas, M.A., and Marn-Diazaraque,
  J.M. (2002). "A Matrix Variate Generalization of the Power Exponential
  Family of Distributions". \emph{Communications in Statistics, Part A -
  Theory and Methods} [Split from: J(CommStat)], 31(12), p. 2167--2182.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{chol}},
  \code{\link{dlaplace}},
  \code{\link{dmvlc}},
  \code{\link{dmvnc}},
  \code{\link{dmvnpc}},
  \code{\link{dnorm}},
  \code{\link{dnormp}},
  \code{\link{dnormv}}, and
  \code{\link{dpe}}.
}
\examples{
library(LaplacesDemon)
n <- 100
k <- 3
x <- matrix(runif(n*k),n,k)
mu <- matrix(runif(n*k),n,k)
Sigma <- diag(k)
U <- chol(Sigma)
dmvpec(x, mu, U, kappa=1)
X <- rmvpec(n, mu, U, kappa=1)
joint.density.plot(X[,1], X[,2], color=TRUE)
}
\keyword{Distribution}