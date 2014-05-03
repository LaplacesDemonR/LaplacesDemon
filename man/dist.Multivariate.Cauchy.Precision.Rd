\name{dist.Multivariate.Cauchy.Precision}
\alias{dmvcp}
\alias{rmvcp}
\title{Multivariate Cauchy Distribution: Precision Parameterization}
\description{
  These functions provide the density and random number generation for
  the multivariate Cauchy distribution. These functions use the
  precision parameterization.
}
\usage{
dmvcp(x, mu, Omega, log=FALSE)
rmvcp(n=1, mu, Omega)
}
\arguments{
  \item{x}{This is either a vector of length \eqn{k} or a matrix with 
       a number of columns, \eqn{k}, equal to the number of columns in
       precision matrix \eqn{\Omega}{Omega}.}
  \item{n}{This is the number of random draws.}
  \item{mu}{This is a numeric vector representing the location parameter,
       \eqn{\mu}{mu} (the mean vector), of the multivariate
       distribution. It must be of length \eqn{k}, as defined above.}
  \item{Omega}{This is a \eqn{k \times k}{k x k} positive-definite
       precision matrix \eqn{\Omega}{Omega}.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
       density is returned.}
}
\details{
  \itemize{
  \item Application: Continuous Multivariate
  \item Density: \deqn{p(\theta) =
    \frac{\Gamma((1+k)/2)}{\Gamma(1/2)1^{k/2}\pi^{k/2}} |\Omega|^{1/2}
    (1 + (\theta-\mu)^T \Omega (\theta-\mu))^{-(1+k)/2}}{p(theta) =
    (Gamma((nu+k)/2) / (Gamma(1/2)*1^(k/2)*pi^(k/2))) * |Omega|^(1/2) *
    (1 + (theta-mu)^T Omega (theta-mu))^(-(1+k)/2)}
  \item Inventor: Unknown (to me, anyway)
  \item Notation 1: \eqn{\theta \sim \mathcal{MC}_k(\mu,
    \Omega^{-1})}{theta ~ MC[k](mu, Omega^(-1))}
  \item Notation 2: \eqn{p(\theta) = \mathcal{MC}_k(\theta | \mu,
    \Omega^{-1})}{p(theta) = MC[k](theta | mu, Omega^(-1))}
  \item Parameter 1: location vector \eqn{\mu}{mu}
  \item Parameter 2: positive-definite \eqn{k \times k}{k x k} precision
  matrix \eqn{\Omega}{Omega}
  \item Mean: \eqn{E(\theta) = \mu}{E(theta) = mu}
  \item Variance: \eqn{var(\theta) = undefined}{var(theta) = undefined}
  \item Mode: \eqn{mode(\theta) = \mu}{mode(theta) = mu}
}

The multivariate Cauchy distribution is a multidimensional extension
of the one-dimensional or univariate Cauchy distribution. A random
vector is considered to be multivariate Cauchy-distributed if every
linear combination of its components has a univariate Cauchy
distribution. The multivariate Cauchy distribution is equivalent to a
multivariate t distribution with 1 degree of freedom.

The Cauchy distribution is known as a pathological distribution because
its mean and variance are undefined, and it does not satisfy the central
limit theorem.

It is usually parameterized with mean and a covariance matrix, or in
Bayesian inference, with mean and a precision matrix, where the
precision matrix is the matrix inverse of the covariance matrix. These
functions provide the precision parameterization for convenience and
familiarity. It is easier to calculate a multivariate Cauchy density
with the precision parameterization, because a matrix inversion can be
avoided.

This distribution has a mean parameter vector \eqn{\mu}{mu} of length
\eqn{k}, and a \eqn{k \times k}{k x k} precision matrix \eqn{\Omega}{Omega},
which must be positive-definite.
}
\value{
  \code{dmvcp} gives the density and 
  \code{rmvcp} generates random deviates.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{dcauchy}},
  \code{\link{dmvc}},
  \code{\link{dmvt}},
  \code{\link{dmvtp}}, and
  \code{\link{dwishart}}.
}
\examples{
library(LaplacesDemon)
x <- seq(-2,4,length=21)
y <- 2*x+10
z <- x+cos(y) 
mu <- c(1,12,2)
Omega <- matrix(c(1,2,0,2,5,0.5,0,0.5,3), 3, 3)
f <- dmvcp(cbind(x,y,z), mu, Omega)

X <- rmvcp(1000, rep(0,2), diag(2))
X <- X[rowSums((X >= quantile(X, probs=0.025)) &
     (X <= quantile(X, probs=0.975)))==2,]
joint.density.plot(X[,1], X[,2], color=TRUE)
}
\keyword{Distribution}


