\name{dist.Multivariate.t.Precision}
\alias{dmvtp}
\alias{rmvtp}
\title{Multivariate t Distribution: Precision Parameterization}
\description{
  These functions provide the density and random number generation for
  the multivariate t distribution, otherwise called the multivariate
  Student distribution. These functions use the precision
  parameterization.
}
\usage{
dmvtp(x, mu, Omega, nu=Inf, log=FALSE)
rmvtp(n=1, mu, Omega, nu=Inf)
}
\arguments{
  \item{x}{This is either a vector of length \eqn{k} or a matrix with 
       a number of columns, \eqn{k}, equal to the number of columns in
       precision matrix \eqn{\Omega}{Omega}.}
  \item{n}{This is the number of random draws.}
  \item{mu}{This is a numeric vector representing the location parameter,
       \eqn{\mu}{mu} (the mean vector), of the multivariate distribution
       (equal to the expected value when \code{df > 1}, otherwise
       represented as \eqn{\nu > 1}{nu > 1}). It must be of length
       \eqn{k}, as defined above.}
  \item{Omega}{This is a \eqn{k \times k}{k x k} positive-definite
       precision matrix \eqn{\Omega}{Omega}.}
  \item{nu}{This is the degrees of freedom \eqn{\nu}{nu}, which must be
       positive.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
       density is returned.}
}
\details{
  \itemize{
  \item Application: Continuous Multivariate
  \item Density: \deqn{p(\theta) =
    \frac{\Gamma((\nu+k)/2)}{\Gamma(\nu/2)\nu^{k/2}\pi^{k/2}}
    |\Omega|^{1/2} (1 + \frac{1}{\nu} (\theta-\mu)^T \Omega
    (\theta-\mu))^{-(\nu+k)/2}}{p(theta) = (Gamma((nu+k)/2) /
    (Gamma(nu/2)*nu^(k/2)*pi^(k/2))) * |Omega|^(1/2) *
    (1 + (1/nu) (theta-mu)^T Omega (theta-mu))^(-(nu+k)/2)}
  \item Inventor: Unknown (to me, anyway)
  \item Notation 1: \eqn{\theta \sim \mathrm{t}_k(\mu,
    \Omega^{-1}, \nu)}{theta ~ t[k](mu, Omega^(-1), nu)}
  \item Notation 2: \eqn{p(\theta) = \mathrm{t}_k(\theta | \mu,
    \Omega^{-1}, \nu)}{p(theta) = t[k](theta | mu, Omega^(-1), \nu)}
  \item Parameter 1: location vector \eqn{\mu}{mu}
  \item Parameter 2: positive-definite \eqn{k \times k}{k x k} precision
  matrix \eqn{\Omega}{Omega}
  \item Parameter 3: degrees of freedom \eqn{\nu > 0}{nu > 0}
  \item Mean: \eqn{E(\theta) = \mu}{E(theta) = mu}, for \eqn{\nu > 1}{nu
    > 1}, otherwise undefined
  \item Variance: \eqn{var(\theta) = \frac{\nu}{\nu - 2}
    \Omega^{-1}}{var(theta) = (nu / (nu - 2))*Omega^(-1)}, for \eqn{\nu
    > 2}{nu> 2}
  \item Mode: \eqn{mode(\theta) = \mu}{mode(theta) = mu}
}

The multivariate t distribution, also called the multivariate Student or
multivariate Student t distribution, is a multidimensional extension of the
one-dimensional or univariate Student t distribution. A random vector is
considered to be multivariate t-distributed if every linear
combination of its components has a univariate Student t-distribution.

It is usually parameterized with mean and a covariance matrix, or in
Bayesian inference, with mean and a precision matrix, where the
precision matrix is the matrix inverse of the covariance matrix. These
functions provide the precision parameterization for convenience and
familiarity. It is easier to calculate a multivariate t density
with the precision parameterization, because a matrix inversion can be
avoided.

This distribution has a mean parameter vector \eqn{\mu}{mu} of length
\eqn{k}, and a \eqn{k \times k}{k x k} precision matrix
\eqn{\Omega}{Omega}, which must be positive-definite. When degrees of
freedom \eqn{\nu=1}{nu=1}, this is the multivariate Cauchy distribution.
}
\value{
  \code{dmvtp} gives the density and 
  \code{rmvtp} generates random deviates.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{dwishart}},
  \code{\link{dmvc}},
  \code{\link{dmvcp}},
  \code{\link{dmvt}},
  \code{\link{dst}},
  \code{\link{dstp}}, and
  \code{\link{dt}}.
}
\examples{
library(LaplacesDemon)
x <- seq(-2,4,length=21)
y <- 2*x+10
z <- x+cos(y) 
mu <- c(1,12,2)
Omega <- matrix(c(1,2,0,2,5,0.5,0,0.5,3), 3, 3)
nu <- 4
f <- dmvtp(cbind(x,y,z), mu, Omega, nu)
X <- rmvtp(1000, c(0,1,2), diag(3), 5)
joint.density.plot(X[,1], X[,2], color=TRUE)
}
\keyword{Distribution}


