\name{dist.Multivariate.t.Precision.Cholesky}
\alias{dmvtpc}
\alias{rmvtpc}
\title{Multivariate t Distribution: Precision-Cholesky Parameterization}
\description{
  These functions provide the density and random number generation for
  the multivariate t distribution, otherwise called the multivariate
  Student distribution. These functions use the precision and Cholesky
  parameterization.
}
\usage{
dmvtpc(x, mu, U, nu=Inf, log=FALSE)
rmvtpc(n=1, mu, U, nu=Inf)
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
  \item{U}{This is a \eqn{k \times k}{k x k} upper-triangular of the
       precision matrix that is Cholesky fator \eqn{\textbf{U}}{U} of
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
avoided. The precision matrix is replaced with an upper-triangular
\eqn{k \times k}{k x k} matrix that is Cholesky factor
\eqn{\textbf{U}}{U}, as per the \code{\link{chol}} function for Cholesky
decomposition.

This distribution has a mean parameter vector \eqn{\mu}{mu} of length
\eqn{k}, and a \eqn{k \times k}{k x k} precision matrix
\eqn{\Omega}{Omega}, which must be positive-definite. When degrees of
freedom \eqn{\nu=1}{nu=1}, this is the multivariate Cauchy distribution.

In practice, \eqn{\textbf{U}}{U} is fully unconstrained for proposals
when its diagonal is log-transformed. The diagonal is exponentiated
after a proposal and before other calculations. Overall, the Cholesky
parameterization is faster than the traditional parameterization.
Compared with \code{dmvtp}, \code{dmvtpc} must additionally
matrix-multiply the Cholesky back to the precision matrix, but it
does not have to check for or correct the precision matrix to
positive-definiteness, which overall is slower. Compared with
\code{rmvtp}, \code{rmvtpc} is faster because the Cholesky decomposition
has already been performed.
}
\value{
  \code{dmvtpc} gives the density and 
  \code{rmvtpc} generates random deviates.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{chol}},
  \code{\link{dwishartc}},
  \code{\link{dmvc}},
  \code{\link{dmvcp}},
  \code{\link{dmvtc}},
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
U <- chol(Omega)
nu <- 4
f <- dmvtpc(cbind(x,y,z), mu, U, nu)
X <- rmvtpc(1000, c(0,1,2), U, 5)
joint.density.plot(X[,1], X[,2], color=TRUE)
}
\keyword{Distribution}


