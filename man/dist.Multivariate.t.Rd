\name{dist.Multivariate.t}
\alias{dmvt}
\alias{rmvt}
\title{Multivariate t Distribution}
\description{
  These functions provide the density and random number generation for
  the multivariate t distribution, otherwise called the multivariate
  Student distribution.
}
\usage{
dmvt(x, mu, S, df=Inf, log=FALSE)
rmvt(n=1, mu, S, df=Inf)
}
\arguments{
  \item{x}{This is either a vector of length \eqn{k} or a matrix with
    a number of columns, \eqn{k}, equal to the number of columns in
    scale matrix \eqn{\textbf{S}}{S}.}
  \item{n}{This is the number of random draws.}
  \item{mu}{This is a numeric vector or matrix representing the location
    parameter,\eqn{\mu}{mu} (the mean vector), of the multivariate
    distribution (equal to the expected value when \code{df > 1},
    otherwise represented as \eqn{\nu > 1}{nu > 1}). When a vector, it
    must be of length \eqn{k}, or must have \eqn{k} columns as a matrix,
    as defined above.}
  \item{S}{This is a \eqn{k \times k}{k x k} positive-definite scale
    matrix \eqn{\textbf{S}}{S}, such that \code{S*df/(df-2)} is the
    variance-covariance matrix when \code{df > 2}. A vector of
    length 1 is also allowed (in this case, \eqn{k=1} is set).}
  \item{df}{This is the degrees of freedom, and is often represented
    with \eqn{\nu}{nu}.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
  \item Application: Continuous Multivariate
  \item Density: \deqn{p(\theta) =
    \frac{\Gamma[(\nu+k)/2]}{\Gamma(\nu/2)\nu^{k/2}\pi^{k/2}|\Sigma|^{1/2}[1
      + (1/\nu)(\theta-\mu)^{\mathrm{T}} \Sigma^{-1}
      (\theta-\mu)]^{(\nu+k)/2}}}{p(theta) = Gamma[(nu+k)/2] /
    {Gamma(nu/2)nu^(k/2)pi^(k/2)|Sigma|^(1/2)[1 +
      (1/nu)(theta-mu)^T*Sigma^(-1)(theta-mu)]^[(nu+k)/2]}}
  \item Inventor: Unknown (to me, anyway)
  \item Notation 1: \eqn{\theta \sim \mathrm{t}_k(\mu,
    \Sigma, \nu)}{theta ~ t[k](mu, Sigma, nu)}
  \item Notation 2: \eqn{p(\theta) = \mathrm{t}_k(\theta | \mu,
    \Sigma, \nu)}{p(theta) = t[k](theta | mu, Sigma, nu)}
  \item Parameter 1: location vector \eqn{\mu}{mu}
  \item Parameter 2: positive-definite \eqn{k \times k}{k x k} scale
    matrix \eqn{\Sigma}{Sigma}
  \item Parameter 3: degrees of freedom \eqn{\nu > 0}{nu > 0} (df in the
    functions)
  \item Mean: \eqn{E(\theta) = \mu}{E(theta) = mu}, for \eqn{\nu > 1}{nu
    > 1}, otherwise undefined
  \item Variance: \eqn{var(\theta) = \frac{\nu}{\nu - 2}
    \Sigma}{var(theta) = (nu / (nu - 2))*Sigma}, for \eqn{\nu > 2}{nu > 2}
  \item Mode: \eqn{mode(\theta) = \mu}{mode(theta) = mu}
}

The multivariate t distribution, also called the multivariate Student or
multivariate Student t distribution, is a multidimensional extension of the
one-dimensional or univariate Student t distribution. A random vector is
considered to be multivariate t-distributed if every linear
combination of its components has a univariate Student t-distribution.
This distribution has a mean parameter vector \eqn{\mu}{mu} of length
\eqn{k}, and a \eqn{k \times k}{k x k} scale matrix \eqn{\textbf{S}}{S},
which must be positive-definite. When degrees of freedom
\eqn{\nu=1}{nu=1}, this is the multivariate Cauchy distribution.
}
\value{
  \code{dmvt} gives the density and 
  \code{rmvt} generates random deviates.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{dinvwishart}},
  \code{\link{dmvc}},
  \code{\link{dmvcp}},
  \code{\link{dmvtp}},
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
S <- matrix(c(1,2,0,2,5,0.5,0,0.5,3), 3, 3)
df <- 4
f <- dmvt(cbind(x,y,z), mu, S, df)
X <- rmvt(1000, c(0,1,2), S, 5)
joint.density.plot(X[,1], X[,2], color=TRUE)
}
\keyword{Distribution}


