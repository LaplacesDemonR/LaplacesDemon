\name{dist.Multivariate.Cauchy}
\alias{dmvc}
\alias{rmvc}
\title{Multivariate Cauchy Distribution}
\description{
  These functions provide the density and random number generation for
  the multivariate Cauchy distribution.
}
\usage{
dmvc(x, mu, S, log=FALSE)
rmvc(n=1, mu, S)
}
\arguments{
  \item{x}{This is either a vector of length \eqn{k} or a matrix with 
       a number of columns, \eqn{k}, equal to the number of columns in
       scale matrix \eqn{\textbf{S}}{S}.}
  \item{n}{This is the number of random draws.}
  \item{mu}{This is a numeric vector representing the location parameter,
       \eqn{\mu}{mu} (the mean vector), of the multivariate distribution
       It must be of length \eqn{k}, as defined above.}
  \item{S}{This is a \eqn{k \times k}{k x k} positive-definite scale
       matrix \eqn{\textbf{S}}{S}.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
       density is returned.}
}
\details{
  \itemize{
  \item Application: Continuous Multivariate
  \item Density: \deqn{p(\theta) =
    \frac{\Gamma[(1+k)/2]}{\Gamma(1/2)1^{k/2}\pi^{k/2}|\Sigma|^{1/2}[1+(\theta-\mu)^{\mathrm{T}}\Sigma^{-1}(\theta-\mu)]^{(1+k)/2}}}{p(theta)
    = Gamma[(1+k)/2] /
    {Gamma(1/2)1^(k/2)pi^(k/2)|Sigma|^(1/2)[1+(theta-mu)^T*Sigma^(-1)(theta-mu)]^[(1+k)/2]}}
  \item Inventor: Unknown (to me, anyway)
  \item Notation 1: \eqn{\theta \sim \mathcal{MC}_k(\mu, \Sigma)}{theta ~
    MC[k](mu, Sigma)}
  \item Notation 2: \eqn{p(\theta) = \mathcal{MC}_k(\theta | \mu,
    \Sigma)}{p(theta) = MC[k](theta | mu, Sigma)}
  \item Parameter 1: location vector \eqn{\mu}{mu}
  \item Parameter 2: positive-definite \eqn{k \times k}{k x k} scale
  matrix \eqn{\Sigma}{Sigma}
  \item Mean: \eqn{E(\theta) = \mu}{E(theta) = mu}
  \item Variance: \eqn{var(\theta) = undefined}{var(theta) = undefined}
  \item Mode: \eqn{mode(\theta) = \mu}{mode(theta) = mu}
}

The multivariate Cauchy distribution is a multidimensional extension of the
one-dimensional or univariate Cauchy distribution. The multivariate
Cauchy distribution is equivalent to a multivariate t distribution with
1 degree of freedom. A random vector is considered to be multivariate
Cauchy-distributed if every linear combination of its components has a
univariate Cauchy distribution.

The Cauchy distribution is known as a pathological distribution because
its mean and variance are undefined, and it does not satisfy the central
limit theorem.
}
\value{
  \code{dmvc} gives the density and 
  \code{rmvc} generates random deviates.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{dcauchy}},
  \code{\link{dinvwishart}},
  \code{\link{dmvcp}},
  \code{\link{dmvt}}, and
  \code{\link{dmvtp}}.
}
\examples{
library(LaplacesDemon)
x <- seq(-2,4,length=21)
y <- 2*x+10
z <- x+cos(y) 
mu <- c(1,12,2)
Sigma <- matrix(c(1,2,0,2,5,0.5,0,0.5,3), 3, 3)
f <- dmvc(cbind(x,y,z), mu, Sigma)

X <- rmvc(1000, rep(0,2), diag(2))
X <- X[rowSums((X >= quantile(X, probs=0.025)) &
     (X <= quantile(X, probs=0.975)))==2,]
joint.density.plot(X[,1], X[,2], color=TRUE)
}
\keyword{Distribution}


