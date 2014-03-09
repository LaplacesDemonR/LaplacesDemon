\name{Precision}
\alias{Cov2Prec}
\alias{Prec2Cov}
\alias{prec2sd}
\alias{prec2var}
\alias{sd2prec}
\alias{sd2var}
\alias{var2prec}
\alias{var2sd}
\title{Precision}
\description{
  Bayesians often use precision rather than variance. These are elementary
  utility functions to facilitate conversions between precision,
  standard deviation, and variance regarding scalars, vectors, and
  matrices, and these functions are designed for those who are new to
  Bayesian inference. The names of these functions consist of two
  different scale parameters, separated by a `2', and capital letters
  refer to matrices while lower case letters refer to scalars and
  vectors. For example, the \code{Prec2Cov} function converts a
  precision matrix to a covariance matrix, while the \code{prec2sd}
  function converts a scalar or vector of precision parameters to
  standard deviation parameters.

  The modern Bayesian use of precision developed because it was more
  straightforward in a normal distribution to estimate precision
  \eqn{\tau}{tau} with a gamma distribution as a conjugate prior, than
  to estimate \eqn{\sigma^2}{sigma^2} with an inverse-gamma distribution
  as a conjugate prior. Today, conjugacy is usually considered to be
  merely a convenience, and in this example, a non-conjugate half-Cauchy
  prior distribution is recommended as a weakly informative prior
  distribution for scale parameters.
}
\usage{
Cov2Prec(Cov)
Prec2Cov(Prec)
prec2sd(prec=1)
prec2var(prec=1)
sd2prec(sd=1)
sd2var(sd=1)
var2prec(var=1)
var2sd(var=1)
}
\arguments{
  \item{Cov}{This is a covariance matrix, usually represented as
    \eqn{\Sigma}{Sigma}.}
  \item{Prec}{This is a precision matrix, usually represented as
    \eqn{\Omega}{Omega}.}
  \item{prec}{This is a precision scalar or vector, usually represented as
    \eqn{\tau}{tau}.}
  \item{sd}{This is a standard deviation scalar or vector, usually
    represented as \eqn{\sigma}{sigma}.}
  \item{var}{This is a variance scalar or vector, usually
    represented as \eqn{\sigma^2}{sigma^2}.}
}
\details{Bayesians often use precision rather than variance, where
  precision is the inverse of the variance. For example, a linear
  regression may be represented equivalently as \eqn{\textbf{y} \sim
    \mathcal{N}(\mu, \sigma^2)}{y ~ N(mu, sigma^2)}, or \eqn{\textbf{y}
    \sim \mathcal{N}(\mu, \tau^{-1})}{y ~ N(mu, tau^(-1))}, where
  \eqn{\sigma^2}{sigma^2} is the variance, and \eqn{\tau}{tau} is the
  precision, which is the inverse of the variance.}
\value{
  \item{Cov2Prec}{
    This returns a precision matrix, \eqn{\Omega}{Omega}, from a
    covariance matrix, \eqn{\Sigma}{Sigma}, where \eqn{\Omega =
    \Sigma^{-1}}{Omega = Sigma^(-1)}.}
  \item{Prec2Cov}{
    This returns a covariance matrix, \eqn{\Sigma}{Sigma}, from a
    precision matrix, \eqn{\Omega}{Omega}, where \eqn{\Sigma =
    \Omega^{-1}}{Sigma = Omega^(-1)}.}
  \item{prec2sd}{
    This returns a standard deviation, \eqn{\sigma}{sigma}, from a
    precision, \eqn{\tau}{tau}, where \eqn{\sigma =
    \sqrt{\tau^{-1}}}{sigma = sqrt(tau^(-1))}.}
  \item{prec2var}{
    This returns a variance, \eqn{\sigma^2}{sigma^2}, from a precision,
    \eqn{\tau}{tau}, where \eqn{\sigma^2 = \tau^{-1}}{sigma^2 =
    tau^(-1)}.}
  \item{sd2prec}{
    This returns a precision, \eqn{\tau}{tau}, from a standard
    deviation, \eqn{\sigma}{sigma}, where \eqn{\tau = \sigma^{-2}}{tau =
    sigma^(-2)}.}
  \item{sd2var}{
    This returns a variance, \eqn{\sigma^2}{sigma^2}, from a standard
    deviation, \eqn{\sigma}{sigma}, where \eqn{\sigma^2 = \sigma
    \sigma}{sigma^2 = sigma x sigma}.}
  \item{var2prec}{
    This returns a precision, \eqn{\tau}{tau}, from a variance,
    \eqn{\sigma^2}{sigma^2}, where \eqn{\tau = \frac{1}{\sigma^2}}{tau = 1
    / sigma^2}.}
  \item{var2sd}{
    This returns a standard deviation, \eqn{\sigma}{sigma}, from a
    variance, \eqn{\sigma^2}{sigma^2}, where \eqn{\sigma =
    \sqrt{\sigma^2}}{sigma = sqrt(sigma^2)}.}
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{\code{\link{Cov2Cor}}}
\examples{
library(LaplacesDemon)
Cov2Prec(matrix(c(1,0.1,0.1,1),2,2))
Prec2Cov(matrix(c(1,0.1,0.1,1),2,2))
prec2sd(0.5)
prec2var(0.5)
sd2prec(1.4142)
sd2var(01.4142)
var2prec(2)
var2sd(2)
}
\keyword{Utility}