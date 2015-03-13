\name{dist.Asymmetric.Multivariate.Laplace}
\alias{daml}
\alias{raml}
\title{Asymmetric Multivariate Laplace Distribution}
\description{
  These functions provide the density and random generation for the
  asymmetric multivariate Laplace distribution with location and skew
  parameter \eqn{\mu}{mu} and covariance \eqn{\Sigma}{Sigma}.
}
\usage{
daml(x, mu, Sigma, log=FALSE)
raml(n, mu, Sigma)
}
\arguments{
  \item{x}{This is a \eqn{N \times K}{N x K} matrix of data, or a
    vector of length \eqn{K}.}
  \item{n}{This is the number of observations, which must be a positive
    integer that has length 1.}
  \item{mu}{This is the location and skew parameter \eqn{\mu}{mu}. This
    may be a \eqn{N \times K}{N x K} matrix, or a vector of length
    \eqn{K}.}
  \item{Sigma}{This is the \eqn{K \times K}{K x K} positive-definite
    covariance matrix \eqn{\Sigma}{Sigma}.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
    \item Application: Continuous Multivariate
    \item Density: \eqn{p(\theta) =
      \frac{2\exp(\theta\Omega\theta)}{(2\pi)^{k/2}|\Sigma|^0.5}
      \frac{\theta\Omega\theta}{2 + \mu\Omega\mu}^{(2-k)/4}
      K_{(2-k)/2}(\sqrt{(2 + \mu\Omega\mu)(\theta\Omega\theta)})}{
      p(theta) = ((2exp(theta Omega theta)) /
      ((2pi)^(k/2)|Sigma|^(0.5)))
      ((theta Omega theta) / (2 + mu Omega mu))^((k-2)/4)
      K[(2-k)/2](sqrt((2 + mu Omega mu)(theta Omega theta)))}
    \item Inventor: Kotz, Kozubowski, and Podgorski (2003)
    \item Notation 1: \eqn{\theta \sim \mathcal{AL}_K(\mu,
      \Sigma)}{theta ~ AL[K](mu, Sigma)}
    \item Notation 2: \eqn{p(\theta) = \mathcal{AL}_K(\theta | \mu,
      \Sigma)}{p(theta) = AL[K](theta | mu, Sigma)}
    \item Parameter 1: location-skew parameter \eqn{\mu}{mu}
    \item Parameter 2: positive-definite covariance matrix
      \eqn{\Sigma}{Sigma}
    \item Mean: Unknown
    \item Variance: Unknown
    \item Mode: \eqn{mode(\theta) = \mu}{mode(theta) = mu}
}

The asymmetric multivariate Laplace distribution of Kotz, Kozubowski,
and Podgorski (2003) is a multivariate extension of the univariate,
asymmetric Laplace distribution. It is parameterized according to
two parameters: location-skew parameter \eqn{\mu}{mu} and positive-definite
covariance matrix \eqn{\Sigma}{Sigma}. Location and skew occur in the
same parameter. When \eqn{\mu=0}{mu=0}, the density is the (symmetric)
multivariate Laplace of Anderson (1992). As each location deviates from
zero, the marginal distribution becomes more skewed. Since location and
skew are combined, it is appropriate for zero-centered variables, such
as a matrix of centered and scaled dependent variables in cluster
analysis, factor analysis, multivariate regression, or multivariate
time-series.

The asymmetric multivariate Laplace distribution is also discussed
earlier in Kozubowski and Podgorski (2001), and is well-suited for
financial modeling via multivariate regression, specifically with
currency exchange rates. Cajigas and Urga (2005) fit residuals in a
multivariate GARCH model with the asymmetric multivariate Laplace
distribution, regarding stocks and bonds. They find that it
"overwhelmingly outperforms" normality.
}
\value{
  \code{daml} gives the density, and
  \code{raml} generates random deviates.
}
\references{
  Anderson, D.N. (1992). "A Multivariate Linnik Distribution".
  \emph{Statistical Probability Letters}, 14, p. 333--336.
  
  Cajigas, J.P. and Urga, G. (2005) "Dynamic Conditional Correlation
  Models with Asymmetric Laplace Innovations". Centre for Economic
  Analysis: Cass Business School.
  
  Kotz, S., Kozubowski, T.J., and Podgorski, K. (2003). "An
  Asymmetric Multivariate Laplace Distribution". Working Paper.
    
  Kozubowski, T.J. and Podgorski, K. (2001). "Asymmetric Laplace
  Laws and Modeling Financial Data". \emph{Mathematical and Computer
  Modelling}, 34, p. 1003--1021.
}
\seealso{\code{\link{dalaplace}} and
  \code{\link{dmvl}}
}
\examples{
library(LaplacesDemon)
x <- daml(c(1,2,3), c(0,1,2), diag(3))
X <- raml(1000, c(0,1,2), diag(3))
joint.density.plot(X[,1], X[,2], color=FALSE)
}
\keyword{Distribution}