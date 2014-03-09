\name{dist.Multivariate.Polya}
\alias{dmvpolya}
\alias{rmvpolya}
\title{Multivariate Polya Distribution}
\description{
  These functions provide the density and random number generation for
  the multivariate Polya distribution.
}
\usage{
dmvpolya(x, alpha, log=FALSE)
rmvpolya(n, alpha)
}
\arguments{
  \item{x}{This is data or parameters in the form of a vector of length
    \eqn{k}.}
  \item{n}{This is the number of random draws to take from the
    distribution.}
  \item{alpha}{This is shape vector \eqn{\alpha}{alpha} with length
    \eqn{k}.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
       density is returned.}
}
\details{
  \itemize{
  \item Application: Discrete Multivariate
  \item Density: \deqn{p(\theta) = \frac{N!}{\prod_k N_k!} \frac{(\sum_k
      \alpha_k - 1)!}{(\sum_k \theta_k + \sum_k \alpha_k - 1)!}
    \frac{\prod (\theta + \alpha - 1)!}{(\alpha - 1)!}}{p(theta) = (N! /
    prod(N[k]!)) * ((sum alpha[k] - 1)! / (sum theta[k] + sum alpha[k] -
    1)!) * prod((theta + alpha - 1)! / (alpha - 1)!)}
  \item Inventor: George Polya (1887-1985)
  \item Notation 1: \eqn{\theta \sim \mathcal{MPO}(\alpha)}{theta ~ MPO(alpha)}
  \item Notation 3: \eqn{p(\theta) = \mathcal{MPO}(\theta |
    \alpha)}{p(theta) = MPO(theta | alpha)}
  \item Parameter 1: shape parameter vector \eqn{\alpha}{alpha}
  \item Mean: \eqn{E(\theta) = }{E(theta) = }
  \item Variance: \eqn{var(\theta) =}{var(theta) = }
  \item Mode: \eqn{mode(\theta) = }{mode(theta) = }
}

The multivariate Polya distribution is named after George Polya
(1887-1985). It is also called the Dirichlet compound multinomial
distribution or the Dirichlet-multinomial distribution. The multivariate
Polya distribution is a compound probability distribution, where a
probability vector \eqn{p} is drawn from a Dirichlet distribution with
parameter vector \eqn{\alpha}{alpha}, and a set of \eqn{N} discrete
samples is drawn from the categorical distribution with probability
vector \eqn{p} and having \eqn{K} discrete categories. The compounding
corresponds to a Polya urn scheme. In document classification, for
example, the distribution is used to represent probabilities over word
counts for different document types. The multivariate Polya distribution
is a multivariate extension of the univariate Beta-binomial distribution.
}
\value{
  \code{dmvpolya} gives the density and \code{rmvpolya} generates random
  deviates.
}
\author{Statisticat, LLC \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{dcat}},
  \code{\link{ddirichlet}}, and
  \code{\link{dmultinom}}.
}
\examples{
library(LaplacesDemon)
dmvpolya(x=1:3, alpha=1:3, log=TRUE)
x <- rmvpolya(1000, c(0.1,0.3,0.6))
}
\keyword{Distribution}