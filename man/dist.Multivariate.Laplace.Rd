\name{dist.Multivariate.Laplace}
\alias{dmvl}
\alias{rmvl}
\title{Multivariate Laplace Distribution}
\description{
  These functions provide the density and random number generation
  for the multivariate Laplace distribution.
}
\usage{
dmvl(x, mu, Sigma, log=FALSE)
rmvl(n, mu, Sigma)
}
\arguments{
  \item{x}{This is data or parameters in the form of a vector of length
       \eqn{k} or a matrix with \eqn{k} columns.}
  \item{n}{This is the number of random draws.}
  \item{mu}{This is mean vector \eqn{\mu}{mu} with length \eqn{k} or
       matrix with \eqn{k} columns.}
  \item{Sigma}{This is the \eqn{k \times k}{k x k} covariance matrix
       \eqn{\Sigma}{Sigma}.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
       density is returned.}
}
\details{
  \itemize{
  \item Application: Continuous Multivariate
  \item Density: \deqn{p(\theta) = \frac{2}{(2\pi)^{k/2}
      |\Sigma|^{1/2}} \frac{(\pi/(2\sqrt{2(\theta - \mu)^T \Sigma^{-1}
      (\theta - \mu)}))^{1/2} \exp(-\sqrt{2(\theta - \mu)^T \Sigma^{-1}
      (\theta - \mu)})}{\sqrt{((\theta - \mu)^T \Sigma^{-1} (\theta -
      \mu) / 2)}^{k/2-1}}}{p(theta) = (2 / ((2*pi)^(k/2) *
    |Sigma|^(1/2)))
    ((sqrt(pi/(2*sqrt(2*(theta-mu)^TSigma^(-1)(theta-mu)))) *
    exp(-sqrt(2*(theta-mu)^TSigma^(-1)(theta-mu)))) /
    sqrt((theta-mu)^TSigma^(-1)(theta-mu)/2)^(k/2-1))}
  \item Inventor: Fang et al. (1990)
  \item Notation 1: \eqn{\theta \sim \mathcal{MVL}(\mu, \Sigma)}{theta ~
    MVL(mu, Sigma)}
  \item Notation 2: \eqn{\theta \sim \mathcal{L}_k(\mu, \Sigma)}{theta ~
    L[k](mu, Sigma)}
  \item Notation 3: \eqn{p(\theta) = \mathcal{MVL}(\theta | \mu,
    \Sigma)}{p(theta) = MVL(theta | mu, Sigma)}
  \item Notation 4: \eqn{p(\theta) = \mathcal{L}_k(\theta | \mu,
    \Sigma)}{p(theta) = L[k](theta | mu, Sigma)}
  \item Parameter 1: location vector \eqn{\mu}{mu}
  \item Parameter 2: positive-definite \eqn{k \times k}{k x k}
  covariance matrix \eqn{\Sigma}{Sigma}
  \item Mean: \eqn{E(\theta) = \mu}{E(theta) = mu}
  \item Variance: \eqn{var(\theta) = \Sigma}{var(theta) = Sigma}
  \item Mode: \eqn{mode(\theta) = \mu}{mode(theta) = mu}
}

The multivariate Laplace distribution is a multidimensional extension of
the one-dimensional or univariate symmetric Laplace distribution. There
are multiple forms of the multivariate Laplace distribution.

The bivariate case was introduced by Ulrich and Chen (1987), and the
first form in larger dimensions may have been Fang et al. (1990), which
requires a Bessel function. Alternatively, multivariate Laplace was soon
introduced as a special case of a multivariate Linnik distribution
(Anderson, 1992), and later as a special case of the multivariate power
exponential distribution (Fernandez et al., 1995; Ernst, 1998). Bayesian
considerations appear in Haro-Lopez and Smith (1999). Wainwright and
Simoncelli (2000) presented multivariate Laplace as a Gaussian scale
mixture. Kotz et al. (2001) present the distribution formally. Here, the
density is calculated with the asymptotic formula for the Bessel
function as presented in Wang et al. (2008).

The multivariate Laplace distribution is an attractive alternative to
the multivariate normal distribution due to its wider tails, and remains
a two-parameter distribution (though alternative three-parameter forms
have been introduced as well), unlike the three-parameter multivariate t
distribution, which is often used as a robust alternative to the
multivariate normal distribution.
}
\value{
  \code{dmvl} gives the density, and
  \code{rmvl} generates random deviates.
}
\references{
  Anderson, D.N. (1992). "A Multivariate Linnik Distribution".
  \emph{Statistical Probability Letters}, 14, p. 333--336.

  Eltoft, T., Kim, T., and Lee, T. (2006). "On the Multivariate Laplace
  Distribution". \emph{IEEE Signal Processing Letters}, 13(5),
  p. 300--303.

  Ernst, M. D. (1998). "A Multivariate Generalized Laplace
  Distribution". \emph{Computational Statistics}, 13, p. 227--232.

  Fang, K.T., Kotz, S., and Ng, K.W. (1990). "Symmetric Multivariate and
  Related Distributions". Monographs on Statistics and Probability, 36,
  Chapman-Hall, London.

  Fernandez, C., Osiewalski, J. and Steel, M.F.J. (1995). "Modeling and
  Inference with v-spherical Distributions". \emph{Journal of the
  American Statistical Association}, 90, p. 1331--1340.
  
  Gomez, E., Gomez-Villegas, M.A., and Marin, J.M. (1998). "A
  Multivariate Generalization of the Power Exponential Family of
  Distributions". \emph{Communications in Statistics-Theory and
  Methods}, 27(3), p. 589--600.

  Haro-Lopez, R.A. and Smith, A.F.M. (1999). "On Robust Bayesian
  Analysis for Location and Scale Parameters". \emph{Journal of
  Multivariate Analysis}, 70, p. 30--56.

  Kotz., S., Kozubowski, T.J., and Podgorski, K. (2001). "The Laplace
  Distribution and Generalizations: A Revisit with Applications to
  Communications, Economics, Engineering, and Finance". Birkhauser:
  Boston, MA.

  Ulrich, G. and Chen, C.C. (1987). "A Bivariate Double Exponential
  Distribution and its Generalization". \emph{ASA Proceedings on
  Statistical Computing}, p. 127--129.

  Wang, D., Zhang, C., and Zhao, X. (2008). "Multivariate Laplace
  Filter: A Heavy-Tailed Model for Target Tracking". \emph{Proceedings
  of the 19th International Conference on Pattern Recognition}: FL.

  Wainwright, M.J. and Simoncelli, E.P. (2000). "Scale Mixtures of
  Gaussians and the Statistics of Natural Images". \emph{Advances in
  Neural Information Processing Systems}, 12, p. 855--861.

}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{daml}},
  \code{\link{dlaplace}},
  \code{\link{dmvn}},
  \code{\link{dmvnp}},
  \code{\link{dmvpe}},
  \code{\link{dmvt}},
  \code{\link{dnorm}},
  \code{\link{dnormp}}, and
  \code{\link{dnormv}}.
}
\examples{
library(LaplacesDemon)
x <- dmvl(c(1,2,3), c(0,1,2), diag(3))
X <- rmvl(1000, c(0,1,2), diag(3))
joint.density.plot(X[,1], X[,2], color=TRUE)
}
\keyword{Distribution}