\name{dist.Dirichlet}
\alias{ddirichlet}
\alias{rdirichlet}
\title{Dirichlet Distribution}
\description{
  This is the density function and random generation from the Dirichlet
  distribution.
}
\usage{
ddirichlet(x, alpha, log=FALSE)
rdirichlet(n, alpha)
}
\arguments{
  \item{x}{This is a vector containing a single deviate or matrix
    containing one random deviate per row. Each vector, or matrix row,
    must sum to 1.}
  \item{n}{This is the number of random deviates to generate.}
  \item{alpha}{This is a vector or matrix of shape parameters.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
    \item Application: Continuous Multivariate
    \item Density: \deqn{p(\theta) = \frac{\gamma(\alpha_1 + \dots +
	\alpha_k)}{\gamma \alpha_1 \dots \gamma \alpha_k}
      \theta^{(\alpha[1]-1)}_1 \dots \theta^{(\alpha[k]-1)}_k, \quad
      \theta_1, \dots, \theta_k > 0, \quad \sum^k_{j=1} \theta_j = 1}{p(theta)
      = (gamma(alpha[1] + ... + alpha[k]) / 
          (gamma(alpha[1]) ... gamma(alpha[k]))) * 
          theta[1]^(alpha[1] - 1) ... theta[k]^(alpha[k] - 1)
          theta[1],...,theta[k] >= 0, 
          (the sum of j=1 to k of) theta[j] = 1}
    \item Inventor: Johann Peter Gustav Lejeune Dirichlet (1805-1859)
    \item Notation 1: \eqn{\theta \sim}{theta ~}
    Dirichlet(\eqn{\alpha_1,\dots,\alpha_k}{alpha[1],..., alpha[k]})
    \item Notation 2: \eqn{p(\theta) =}{p(theta) =} Dirichlet(\eqn{\theta
      | \alpha_1,\dots,\alpha_k}{theta | alpha[1],...,alpha[k]})
    \item Notation 3: \eqn{\theta \sim
      \mathcal{DIR}(\alpha_1,\dots,\alpha_k)}{theta ~ Dir(alpha[1],...,
  alpha[k])}
    \item Notation 4: \eqn{p(\theta) = \mathcal{DIR}(\theta |
      \alpha_1,\dots,\alpha_k)}{p(theta) = Dir(theta |
      alpha[1],...,alpha[k])}
    \item Parameter: 'prior sample sizes' \eqn{\alpha_j > 0, \alpha_0 =
      \sum^k_{j=1} \alpha_j}{alpha[j] > 0, alpha[0] = (the sum from j=1
      to k of) alpha[j]}
    \item Mean: \eqn{E(\theta_j) =
      \frac{\alpha_j}{\alpha_0}}{E(theta[j]) = alpha[j] / alpha[0]}
    \item Variance: \eqn{var(\theta_j) = \frac{\alpha_j (\alpha_0 -
	\alpha_j)}{\alpha^2_0 (\alpha_0 + 1)}}{var(theta[j]) = (alpha[j]
      * (alpha[0] - alpha[j])) / (alpha[0]^2 * (alpha[0]+ 1))}
    \item Covariance: \eqn{cov(\theta_i, \theta_j) = - \frac{\alpha_i
	\alpha_j}{\alpha^2_0 (\alpha_0 + 1)}}{cov(theta[i], theta[j]) =
      - ((alpha[i]*alpha[j]) / (alpha[0]^2 * (alpha[0] + 1)))}
    \item Mode: \eqn{mode(\theta_j) = \frac{\alpha_j - 1}{\alpha_0 -
	k}}{mode(theta[j]) = (alpha[j] - 1) / (alpha[0] - k)}
}

The Dirichlet distribution is the multivariate generalization of the
univariate beta distribution. Its probability density function returns
the belief that the probabilities of \eqn{k} rival events are
\eqn{\theta_j}{theta[j]} given that each event has been observed
\eqn{\alpha_j - 1}{alpha[j] - 1} times.

The Dirichlet distribution is commonly used as a prior distribution in
Bayesian inference. The Dirichlet distribution is the conjugate prior
distribution for the parameters of the categorical and multinomial
distributions.

A very common special case is the symmetric Dirichlet distribution,
where all of the elements in parameter vector \eqn{\alpha}{alpha} have
the same value. Symmetric Dirichlet distributions are often used as 
vague or weakly informative Dirichlet prior distributions, so that one
component is not favored over another. The single value that is entered
into all elements of \eqn{\alpha}{alpha} is called the concentration
parameter.
}
\value{
  \code{ddirichlet} gives the density and
  \code{rdirichlet} generates random deviates.
}
\seealso{
  \code{\link{dbeta}},
  \code{\link{dcat}},
  \code{\link{dmvpolya}},
  \code{\link{dmultinom}}, and
  \code{\link{TransitionMatrix}}.
}
\examples{
library(LaplacesDemon)
x <- ddirichlet(c(.1,.3,.6), c(1,1,1))
x <- rdirichlet(10, c(1,1,1))
}
\keyword{Distribution}
