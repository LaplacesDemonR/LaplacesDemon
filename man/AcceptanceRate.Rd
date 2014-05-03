\name{AcceptanceRate}
\alias{AcceptanceRate}
\title{Acceptance Rate}
\description{
  The \code{Acceptance.Rate} function calculates the acceptance rate per
  chain from a matrix of posterior MCMC samples.
}
\usage{
AcceptanceRate(x)
}
\arguments{
  \item{x}{This required argument accepts a \eqn{S \times J}{S x J}
    numeric matrix of \eqn{S} posterior samples for \eqn{J} variables,
    such as \code{Posterior1} or \code{Posterior2} from an object of
    class \code{demonoid}.}
}
\details{
  The acceptance rate of an MCMC algorithm is the percentage of
  iterations in which the proposals were accepted.

  Optimal Acceptance Rates\cr
  The optimal acceptance rate varies with the number of parameters and
  by algorithm. Algorithms with componentwise Gaussian proposals have an
  optimal acceptance rate of 0.44, regardless of the number of
  parameters. Algorithms that update with multivariate Gaussian
  proposals tend to have an optimal acceptance rate that ranges from
  0.44 for one parameter (one IID Gaussian target distribution) to 0.234
  for an infinite number of parameters (IID Gaussian target
  distributions), and 0.234 is approached quickly as the number of
  parameters increases. The AHMC, HMC, and THMC algorithms have an
  optimal acceptance rate of 0.67, except with the algorithm
  specification \code{L=1}, where the optimal acceptance rate is
  0.574. The target acceptance rate is specified in HMCDA and NUTS, and
  the recommended rate is 0.65 and 0.60 respectively. Some algorithms
  have an acceptance rate of 1, such as AGG, ESS, GG, GS (MISS only),
  SGLD, or Slice.

  Global and Local Acceptance Rates\cr
  \code{\link{LaplacesDemon}} reports the global acceptance rate for the
  un-thinned chains. However, componentwise algorithms make a proposal
  per parameter, and therefore have a local acceptance rate for each
  parameter. Since only the global acceptance rate is reported, the
  \code{AcceptanceRate} function may be used to calculate the local
  acceptance rates from a matrix of un-thinned posterior samples.

  Thinning\cr
  Thinned samples tend to have higher local acceptance rates than
  un-thinned samples. With enough samples and enough thinning, local
  acceptance rates approach 1. Local acceptance rates do not need to
  approach the optimal acceptance rates above. Conversely, local
  acceptance rates do not need to approach 1, because too much
  information may possibly be discarded by thinning. For more
  information on thinning, see the \code{\link{Thin}} function.

  Diagnostics\cr
  The \code{AcceptanceRate} function may be used to calculate local
  acceptance rates on a matrix of thinned or un-thinned samples. Any
  chain with a local acceptance rate that is an outlier may be studied
  for reasons that may cause the outlier. A local acceptance rate
  outlier does not violate theory and is often acceptable, but may
  indicate a potential problem. Only some of the many potential problems
  include: identifiability, model misspecification, multicollinearity,
  multimodality, choice of prior distributions, or becoming trapped in a
  low-probability space. The solution to local acceptance rate outliers
  tends to be either changing the MCMC algorithm or re-specifying the
  model or priors. For example, an MCMC algorithm that makes
  multivariate Gaussian proposals for a large number of parameters may
  have low global and local acceptance rates when far from the target
  distributions.
}
\value{
  The \code{AcceptanceRate} function returns a vector of acceptance
  rates, one for each chain.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{LaplacesDemon}},
  \code{\link{MISS}},
  \code{\link{PosteriorChecks}}, and
  \code{\link{Thin}}.
}
\examples{
library(LaplacesDemon)
AcceptanceRate(matrix(rnorm(5000),1000,5))
}
\keyword{Diagnostic}
\keyword{MCMC}
\keyword{Utility}
