\name{ESS}
\alias{ESS}
\title{Effective Sample Size due to Autocorrelation}
\description{
  This function may be used to estimate the effective sample size (ESS)
  (not to be confused with Elliptical Slice Sampling) of a continuous
  target distribution, where the sample size is reduced by
  autocorrelation. ESS is a measure of how well each continuous chain is
  mixing.

  ESS is a univariate function that is often applied to each continuous,
  marginal posterior distribution. A multivariate form is not
  included. By chance alone due to multiple independent tests, 5\% of
  the continuous parameters may indicate that ESS is below a user
  threshold of acceptability, such as 100, even when above the
  threshold. Assessing convergence is difficult.
}
\usage{ESS(x)}
\arguments{
  \item{x}{This required argument is a vector or matrix of posterior
    samples.}
}
\details{
  Effective Sample Size (ESS) was recommended by Radford Neal in the
  panel discussion of Kass et al. (1998). When a continuous, marginal
  posterior distribution is sampled with a Markov chain Monte Carlo
  (MCMC) algorithm, there is usually autocorrelation present in the
  samples. More autocorrelation is associated with less posterior
  sampled information, because the information in the samples is
  autocorrelated, or put another way, successive samples are not
  independent from earlier samples. This reduces the effective sample
  size of, and precision in representing, the continuous, marginal
  posterior distribution. \code{ESS} is one of the criteria in the
  \code{\link{Consort}} function, where stopping the MCMC updates is
  not recommended until \code{ESS} \eqn{\ge 100}. Although the need
  for precision of each modeler differs with each model, it is often
  a good goal to obtain \code{ESS} \eqn{= 1000}.

  \code{ESS} is related to the integrated autocorrelation time (see
  \code{\link{IAT}} for more information).

  ESS is usually defined as

  \deqn{\mathrm{ESS}(\theta) = \frac{S}{1 + 2 \sum^{\infty}_{k=1} \rho_k
    (\theta)},}{ESS(theta) = S / (1 + 2 sum[k] rho[k] (theta)),}

  where \eqn{S}{S} is the number of posterior samples,
  \eqn{\rho_k}{rho[k]} is the autocorrelation at lag \eqn{k}{k}, and
  \eqn{\theta}{theta} is the vector of marginal posterior samples. The
  infinite sum is often truncated at lag \eqn{k}{k} when
  \eqn{\rho_k (\theta) < 0.05}{rho[k](theta) < 0.05}. Just as with the
  \code{effectiveSize} function in the \code{coda} package, the
  \code{AIC} argument in the \code{ar} function is used to estimate the
  order.

  ESS is a measure of how well each continuous chain is mixing, and a
  continuous chain mixes better when in the target distribution. This
  does not imply that a poorly mixing chain still searching for its
  target distribution will suddenly mix well after finding it, though
  mixing should improve. A poorly mixing continuous chain does not
  necessarily indicate problems. A smaller ESS is often due to
  correlated parameters, and is commonly found with scale
  parameters. Posterior correlation may be obtained from the
  \code{\link{PosteriorChecks}} function, and plotted with the 
  \code{\link{plotMatrix}} function. Common remedies for poor mixing
  include re-parameterizing the model or trying a different MCMC
  algorithm that better handles correlated parameters. Slow mixing is
  indicative of an inefficiency in which a continuous chain takes longer
  to find its target distribution, and once found, takes longer to
  explore it. Therefore, slow mixing results in a longer required
  run-time to find and adequately represent the continuous target
  distribution, and increases the chance that the user may make
  inferences from a less than adequate representation of the continuous
  target distribution.

  There are many methods of re-parameterization to improve mixing. It
  is helpful when predictors are centered and scaled, such as with the
  \code{\link{CenterScale}} function. Parameters for predictors are
  often assigned prior distributions that are independent per parameter,
  in which case an exchangeable prior distribution or a multivariate
  prior distribution may help. If a parameter with poor mixing is
  bounded with the \code{\link{interval}} function, then
  transforming it to the real line (such as with a log transformation
  for a scale parameter) is often helpful, since constraining a
  parameter to an interval often reduces ESS. Another method is to
  re-parameterize so that one or more latent variables represent the
  process that results in slow mixing. Such re-parameterization uses
  data augmentation.

  This is numerically the same as the \code{effectiveSize} function in
  the \code{coda} package, but programmed to accept a simple vector or
  matrix so it does not require an \code{mcmc} or \code{mcmc.list}
  object, and the result is bound to be less than or equal to the
  original number of samples.
}
\value{
  A vector is returned, and each element is the effective sample size
  (ESS) for a corresponding column of \code{x}, after autocorrelation has
  been taken into account.
}
\references{
  Kass, R.E., Carlin, B.P., Gelman, A., and Neal, R. (1998). "Markov
  Chain Monte Carlo in Practice: A Roundtable Discussion". \emph{The
  American Statistician}, 52, p. 93--100. 
}
\seealso{
  \code{\link{CenterScale}},
  \code{\link{Consort}},
  \code{\link{IAT}},
  \code{\link{interval}},
  \code{\link{LaplacesDemon}},
  \code{\link{plotMatrix}}, and
  \code{\link{PosteriorChecks}}.
}
\keyword{Diagnostic}
\keyword{MCMC}