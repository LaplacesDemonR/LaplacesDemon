\name{Gelman.Diagnostic}
\alias{Gelman.Diagnostic}
\title{Gelman and Rubin's MCMC Convergence Diagnostic}
\description{
  Gelman and Rubin (1992) proposed a general approach to monitoring
  convergence of MCMC output in which \eqn{m > 1} parallel chains are
  updated with initial values that are overdispersed relative to each
  target distribution, which must be normally distributed. Convergence
  is diagnosed when the chains have `forgotten' their initial values,
  and the output from all chains is indistinguishable. The
  \code{Gelman.Diagnostic} function makes a comparison of within-chain
  and between-chain variances, and is similar to a classical analysis of
  variance. A large deviation between these two variances indicates
  non-convergence.

  This diagnostic is popular as a stopping rule, though it requires
  parallel chains. The \code{\link{LaplacesDemon.hpc}} function is an
  extension of \code{\link{LaplacesDemon}} to enable parallel chains.
  As an alternative, the popular single-chain stopping rule is based on
  \code{\link{MCSE}}.
}
\usage{Gelman.Diagnostic(x, confidence=0.95, transform=FALSE)}
\arguments{
  \item{x}{This required argument accepts an object of class
    \code{demonoid.hpc}, or a list of multiple objects of class
    \code{demonoid}, where the number of components in the list 
    is the number of chains.}
  \item{confidence}{This is the coverage probability of the confidence
    interval for the potential scale reduction factor (PSRF).}
  \item{transform}{Logical. If \code{TRUE}, then marginal posterior
    distributions in \code{x} may be transformed to improve the
    normality of the distribution, which is assumed. A log-transform is
    applied to marginal posterior distributions in the interval \eqn{(0,
      \infty]}{(0, Inf]}, or a logit-transform is applied to marginal
    posterior distributions in the interval \eqn{(0,1)}.}
}
\value{
  A list is returned with the following components:
  
  \item{PSRF}{This is a list containing the point-estimates of the
    potential scale reduction factor (labelled \code{Point Est.}) and
    the associated upper confidence limits (labelled \code{Upper C.I.}).}
  \item{MPSRF}{This is the point-estimate of the multivariate potential
    scale reduction factor.}
}
\details{
  To use the \code{Gelman.Diagnostic} function, the user must first have
  multiple MCMC chains for the same model, and three chains is usually
  sufficient. The easiest way to obtain multiple chains is with the
  \code{\link{LaplacesDemon.hpc}} function.

  Although the \code{\link{LaplacesDemon}} function does not
  simultaneously update multiple MCMC chains, it is easy enough to
  obtain multiple chains, and if the computer has multiple processors
  (which is common), then multiple chains may be obtained simultaneously
  as follows. The model file may be opened in separate, concurrent R
  sessions, and it is recommended that a maximum number of sessions is
  equal to the number of processors, minus one. Each session constitutes
  its own chain, and the code is identical, except the initial values
  should be randomized with the \code{\link{GIV}} function so the chains
  begin in different places. The resulting object of class
  \code{demonoid} for each chain is saved, all objects are read into one
  session, put into a list, and passed to the \code{Gelman.Diagnostic}
  function.

  Initial values must be overdispersed with respect to each target
  distribution, though these distributions are unknown in the
  beginning. Since the \code{Gelman.Diagnostic} function relies heavily
  on overdispersion with respect to the target distribution, the user
  should consider using MCMC twice, first to estimate the target
  distributions, and secondly to overdisperse initial values with
  respect to them. This may help identify multimodal target
  distributions. If multiple modes are found, it remain possible that
  more modes exist. When multiple modes are found, and if chains are
  combined with the \code{\link{Combine}} function, each mode is
  probably not represented in a proportion correct to the distribution.
    
  The `potential scale reduction factor' (PSRF) is an estimated factor
  by which the scale of the current distribution for the target
  distribution might be reduced if the simulations were continued for
  an infinite number of iterations. Each PSRF declines to 1 as the
  number of iterations approaches infinity. PSRF is also often
  represented as R-hat. PSRF is calculated for each marginal posterior
  distribution in \code{x}, together with upper and lower confidence
  limits. Approximate convergence is diagnosed when the upper limit is
  close to 1. The recommended proximity of each PSRF to 1 varies with
  each problem, but a general goal is to achieve PSRF < 1.1. PSRF is an
  estimate of how much narrower the posterior might become with an
  infinite number of iterations. When PSRF = 1.1, for example, it may be
  interpreted as a potential reduction of 10\% in posterior interval
  width, given infinite iterations. The multivariate form bounds above
  the potential scale reduction factor for any linear combination of the
  (possibly transformed) variables.
  
  The confidence limits are based on the assumption that the
  target distribution is stationary and normally distributed. The
  \code{transform} argument may be used to improve the normal
  approximation.

  A large PSRF indicates that the between-chain variance is
  substantially greater than the within-chain variance, so that longer
  simulation is needed. If a PSRF is close to 1, then the associated
  chains are likely to have converged to one target distribution. A
  large PSRF (perhaps generally when a PSRF > 1.2) indicates convergence
  failure, and can indicate the presence of a multimodal marginal
  posterior distribution in which different chains may have converged
  to different local modes (see \code{\link{is.multimodal}}), or the
  need to update the associated chains longer, because burn-in (see
  \code{\link{burnin}}) has yet to be completed.
  
  The \code{Gelman.Diagnostic} is essentially the same as the
  \code{gelman.diag} function in the \code{coda} package, but here it is 
  programmed to work with objects of class \code{demonoid}.

  There are two ways to estimate the variance of the stationary
  distribution: the mean of the empirical variance within each chain,
  \eqn{W}, and the empirical variance from all chains combined, which
  can be expressed as

  \deqn{ \widehat{\sigma}^2 = \frac{(n-1) W}{n} +
    \frac{B}{n}}{sigma.hat^2 = (n-1)W/n + B/n} 

  where \eqn{n} is the number of iterations and \eqn{B/n} is the
  empirical between-chain variance.

  If the chains have converged, then both estimates are
  unbiased. Otherwise the first method will \emph{underestimate} the
  variance, since the individual chains have not had time to range all
  over the stationary distribution, and the second method will
  \emph{overestimate} the variance, since the initial values were
  chosen to be overdispersed (and this assumes the target distribution
  is known, see above).

  This convergence diagnostic is based on the assumption that each
  target distribution is normal. A Bayesian probability interval (see
  \code{\link{p.interval}}) can be constructed using a t-distribution
  with mean
  
  \deqn{\widehat{\mu}=\mbox{Sample mean of all chains combined,}}{mu.hat
    = Sample mean of all chains combined,}

  variance

  \deqn{\widehat{V} = \widehat{\sigma}^2 + \frac{B}{mn},}{V.hat =
    sigma.hat2 + B/(mn),}

  and degrees of freedom estimated by the method of moments

  \deqn{d = \frac{2\widehat{V}^2}{\mbox{Var}(\widehat{V})}}{d = 2*V.hat^2 /
    Var(V.hat)}

  Use of the t-distribution accounts for the fact that the mean and
  variance of the posterior distribution are estimated. The convergence
  diagnostic itself is

  \deqn{R=\sqrt{\frac{(d+3) \widehat{V}}{(d+1)W}}}{R=sqrt((d+3) V.hat /((d+1)W)}

  Values substantially above 1 indicate lack of convergence. If the
  chains have not converged, then Bayesian probability intervals based
  on the t-distribution are too wide, and have the potential to shrink
  by this factor if the MCMC run is continued.

  The multivariate version of Gelman and Rubin's diagnostic was proposed
  by Brooks and Gelman (1998). Unlike the univariate proportional scale
  reduction factor, the multivariate version does not include an
  adjustment for the estimated number of degrees of freedom.
}
\references{
  Brooks, S.P. and Gelman, A. (1998). "General Methods for Monitoring
  Convergence of Iterative Simulations". \emph{Journal of Computational
  and Graphical Statistics}, 7, p. 434--455.

  Gelman, A. and Rubin, D.B. (1992). "Inference from Iterative Simulation
  using Multiple Sequences". \emph{Statistical Science}, 7, p. 457--511.
}
\seealso{
  \code{\link{Combine}},
  \code{\link{GIV}},
  \code{\link{is.multimodal}},
  \code{\link{LaplacesDemon}},
  \code{\link{LaplacesDemon.hpc}},
  \code{\link{MCSE}}, and
  \code{\link{p.interval}}.
}
\examples{
#library(LaplacesDemon)
###After updating multiple chains with LaplacesDemon.hpc, do:
#Gelman.Diagnostic(Fit)
}
\keyword{Diagnostic}
\keyword{MCMC}
\keyword{Multiple Chains}
\keyword{Parallel Chains}
