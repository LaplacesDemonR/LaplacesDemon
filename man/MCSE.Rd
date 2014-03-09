\name{MCSE}
\alias{MCSE}
\alias{MCSS}
\title{Monte Carlo Standard Error}
\description{
  Monte Carlo Standard Error (MCSE) is an estimate of the inaccuracy of
  Monte Carlo samples, usually regarding the expectation of posterior
  samples, \eqn{\mathrm{E}(\theta)}{E(theta)}, from Monte Carlo or
  Markov chain Monte Carlo (MCMC) algorithms, such as with the
  \code{\link{LaplacesDemon}} or \code{\link{LaplacesDemon.hpc}}
  functions. MCSE approaches zero as the number of independent posterior
  samples approaches infinity. MCSE is essentially a standard deviation
  around the posterior mean of the samples,
  \eqn{\mathrm{E}(\theta)}{E(theta)}, due to uncertainty associated with
  using an MCMC algorithm, or Monte Carlo methods in general.

  The acceptable size of the MCSE depends on the acceptable uncertainty
  associated around the marginal posterior mean,
  \eqn{\mathrm{E}(\theta)}{E(theta)}, and the goal of inference. It has
  been argued that MCSE is generally unimportant when the goal of
  inference is \eqn{\theta}{theta} rather than
  \eqn{\mathrm{E}(\theta)}{E(theta)} (Gelman et al., 2004, p. 277), and
  that a sufficient \code{\link{ESS}} is more important. Others perceive
  MCSE to be a vital part of reporting any Bayesian model, and as a
  stopping rule (Flegal et al., 2008).

  In \code{\link{LaplacesDemon}}, MCSE is part of the posterior
  summaries because it is easy to estimate, and Laplace's Demon prefers
  to continue updating until each MCSE is less than 6.27\% of its
  associated marginal posterior standard deviation (for more information
  on this stopping rule, see the \code{\link{Consort}} function), since
  MCSE has been demonstrated to be an excellent stopping rule.

  Acceptable error may be specified, if known, in the \code{MCSS}
  (Monte Carlo Sample Size) function to estimate the required number of
  posterior samples.

  \code{MCSE} is a univariate function that is often applied to each
  marginal posterior distribution. A multivariate form is not
  included. By chance alone due to multiple independent tests, 5\% of
  the parameters should indicate unacceptable MSCEs, even when
  acceptable. Assessing convergence is difficult.
}
\usage{
MCSE(x, method="IMPS", batch.size="sqrt", warn=FALSE)
MCSS(x, a)
}
\arguments{
  \item{x}{This is a vector of posterior samples for which MCSE or MCSS
    will be estimated.}
  \item{a}{This is a scalar argument of acceptable error for the mean of
    \code{x}, and \code{a} must be positive. As acceptable error
    decreases, the required number of samples increases.}
  \item{method}{This is an optional argument for the method of MCSE
    estimation, and defaults to Geyer's \code{"IMPS"} method. Optional
    methods include \code{"sample.variance"} and \code{"batch.mean"}.
    Note that \code{"batch.mean"} is recommended only when the number of
    posterior samples is at least 1,000.}
  \item{batch.size}{This is an optional argument that corresponds only
    with \code{method="batch.means"}, and determines either the size of
    the batches (accepting a numerical argument) or the method of
    creating the size of batches, which is either \code{"sqrt"} or
    \code{"cuberoot"}, and refers to the length of \code{x}. The default
    argument is \code{"sqrt"}.}
  \item{warn}{Logical. If \code{warn=TRUE}, then a warning is provided
    with \code{method="batch.means"} whenever posterior sample size is
    less than 1,000, or a warning is produced when more autcovariance
    is recommended with \code{method="IMPS"}.}
}
\details{
  The default method for estimating MCSE is Geyer's Initial Monotone
  Positive Sequence (IMPS) estimator (Geyer, 1992), which takes the
  asymptotic variance into account and is time-series based. This method
  goes by other names, such as Initial Positive Sequence (IPS).

  The simplest method for estimating MCSE is to modify the formula for
  standard error, \eqn{\sigma(\textbf{x}) / \sqrt{N}}{sigma(x) /
  sqrt(N)}, to account for non-independence in the sequence
  \eqn{\textbf{x}}{x} of posterior samples. Non-independence is
  estimated with the \code{ESS} function for Effective Sample Size (see
  the \code{\link{ESS}} function for more details), where \eqn{M =
  ESS(\textbf{x})}{m = ESS(x)}, and MCSE is
  \eqn{\sigma(\textbf{x}) / \sqrt{M}}{sigma(x) / sqrt(M)}. Although this
  is the fastest and easiest method of estimation, it does not
  incorporate an estimate of the asymptotic variance of
  \eqn{\textbf{x}}{x}.

  The batch means method (Jones et al., 2006; Flegal et al., 2008)
  separates elements of \eqn{\textbf{x}}{x} into batches and estimates
  MCSE as a function of multiple batches. This method is excellent, but
  is not recommended when the number of posterior samples is less than
  1,000. These journal articles also assert that MCSE is a better
  stopping rule than MCMC convergence diagnostics.

  The \code{MCSS} function estimates the required number of posterior
  samples, given the user-specified acceptable error, posterior samples
  \code{x}, and the observed variance (rather than asymptotic
  variance). Due to the observed variance, this is a rough estimate.
}
\references{
  Flegal, J.M., Haran, M., and Jones, G.L. (2008). "Markov chain Monte
  Carlo: Can We Trust the Third Significant Figure?". \emph{Statistical
  Science}, 23, p. 250--260.

  Gelman, A., Carlin, J., Stern, H., and Rubin, D. (2004). "Bayesian
  Data Analysis, Texts in Statistical Science, 2nd ed.". Chapman and
  Hall, London.

  Geyer, C.J. (1992). "Practical Markov Chain Monte Carlo".
  \emph{Statistical Science}, 7, 4, p. 473--483.
  
  Jones, G.L., Haran, M., Caffo, B.S., and Neath, R. (2006). "Fixed-Width
  Output Analysis for Markov chain Monte Carlo". \emph{Journal of the
  American Statistical Association}, 101(1), p. 1537--1547.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{Consort}},
  \code{\link{ESS}},
  \code{\link{LaplacesDemon}}, and
  \code{\link{LaplacesDemon.hpc}}.
}
\examples{
library(LaplacesDemon)
x <- rnorm(1000)
MCSE(x)
MCSE(x, method="batch.means")
MCSE(x, method="sample.variance")
MCSS(x, a=0.01)
}
\keyword{Diagnostic}
\keyword{MCMC}
\keyword{MCSE}
\keyword{MCSS}
\keyword{Utility}
