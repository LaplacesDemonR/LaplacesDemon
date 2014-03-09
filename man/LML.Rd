\name{LML}
\alias{LML}
\title{Logarithm of the Marginal Likelihood}
\description{
  This function approximates the logarithm of the marginal likelihood
  (LML), where the marginal likelihood is also called the integrated
  likelihood or the prior predictive distribution of \eqn{\textbf{y}}{y}
  in Bayesian inference. The marginal likelihood is

  \deqn{p(\textbf{y}) = \int p(\textbf{y} | \Theta)p(\Theta) d\Theta}{p(y)
    = integral p(y | Theta)p(Theta) d Theta}

  The prior predictive distribution indicates what \eqn{\textbf{y}}{y}
  should look like, given the model, before \eqn{\textbf{y}}{y} has been
  observed. The presence of the marginal likelihood of
  \eqn{\textbf{y}}{y} normalizes the joint posterior distribution,
  \eqn{p(\Theta|\textbf{y})}{p(Theta|y)}, ensuring it is a proper
  distribution and integrates to one (see \code{\link{is.proper}}). The
  marginal likelihood is the denominator of Bayes' theorem, and is often
  omitted, serving as a constant of proportionality. Several methods of
  approximation are available.
}
\usage{
LML(Model=NULL, Data=NULL, Modes=NULL, theta=NULL, LL=NULL, Covar=NULL,
method="NSIS")
}
\arguments{
  \item{Model}{This is the model specification for the model that was
    updated either in \code{\link{IterativeQuadrature}},
    \code{\link{LaplaceApproximation}}, \code{\link{LaplacesDemon}},
    \code{\link{LaplacesDemon.hpc}}, or \code{\link{VariationalBayes}}.
    This argument is used only with the \code{LME} method.}
  \item{Data}{This is the list of data passed to the model
    specification. This argument is used only with the \code{LME}
    method.}
  \item{Modes}{This is a vector of the posterior modes (or medians, in
    the case of MCMC). This argument is used only with the \code{GD} or
    \code{LME} methods.}
  \item{theta}{This is a matrix of posterior samples (parameters only),
    and is specified only with the \code{GD}, \code{HME}, or
    \code{NSIS} methods.}
  \item{LL}{This is a vector of MCMC samples of the log-likelihood, and
    is specified only with the \code{GD}, code{HME}, or \code{NSIS}
    methods.}
  \item{Covar}{This argument accepts the covariance matrix of the
    posterior modes, and is used only with the \code{GD} or \code{LME}
    methods.}
  \item{method}{The method may be \code{"GD"}, \code{"HME"},
    \code{"LME"}, or \code{"NSIS"}, and defaults to \code{"NSIS"}.
    \code{"GD"} uses the Gelfand-Dey estimator, \code{"HME"} uses the
    Harmonic Mean Estimator, \code{"LME"} uses the Laplace-Metropolis
    Estimator, and \code{"NSIS"} uses nonparametric self-normalized
    importance sampling (NSIS).}
}
\details{
  Generally, a user of \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}}, \code{\link{LaplacesDemon.hpc}},
  \code{\link{PMC}}, or \code{\link{VariationalBayes}} does not need to
  use the \code{LML} function, because these methods already include
  it. However, \code{LML} may be called by the user, should the user
  desire to estimate the logarithm of the marginal likelihood with a
  different method, or with non-stationary chains. The
  \code{\link{LaplacesDemon}} and \code{\link{LaplacesDemon.hpc}}
  functions only call \code{LML} when all parameters are stationary, and
  only with non-adaptive algorithms.

  The \code{GD} method, where GD stands for Gelfand-Dey (1994), is a
  modification of the harmonic mean estimator (HME) that results in a
  more stable estimator of the logarithm of the marginal
  likelihood. This method is unbiased, simulation-consistent, and
  usually satisfies the Gaussian central limit theorem.

  The \code{HME} method, where HME stands for harmonic mean estimator,
  of Newton-Raftery (1994) is the easiest, and therefore fastest,
  estimation of the logarithm of the marginal likelihood. However, it is
  an unreliable estimator and should be avoided, because small
  likelihood values can overly influence the estimator, variance is
  often infinite, and the Gaussian central limit theorem is usually not
  satisfied. It is included here for completeness. There is not a
  function in this package that uses this method by default. Given
  \eqn{N} samples, the estimator is
  \eqn{1/[\frac{1}{N} \sum_N \exp(-LL)]}{1 / [1/N sum_N exp(-LL)]}.

  The \code{LME} method uses the Laplace-Metropolis Estimator (LME), in
  which the estimation of the Hessian matrix is approximated
  numerically. It is the slowest method here, though it returns an
  estimate in more cases than the other methods. The supplied
  \code{Model} specification must be executed a number of times equal
  to \eqn{k^2 \times 4}{k^2 x 2}, where \eqn{k} is the number of
  parameters. In large dimensions, this is very slow. The
  Laplace-Metropolis Estimator is inappropriate with hierarchical
  models. The \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}}, and \code{\link{VariationalBayes}}
  functions use \code{LME} when it has converged and \code{sir=FALSE},
  in which case it uses the posterior means or modes, and is itself
  Laplace Approximation.

  The Laplace-Metropolis Estimator (LME) is the logarithmic form of
  equation 4 in Lewis and Raftery (1997). In a non-hierarchical model,
  the marginal likelihood may easily be approximated with the
  Laplace-Metropolis Estimator for model \eqn{m} as
  
  \deqn{p(\textbf{y}|m) =
  (2\pi)^{d_m/2}|\Sigma_m|^{1/2}p(\textbf{y}|\Theta_m,m)p(\Theta_m|m)}{p(y|m)
  = (2*pi)^(d_m/2) |Sigma_m|^(1/2) p(y|Theta_m, m)p(Theta_m|m)}

  where \eqn{d} is the number of parameters and \eqn{\Sigma}{Sigma} is
  the inverse of the negative of the approximated Hessian matrix of
  second derivatives.
  
  As a rough estimate of Kass and Raftery (1995), LME is worrisome when
  the sample size of the data is less than five times the number of
  parameters, and LME should be adequate in most problems when the
  sample size of the data exceeds twenty times the number of parameters
  (p. 778).

  The \code{NSIS} method is essentially the \code{MarginalLikelihood}
  function in the \code{MargLikArrogance} package. After \code{HME},
  this is the fastest method available here. The
  \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}}, and \code{\link{VariationalBayes}}
  functions use \code{NSIS} when converged and \code{sir=TRUE}. The
  \code{\link{LaplacesDemon}}, \code{\link{LaplacesDemon.hpc}}, and
  \code{\link{PMC}} functions use \code{NSIS}. At least 301 stationary
  samples are required, and the number of parameters cannot exceed half
  the number of stationary samples.
}
\value{
  \code{LML} returns a list with two components:
  \item{LML}{
    This is an approximation of the logarithm of the marginal
    likelihood (LML), which is notoriously difficult to estimate. For this
    reason, several methods are provided. The marginal likelihood is
    useful when comparing models, such as with Bayes factors in the
    \code{\link{BayesFactor}} function. When the method fails, \code{NA}
    is returned, and it is most likely that the joint posterior is
  improper (see \code{\link{is.proper}}).}
  \item{VarCov}{
    This is a variance-covariance matrix, and is the negative inverse of
    the Hessian matrix, if estimated. The \code{GD}, \code{HME}, and
    \code{NSIS} methods do not estimate \code{VarCov}, and return
    \code{NA}.}
}
\references{
  Gelfand, A.E. and Dey, D.K. (1994). "Bayesian Model Choice:
  Asymptotics and Exact Calculations". \emph{Journal of the Royal
  Statistical Society}, Series B 56, p. 501--514.

  Kass, R.E. and Raftery, A.E. (1995). "Bayes Factors". \emph{Journal
  of the American Statistical Association}, 90(430), p. 773--795.

  Lewis, S.M. and Raftery, A.E. (1997). "Estimating Bayes Factors via
  Posterior Simulation with the Laplace-Metropolis Estimator".
  \emph{Journal of the American Statistical Association}, 92,
  p. 648--655.

  Newton, M.A. and Raftery, A.E. (1994). "Approximate Bayesian
  Inference by the Weighted Likelihood Bootstrap". \emph{Journal of the
  Royal Statistical Society}, Series B 3, p. 3--48.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{BayesFactor}},
  \code{\link{is.proper}},
  \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{LaplacesDemon.hpc}},
  \code{\link{PMC}}, and
  \code{\link{VariationalBayes}}.
}
\examples{
### If a model object were created and called Fit, then:
#
### Applying HME to an object of class demonoid or pmc:
#LML(LL=Fit$Deviance*(-1/2), method="HME")
#
### Applying LME to an object of class demonoid:
#LML(Model, MyData, Modes=apply(Fit$Posterior1, 2, median), method="LME")
#
### Applying NSIS to an object of class demonoid
#LML(theta=Fit$Posterior1, LL=Fit$Deviance*-(1/2), method="NSIS")
#
### Applying LME to an object of class iterquad:
#LML(Model, MyData, Modes=Fit$Summary1[,1], method="LME")
#
### Applying LME to an object of class laplace:
#LML(Model, MyData, Modes=Fit$Summary1[,1], method="LME")
#
### Applying LME to an object of class vb:
#LML(Model, MyData, Modes=Fit$Summary1[,1], method="LME")
}
\keyword{Harmonic Mean Estimator}
\keyword{Integrated Likelihood}
\keyword{Laplace-Metropolis Estimator}
\keyword{Marginal Likelihood}
\keyword{Nonparametric Self-Normalized Importance Sampling}
\keyword{Prior Predictive Distribution}
