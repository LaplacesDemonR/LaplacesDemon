\name{PosteriorChecks}
\alias{PosteriorChecks}
\title{Posterior Checks}
\description{
  Not to be confused with posterior predictive checks, this function
  provides additional information about the marginal posterior
  distributions of continuous parameters, such as the probability that
  each posterior coefficient of the parameters (referred to generically
  as \eqn{\theta}{theta}), is greater than zero
  [\eqn{p(\theta > 0)}{p(theta > 0)}], the estimated number of modes,
  the kurtosis and skewness of the posterior distributions, the burn-in
  of each chain (for MCMC only), integrated autocorrelation time,
  independent samples per minute, and acceptance rate. A posterior
  correlation matrix is provided only for objects of class
  \code{demonoid} or \code{pmc}.

  For discrete parameters, see the \code{\link{Hangartner.Diagnostic}}.
}
\usage{
PosteriorChecks(x, Parms)
}
\arguments{
  \item{x}{This required argument accepts an object of class
    \code{demonoid}, \code{iterquad}, \code{laplace}, \code{pmc}, or
    \code{vb}.}
  \item{Parms}{This argument accepts a vector of quoted strings to be
    matched for selecting parameters. This argument defaults to
    \code{NULL} and selects every parameter. Each quoted string is
    matched to one or more parameter names with the \code{grep}
    function. For example, if the user specifies \code{Parms=c("eta",
      "tau")}, and if the parameter names are beta[1], beta[2], eta[1],
    eta[2], and tau, then all parameters will be selected, because the
    string \code{eta} is within \code{beta}. Since \code{grep} is used,
    string matching uses regular expressions, so beware of
    meta-characters, though these are acceptable: ".", "[", and "]".}
}
\details{\code{PosteriorChecks} is a supplemental function that returns
  a list with two components. Following is a summary of popular uses of
  the \code{PosteriorChecks} function.

  First (and only for MCMC users), the user may be considering the
  current MCMC algorithm versus others. In this case, the
  \code{PosteriorChecks} function is often used to find the two MCMC
  chains with the highest \code{\link{IAT}}, and these chains are
  studied for non-randomness with a joint trace plot, via the
  \code{\link{joint.density.plot}} function. The best algorithm has the
  chains with the highest independent samples per minute (ISM).

  Posterior correlation may be studied between model updates as well as
  after a model seems to have converged. While frequentists consider
  multicollinear predictor variables, Bayesians tend to consider
  posterior correlation of the parameters. Models with multicollinear
  parameters take more iterations to converge. Hierarchical models often
  have high posterior correlations. Posterior correlation often
  contributes to a lower effective sample size (\code{\link{ESS}}).
  Common remedies include transforming the predictors,
  re-parameterization to reduce posterior correlation, using WIPs
  (Weakly-Informative Priors), or selecting a different numerical
  approximation algorithm. An example of re-parameterization is to
  constrain related parameters to sum to zero. Another approach is to
  specify the parameters according to a multivariate distribution that
  is assisted by estimating a covariance matrix. Some algorithms are
  more robust to posterior correlation than others. For example,
  posterior correlation should generally be less problematic for twalk
  than AMWG in \code{\link{LaplacesDemon}}. Posterior correlation may be
  plotted with the \code{\link{plotMatrix}} function, and may be useful
  for blocking parameters. For more information on blockwise sampling,
  see the \code{\link{Blocks}} function.

  After a user is convinced of the applicability of the current MCMC
  algorithm, and that the chains have converged, \code{PosteriorChecks}
  is often used to identify multimodal marginal posterior distributions
  for further study or model re-specification.

  Although many marginal posterior distributions appear normally
  distributed, there is no such assumption. Nonetheless, a marginal
  posterior distribution tends to be distributed the same as its prior
  distribution. If a parameter has a prior specified with a Laplace
  distribution, then the marginal posterior distribution tends also to
  be Laplace-distributed. In the common case of normality, kurtosis
  and skewness may be used to identify discrepancies between the prior
  and posterior, and perhaps this should be called a `prior-posterior
  check'.

  Lastly, parameter importance may be considered, in which case it is
  recommended to be considered simultaneously with variable importance
  from the \code{\link{Importance}} function.
}
\value{
  \code{PosteriorChecks} returns an object of class
  \code{posteriorchecks} that is a list with the following components:
  \item{Posterior.Correlation}{
    This is a correlation matrix of the parameters selected with the
    \code{Parms} argument. This component is returned as \code{NA} for
    objects of classes \code{"laplace"} or \code{"vb"}.}
  \item{Posterior.Summary}{This is a matrix in which each row is a
    parameter and there are eight columns: p(theta > 0), N.Modes,
    Kurtosis, Skewness, Burn-In, IAT, ISM, and AR. The first column,
    p(theta > 0), indicates parameter importance by reporting how much
    of the distribution is greater than zero. An important parameter
    distribution will have a result at least as extreme as 0.025 or
    0.975, and an unimportant parameter distribution is centered at
    0.5. This is not the importance of the associated variable relative
    to how well the model fits the data. For variable importance, see
    the \code{\link{Importance}} function. The second column, N.Modes,
    is the number of modes, estimated with the \code{\link{Modes}}
    function. Kurtosis and skewness are useful posterior checks that may
    suggest that a posterior distribution is non-normal or does not fit
    well with a distributional assumption, assuming a distributional
    assumption exists, which it may not. The burn-in is estimated for
    each chain (only for objects of class \code{demonoid} with the
    \code{\link{burnin}} function. The integrated autocorrelation
    time is estimated with \code{\link{IAT}}. The number of independent
    samples per minute (ISM) is calculated for objects of class
    \code{"demonoid"} as \code{\link{ESS}} divided by minutes. Lastly,
    the local acceptance rate of each MCMC chain is calculated with the
    \code{\link{AcceptanceRate}} function, and is set to 1 for objects
    of class \code{iterquad}, \code{laplace}, \code{pmc}, or \code{vb}.}
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{AcceptanceRate}},
  \code{\link{Blocks}},
  \code{\link{burnin}},
  \code{\link{ESS}},
  \code{\link{Hangartner.Diagnostic}},
  \code{\link{joint.density.plot}},
  \code{\link{IAT}},
  \code{\link{Importance}},
  \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{Modes}},
  \code{\link{plotMatrix}},
  \code{\link{PMC}}, and
  \code{\link{VariationalBayes}}.
}
\examples{### See the LaplacesDemon function for an example.}
\keyword{Diagnostic}
\keyword{MCMC}
\keyword{Monte Carlo}
\keyword{Multicollinearity}
\keyword{Posterior Correlation}
\keyword{Utility}