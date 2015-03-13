\name{Consort}
\alias{Consort}
\title{Consort with Laplace's Demon}
\description{
  This may be used to consort with Laplace's Demon regarding an object
  of class \code{demonoid}. Laplace's Demon will offer suggestions.
}
\usage{Consort(object)}
\arguments{
  \item{object}{This required argument is an object of class
    \code{demonoid}. For more information, see the
    \code{\link{LaplacesDemon}} function.}
}
\details{
  First, \code{Consort} calls \code{print.demonoid}, which prints most
  of the components to the screen from the supplied object of class
  \code{demonoid}.

  Second, Laplace's Demon considers a combination of five conditions
  when making the largest part of its suggestion. These conditions are:
  the algorithm, acceptance rate, MCSE, ESS, and stationarity. Other
  things are considered as well, such as the recommended thinning value
  is used to suggest a new number of iterations, how fast the algorithm
  is expected to be, and if the condition of diminishing adaptation
  (also called the vanishing adaptation condition) was met (for an
  adaptive algorithm). Diminishing adaptation occurs only when the
  absolute value of the proposed variances trends downward (toward zero)
  over the course of all adaptations. When an algorithm is adaptive and
  it does not have diminishing adaptations, the \code{Consort} function
  will suggest a different adaptive algorithm. The \code{Periodicity}
  argument is suggested to be set equal to the value of
  \code{Rec.Thinning}.

  Appeasement applies only when all parameters are continuous.The
  \code{\link{Hangartner.Diagnostic}} should be considered for discrete
  parameters.

  Appeasement Conditions
  \itemize{
    \item Algorithm: The final algorithm must be non-adaptive, so that
    the Markov property holds. This is conservative. A user may have
    an adaptive (non-final) algorithm in which adaptations in the latest
    update are stationary, or no longer diminishing. Laplace's Demon is
    unaware of previous updates, and conservatively interprets this as
    failing to meet the condition of diminishing adaptation, when the
    output may be satisfactory. On the other hand, if the adaptive
    algorithm has essentially stopped adapting, and if there is a
    non-adaptive version, then the user should consider switching to
    the non-adaptive algorithm. User discretion is advised.
    \item Acceptance Rate: The acceptance rate is considered
    satisfactory if it is within the interval [15\%,50\%] for most
    algorithms. Some algorithms have different recommended intervals.
    \item MCSE: The Monte Carlo Standard Error (MCSE) is considered
    satisfactory for each target distribution if it is less than 6.27\%
    of the standard deviation of the target distribution. This allows
    the true mean to be within 5\% of the area under a Gaussian
    distribution around the estimated mean. The \code{\link{MCSE}}
    function is used. Toft et al. (2007) propose a stricter criterion of
    5\%. The criterion of 6.27\% for this stopping rule is arbitrary,
    and may be too lenient or strict, depending on the needs of the
    user. Nonetheless, it has performed well, and this type of stopping
    rule has been observed to perform better than MCMC convergence
    diagnostics (Flegal et al., 2008).
    \item ESS: The effective sample size (ESS) is considered
    satisfactory for each target distribution if it is at least 100,
    which is usually enough to describe 95\% probability intervals (see
    \code{\link{p.interval}} and \code{\link{LPL.interval}} for more
    information). The \code{\link{ESS}} function is used. When this
    criterion is unmet, the name of the worst mixing chain in Summary1
    appears.
    \item Stationarity: Each target distribution is considered
    satisfactory if it is estimated to be stationary with the
    \code{\link{BMK.Diagnostic}} function.
}
  Bear in mind that the MCSE, ESS, and stationarity criteria are all
  univariate measures applied to each marginal posterior
  distribution. Multivariate forms are not included. By chance alone due
  to multiple independent tests, 5\% of these diagnostics should
  indicate non-convergence when 'convergence' exists. In contrast, even
  one non-convergent nuisance parameter is associated with
  non-convergence in all other parameters. Assessing convergence is
  difficult.

  If all five conditions are satisfactory, then Laplace's Demon is
  appeased. Otherwise, Laplace's Demon will suggest and supply R
  code that is ready to be copy/pasted and executed.

  To visualize the MCSE-based stopping rule, run the following code:

  \code{x <- seq(from=-3, to=3, by=0.1);}
  \code{plot(x, dnorm(x,0,1), type="l");}
  \code{abline(v=-0.0627); abline(v=0.0627);}
  \code{abline(v=2*-0.0627, col="red"); abline(v=2*0.0627, col="red")}

  The black vertical lines show the standard error, and the red vertical
  lines show the 95\% interval.

  If the user has an object of class \code{demonoid.hpc}, then the
  \code{Consort} function may be still be applied, but a particular
  chain in the object must be specified as a component in a list. For
  example, with an object called \code{Fit} and a goal of consorting
  over the second chain, the code would be: \code{Consort(Fit[[2]])}.

  The Demonic Suggestion is usually very helpful, but should not be
  followed blindly. Do not let it replace critical thinking. For
  example, \code{Consort} may find that diminishing adaptation is unmet,
  and recommend a different algorithm. However, the user may be
  convinced that the current algorithm is best, and believe instead that
  MCMC found a local solution, and is leaving it to find the global
  solution, in which case adaptations may increase again. Diminishing
  adaptation may have occurred in a previous run, and is not found in
  the current run because adaptation is essentially finished. If either
  of these is true, then it may be best to ignore the newly suggested
  algorithm, and continue with the current algorithm. The suggested code
  may be helpful, but it is merely a suggestion.

  If achieving the appeasement of Laplace's Demon is difficult, consider
  ignoring the MCSE criterion and terminate when all other criteria have
  been met, placing special emphasis on ESS.
}
\references{
  Flegal, J.M., Haran, M., and Jones, G.L. (2008). "Markov chain Monte
  Carlo: Can We Trust the Third Significant Figure?". \emph{Statistical
  Science}, 23, p. 250--260.

  Toft, N., Innocent, G., Gettinby, G., and Reid, S. (2007). "Assessing
  the Convergence of Markov Chain Monte Carlo Methods: An Example from
  Evaluation of Diagnostic Tests in Absence of a Gold Standard".
  \emph{Preventive Veterinary Medicine}, 79, p. 244--256.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{BMK.Diagnostic}},
  \code{\link{ESS}},
  \code{\link{Hangartner.Diagnostic}},
  \code{\link{LaplacesDemon}},
  \code{\link{LaplacesDemon.hpc}},
  \code{\link{LPL.interval}},
  \code{\link{MCSE}}, and
  \code{\link{p.interval}}.
}
\keyword{Diminishing Adaptation}
\keyword{MCMC}
