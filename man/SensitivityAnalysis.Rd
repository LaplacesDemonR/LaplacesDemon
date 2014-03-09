\name{SensitivityAnalysis}
\alias{SensitivityAnalysis}
\title{Sensitivity Analysis}
\description{
  This function performs an elementary sensitivity analysis for two
  models regarding marginal posterior distributions and posterior
  inferences.
}
\usage{
SensitivityAnalysis(Fit1, Fit2, Pred1, Pred2)
}
\arguments{
  \item{Fit1}{This argument accepts an object of class \code{demonoid},
    \code{iterquad}, \code{laplace}, \code{pmc}, or \code{vb}.}
  \item{Fit2}{This argument accepts an object of class \code{demonoid},
    \code{iterquad}, \code{laplace}, \code{pmc}, or \code{vb}.}
  \item{Pred1}{This argument accepts an object of class
    \code{demonoid.ppc}, \code{iterquad.ppc}, \code{laplace.ppc},
    \code{pmc.ppc}, or \code{vb.ppc}.}
  \item{Pred2}{This argument accepts an object of class
    \code{demonoid.ppc}, \code{iterquad.ppc}, \code{laplace.ppc},
    \code{pmc.ppc}, or \code{vb.ppc}.}
}
\details{
  Sensitivity analysis is concerned with the influence from changes to the
  inputs of a model on the output. Comparing differences resulting from
  different prior distributions is the most common application of
  sensitivity analysis, though results from different likelihoods may be
  compared as well. The outputs of interest are the marginal posterior
  distributions and posterior inferences.

  There are many more methods of conducting a sensitivity analysis than
  exist in the \code{SensitivityAnalysis} function. For more
  information, see Oakley and O'Hagan (2004). The \code{\link{SIR}}
  function is useful for approximating changes in the posterior due to
  small changes in prior distributions.

  The \code{SensitivityAnalysis} function compares marginal posterior
  distributions and posterior predictive distributions. Specifically,
  it calculates the probability that each distribution in \code{Fit1}
  and \code{Pred1} is greater than the associated distribution in
  \code{Fit2} and \code{Pred2}, and returns a variance ratio of each
  pair of distributions. If the probability is \eqn{0.5} that a
  distribution is greater than another, or if the variance ratio is
  \eqn{1}, then no difference is found due to the inputs.

  Additional comparisons and methods are currently outside the scope of
  the \code{SensitivityAnalysis} function. The \code{\link{BayesFactor}}
  function may also be considered, as well as comparing posterior
  predictive checks resulting from \code{\link{summary.demonoid.ppc}},
  \code{\link{summary.iterquad.ppc}},
  \code{\link{summary.laplace.ppc}}, \code{\link{summary.pmc.ppc}}, or
  \code{\link{summary.vb.ppc}}.

  Regarding marginal posterior distributions, the
  \code{SensitivityAnalysis} function compares only distributions with
  identical parameter names. For example, suppose a statistician
  conducts a sensitivity analysis to study differences resulting from
  two prior distributions: a normal distribution and a Student t
  distribution. These distributions have two and three parameters,
  respectively. The statistician has named the parameters \code{beta}
  and \code{sigma} for the normal distribution, while for the Student
  t distribution, the parameters are named \code{beta}, \code{sigma},
  and \code{nu}. In this case, the \code{SensitivityAnalysis} function
  compares the marginal posterior distributions for \code{beta} and
  \code{sigma}, though \code{nu} is ignored because it is not in both
  models. If the statistician does not want certain parameters compared,
  then differing parameter names should be assigned.

  Robust Bayesian analysis is a very similar topic, and often called
  simply Bayesian sensitivity analysis. In robust Bayesian analysis, the
  robustness of answers from a Bayesian analysis to uncertainty about
  the precise details of the analysis is studied. An answer is
  considered robust if it does not depend sensitively on the assumptions
  and inputs on which it is based. Robust Bayes methods acknowledge that
  it is sometimes very difficult to come up with precise distributions
  to be used as priors. Likewise the appropriate likelihood function
  that should be used for a particular problem may also be in doubt. In
  a robust Bayesian analysis, a standard Bayesian analysis is applied to
  all possible combinations of prior distributions and likelihood
  functions selected from classes of priors and likelihoods considered
  empirically plausible by the statistician.
}
\value{
  This function returns a list with the following components:
  \item{Posterior}{This is a \eqn{J \times 2}{J x 2} matrix of \eqn{J}
    marginal posterior distributions. Column names are "p(Fit1 > Fit2)"
    and "var(Fit1) / var(Fit2)".}
  \item{Post.Pred.Dist}{This is a \eqn{N \times 2}{N x 2} matrix of \eqn{N}
    posterior predictive distributions. Column names are
    "p(Pred1 > Pred2)" and "var(Pred1) / var(Pred2)".}
}
\references{
  Berger, J.O. (1984). "The Robust Bayesian Viewpoint (with
  discussion)". In J. B. Kadane, editor, Robustness of Bayesian
  Analyses, p. 63--144. North-Holland, Amsterdam.

  Berger, J.O. (1985). "Statistical Decision Theory and Bayesian
  Analysis". Springer-Verlag, New York.

  Berger, J.O. (1994). "An Overview of Robust Bayesian Analysis
  (with discussion)". Test, 3, p. 5--124.
  
  Oakley, J. and O'Hagan, A. (2004). "Probabilistic Sensitivity Analysis
  of Complex Models: a Bayesian Approach". \emph{Journal of the Royal
  Statistical Society, Series B}, 66, p. 751--769.

  Weiss, R. (1995). "An Approach to Bayesian Sensitivity Analysis".
  \emph{Journal of the Royal Statistical Society, Series B}, 58,
  p. 739--750.
}
\author{Statisticat, LLC \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{BayesFactor}},
  \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{PMC}},
  \code{\link{predict.demonoid}},
  \code{\link{predict.iterquad}},
  \code{\link{predict.laplace}},
  \code{\link{predict.pmc}},
  \code{\link{SIR}},
  \code{\link{summary.demonoid.ppc}},
  \code{\link{summary.iterquad.ppc}},
  \code{\link{summary.laplace.ppc}},
  \code{\link{summary.pmc.ppc}}, and
  \code{\link{VariationalBayes}}.
}
\examples{
#sa <- SensitivityAnalysis(Fit1, Fit2, Pred1, Pred2)
#sa
}
\keyword{Sensitivity}
