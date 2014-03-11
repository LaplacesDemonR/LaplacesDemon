\name{Importance}
\alias{Importance}
\title{Variable Importance}
\description{
  The \code{Importance} function considers variable importance (or
  predictor importance) to be the effect that the variable has on
  replicates \eqn{\textbf{y}^{rep}}{y^rep} (or
  \eqn{\textbf{Y}^{rep}}{Y^rep}) when the variable is removed from the
  model by setting it equal to zero. Here, variable importance is
  considered in terms of the comparison of posterior predictive
  checks. This may be considered to be a form of sensitivity analysis,
  and can be useful for model revision, variable selection, and model
  interpretation.

  Currently, this function only tests the variable importance of design
  matrix \eqn{\textbf{X}}{X}.
}
\usage{
Importance(object, Model, Data, Categorical=FALSE, Discrep, d=0, CPUs=1,
  Type="PSOCK")
}
\arguments{
  \item{object}{An object of class \code{demonoid}, \code{iterquad},
    \code{laplace}, \code{pmc}, or \code{vb} is required.}
  \item{Model}{The model specification function is required.}
  \item{Data}{A data set in a list is required. The dependent variable
    is required to be named either \code{y} or \code{Y}. The
    \code{Importance} function will sequentially remove each column
    vector in \code{X}, so \code{X} is required to be in data set
    \code{Data}.}
  \item{Categorical}{Logical. If \code{TRUE}, then \code{y} and
    \code{yhat} are considered to be categorical (such as y=0 or y=1),
    rather than continuous. This defaults to \code{FALSE}.}
  \item{Discrep}{This optional argument allows a discrepancy statistic to
    be included. For more information on discrepancy statistics, see
    \code{\link{summary.demonoid.ppc}}.}
  \item{d}{This is an optional integer to be used with the
    \code{Discrep} argument above, and it defaults to \code{d=0}. For
    more information on discrepancy, see
    \code{\link{summary.demonoid.ppc}}.}
  \item{CPUs}{This argument accepts an integer that specifies the number
    of central processing units (CPUs) of the multicore computer or
    computer cluster. This argument defaults to \code{CPUs=1}, in which
    parallel processing does not occur.}
  \item{Type}{This argument specifies the type of parallel processing to
    perform, accepting either \code{Type="PSOCK"} or
    \code{Type="MPI"}.}
}
\details{Variable importance is defined here as the impact of each
  variable (predictor, or column vector) in design matrix
  \eqn{\textbf{X}}{X} on \eqn{\textbf{y}^{rep}}{y^rep} (or
  \eqn{\textbf{Y}^{rep}}{Y^rep}), when the variable is removed.

  First, the full model is predicted with the
  \code{\link{predict.demonoid}}, \code{\link{predict.iterquad}},
  \code{\link{predict.laplace}}, \code{\link{predict.pmc}}, or
  \code{\link{predict.vb}} function, and summarized with the
  \code{\link{summary.demonoid.ppc}},
  \code{\link{summary.iterquad.ppc}}, \code{\link{summary.laplace.ppc}},
  \code{\link{summary.pmc.ppc}}, or \code{\link{summary.vb.ppc}}
  function, respectively. The results are stored in the first row of the
  output. Each successive row in the output corresponds to the
  application of \code{predict} and \code{summary} functions, but with
  each variable in design matrix \eqn{\textbf{X}}{X} being set to zero
  and effectively removed. The results show the impact of sequentially
  removing each predictor.

  The criterion for variable importance may differ from model to
  model. As a default, BPIC is recommended. The Bayesian Predictive
  Information Criterion (BPIC) was introduced by Ando (2007). BPIC is a
  variation of the Deviance Information Criterion (DIC) that has been
  modified for predictive distributions. For more information on DIC
  (Spiegelhalter et al., 2002), see the accompanying vignette entitled
  "Bayesian Inference". \eqn{BPIC = Dbar + 2pD}.

  With BPIC, variable importance has a positive relationship, such that
  larger values indicate a more important variable, because removing
  that variable resulted in a worse fit to the data. The best model
  has the lowest BPIC.

  In a model in which the dependent variable is not categorical, it is
  also recommended to consider the L-criterion (Laud and Ibrahim, 1995),
  provided that sample size is small enough that it does not result in
  \code{Inf}. For more information on the L-criterion, see the
  accompanying vignette entitled "Bayesian Inference".

  With the L-criterion, variable importance has a positive relationship, 
  such that larger values indicate a more important variable, because
  removing that variable resulted in a worse fit to the data. Ibrahim
  (1995) recommended considering the model with the lowest
  L-criterion, say as \eqn{L_1}{L[1]}, and the model with the closest
  L-criterion, say as \eqn{L_2}{L[2]}, and creating a comparison score
  as \eqn{\phi = (L_2-L_1)/S_L}{phi = (L[2]-L[1])/S[L]}, where
  \code{S.L} is from the \eqn{L_1}{L[1]} model. If the comparison score,
  \eqn{\phi}{phi} is less than 2, then \eqn{L_2}{L[2]} is within 2
  standard deviations of \eqn{L_1}{L[1]}, and is the recommended cut-off
  for model choice.

  The \code{Importance} function may suggest that a model fits the data
  better with a variable removed. In which case, the user may
  choose to leave the variable in the model (perhaps the model is
  misspecified without the variable), investigate and possibly
  re-specify the relationship between the independent and dependent
  variable(s), or remove the variable and update the model again.

  In contrast to variable importance, the \code{\link{PosteriorChecks}}
  function calculates parameter importance, which is the probability
  that each parameter's marginal posterior distribution is greater than
  zero, where an important parameter does not include zero in its
  probability interval (see \code{\link{p.interval}}). Parameter
  importance and variable importance may disagree, and both should be
  studied.

  The \code{Importance} function tends to indicate that a model fits the
  data better when variables are removed that have parameters with
  marginal posterior distributions that include 0 in the 95\%
  probability interval (variables associated with lower parameter
  importance).

  Often, in complicated models, it is difficult to assess variable
  importance by examining the marginal posterior distribution of the
  associated parameter(s). Consider polynomial regression, in which each
  variable may have multiple parameters.

  The information provided by the \code{Importance} function may be used
  for model revision, or reporting the relative importance of variables.

  The \code{\link{plot.importance}} function is available to plot the
  output of the \code{Importance} function according to BPIC, predictive
  concordance (Gelfand, 1996), the selected discrepancy statistic
  (Gelman et al., 1996), or the L-criterion.

  Parallel processing may be performed when the user specifies
  \code{CPUs} to be greater than one, implying that the specified number
  of CPUs exists and is available. Parallelization may be performed on a
  multicore computer or a computer cluster. Either a Simple Network of
  Workstations (SNOW) or Message Passing Interface is used (MPI). With
  small data sets and few samples, parallel processing may be slower,
  due to computer network communication. With larger data sets and more
  samples, the user should experience a faster run-time.
}
\value{
  \code{Importance} returns an object of class \code{importance}, which
  is a matrix with a number of rows equal to the number of columns in
  design matrix \eqn{\textbf{X}}{X} + 1 (including the full model), and
  4 columns, which are BPIC, Concordance (or Mean.Lift if categorical),
  Discrep, and L-criterion. Each row represents a model with a predictor
  in \eqn{\textbf{X}}{X} removed (except for the first row, which is the
  full model), and the resulting posterior predictive checks. For
  non-categorical dependent variables, an attribute is returned with the
  object, and the attribute is a vector of \code{S.L}, the calibration
  number of the L-criterion.
}
\references{
  Ando, T. (2007). "Bayesian Predictive Information Criterion for
  the Evaluation of Hierarchical Bayesian and Empirical Bayes Models".
  \emph{Biometrika}, 94(2), p. 443--458.

  Gelfand, A. (1996). "Model Determination Using Sampling Based
  Methods". In Gilks, W., Richardson, S., Spiegehalter, D., Chapter 9 in
  Markov Chain Monte Carlo in Practice. Chapman and Hall: Boca Raton, FL.

  Laud, P.W. and Ibrahim, J.G. (1995). "Predictive Model
  Selection". \emph{Journal of the Royal Statistical Society}, B 57,
  p. 247--262.

  Spiegelhalter, D.J., Best, N.G., Carlin, B.P., and van der Linde, A.
  (2002). "Bayesian Measures of Model Complexity and Fit (with
  Discussion)". \emph{Journal of the Royal Statistical Society}, B 64,
  p. 583--639.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{is.importance}},
  \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{PMC}},
  \code{\link{plot.importance}},
  \code{\link{PosteriorChecks}},
  \code{\link{p.interval}},
  \code{\link{predict.demonoid}},
  \code{\link{predict.iterquad}},
  \code{\link{predict.laplace}},
  \code{\link{predict.pmc}},
  \code{\link{predict.vb}},
  \code{\link{summary.demonoid.ppc}},
  \code{\link{summary.iterquad.ppc}},
  \code{\link{summary.laplace.ppc}},
  \code{\link{summary.pmc.ppc}},
  \code{\link{summary.vb.ppc}}, and
  \code{\link{VariationalBayes}}.
}
\examples{
#First, update the model with the LaplacesDemon function, such as
#the example with linear regression, creating an object called Fit.
#Then
#Importance(Fit, Model, MyData, Discrep="Chi-Square", CPUs=1)
}
\keyword{BPIC}
\keyword{High Performance Computing}
\keyword{L-criterion}
\keyword{Model Selection}
\keyword{Posterior Predictive Checks}