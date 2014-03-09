\name{Validate}
\alias{Validate}
\title{Holdout Validation}
\description{
  This function performs holdout validation on an object of class
  \code{demonoid} or \code{pmc}, given both a modeled and validation
  data set.
}
\usage{
Validate(object, Model, Data, plot=FALSE, PDF=FALSE)
}
\arguments{
  \item{object}{This is an object of class \code{demonoid} or
    \code{pmc}.}
  \item{Model}{This is a model specification function for
    \code{\link{LaplacesDemon}} or \code{\link{PMC}}.}
  \item{Data}{This is a list that contains two lists of data, as
    specified for \code{\link{LaplacesDemon}}. The first component in
    the list is the list of modeled data, and the second component in
    the list is the list of validation data.}
  \item{plot}{Logical. When \code{plot=TRUE}, two plots are
    displayed. The upper plot shows the density of the modeled deviance
    in black and the density of the validation deviance in red. The
    lower plot shows the density of the change in deviance in gray. The
    \code{plot} argument defaults to \code{FALSE}.}
  \item{PDF}{Logical. When \code{PDF=TRUE} (and \code{plot=TRUE}), the
    plot is saved as a .pdf file. The \code{PDF} argument defaults to
    \code{FALSE}.}
}
\details{
  There are numerous ways to validate a model. In this context,
  validation means to assess the predictive performance of a model on
  out-of-sample data. If reasonable, leave-one-out cross-validation
  (LOOCV) via the conditional predictive ordinate (CPO) should be
  considered when using \code{\link{LaplacesDemon}} or
  \code{\link{PMC}}. For more information on CPO, see the accompanying
  vignettes entitled "Bayesian Inference" and "Examples". CPO is
  unavailable when using \code{\link{LaplaceApproximation}} or
  \code{\link{VariationalBayes}}.

  For \code{\link{LaplaceApproximation}} or
  \code{\link{VariationalBayes}}, it is recommended that the user
  perform holdout validation by comparing posterior predictive checks,
  comparing the differences in the specified discrepancy measure.

  When LOOCV is unreasonable, popular alternatives include k-fold
  cross-validation and holdout validation. Although k-fold
  cross-validation is not performed explicitly here, the user may
  accomplish it with some effort. Of these methods, holdout validation
  includes the most bias, but is the most common in applied use, since
  only one model is fitted, rather than \eqn{k-1} models in k-fold
  cross-validation. The \code{Validate} function performs holdout
  validation.

  For holdout validation, the observed data is sampled randomly into two
  data sets of approximately equal size, or three data sets that
  consists of two data sets of approximately equal size and a remainder
  data set. Of the two data sets approximately equal in size, one is
  called the modeled (or training) data set, and the other is called the
  validation (or test) data set. The modeled data set is used when
  updating the model. After the model is updated, both data sets are
  predicted in the \code{Validate} function, given the model. Predictive
  loss is estimated for the validation data set, relative to the modeled
  data set.

  Predictive loss is associated with overfitting, differences between
  the model and validation data set, or model misspecification. Bayesian
  inference is reputed to be much more robust to overfitting than
  frequentist inference.

  There are many ways to measure predictive loss, and within each
  approach, there are usually numerous possible loss functions. The
  log-likelihood of the model is a popular approximate utility function,
  and consequently, the deviance of the model is a popular loss
  function.

  A vector of model-level (rather than record-level) deviance
  samples is returned with each object of class \code{demonoid} or
  \code{pmc}. The \code{Validate} function obtains this vector for each
  data set, and then calculates the Bayesian Predictive Information
  Criterion (BPIC), as per Ando (2007). BPIC is a variation of the
  Deviance Information Criterion (DIC) that has been modified for
  predictive distributions. For more information on DIC (Spiegelhalter
  et al., 2002), see the accompanying vignette entitled "Bayesian
  Inference". The goal is to minimize BPIC.

  When DIC is applied after the model, such as with a predictive
  distribution, it is positively biased, or too small. The bias is due
  to the same data \eqn{\textbf{y}}{y} being used both to construct the
  posterior distributions and to evaluate pD, the penalty term for model
  complexity. For example, for validation data set
  \eqn{\textbf{y}_{new}}{ynew}, BPIC is:

  \deqn{BPIC = -2\mathrm{log}[p(\textbf{y}_{new}|\textbf{y},\Theta)] +
    2pD}{BPIC = -2log[p(ynew|y,Theta)] + 2pD}

  When \code{plot=TRUE}, the distributions of the modeled and validation
  deviances are plotted above, and the lower plot is the modeled
  deviance subtracted from the validation deviance. When positive, this
  distribution of the change in deviance is the loss in predictive
  deviance associated with moving from the modeled data set to the
  validation data set.

  After using the \code{Validate} function, the user is encouraged to
  perform posterior predictive checks on each data set via the
  \code{\link{summary.demonoid.ppc}} or \code{\link{summary.pmc.ppc}}
  function.
}
\value{
  This function returns a list with three components. The first two
  components are also lists. Each list consists of \code{y},
  \code{yhat}, and \code{Deviance}. The third component is a matrix that
  reports the expected deviance, pD, and BPIC. The object is of class
  \code{demonoid.val} for \code{\link{LaplacesDemon}}, or \code{pmc.val}
  when associated with \code{\link{PMC}}.
}
\references{
  Ando, T. (2007). "Bayesian Predictive Information Criterion for
  the Evaluation of Hierarchical Bayesian and Empirical Bayes Models".
  \emph{Biometrika}, 94(2), p. 443--458.
  
  Spiegelhalter, D.J., Best, N.G., Carlin, B.P., and van der Linde, A.
  (2002). "Bayesian Measures of Model Complexity and Fit (with
  Discussion)". \emph{Journal of the Royal Statistical Society}, B 64,
  p. 583--639.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{PMC}}, and
  \code{\link{VariationalBayes}}.
}
\examples{
library(LaplacesDemon)
#Given an object called Fit of class demonoid, a Model specification,
#and a modeled data set (MyData.M) and validation data set (MyData.V):
#Validate(Fit, Model, Data=list(MyData.M=MyData.M, MyData.V=MyData.V))
}
\keyword{BPIC}
\keyword{Model Selection}
\keyword{Posterior Predictive Checks}
\keyword{Utility}
