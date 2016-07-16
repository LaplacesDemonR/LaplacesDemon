\name{summary.iterquad.ppc}
\alias{summary.iterquad.ppc}
\title{Posterior Predictive Check Summary}
\description{
  This may be used to summarize either new, unobserved instances of
  \eqn{\textbf{y}}{y} (called \eqn{\textbf{y}^{new}}{y[new]}) or
  replicates of \eqn{\textbf{y}}{y} (called
  \eqn{\textbf{y}^{rep}}{y[rep]}). Either \eqn{\textbf{y}^{new}}{y[new]} or
  \eqn{\textbf{y}^{rep}}{y[rep]} is summarized, depending on
  \code{\link{predict.iterquad}}.
}
\usage{\method{summary}{iterquad.ppc}(object, Categorical, Rows,
     Discrep, d, Quiet, \dots)}
\arguments{
  \item{object}{An object of class \code{iterquad.ppc} is required.}
  \item{Categorical}{Logical. If \code{TRUE}, then \code{y} and
    \code{yhat} are considered to be categorical (such as y=0 or y=1),
    rather than continuous.}
  \item{Rows}{An optional vector of row numbers, for example
    \code{c(1:10)}. All rows will be estimated, but only these rows will
    appear in the summary.}
  \item{Discrep}{A character string indicating a discrepancy
    test. \code{Discrep} defaults to \code{NULL}. Valid character
    strings when \code{y} is continuous are: \code{"Chi-Square"},
    \code{"Chi-Square2"}, \code{"Kurtosis"}, \code{"L-criterion"},
    \code{"MASE"}, \code{"MSE"}, \code{"PPL"}, \code{"Quadratic Loss"},
    \code{"Quadratic Utility"}, \code{"RMSE"}, \code{"Skewness"},
    \code{"max(yhat[i,]) > max(y)"}, \code{"mean(yhat[i,]) > mean(y)"},
    \code{"mean(yhat[i,] > d)"}, \code{"mean(yhat[i,] > mean(y))"},
    \code{"min(yhat[i,]) < min(y)"}, \code{"round(yhat[i,]) = d"}, and
    \code{"sd(yhat[i,]) > sd(y)"}. Valid character strings when \code{y}
    is categorical are: \code{"p(yhat[i,] != y[i])"}. Kurtosis and
    skewness are not discrepancies, but are included here for convenience.}
  \item{d}{This is an optional integer to be used with the
    \code{Discrep} argument above, and it defaults to \code{d=0}.}
  \item{Quiet}{This logical argument defaults to \code{FALSE} and will
    print results to the console. When \code{TRUE}, results are not
    printed.}
  \item{\dots}{Additional arguments are unused.}
}
\details{
  This function summarizes an object of class \code{iterquad.ppc}, which
  consists of posterior predictive checks on either
  \eqn{\textbf{y}^{new}}{y[new]} or \eqn{\textbf{y}^{rep}}{y[rep]},
  depending respectively on whether unobserved instances of
  \eqn{\textbf{y}}{y} or the model sample of \eqn{\textbf{y}}{y} was
  used in the \code{\link{predict.iterquad}} function. The deviance and
  monitored variables are also summarized.
  
  The purpose of a posterior predictive check is to assess how well (or
  poorly) the model fits the data, or to assess discrepancies between
  the model and the data. For more information on posterior predictive
  checks, see
  \url{https://web.archive.org/web/20150215050702/http://www.bayesian-inference.com/posteriorpredictivechecks}.

  When \eqn{\textbf{y}}{y} is continuous and known, this function
  estimates the predictive concordance between \eqn{\textbf{y}}{y} and
  \eqn{\textbf{y}^{rep}}{y[rep]} as per Gelfand (1996), and the
  predictive quantile (PQ), which is for record-level outlier detection
  used to calculate Gelfand's predictive concordance.

  When \eqn{\textbf{y}}{y} is categorical and known, this function
  estimates the record-level lift, which is
  \code{p(yhat[i,] = y[i]) / [p(y = j) / n]}, or
  the number of correctly predicted samples over the rate of that
  category of \eqn{\textbf{y}}{y} in vector \eqn{\textbf{y}}{y}.

  A discrepancy measure is an approach to studying discrepancies between
  the model and data (Gelman et al., 1996). Below is a list of
  discrepancy measures, followed by a brief introduction to discrepancy
  analysis:

  \itemize{
    \item The \code{"Chi-Square"} discrepancy measure is the chi-square
  goodness-of-fit test that is recommended by Gelman. For each record
  i=1:N, this returns (y[i] - E(y[i]))^2 / var(yhat[i,]).
    \item The \code{"Chi-Square2"} discrepancy measure returns the
  following for each record: Pr(chisq.rep[i,] > chisq.obs[i,]), where
  chisq.obs[i,] <- (y[i] - E(y[i]))^2 / E(y[i]), and chisq.rep[i,] <-
  (yhat[i,] - E(yhat[i,]))^2 / E(yhat[i,]), and the overall
  discrepancy is the percent of records that were outside of the 95\%
  quantile-based probability interval (see \code{\link{p.interval}}).
    \item The \code{"Kurtosis"} discrepancy measure returns the kurtosis
  of \eqn{\textbf{y}^{rep}}{y[rep]} for each record, and the discrepancy
  statistic is the mean for all records. This does not measure
  discrepancies between the model and data, and is useful for finding
  kurtotic replicate distributions.
    \item The \code{"L-criterion"} discrepancy measure of Laud and Ibrahim
  (1995) provides the record-level combination of two components (see
  below), and the discrepancy statistic is the sum, \code{L}, as well as
  a calibration number, \code{S.L}. For more information on the
  L-criterion, see the accompanying vignette entitled "Bayesian
  Inference".
    \item The \code{"MASE"} (Mean Absolute Scaled Error) is a
  discrepancy measure for the accuracy of time-series forecasts,
  estimated as \code{(|y - yhat|) / mean(abs(diff(y)))}. The discrepancy
  statistic is the mean of the record-level values.
    \item The \code{"MSE"} (Mean Squared Error) discrepancy measure
  provides the MSE for each record across all replicates, and the
  discrepancy statistic is the mean of the record-level MSEs. MSE and
  quadratic loss are identical.
    \item The \code{"PPL"} (Posterior Predictive Loss) discrepancy
  measure of Gelfand and Ghosh (1998) provides the record-level
  combination of two components: one involves the predictive variance
  and the other includes the accuracy of the means of the predictive
  distribution. The \code{d=0} argument applies the following weight to
  the accuracy component, which is then added to the variance component:
  \eqn{d/(d+1)}. For \eqn{\textbf{y}^{new}}{y[new]}, use \eqn{d=0}. For
  \eqn{\textbf{y}^{rep}}{y[rep]} and model comparison, \eqn{d} is
  commonly set to 1, 10, or 100000. Larger values of \eqn{d} put more
  stress on fit and downgrade the precision of the estimates.
    \item The \code{"Quadratic Loss"} discrepancy measure provides the
  mean quadratic loss for each record across all replicates, and the
  discrepancy statistic is the mean of the record-level mean quadratic
  losses. Quadratic loss and MSE are identical, and quadratic loss is
  the negative of quadratic utility.
    \item The \code{"Quadratic Utility"} discrepancy measure provides
  the mean quadratic utility for each record across all replicates, and
  the discrepancy statistic is the mean of the record-level mean
  quadratic utilities. Quadratic utility is the negative of quadratic
  loss.
    \item The \code{"RMSE"} (Root Mean Squared Error) discrepancy
  measure provides the RMSE for each record across all replicates, and
  the discrepancy statistic is the mean of the record-level RMSEs.
    \item The \code{"Skewness"} discrepancy measure returns the skewness
  of \eqn{\textbf{y}^{rep}}{y[rep]} for each record, and the discrepancy
  statistic is the mean for all records. This does not measure
  discrepancies between the model and data, and is useful for finding
  skewed replicate distributions.
    \item The \code{"max(yhat[i,]) > max(y)"} discrepancy measure
  returns a record-level indicator when a record's maximum
  \eqn{\textbf{y}^{rep}_i}{y[rep][i]} exceeds the maximum of
  \eqn{\textbf{y}}{y}. The discrepancy statistic is the mean of the
  record-level indicators, reporting the proportion of records with
  replications that exceed the maximum of \eqn{\textbf{y}}{y}.
    \item The \code{"mean(yhat[i,]) > mean(y)"} discrepancy measure
    returns a record-level indicator when the mean of a record's
  \eqn{\textbf{y}^{rep}_i}{y[rep][i]} is greater than the mean of
  \eqn{\textbf{y}}{y}. The discrepancy statistic is the mean of the
  record-level indicators, reporting the proportion of records with
  mean replications that exceed the mean of \eqn{\textbf{y}}{y}.
    \item The \code{"mean(yhat[i,] > d)"} discrepancy measure returns a
  record-level proportion of \eqn{\textbf{y}^{rep}_i}{y[rep][i]} that
  exceeds a specified value, \code{d}. The discrepancy statistic is the
  mean of the record-level proportions.
    \item The \code{"mean(yhat[i,] > mean(y))"} discrepancy measure
  returns a record-level proportion of
  \eqn{\textbf{y}^{rep}_i}{y[rep][i]} that exceeds the mean of
  \eqn{\textbf{y}}{y}. The discrepancy statistic is the mean of the
  record-level proportions.
    \item The \code{"min(yhat[i,]) < min(y)"} discrepancy measure
  returns a record-level indicator when a record's minimum
  \eqn{\textbf{y}^{rep}_i}{y[rep][i]} is less than the minimum of
  \eqn{\textbf{y}}{y}. The discrepancy statistic is the mean of the
  record-level indicators, reporting the proportion of records with
  replications less than the minimum of \eqn{\textbf{y}}{y}.
    \item The \code{"round(yhat[i,]) = d"} discrepancy measure returns a
  record-level proportion of \eqn{\textbf{y}^{rep}_i}{y[rep][i]} that,
  when rounded, is equal to a specified discrete value, \code{d}. The
  discrepancy statistic is the mean of the record-level proportions.
    \item The \code{"sd(yhat[i,]) > sd(y)"} discrepancy measure returns a
  record-level indicator when the standard deviation of replicates is
  larger than the standard deviation of all of \eqn{\textbf{y}}{y}. The
  discrepancy statistic is the mean of the record-level indicators,
  reporting the proportion of records with larger standard deviations
  than \eqn{\textbf{y}}{y}.
    \item The \code{"p(yhat[i,] != y[i])"} discrepancy measure returns
  the record-level probability that \eqn{\textbf{y}^{rep}_i}{y[rep][i]}
  is not equal to \eqn{\textbf{y}}{y}. This is valid when
  \eqn{\textbf{y}}{y} is categorical and \code{yhat} is the predicted
  category. The probability is the proportion of replicates.
}
  After observing a discrepancy statistic, the user attempts to improve
  the model by revising the model to account for discrepancies between
  data and the current model. This approach to model revision relies on
  an analysis of the discrepancy statistic. Given a discrepancy measure
  that is based on model fit, such as the L-criterion, the user may
  correlate the record-level discrepancy statistics with the dependent
  variable, independent variables, and interactions of independent
  variables. The discrepancy statistic should not correlate with the
  dependent and independent variables. Interaction variables may be
  useful for exploring new relationships that are not in the current
  model. Alternatively, a decision tree may be applied to the
  record-level discrepancy statistics, given the independent variables,
  in an effort to find relationships in the data that may be helpful
  in the model. Model revision may involve the addition of a finite
  mixture component to account for outliers in discrepancy, or
  specifying the model with a distribution that is more robust to
  outliers. There are too many suggestions to include here, and
  discrepancy analysis varies by model.
}
\value{
  This function returns a list with the following components:
  \item{BPIC}{The Bayesian Predictive Information Criterion (BPIC) was
    introduced by Ando (2007). BPIC is a variation of the Deviance
    Information Criterion (DIC) that has been modified for predictive
    distributions. For more information on DIC (Spiegelhalter
    et al., 2002), see the accompanying vignette entitled "Bayesian
    Inference". \eqn{BPIC = Dbar + 2pD}. The goal is to minimize BPIC.}
  \item{Concordance}{This is the percentage of the records of y that are
    within the 95\% quantile-based probability interval (see
    \code{\link{p.interval}}) of \eqn{\textbf{y}^{rep}}{y[rep]}.
    Gelfand's suggested goal is to achieve 95\% predictive concordance.
    Lower percentages indicate too many outliers and a poor fit of the
    model to the data, and higher percentages may suggest overfitting.
    Concordance occurs only when \eqn{\textbf{y}}{y} is continuous.}
  \item{Mean Lift}{This is the mean of the record-level lifts, and
    occurs only when \eqn{\textbf{y}}{y} is specified as categorical
    with \code{Categorical=TRUE}.}
  \item{Discrepancy.Statistic}{This is only reported if the
    \code{Discrep} argument receives a valid discrepancy measure as
    listed above. The \code{Discrep} applies to each record of
    \eqn{\textbf{y}}{y}, and the \code{Discrepancy.Statistic} reports
    the results of the discrepancy measure on the entire data set. For
    example, if \code{Discrep="min(yhat[i,]) < min(y)"}, then the
    overall result is the proportion of records in which the minimum
    sample of yhat was less than the overall minimum
    \eqn{\textbf{y}}{y}. This is \code{Pr(min(yhat[i,]) < min(y) | y,
      Theta)}, where \code{Theta} is the parameter set.}
  \item{L-criterion}{The L-criterion (Laud and Ibrahim, 1995) was
    developed for model and variable selection. It is a sum of two
    components: one involves the predictive variance and the other
    includes the accuracy of the means of the predictive
    distribution. The L-criterion measures model performance with a
    combination of how close its predictions are to the observed data
    and variability of the predictions. Better models have smaller
    values of \code{L}. \code{L} is measured in the same units as
    the response variable, and measures how close the data vector
    \eqn{\textbf{y}}{y} is to the predictive distribution. In addition
    to the value of \code{L}, there is a value for \code{S.L}, which is
    the calibration number of \code{L}, and is useful in determining how
    much of a decrease is necessary between models to be noteworthy.}
  \item{Monitor}{This is a \eqn{N \times 5}{N x 5} matrix, where \eqn{N}
    is the number of monitored variables and there are 5 columns, as
    follows: Mean, SD, LB (the 2.5\% quantile), Median, and UB (the
    97.5\% quantile).}
  \item{Summary}{When \eqn{\textbf{y}}{y} is continuous, this is a
    \eqn{N \times 8}{N x 8} matrix, where \eqn{N} is the number of
    records of \eqn{\textbf{y}}{y} and there are 8 columns, as follows:
    y, Mean, SD, LB (the 2.5\% quantile), Median, UB (the 97.5\%
    quantile), PQ (the predictive quantile, which is
    \eqn{Pr(\textbf{y}^{rep} \ge \textbf{y})}{Pr(y[rep] >= y)}), and
    Test, which shows the record-level result of a test, if
    specified. When \eqn{\textbf{y}}{y} is categorical, this matrix has
    a number of columns equal to the number of categories of
    \eqn{\textbf{y}}{y} plus 3, also including \code{y}, \code{Lift},
    and \code{Discrep}.}
}
\references{
  Ando, T. (2007). "Bayesian Predictive Information Criterion for
  the Evaluation of Hierarchical Bayesian and Empirical Bayes Models".
  \emph{Biometrika}, 94(2), p. 443--458.
  
  Gelfand, A. (1996). "Model Determination Using Sampling Based
  Methods". In Gilks, W., Richardson, S., Spiegehalter, D., Chapter 9 in
  Markov Chain Monte Carlo in Practice. Chapman and Hall: Boca Raton, FL.

  Gelfand, A. and Ghosh, S. (1998). "Model Choice: A Minimum Posterior
  Predictive Loss Approach". \emph{Biometrika}, 85, p. 1--11.

  Gelman, A., Meng, X.L., and Stern H. (1996). "Posterior Predictive
  Assessment of Model Fitness via Realized Discrepancies".
  \emph{Statistica Sinica}, 6, p. 733--807.

  Laud, P.W. and Ibrahim, J.G. (1995). "Predictive Model Selection".
  \emph{Journal of the Royal Statistical Society}, B 57, p. 247--262.

  Spiegelhalter, D.J., Best, N.G., Carlin, B.P., and van der Linde, A.
  (2002). "Bayesian Measures of Model Complexity and Fit (with
  Discussion)". \emph{Journal of the Royal Statistical Society}, B 64,
  p. 583--639.
}
\author{Statisticat, LLC.}
\seealso{
  \code{\link{IterativeQuadrature}},
  \code{\link{predict.iterquad}}, and
  \code{\link{p.interval}}.
}
\examples{### See the IterativeQuadrature function for an example.}
\keyword{BPIC}
\keyword{Chi-Square}
\keyword{Kurtosis}
\keyword{L-criterion}
\keyword{MASE}
\keyword{MSE}
\keyword{Posterior Predictive Checks}
\keyword{Posterior Predictive Loss}
\keyword{Quadratic Loss}
\keyword{Quadratic Utility}
\keyword{RMSE}
\keyword{Skewness}
\keyword{summary}
