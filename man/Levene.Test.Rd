\name{Levene.Test}
\alias{Levene.Test}
\title{Levene's Test}
\description{
  The \code{Levene.Test} function is a Bayesian form of Levene's test
  (Levene, 1960) of equality of variances.
}
\usage{
Levene.Test(x, Method="U", G=NULL, Data=NULL)
}
\arguments{
  \item{x}{This required argument must be an object of class
    \code{demonoid.ppc}, \code{iterquad.ppc}, \code{laplace.ppc},
    \code{pmc.ppc}, or \code{vb.ppc}.}
  \item{Method}{The method defaults to \code{U} for a univariate
    dependent variable (DV), y. When the DV is multivariate,
    \code{Method="C"} applies Levene's test to each column associated
    in Y. When \code{Method="R"}, Levene's test is applied to each
    row associated in Y.}
  \item{G}{This argument defaults to \code{NULL}, or is required to
    have the same dimensions as the DV. For example, if the DV is
    univariate, then G must have a length equal to y, which is usually
    represented with a length of N for the number of records or T for
    the number of time-periods. If the DV is multivariate, then
    \code{G} must be a matrix, like Y, and have the same number of
    rows and columns. The purpose of the \code{G} argument is to allow
    the user to specify each element of y or Y to be in a particular
    group, so the variance of the groups can be tested. As such, each
    element of \code{G} must consist of an integer from one to the
    number of groups desired to be tested. The reason this test allows
    this degree of specificity is so the user can specify groups, such
    as according to covariate levels. By default, 4 groups are
    specified, where the first quarter of the records are group 1 and
    the last quarter of the records are group 4.}
  \item{Data}{This argument is required when the DV is multivariate,
    hence when \code{Method="C"} or \code{Method="R"}. The DV is
    required to be named Y.}
}
\details{
  This function is a Bayesian form of Levene's test. Levene's test is
  used to assess the probability of the equality of residual variances
  in different groups. When residual variance does not differ by group,
  it is often called homoscedastic (or homoskedastic) residual variance.
  Homoskedastic residual variance is a common assumption. An advantage
  of Levene's test to other tests of homoskedastic residual variance is
  that Levene's test does not require normality of the residuals.

  The \code{Levene.Test} function estimates the test statistic,
  \eqn{W}{W}, as per Levene's test. This Bayesian form, however,
  estimates \eqn{W}{W} from the observed residuals as
  \eqn{W^{obs}}{W.obs}, and \eqn{W}{W} from residuals that are
  replicated from a homoskedastic process as \eqn{W^{rep}}{W.rep}.
  Further, \eqn{W^{obs}}{W.obs} and \eqn{W^{rep}}{W.rep} are
  estimated for each posterior sample. Finally, the probability that
  the distribution of \eqn{W^{obs}}{W.obs} is greater than the
  distribution of \eqn{W^{rep}}{W.rep} is reported (see below).
}
\value{
  The \code{Levene.Test} function returns a plot (or for multivariate Y,
  a series of plots), and a vector with a length equal to the number of
  Levene's tests conducted.

  One plot is produced per univariate application of Levene's test. Each
  plot shows the test statistic W, both from the observed process
  (W.obs as a black density) and the replicated process (W.rep as a red
  line). The mean of W.obs is reported, along with its 95\% quantile-based
  probability interval (see \code{\link{p.interval}}), the probability
  \eqn{p(W^{obs} > W^{rep})}{p(W.obs > W.rep)}, and the indicated
  results, either homoskedastic or heteroskedastic.

  Each element of the returned vector is the probability
  \eqn{p(W^{obs} > W^{rep})}{p(W.obs > W.rep)}. When the probability
  is \eqn{p < 0.025}{p < 0.025} or \eqn{p > 0.975}{p > 0.975},
  heteroskedastic variance is indicated. Otherwise, the variances of
  the groups are assumed not to differ effectively.
}
\references{
  Levene, H. (1960). "Robust Tests for Equality of Variances". In
  I. Olkins, S. G. Ghurye, W. Hoeffding, W. G. Madow, & H. B. Mann
  (Eds.), \emph{Contributions to Probability and Statistics},
  p. 278--292. Stanford University Press: Stanford, CA.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{PMC}},
  \code{\link{p.interval}}, and
  \code{\link{VariationalBayes}}.
}
\examples{
#First, update the model with IterativeQuadrature, LaplaceApproximation,
#  LaplacesDemon, PMC, or VariationalBayes.
#Then, use the predict function, creating, say, object Pred.
#Finally:
#Levene.Test(Pred)
}
\keyword{Homoskedasticity}
\keyword{Utility}
