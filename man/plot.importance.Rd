\name{plot.importance}
\alias{plot.importance}
\title{Plot Variable Importance}
\description{
  This may be used to plot variable importance with BPIC, predictive
  concordance, a discrepancy statistic, or the L-criterion regarding an
  object of class \code{importance}.
}
\usage{\method{plot}{importance}(x, Style="BPIC", \dots)}
\arguments{
  \item{x}{This required argument is an object of class
    \code{importance}.}
  \item{Style}{When \code{Style="BPIC"}, BPIC is shown, and \code{BPIC}
    is the default. Otherwise, predictive concordance is plotted when
    \code{Style="Concordance"}, a discrepancy statistic is plotted when
    \code{Style="Discrep"}, or the L-criterion is plotted when
    \code{Style="L-criterion"}.}
  \item{\dots}{Additional arguments are unused.}
}
\details{
  The x-axis is either BPIC (Ando, 2007), predictive concordance
  (Gelfand, 1996), a discrepancy statistic (Gelman et al., 1996), or the
  L-criterion (Laud and Ibrahim, 1995) of the \code{\link{Importance}}
  function (depending on the \code{Style} argument), and variables are
  on the y-axis. A more important variable is associated with a dot that
  is plotted farther to the right. For more information on variable
  importance, see the \code{\link{Importance}} function.
}
\references{
  Ando, T. (2007). "Bayesian Predictive Information Criterion for
  the Evaluation of Hierarchical Bayesian and Empirical Bayes Models".
  \emph{Biometrika}, 94(2), p. 443--458.

  Gelfand, A. (1996). "Model Determination Using Sampling Based
  Methods". In Gilks, W., Richardson, S., Spiegehalter, D., Chapter 9 in
  Markov Chain Monte Carlo in Practice. Chapman and Hall: Boca Raton, FL.

  Gelman, A., Meng, X.L., and Stern H. (1996). "Posterior Predictive
  Assessment of Model Fitness via Realized Discrepancies".
  \emph{Statistica Sinica}, 6, p. 733--807.
  
  Laud, P.W. and Ibrahim, J.G. (1995). "Predictive Model
  Selection". \emph{Journal of the Royal Statistical Society}, B 57,
  p. 247--262.
}
\author{Statisticat, LLC \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{Importance}}}
\keyword{Plot}