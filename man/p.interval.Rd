\name{p.interval}
\alias{p.interval}
\title{Probability Interval}
\description{
  This function returns one or more probability intervals of posterior
  samples.
}
\usage{
p.interval(obj, HPD=TRUE, MM=TRUE, prob=0.95, plot=FALSE, PDF=FALSE, \dots)
}
\arguments{
  \item{obj}{This can be either a vector or matrix of posterior samples,
    or an object of class \code{demonoid}, \code{iterquad},
    \code{laplace}, \code{pmc}, or \code{vb}. If it is an object of
    class \code{demonoid}, then it will use only stationary posterior
    samples and monitored target distributions (automatically discarding
    the burn-in; if stationarity does not exist, then it will use all
    samples).}
  \item{HPD}{Logical. This argument defaults to \code{TRUE}, in which
    case one or more Highest Posterior Density (HPD) intervals is
    returned. When \code{FALSE}, one or more quantile-based probability
    intervals is returned.}
  \item{MM}{Logical. This argument defaults to \code{TRUE}, in which
    case each column vector is checked for multimodality, and if found,
    the multimodal form of a Highest Posterior Density (HPD) interval is
    additionally estimated, even when \code{HPD=FALSE}.}
  \item{prob}{This is a numeric scalar in the interval (0,1) giving the
    target probability interval, and defaults to 0.95, representing a
    95\% probability interval. A 95\% probability interval, for example,
    is an interval that contains 95\% of a posterior probability
    distribution.}
  \item{plot}{Logical. When \code{plot=TRUE}, each kernel density is
    plotted and shaded gray, and the area under the curve within the
    probability interval is shaded black. If the kernel density is
    considered to be multimodal, then up to three intervals are shaded
    black. A vertical, red, dotted line is added at zero. The
    \code{plot} argument defaults to \code{FALSE}.}
  \item{PDF}{Logical. When \code{PDF=TRUE}, and only when
    \code{plot=TRUE}, plots are saved as a .pdf file in the working
    directory.}
  \item{\dots}{Additional arguments are unused.}
}
\details{
  A probability interval, also called a credible interval or Bayesian
  confidence interval, is an interval in the domain of a posterior
  probability distribution. When generalized to multivariate forms, it
  is called a probability region (or credible region), though some
  sources refer to a probability region (or credible region) as the
  area within the probability interval. Bivariate probability regions
  may be plotted with the \code{\link{joint.pr.plot}} function.

  The \code{p.interval} function may return different probability
  intervals: a quantile-based probability interval, a unimodal
  Highest Posterior Density (HPD) interval, and multimodal HPD
  intervals. Another type of probability interval is the Lowest
  Posterior Loss (LPL) interval, which is calculated with the
  \code{\link{LPL.interval}} function.

  The quantile-based probability interval is used most commonly,
  possibly because it is simple, the fastest to calculate, invariant
  under transformation, and more closely resembles the frequentist
  confidence interval. The lower and upper bounds of the
  quantile-based probability interval are calculated with the
  \code{quantile} function. A 95\% quantile-based probability interval
  reports the values of the posterior probability distribution that
  indicate the 2.5\% and 97.5\% quantiles, which contain the central
  95\% of the distribution. The quantile-based probability interval is
  centered around the median and has equal-sized tails.
  
  The HPD (highest posterior density) interval is identical to
  the quantile-based probability interval when the posterior probability
  distribution is unimodal and symmetric. Otherwise, the HPD interval
  is the smallest interval, because it is estimated as the interval
  that contains the highest posterior density. Unlike the quantile-based
  probability interval, the HPD interval could be one-tailed or
  two-tailed, whichever is more appropriate. However, unlike the
  quantile-based interval, the HPD interval is not invariant to
  reparameterization (Bernardo, 2005).

  The unimodal HPD interval is estimated from the empirical CDF of the
  sample for each parameter (or deviance or monitored variable) as the
  shortest interval for which the difference in the ECDF values of the
  end-points is the user-specified probability width. This assumes the
  distribution is not severely multimodal.

  As an example, imagine an exponential posterior distribution. A
  quantile-based probability interval would report the highest density
  region near zero to be outside of its interval. In contrast, the
  unimodal HPD interval is recommended for such skewed posterior
  distributions.

  When \code{MM=TRUE}, the \code{\link{is.multimodal}} function is
  applied to each column vector after the unimodal interval (either
  quantile-based or HPD) is estimated. If multimodality is found, then
  multimodal HPD intervals are estimated with kernel density and
  printed to the screen as a character string. The original unimodal
  intervals are returned in the output matrix, because the matrix is
  constrained to have a uniform number of columns per row, and because
  multimodal HPD intervals may be disjoint.

  Disjoint multimodal HPD intervals have multiple intervals for one
  posterior probability distribution. An example may be when there is a
  bimodal, Gaussian distribution with means -10 and 10, variances of 1
  and 1, and a 95\% probability interval is specified. In this case,
  there is not enough density between these two distant modes to have
  only one probability interval.

  The user should also consider \code{\link{LPL.interval}}, since it is
  invariant to reparameterization like the quantile-based probability
  interval, but could be one- or two-tailed, whichever is more
  appropriate, like the HPD interval. A comparison of the quantile-based
  probability interval, HPD interval, and LPL interval is available
  here: \url{https://web.archive.org/web/20150214090353/http://www.bayesian-inference.com/credible}.
}
\value{
  A matrix is returned with rows corresponding to the parameters (or
  deviance or monitored variables), and columns \code{"Lower"} and
  \code{"Upper"}. The elements of the matrix are the unimodal
  probability intervals. The attribute \code{"Probability"} is the
  user-selected probability width. If \code{MM=TRUE} and multimodal
  posterior distributions are found, then multimodal HPD intervals are
  printed to the screen in a character string.
}
\references{
  Bernardo, J.M. (2005). "Intrinsic Credible Regions: An Objective
  Bayesian Approach to Interval Estimation". \emph{Sociedad de
  Estadistica e Investigacion Operativa}, 14(2), p. 317--384.
}
\author{Statisticat, LLC}
\seealso{
  \code{\link{is.multimodal}},
  \code{\link{IterativeQuadrature}},
  \code{\link{joint.pr.plot}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{LPL.interval}},
  \code{\link{PMC}}, and
  \code{\link{VariationalBayes}}.
}
\examples{
##First, update the model with the LaplacesDemon function.
##Then
#p.interval(Fit, HPD=TRUE, MM=TRUE, prob=0.95)
}
\keyword{Credible Interval}
\keyword{Credible Region}
\keyword{Credible Set}
\keyword{HPD}
\keyword{HPDI}
\keyword{Probability Interval}
\keyword{Probability Region}
\keyword{Probability Set}
