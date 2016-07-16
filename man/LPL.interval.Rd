\name{LPL.interval}
\alias{LPL.interval}
\title{Lowest Posterior Loss Interval}
\description{
  This function returns the Lowest Posterior Loss (LPL) interval for one
  parameter, given samples from the density of its prior distribution
  and samples of the posterior distribution.
}
\usage{
LPL.interval(Prior, Posterior, prob=0.95, plot=FALSE, PDF=FALSE)
}
\arguments{
  \item{Prior}{This is a vector of samples of the prior density.}
  \item{Posterior}{This is a vector of posterior samples.}
  \item{prob}{This is a numeric scalar in the interval (0,1) giving the
    Lowest Posterior Loss (LPL) interval, and defaults to 0.95,
    representing a 95\% LPL interval.}
  \item{plot}{Logical. When \code{plot=TRUE}, two plots are
    produced. The top plot shows the expected posterior loss. The LPL
    region is shaded black, and the area outside the region is gray. The
    bottom plot shows LPL interval of \eqn{\theta}{theta} on the kernel
    density of \eqn{\theta}{theta}. Again, the LPL region is shaded
    black, and the outside area is gray. A vertical, red, dotted line is
    added at zero for both plots. The \code{plot} argument defaults to
    \code{FALSE}. The plot treats the distribution as if it were
    unimodal; disjoint regions are not estimated here. If multimodality
    should result in disjoint regions, then consider using HPD intervals
    in the \code{\link{p.interval}} function.}
  \item{PDF}{Logical. When \code{PDF=TRUE}, and only when
    \code{plot=TRUE}, plots are saved as a .pdf file in the working
    directory.}
}
\details{
  The Lowest Posterior Loss (LPL) interval (Bernardo, 2005), or LPLI, is
  a probability interval based on intrinsic discrepancy loss between
  prior and posterior distributions. The expected posterior loss
  is the loss associated with using a particular value
  \eqn{\theta_i \in \theta}{theta[i] in theta} of the parameter as the
  unknown true value of \eqn{\theta}{theta} (Bernardo, 2005). Parameter
  values with smaller expected posterior loss should always be
  preferred. The LPL interval includes a region in which all parameter
  values have smaller expected posterior loss than those outside the
  region.
 
  Although any loss function could be used, the loss function should be
  invariant under reparameterization. Any intrinsic loss function is
  invariant under reparameterization, but not necessarily invariant
  under one-to-one transformations of data \eqn{\textbf{x}}{x}. When a
  loss function is also invariant under one-to-one transformations, it
  is usually also invariant when reduced to a sufficient statistic. Only
  an intrinsic loss function that is invariant when reduced to a
  sufficient statistic should be considered.

  The intrinsic discrepancy loss is easily a superior loss function to
  the overused quadratic loss function, and is more appropriate than
  other popular measures, such as Hellinger distance, Kullback-Leibler
  divergence (\code{\link{KLD}}), and Jeffreys logarithmic
  divergence. The intrinsic discrepancy loss is also an
  information-theory related divergence measure. Intrinsic discrepancy
  loss is a symmetric, non-negative loss function, and is a continuous,
  convex function. Intrinsic discrepancy loss was introduced
  by Bernardo and Rueda (2002) in a different context: hypothesis
  testing. Formally, it is:

  \deqn{\delta f(p_2,p_1) = min[\kappa(p_2 | p_1), \kappa(p_1 |
    p_2)]}{delta f(p[2],p[1]) = min[kappa(p[2] | p[1]), kappa(p[1] | p[2])]}

  where \eqn{\delta}{delta} is the discrepancy, \eqn{\kappa}{kappa} is
  the \code{\link{KLD}}, and \eqn{p_1}{p[1]} and \eqn{p_2}{p[2]} are the
  probability distributions. The intrinsic discrepancy loss is the loss
  function, and the expected posterior loss is the mean of the directed
  divergences.
  
  The LPL interval is also called an intrinsic credible interval or
  intrinsic probability interval, and the area inside the interval is
  often called an intrinsic credible region or intrinsic probability
  region.

  In practice, whether a reference prior or weakly informative prior
  (WIP) is used, the LPL interval is usually very close to the HPD
  interval, though the posterior losses may be noticeably different. If
  LPL used a zero-one loss function, then the HPD interval would be
  produced. An advantage of the LPL interval over HPD interval (see
  \code{\link{p.interval}}) is that the LPL interval is invariant to
  reparameterization. This is due to the invariant reparameterization
  property of reference priors. The quantile-based probability interval
  is also invariant to reparameterization. The LPL interval enjoys the
  same advantage as the HPD interval does over the quantile-based
  probability interval: it does not produce equal tails when
  inappropriate.

  Compared with probability intervals, the LPL interval is slightly less
  convenient to calculate. Although the prior distribution is specified
  within the \code{Model} specification function, the user must specify
  it for the \code{LPL.interval} function as well. A comparison of the
  quantile-based probability interval, HPD interval, and LPL interval is
  available here: \url{https://web.archive.org/web/20150214090353/http://www.bayesian-inference.com/credible}.
}
\value{
  A matrix is returned with one row and two columns. The row represents
  the parameter and the column names are \code{"Lower"} and
  \code{"Upper"}. The elements of the matrix are the lower and upper
  bounds of the LPL interval.
}
\references{
  Bernardo, J.M. (2005). "Intrinsic Credible Regions: An Objective
  Bayesian Approach to Interval Estimation". \emph{Sociedad de
  Estadistica e Investigacion Operativa}, 14(2), p. 317--384.

  Bernardo, J.M. and Rueda, R. (2002). "Bayesian Hypothesis Testing: A
  Reference Approach". \emph{International Statistical Review}, 70,
  p. 351--372.
}
\author{Statisticat, LLC.}
\seealso{
  \code{\link{KLD}},
  \code{\link{p.interval}},
  \code{\link{LaplacesDemon}}, and
  \code{\link{PMC}}.}
\examples{
library(LaplacesDemon)
#Although LPL is intended to be applied to output from LaplacesDemon or
#PMC, here is an example in which p(theta) ~ N(0,100), and
#p(theta | y) ~ N(1,10), given 1000 samples.
theta <- rnorm(1000,1,10)
LPL.interval(Prior=dnorm(theta,0,100^2), Posterior=theta, prob=0.95,
  plot=TRUE)
#A more practical example follows, but it assumes a model has been
#updated with LaplacesDemon or PMC, the output object is called Fit, and
#that the prior for the third parameter is normally distributed with
#mean 0 and variance 100: 
#temp <- Fit$Posterior2[,3]
#names(temp) <- colnames(Fit$Posterior2)[3]
#LPL.interval(Prior=dnorm(temp,0,100^2), Posterior=temp, prob=0.95,
#     plot=TRUE, PDF=FALSE)
}
\keyword{Credible Interval}
\keyword{Credible Region}
\keyword{Credible Set}
\keyword{LPL}
\keyword{LPLI}
\keyword{Probability Interval}
\keyword{Probability Region}
\keyword{Probability Set}
