\name{CenterScale}
\alias{CenterScale}
\title{Centering and Scaling}
\description{
  This function either centers and scales a continuous variable and
  provides options for binary variables, or returns an untransformed
  variable from a centered and scaled variable.
}
\usage{
CenterScale(x, Binary="none", Inverse=FALSE, mu, sigma, Range, Min)
}
\arguments{
  \item{x}{This is a vector to be centered and scaled, or to be
    untransformed if \code{Inverse=TRUE}.}
  \item{Binary}{This argument indicates how binary variables will be
    treated, and defaults to \code{"none"}, which keeps the original
    scale, or transforms the variable to the 0-1 range, if not already
    there. With \code{"center"}, it will center the binary variable by
    subtracting the mean. With \code{"center0"}, it centers the binary
    variable at zero, recoding a 0 to -0.5, and a 1 to 0.5. Finally,
    \code{"centerscale"} will center and scale the binary variable,
    subtracting the mean and dividing by two standard deviations.}
  \item{Inverse}{Logical. If \code{TRUE}, then a centered and scaled
    variable \code{x} will be transformed to its original, un-centered
    and un-scaled state. This defaults to \code{FALSE}.}
  \item{mu, sigma, Range, Min}{These arguments are required only when
    \code{Inverse=TRUE}, where \code{mu} is the mean, \code{sigma} is
    the standard deviation, \code{Range} is the range, and \code{Min} is
    the minimum of the original \code{x}. \code{Range} and \code{Min}
    are used only when \code{Binary="none"} or \code{Binary="center0"}.}
}
\details{
  Gelman (2008) recommends centering and scaling continuous predictors
  to facilitate MCMC convergence and enable comparisons between
  coefficients of centered and scaled continuous predictors with
  coefficients of untransformed binary predictors. A continuous
  predictor is centered and scaled as follows: \code{x.cs <- (x -
    mean(x)) / (2*sd(x))}. This is an improvement over the usual
  practice of standardizing predictors, which is \code{x.z <- (x -
    mean(x)) / sd(x)}, where coefficients cannot be validly compared
  between binary and continuous predictors.

  In MCMC, such as in \code{\link{LaplacesDemon}}, a centered and
  scaled predictor often results in a higher effective sample size
  (\code{\link{ESS}}), and therefore the chain mixes better. Centering
  and scaling is a method of re-parameterization to improve mixing.

  Griffin and Brown (2013) also assert that the user may not want to
  scale predictors that are measured on the same scale, since scaling
  in this case may increase noisy, low signals. In this case, centering
  (without scaling) is recommended. To center a predictor, subtract its
  mean.

}
\value{
  The \code{CenterScale} function returns a centered and scaled vector,
  or the untransformed vector.
}
\references{
  Gelman, A. (2008). "Scaling Regression Inputs by Dividing by Two Standard
  Devations". \emph{Statistics in Medicine}, 27, p. 2865--2873.

  Griffin, J.E. and Brown, P.J. (2013) "Some Priors for Sparse
  Regression Modelling". \emph{Bayesian Analysis}, 8(3), p. 691--702.
}
\seealso{
  \code{\link{ESS}},
  \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}}, and
  \code{\link{PMC}}.
}
\examples{
### See the LaplacesDemon function for an example in use.
library(LaplacesDemon)
x <- rnorm(100,10,1)
x.cs <- CenterScale(x)
x.orig <- CenterScale(x.cs, Inverse=TRUE, mu=mean(x), sigma=sd(x))
}
\keyword{Center}
\keyword{Centering}
\keyword{Predictor}
\keyword{Scale}
\keyword{Scaling}
\keyword{Transformation}