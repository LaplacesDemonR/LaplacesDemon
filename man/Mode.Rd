\name{Mode}
\alias{is.amodal}
\alias{is.bimodal}
\alias{is.multimodal}
\alias{is.trimodal}
\alias{is.unimodal}
\alias{Mode}
\alias{Modes}
\title{The Mode(s) of a Vector}
\description{
  The mode is a measure of central tendency. It is the value that occurs
  most frequently, or in a continuous probability distribution, it is
  the value with the most density. A distribution may have no modes
  (such as with a constant, or in a uniform distribution when no value
  occurs more frequently than any other), or one or more modes.
}
\usage{
is.amodal(x, min.size=0.1)
is.bimodal(x, min.size=0.1)
is.multimodal(x, min.size=0.1)
is.trimodal(x, min.size=0.1)
is.unimodal(x, min.size=0.1)
Mode(x)
Modes(x, min.size=0.1)
}
\arguments{
  \item{x}{This is a vector in which a mode (or modes) will be sought.}
  \item{min.size}{This is the minimum size that can be considered a
    mode, where size means the proportion of the distribution between
    areas of increasing kernel density estimates.}
}
\details{
  The \code{is.amodal} function is a logical test of whether or not
  \code{x} has a mode. If \code{x} has a mode, then \code{TRUE} is
  returned, otherwise \code{FALSE}.

  The \code{is.bimodal} function is a logical test of whether or not
  \code{x} has two modes. If \code{x} has two modes, then \code{TRUE}
  is returned, otherwise \code{FALSE}.

  The \code{is.multimodal} function is a logical test of whether or not
  \code{x} has multiple modes. If \code{x} has multiple modes, then
  \code{TRUE} is returned, otherwise \code{FALSE}.

  The \code{is.trimodal} function is a logical test of whether or not
  \code{x} has three modes. If \code{x} has three modes, then \code{TRUE}
  is returned, otherwise \code{FALSE}.

  The \code{is.unimodal} function is a logical test of whether or not
  \code{x} has one mode. If \code{x} has one mode, then \code{TRUE}
  is returned, otherwise \code{FALSE}.
  
  The \code{Mode} function returns the most frequent value when \code{x}
  is discrete. If \code{x} is a constant, then it is considered amodal,
  and \code{NA} is returned. If multiple modes exist, this function
  returns only the mode with the highest density, or if two or more
  modes have the same density, then it returns the first mode found.
  Otherwise, the \code{Mode} function returns the value of \code{x}
  associated with the highest kernel density estimate, or the first
  one found if multiple modes have the same density.

  The \code{Modes} function is a simple, deterministic function that
  differences the kernel density of \code{x} and reports a number of
  modes equal to half the number of changes in direction, although the
  \code{min.size} function can be used to reduce the number of modes
  returned, and defaults to 0.1, eliminating modes that do not have at
  least 10\% of the distributional area. The \code{Modes} function
  returns a list with three components: \code{modes}, \code{modes.dens},
  and \code{size}. The elements in each component are ordered according
  to the decreasing density of the modes. The \code{modes} component is
  a vector of the values of \code{x} associated with the modes. The
  \code{modes.dens} component is a vector of the kernel density
  estimates at the modes. The \code{size} component is a vector of the
  proportion of area underneath each mode.

  The \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}}, and \code{\link{VariationalBayes}}
  functions characterize the marginal posterior distributions by
  posterior modes (means) and variance. A related topic is MAP or
  maximum \emph{a posteriori} estimation.

  Otherwise, the results of Bayesian inference tend to report the
  posterior mean or median, along with probability intervals (see
  \code{\link{p.interval}} and \code{\link{LPL.interval}}), rather than
  posterior modes. In many types of models, such as mixture models, the
  posterior may be multimodal. In such a case, the usual recommendation
  is to choose the highest mode if feasible and possible. However, the
  highest mode may be uncharacteristic of the majority of the posterior.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{LPL.interval}},
  \code{\link{p.interval}}, and
  \code{\link{VariationalBayes}}.
}
\examples{
library(LaplacesDemon)
### Below are distributions with different numbers of modes.
x <- c(1,1) #Amodal
x <- c(1,2,2,2,3) #Unimodal
x <- c(1,2) #Bimodal
x <- c(1,3,3,3,3,4,4,4,4,4) #min.size affects the answer
x <- c(1,1,3,3,3,3,4,4,4,4,4) #Trimodal

### And for each of the above, the functions below may be applied.
Mode(x)
Modes(x)
is.amodal(x)
is.bimodal(x)
is.multimodal(x)
is.trimodal(x)
is.unimodal(x)
}
\keyword{Mode}
\keyword{Utility}