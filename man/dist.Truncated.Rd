\name{dist.Truncated}
\alias{dtrunc}
\alias{extrunc}
\alias{ptrunc}
\alias{qtrunc}
\alias{rtrunc}
\alias{vartrunc}
\title{Truncated Distributions}
\description{
  Density, distribution function, quantile function and random
  generation for truncated distributions.
}
\usage{
dtrunc(x, spec, a=-Inf, b=Inf, log=FALSE, ...)
extrunc(spec, a=-Inf, b=Inf, ...)
ptrunc(x, spec, a=-Inf, b=Inf, ...)
qtrunc(p, spec, a=-Inf, b=Inf, ...)
rtrunc(n, spec, a=-Inf, b=Inf, ...)
vartrunc(spec, a=-Inf, b=Inf, ...)
}
\arguments{
  \item{n}{This is a the number of random draws for \code{rtrunc}.}
  \item{p}{This is a vector of probabilities.}
  \item{x}{This is a vector to be evaluated.}
  \item{spec}{The base name of a probability distribution is
    specified here. For example, to estimate the density of a
    truncated normal distribution, enter \code{norm}.}
  \item{a}{This is the lower bound of truncation, which defaults
    to negative infinity.}
  \item{b}{This is the upper bound of truncation, which defaults
    to infinity.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
  \item{\dots}{Additional arguments pertain to the probability
    distribution specified in the \code{spec} argument.}
}
\details{
  A truncated distribution is a conditional distribution that results
  from a priori restricting the domain of some other probability
  distribution. More than merely preventing values outside of truncated
  bounds, a proper truncated distribution integrates to one within the
  truncated bounds. For more information on propriety, see
  \code{\link{is.proper}}. In contrast to a truncated distribution, a
  censored distribution occurs when the probability distribution is
  still allowed outside of a pre-specified range. Here, distributions
  are truncated to the interval \eqn{[a,b]}, such as \eqn{p(\theta) \in
  [a,b]}{p(theta) in [a,b]}.

  The \code{dtrunc} function is often used in conjunction with the
  \code{\link{interval}} function to truncate prior probability
  distributions in the model specification function for use with these
  numerical approximation functions: \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}}, and \code{\link{PMC}}.

  The R code of Nadarajah and Kotz (2006) has been modified to work with
  log-densities.
}
\value{
  \code{dtrunc} gives the density,
  \code{extrunc} gives the expectation,
  \code{ptrunc} gives the distribution function,
  \code{qtrunc} gives the quantile function,
  \code{rtrunc} generates random deviates, and
  \code{vartrunc} gives the variance of the truncated distribution.
}
\references{
  Nadarajah, S. and Kotz, S. (2006). "R Programs for Computing Truncated
  Distributions". \emph{Journal of Statistical Software}, 16,
  Code Snippet 2, p. 1--8.
}
\seealso{
  \code{\link{interval}},
  \code{\link{is.proper}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}}, and
  \code{\link{PMC}}.
}
\examples{
library(LaplacesDemon)
x <- seq(-0.5, 0.5, by = 0.1)
y <- dtrunc(x, "norm", a=-0.5, b=0.5, mean=0, sd=2)
}
\keyword{Distribution}
