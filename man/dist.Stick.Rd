\name{dist.Stick}
\alias{dStick}
\alias{rStick}
\title{Truncated Stick-Breaking Prior Distribution}
\description{
  These functions provide the density and random number generation of
  the original, truncated stick-breaking (TSB) prior distribution given
  \eqn{\theta}{theta} and \eqn{\gamma}{gamma}, as per Ishwaran and James
  (2001).
}
\usage{
dStick(theta, gamma, log=FALSE)
rStick(M, gamma)
}
\arguments{
  \item{M}{This accepts an integer that is equal to one less than the
    number of truncated number of possible mixture components
    (\eqn{M=1}). Unlike most random deviate functions, this is not the
    number of random deviates to return.}
  \item{theta}{This is \eqn{\theta}{theta}, a vector of length
    \eqn{M-1}, where \eqn{M} is the truncated number of possible mixture
    components.}
  \item{gamma}{This is \eqn{\gamma}{gamma}, a scalar, and is usually
    gamma-distributed.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
  \item Application: Discrete Multivariate
  \item Density: \eqn{p(\pi) = \frac{(1-\theta)^{\beta-1}}{\mathrm{B}(1,\beta)}}{p(pi) = ((1-theta)^(beta-1))/(B(1,\beta))}
  \item Inventor: Sethuraman, J. (1994)
  \item Notation 1: \eqn{\pi \sim
    \mathrm{Stick}(\theta,\gamma)}{pi ~ Stick(theta, gamma)}
  \item Notation 2: \eqn{\pi \sim
    \mathrm{GEM}(\theta,\gamma)}{pi ~ Stick(theta, gamma)}
  \item Notation 3: \eqn{p(\pi) = \mathrm{Stick}(\pi | \theta,
    \gamma)}{p(pi) = Stick(pi | theta, gamma)}
  \item Notation 4: \eqn{p(\pi) = \mathrm{GEM}(\pi | \theta,
    \gamma)}{p(pi) = GEM(pi | theta, gamma)}
  \item Parameter 1: shape parameter \eqn{\theta \in (0,1)}{theta in (0,1)}
  \item Parameter 2: shape parameter \eqn{\gamma > 0}{gamma > 0}
  \item Mean: \eqn{E(\pi) = \frac{1}{1+\gamma}}{E(pi) = 1/(1+gamma)}
  \item Variance: \eqn{var(\pi) = \frac{\gamma}{(1+\gamma)^2
      (\gamma+2)}}{var(pi) = gamma / ((1+gamma)^2 (gamma+2))}
  \item Mode: \eqn{mode(\pi) = 0}{mode(pi) = 0}
}
The original truncated stick-breaking (TSB) prior distribution assigns
each \eqn{\theta}{theta} to be beta-distributed with parameters
\eqn{\alpha=1}{alpha=1} and \eqn{\beta=\gamma}{beta=gamma} (Ishwaran and
James, 2001). This distribution is commonly used in truncated Dirichlet
processes (TDPs).
}
\value{
  \code{dStick} gives the density and
  \code{rStick} generates a random deviate vector of length \eqn{M}.
}
\references{
  Ishwaran, H. and James, L. (2001). "Gibbs Sampling Methods for Stick
  Breaking Priors". \emph{Journal of the American Statistical
  Association}, 96(453), p. 161--173.
  
  Sethuraman, J. (1994). "A Constructive Definition of Dirichlet
  Priors". \emph{Statistica Sinica}, 4, p. 639--650.
}
\seealso{
  \code{\link{ddirichlet}},
  \code{\link{dmvpolya}}, and
  \code{\link{Stick}}.
}
\examples{
library(LaplacesDemon)
dStick(runif(4), 0.1)
rStick(4, 0.1)
}
\keyword{Distribution}
