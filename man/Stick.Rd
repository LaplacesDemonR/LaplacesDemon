\name{Stick}
\alias{Stick}
\title{Truncated Stick-Breaking}
\description{
  The \code{Stick} function provides the utility of truncated
  stick-breaking regarding the vector
  \eqn{\theta}{theta}. Stick-breaking is commonly referred to as a
  stick-breaking process, and is used often in a Dirichlet
  process (Sethuraman, 1994). It is commonly associated with
  infinite-dimensional mixtures, but in practice, the `infinite' number
  is truncated to a finite number, since it is impossible to estimate an
  infinite number of parameters (Ishwaran and James, 2001).
}
\usage{
Stick(theta)
}
\arguments{
  \item{theta}{This required argument, \eqn{\theta}{theta} is a vector
    of length \eqn{(M-1)}{M-1} regarding \eqn{M}{M} mixture components.}
}
\details{
  The Dirichlet process (DP) is a stochastic process used in Bayesian
  nonparametric modeling, most commonly in DP mixture models, otherwise
  known as infinite mixture models. A DP is a distribution over
  distributions. Each draw from a DP is itself a discrete
  distribution. A DP is an infinite-dimensional generalization of
  Dirichlet distributions. It is called a DP because it has
  Dirichlet-distributed, finite-dimensional, marginal distributions,
  just as the Gaussian process has Gaussian-distributed,
  finite-dimensional, marginal distributions. Distributions drawn from a
  DP cannot be described using a finite number of parameters, thus the
  classification as a nonparametric model. The truncated stick-breaking
  (TSB) process is associated with a truncated Dirichlet process (TDP).

  An example of a TSB process is cluster analysis, where the number of
  clusters is unknown and treated as mixture components. In such a
  model, the TSB process calculates probability vector \eqn{\pi}{pi}
  from \eqn{\theta}{theta}, given a user-specified maximum number of
  clusters to explore as \eqn{C}{C}, where \eqn{C}{C} is the length of
  \eqn{\theta + 1}{theta + 1}. Vector \eqn{\pi}{pi} is assigned a TSB
  prior distribution (for more information, see \code{\link{dStick}}).

  Elsewhere, each element of \eqn{\theta}{theta} is constrained to the
  interval (0,1), and the original TSB form is beta-distributed with the
  \eqn{\alpha}{alpha} parameter of the beta distribution constrained
  to 1 (Ishwaran and James, 2001). The \eqn{\beta}{beta} hyperparameter
  in the beta distribution is usually gamma-distributed.

  A larger value for a given \eqn{\theta_m}{theta[m]} is associated
  with a higher probability of the associated mixture component,
  however, the proportion changes according to the position of the
  element in the \eqn{\theta}{theta} vector.

  A variety of stick-breaking processes exist. For example, rather than
  each \eqn{\theta}{theta} being beta-distributed, there have been other
  forms introduced such as logistic and probit, among others.
}
\value{
  The \code{Stick} function returns a probability vector wherein each
  element relates to a mixture component.
}
\references{
  Ishwaran, H. and James, L. (2001). "Gibbs Sampling Methods for Stick
  Breaking Priors". \emph{Journal of the American Statistical
  Association}, 96(453), p. 161--173.
  
  Sethuraman, J. (1994). "A Constructive Definition of Dirichlet
  Priors". \emph{Statistica Sinica}, 4, p. 639--650.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{ddirichlet}},
  \code{\link{dmvpolya}}, and 
  \code{\link{dStick}}.
}
\keyword{Stick-Breaking Process}
\keyword{Truncated Stick-Breaking Process}
\keyword{Utility}
