\name{logit}
\alias{invlogit}
\alias{logit}
\title{The logit and inverse-logit functions}
\description{
  The logit and inverse-logit (also called the logistic function) are
  provided.
}
\usage{
invlogit(x)
logit(p)
}
\arguments{
  \item{x}{This object contains real values that will be transformed to
    the interval [0,1].}
  \item{p}{This object contains of probabilities p in the interval [0,1]
    that will be transformed to the real line.}
}
\details{
  The \code{logit} function is the inverse of the sigmoid or logistic
  function, and transforms a continuous value (usually probability
  \eqn{p}) in the interval [0,1] to the real line (where it is usually
  the logarithm of the odds). The \code{logit} function is \eqn{\log(p /
    (1-p))}{log(p / (1 - p))}.
  
  The \code{invlogit} function (called either the inverse logit or the
  logistic function) transforms a real number (usually the logarithm of
  the odds) to a value (usually probability \eqn{p}) in the interval
  [0,1]. The \code{invlogit} function is \eqn{\frac{1}{1 + \exp(-x)}}{1
    / (1 + exp(-x))}.

  If \eqn{p} is a probability, then \eqn{\frac{p}{1-p}}{p/(1-p)} is the
  corresponding odds, while the \code{logit} of \eqn{p} is the logarithm
  of the odds. The difference between the logits of two probabilities is
  the logarithm of the odds ratio. The derivative of probability \eqn{p}
  in a logistic function (such as \code{invlogit}) is: \eqn{\frac{d}{dx}
    = p(1-p)}{(d / dx) = p * (1 - p)}.

  In the LaplacesDemon package, it is common to re-parameterize a model
  so that a parameter that should be in an interval can be updated from
  the real line by using the \code{logit} and \code{invlogit} functions,
  though the \code{\link{interval}} function provides an
  alternative. For example, consider a parameter \eqn{\theta}{theta}
  that must be in the interval [0,1]. The algorithms in
  \code{\link{IterativeQuadrature}}, \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}}, \code{\link{PMC}}, and
  \code{\link{VariationalBayes}} are unaware of the desired interval,
  and may attempt \eqn{\theta}{theta} outside of this interval. One
  solution is to have the algorithms update \code{logit(theta)} rather
  than \code{theta}. After \code{logit(theta)} is manipulated by the
  algorithm, it is transformed via \code{invlogit(theta)} in the model
  specification function, where \eqn{\theta \in [0,1]}{theta in [0,1]}.
}
\value{
  \code{invlogit} returns probability \code{p}, and
  \code{logit} returns \code{x}.
}
\seealso{
  \code{\link{interval}},
  \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{plogis}},
  \code{\link{PMC}},
  \code{\link{qlogis}}, and
  \code{\link{VariationalBayes}}.
}
\examples{
library(LaplacesDemon)
x <- -5:5
p <- invlogit(x)
x <- logit(p)
}
\keyword{Inverse Logit}
\keyword{Link Function}
\keyword{Logistic Function}
\keyword{Logit}
\keyword{Transformation}
