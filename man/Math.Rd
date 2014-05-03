\name{Math}
\alias{GaussHermiteQuadRule}
\alias{Hermite}
\alias{logadd}
\alias{partial}
\title{Math Utility Functions}
\description{
  These are utility functions for math.
}
\usage{
GaussHermiteQuadRule(N)
Hermite(x, N, prob=TRUE)
logadd(x, add=TRUE)
partial(Model, parm, Data, Interval=1e-6, Method="simple")
}
\arguments{
  \item{N}{This required argument accepts a positive integer that
    indicates the number of nodes.}
  \item{x}{This is a numeric vector.}
  \item{add}{Logical. This defaults to \code{TRUE}, in which case
    \eqn{\log(x+y)}{log(x+y)} is performed. Otherwise,
    \eqn{\log(x-y)}{log(x-y)} is performed.}
  \item{Model}{This is a model specification function. For more
    information, see \code{\link{LaplacesDemon}}.}
  \item{parm}{This is a vector parameters.}
  \item{prob}{Logical. This defaults to \code{TRUE}, which uses the
    probabilist's kernel for the Hermite polynomial. Otherwise,
    \code{FALSE} uses the physicist's kernel.}
  \item{Data}{This is a list of data. For more information, see
    \code{\link{LaplacesDemon}}.}
  \item{Interval}{This is the interval of numeric differencing.}
  \item{Method}{This accepts a quoted string, and defaults to
    "simple", which is finite-differencing. Alternatively
    \code{Method="Richardson"} uses Richardson extrapolation, which
    is more accurate, but takes longer to calculate. Another method
    called automatic differentiation is currently unsupported, but
    is even more accurate, and takes even longer to calculate.}
}
\details{
  The \code{GaussHermiteQuadRule} function returns nodes and weights for
  univariate Gauss-Hermite quadrature. The nodes and weights are
  obtained from a tridiagonal eigenvalue problem. Weights are calculated
  from the physicist's (rather than the probabilist's) kernel. This has
  been adapted from the GaussHermite function in the pracma package. The
  \code{\link{GaussHermiteCubeRule}} function is a multivariate version.
  This is used in the \code{\link{IterativeQuadrature}} function.

  The \code{Hermite} function evaluates a Hermite polynomial of degree
  \eqn{N} at \eqn{x}, using either the probabilist's (\code{prob=TRUE})
  or physicist's (\code{prob=FALSE}) kernel. This function was adapted
  from the \code{hermite} function in package EQL.
  
  The \code{logadd} function performs addition (or subtraction) when the
  terms are logarithmic. The equations are:

  \deqn{\log(x+y) = \log(x) + \log(1 + \exp(\log(y) - \log(x)))}{log(x+y) =
    log(x) + log(1 + exp(log(y) - log(x)))}
  \deqn{\log(x-y) = \log(x) + \log(1 - \exp(\log(y) - \log(x)))}{log(x-y) =
    log(x) + log(1 - exp(log(y) - log(x)))}
  
  The \code{partial} function estimates partial derivatives of
  parameters in a model specification with data, using either
  forward finite-differencing or Richardson extrapolation. In calculus,
  a partial derivative of a function of several variables is its
  derivative with respect to one of those variables, with the others
  held constant. Related functions include \code{Jacobian} which returns
  a matrix of first-order partial derivatives, and \code{Hessian}, which
  returns a matrix of second-order partial derivatives of the model
  specification function with respect to its parameters. The
  \code{partial} function is not intended to be called by the user, but
  is used by other functions. This is essentially the \code{grad}
  function in the numDeriv package, but defaulting to forward
  finite-differencing with a smaller interval.
}
\value{
  \code{logadd} returns the result of \eqn{\log(x+y)}{log(x+y)} or
  \eqn{\log(x-y)}{log(x-y)}.
  
  \code{partial} returns a vector of partial derivatives.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{GaussHermiteCubeRule}},
  \code{\link{Hessian}},
  \code{\link{IterativeQuadrature}},
  \code{\link{Jacobian}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}}, and
  \code{\link{VariationalBayes}}.
}
\keyword{Math}