\name{Matrices}
\alias{as.indicator.matrix}
\alias{as.inverse}
\alias{as.parm.matrix}
\alias{as.positive.definite}
\alias{as.positive.semidefinite}
\alias{as.symmetric.matrix}
\alias{Cov2Cor}
\alias{CovEstim}
\alias{GaussHermiteCubeRule}
\alias{Hessian}
\alias{is.positive.definite}
\alias{is.positive.semidefinite}
\alias{is.square.matrix}
\alias{is.symmetric.matrix}
\alias{Jacobian}
\alias{logdet}
\alias{lower.triangle}
\alias{read.matrix}
\alias{SparseGrid}
\alias{TransitionMatrix}
\alias{tr}
\alias{upper.triangle}
\title{Matrix Utility Functions}
\description{
  These are utility functions for working with matrices.
}
\usage{
as.indicator.matrix(x)
as.inverse(x)
as.parm.matrix(x, k, parm, Data, a=-Inf, b=Inf, restrict=FALSE, chol=FALSE)
as.positive.definite(x)
as.positive.semidefinite(x)
as.symmetric.matrix(x, k=NULL)
is.positive.definite(x)
is.positive.semidefinite(x)
is.square.matrix(x)
is.symmetric.matrix(x)
Cov2Cor(Sigma)
CovEstim(Model, parm, Data, Method="Hessian")
GaussHermiteCubeRule(N, dims, rule)
Hessian(Model, parm, Data, Interval=1e-6, Method="Richardson")
Jacobian(Model, parm, Data, Interval=1e-6, Method="simple")
logdet(x)
lower.triangle(x, diag=FALSE)
read.matrix(file, header=FALSE, sep=",", nrow=0, samples=0, size=0, na.rm=FALSE)
SparseGrid(J, K)
TransitionMatrix(theta.y=NULL, y.theta=NULL, p.theta=NULL)
tr(x)
upper.triangle(x, diag=FALSE)
}
\arguments{
  \item{N}{This required argument accepts a positive integer that
    indicates the number of nodes.}
  \item{x}{This is a matrix (though \code{as.symmetric.matrix} also
    accepts vectors).}
  \item{J}{This required argument indicates the dimension of the
    integral and accepts a positive integer.}
  \item{k}{For \code{as.parm.matrix}, this is a required argument,
    indicating the dimension of the matrix. For
    \code{as.symmetric.matrix}, this is an optional argument that
    specifies the dimension of the symmetric matrix. This applies
    only when \code{x} is a vector. It defaults to \code{NULL},
    in which case it calculates \code{k <- (-1 + sqrt(1 + 8 *
      length(x)))/ 2}.}
  \item{K}{This required argument indicates the accuracy and accepts a
    positive integer. Larger values result in many more integration
    nodes.}
  \item{diag}{Logical. If \code{TRUE}, then the elements in the main
    diagonal are also returned.}
  \item{dims}{This required argument indicates the dimension of the
    integral and accepts a positive integer.}
  \item{Sigma}{This is a covariance matrix, \eqn{\Sigma}{Sigma}, and may
    be entered either as a matrix or vector.}
  \item{Model}{This is a model specification function. For more
    information, see \code{\link{LaplacesDemon}}.}
  \item{parm}{This is a vector of parameters passed to the model
    specification.}
  \item{Data}{This is the list of data passed to the model
    specification. For more information, see \code{\link{LaplacesDemon}}.}
  \item{a,b}{These optional arguments allow the elements of \code{x} to
    be bound to the interval \eqn{[a,b]}. For example, elements of a
    correlation matrix are in the interval \eqn{[-1,1]}.}
  \item{restrict}{Logical. If \code{TRUE}, then \code{x[1,1]} is
    restricted to 1. This is useful in multinomial probit, for
    example. The variable, \code{LaplacesDemonMatrix}, is
    created in a new environment, \code{LDEnv} so \code{as.parm.matrix}
    can keep track of changes from iteration to iteration.}
  \item{rule}{This is an optional argument that accepts a univariate
    Gauss-Hermite quadrature rule. Usually, this argument is left empty.
    A rule may be supplied that differs from the traditional rule, such
    as when constraints have been observed, and one or more nodes and
    weights were adjusted.}
  \item{chol}{Logical. If \code{TRUE}, then x is an upper-triangular
    matrix.}
  \item{file}{This is the name of the file from which the numeric data
    matrix will be imported or read.}
  \item{header}{Logical. When \code{TRUE}, the first row of the file
    must contain names of the columns, and will be converted to the
    column names of the numeric matrix. When \code{FALSE}, the first row
    of the file contains data, not column names.}
  \item{Interval}{This accepts a small scalar number for precision.}
  \item{Method}{This accepts a quoted string. For \code{Hessian}, it
    defaults to \code{Method="Richardson"}, which uses Richardson
    extrapolation. For \code{Jacobian}, it defaults to
    \code{Method="simple"}, which uses finite-differencing. Richardson
    Richardson extrapolation is more accurate, but slower to calculate.
    Since error due to finite-differencing propagates through to higher
    derivatives, finite-differencing should not be used when
    approximating a Hessian matrix. Another method called automatic
    differentiation is not currently available here, but should be more
    accurate, though even slower to calculate. Another popular
    alternative is to use the \code{\link{BayesianBootstrap}} on the
    data. For \code{CovEstim}, this accepts \code{Method="Hessian"},
    \code{Method="Identity"} (which simply assigns an identity matrix),
    \code{Method="OPG"} (which calculates the sum of outer products of
    record-level gradients), or \code{Method="Sandwich"}, which is the
    sandwich estimator and combines the Hessian and OPG estimates.}
  \item{nrow}{This is the number of rows of the numeric matrix, and
    defaults to \code{nrow=0}. If the number is known, the function will
    perform noticeably faster when it does not have to check.}
  \item{p.theta}{This accepts a matrix of prior probabilities for a
    transition matrix, and defaults to \code{NULL}. If used, each row
    must sum to 1.}
  \item{samples}{This is the number of samples to take from the numeric
    matrix. When \code{samples=0}, sampling is not performed and the
    entire matrix is returned.}
  \item{sep}{This argument indicates a character with which it will
    separate fields when creating column vectors. For example, a
    read a comma-separated file (.csv), use \code{sep=","}.}
  \item{size}{This is the batch size to be used only when reading a
    numeric matrix that is larger than the available computer memory
    (RAM), and only when \code{samples} is greater than zero. Sampling
    of a big data matrix is performed by first determining the records
    to keep, and then reading batches, one by one, and keeping the
    matching records.}
  \item{theta.y}{This accepts a vector of posterior samples of a
    discrete Markov chain, and defaults to \code{NULL}. If used, the
    order of the samples affects the transition probability.}
  \item{na.rm}{Logical. When \code{TRUE}, rows with missing values are
    removed from the matrix after it is read. Rather than removing
    missing values, the user may consider imputing missing values inside
    the model, or before the model with the \code{\link{MISS}}
    function. Examples of within-model imputation may be found in the
    accompanying "Examples" vignette.}
  \item{y.theta}{This accepts a vector of data that are samples of a
    discrete distribution, and defaults to \code{NULL}. If used, the
    order of the samples affects the transition probability.}
}
\details{
  The \code{as.indicator.matrix} function creates an indicator matrix
  from a vector. This function is useful for converting a discrete
  vector into a matrix in which each column represents one of the
  discrete values, and each occurence of that value in the related
  column is indicated by a one, and is otherwise filled with
  zeroes. This function is similar to the \code{class.ind} function in
  the nnet package.

  The \code{as.inverse} function returns the matrix inverse of
  \code{x}. The \code{solve} function in base R also returns the matrix
  inverse, but \code{solve} can return a matrix that is not symmetric,
  and can fail due to singularities. The \code{as.inverse} function
  tries to use the \code{solve} function to return a matrix inverse, and
  when it fails due to a singularity, \code{as.inverse} uses eigenvalue
  decomposition (in which eigenvalues below a tolerance are replaced
  with the tolerance), and coerces the result to a symmetric
  matrix. This is similar to the \code{solvcov} function in the fpc
  package.

  The \code{as.parm.matrix} function prepares a correlation, covariance,
  or precision matrix in two important ways. First,
  \code{as.parm.matrix} obtains the parameters for the matrix specified
  in the \code{x} argument by matching the name of the matrix in the
  \code{x} argument with any parameters in \code{parm}, given the
  parameter names in the \code{Data} listed in \code{parm.names}. These
  obtained parameters are organized into a matrix as the elements of the
  upper-triangular, including the diagonal. A copy is made, without the
  diagonal, and the lower-triangular is filled in, completing the
  matrix. Second, \code{as.parm.matrix} checks for
  positive-definiteness. If matrix \code{x} is positive-definite, then
  the matrix is stored as a variable called \code{LaplacesDemonMatrix}
  in a new environment called \code{LDEnv}. If matrix \code{x} is not
  positive-definite, then \code{LaplacesDemonMatrix} in \code{LDEnv} is
  sought as a replacement. If this variable exists, then it is used to
  replace the matrix. If not, then the matrix is replaced with an
  identity matrix. Back in the model specification, after using
  \code{as.parm.matrix}, it is recommended that the user also pass the
  resulting matrix back into the \code{parm} vector, so the sampler or
  algorithm knows that the elements of the matrix have changed.

  The \code{as.positive.definite} function returns the nearest
  positive-definite matrix for a matrix that is square and symmetric
  (Higham, 2002). This version is intended only for covariance and
  precision matrices, and has been optimized for speed. A more
  extensible function is \code{nearPD} in the matrixcalc package, which
  is also able to work with correlation matrices, and matrices that are
  asymmetric.

  The \code{as.positive.semidefinite} function iteratively seeks to
  return a square, symmetric matrix that is at least
  positive-semidefinite, by replacing each negative eigenvalue and
  calculating its projection. This is intended only for covariance and
  precision matrices. A similar function is \code{makePsd} in the RTAQ
  package, though it is not iterative, and returns matrices that fail a
  logical check with \code{is.positive.semidefinite}.

  The \code{as.symmetric.matrix} function accepts either a vector or
  matrix, and returns a symmetric matrix. In the case of a vector, it
  can be either all elements of the matrix, or the lower triangular. In
  the case of a \code{x} being entered as a matrix, this function
  tolerates non-finite values in one triangle (say, the lower), as long
  as the corresponding element is finite in the other (say, the upper)
  triangle.

  The \code{Cov2Cor} function converts a covariance matrix into a
  correlation matrix, and accepts the covariance matrix either in matrix
  or vector form. This function may be useful inside a model
  specification and also with converting posterior draws of the elements
  of a covariance matrix to a correlation matrix. \code{Cov2Cor} is an
  expanded form of the \code{cov2cor} function in the \code{stats}
  package, where \code{Cov2Cor} is also able to accept and return a
  vectorized matrix.

  The \code{CovEstim} function estimates a covariance matrix with one of
  several methods. This is mainly used by
  \code{\link{LaplaceApproximation}}, where the \code{parm} argument
  receives the posterior modes. See the \code{CovEst} argument for
  more details.

  The \code{GaussHermiteCubeRule} function returns a matrix of nodes and
  a vector of weights for a \code{dims}-dimensional integral given
  \eqn{N} univariate nodes. The number of multivariate nodes will differ
  from the number of univariate nodes. This function is for use with
  multivariate quadrature, often called cubature. This has been adapted
  from the \code{multiquad} function in the NominalLogisticBiplot
  package. The \code{\link{GaussHermiteQuadRule}} function is a
  univariate version. A customized univariate \code{rule} may be
  supplied when constraints necessitate that one or more nodes and
  weights had to be altered.

  The \code{Hessian} returns a symmetric, Hessian matrix, which is a
  matrix of second partial derivatives. The estimation of the Hessian
  matrix is approximated numerically using Richardson extrapolation by
  default. This is a slow function. This function is not intended to be
  called by the user, but is made available here. This is essentially
  the \code{hessian} function from the numDeriv package, adapted to
  Laplace's Demon.
  
  The \code{is.positive.definite} function is a logical test of whether
  or not a matrix is positive-definite. A \eqn{k \times k}{k x k}
  symmetric matrix \eqn{\textbf{X}}{X} is positive-definite if all of
  its eigenvalues are positive (\eqn{\lambda_i > 0, i \in k}{lambda[i] >
  0, i in k}). All main-diagonal elements must be positive. The
  determinant of a positive-definite matrix is always positive, so a
  positive-definite matrix is always nonsingular. Non-symmetric,
  positive-definite matrices exist, but are not considered here.

  The \code{is.positive.semidefinite} function is a logical test of
  whether or not a matrix is positive-semidefinite. A \eqn{k x k}
  symmetric matrix \eqn{\textbf{X}}{X} is positive-semidefinite if all
  of its eigenvalues are non-negative (\eqn{\lambda_i \ge 0, i \in
  k}{lambda[i] >= 0, i ink}).

  The \code{is.square.matrix} function is a logical test of whether or
  not a matrix is square. A square matrix is a matrix with the same
  number of rows and columns, and is usually represented as a \eqn{k
  \times k}{k x k} matrix \eqn{\textbf{X}}{X}.

  The \code{is.symmetric.matrix} function is a logical test of whether
  or not a matrix is symmetric. A symmetric matrix is a square matrix
  that is equal to its transpose, \eqn{\textbf{X} = \textbf{X}^T}{X =
  t(X)}. For example, where \eqn{i} indexes rows and \eqn{j} indexes
  columns, \eqn{\textbf{X}_{i,j} = \textbf{X}_{j,i}}{X[i,j] =
  X[j,i]}. This differs from the \code{isSymmetric} function in base R
  that is inexact, using \code{all.equal}.

  The \code{Jacobian} function estimates the Jacobian matrix, which is
  a matrix of all first-order partial derivatives of the \code{Model}.
  The Jacobian matrix is estimated by default with forward
  finite-differencing, or optionally with Richardson extrapolation. This
  function is not intended to be called by the user, but is made
  available here. This is essentially the \code{jacobian} function from
  the numDeriv package, adapted to LaplacesDemon.

  The \code{logdet} function returns the logarithm of the determinant of
  a positive-definite matrix via the Cholesky decomposition. The
  determinant is a value associated with a square matrix, and was used
  historically to \emph{determine} if a system of linear equations has a
  unique solution. The term \emph{determinant} was introduced by Gauss,
  where Laplace referred to it as the resultant. When the determinant is
  zero, the matrix is singular and non-invertible; there are either no
  solutions or many solutions. A unique solution exists when the
  determinant is non-zero. The \code{det} function in base R works well
  for small matrices, but can return erroneously return zero in larger
  matrices. It is better to work with the log-determinant.

  The \code{lower.triangle} function returns a vector of the lower
  triangular elements of a matrix, and the diagonal is included when
  \code{diag=TRUE}.

  The \code{read.matrix} function is provided here as one of many
  convenient ways to read a numeric matrix into R. The most common
  method of storing data in R is the data frame, because it is
  versatile. For example, a data frame may contain character, factor,
  and numeric variables together. For iterative estimation, common in
  Bayesian inference, the data frame is much slower than the numeric
  matrix. For this reason, the LaplacesDemon package does not use data
  frames, and has not traditionally accepted character or factor
  data. The \code{read.matrix} function returns either an entire numeric
  matrix, or row-wise samples from a numeric matrix. Samples may be
  taken from a matrix that is too large for available computer memory
  (RAM), such as with big data.

  The \code{SparseGrid} function returns a sparse grid for a
  \eqn{J}-dimensional integral with accuracy \eqn{K}, given
  Gauss-Hermite quadrature rules. A grid of order eqn{K} provides an
  exact result for a polynomial of total order of \eqn{2K - 1} or less.
  \code{SparseGrid} returns a matrix of nodes and a vector of weights.
  A sparse grid is more efficient than the full grid in the
  \code{GaussHermiteCubeRule} function. This has been adapted from the
  SparseGrid package.

  The \code{TransitionMatrix} function has several uses. A user may
  supply a vector of marginal posterior samples of a discrete Markov
  chain as \code{theta.y}, and an observed posterior transition matrix
  is returned. Otherwise, a user may supply data (\code{y.theta}) and/or
  a prior (\code{p.theta}), in which case a posterior transition matrix
  is returned. A common row-wise prior is the dirichlet distribution.
  Transition probabilities are from row element to column element.

  The \code{tr} function returns the trace of a matrix. The trace of a
  matrix is the sum of the elements in the main diagonal of a square
  matrix. For example, the trace of a \eqn{k \times k}{k x k} matrix
  \eqn{\textbf{X}}{X}, is \eqn{\sum_{k=1} \textbf{X}_{k,k}}{sum(k=1)
    X[k,k]}.

  The \code{upper.triangle} function returns a vector of the lower
  triangular elements of a matrix, and the diagonal is included when
  \code{diag=TRUE}.
}
\references{
  Higham, N.J. (2002). "Computing the Nearest Correlation Matrix - a
  Problem from Finance". \emph{IMA Journal of Numerical Analysis}, 22,
  p. 329--343.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{BayesianBootstrap}},
  \code{\link{Cov2Prec}},
  \code{\link{cov2cor}},
  \code{\link{ddirichlet}},
  \code{\link{GaussHermiteQuadRule}},
  \code{\link{isSymmetric}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{lower.tri}},
  \code{\link{MISS}},
  \code{\link{Prec2Cov}},
  \code{\link{solve}}, and
  \code{\link{upper.tri}}.
}
\keyword{Array}
\keyword{High Performance Computing}
\keyword{Math}
\keyword{Matrix}
\keyword{Utility}