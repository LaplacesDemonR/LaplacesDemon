\name{SIR}
\alias{SIR}
\title{Sampling Importance Resampling}
\description{
  The \code{SIR} function performs Sampling Importance Resampling, also
  called Sequential Importance Resampling, and uses a multivariate normal
  proposal density.
}
\usage{
SIR(Model, Data, mu, Sigma, n=1000, CPUs=1, Type="PSOCK")
}
\arguments{
  \item{Model}{This is a model specification function. For more
    information, see \code{\link{LaplaceApproximation}}.}
  \item{Data}{This is a list of data. For more information, see
    \code{\link{LaplaceApproximation}}.}
  \item{mu}{This is a mean vector, \eqn{\mu}{mu}, for a multivariate
    normal distribution, and is usually the posterior means from an
    object of class \code{iterquad} (from
    \code{\link{IterativeQuadrature}}) or class \code{vb} (from
    \code{\link{VariationalBayes}}), or the posterior modes from an
    object of class \code{laplace} (from 
    \code{\link{LaplaceApproximation}}).}
  \item{Sigma}{This is a covariance matrix, \eqn{\Sigma}{Sigma}, for a
    multivariate normal distribution, and is usually the \code{Covar}
    component of an object of class \code{iterquad}, \code{laplace}, or
    \code{vb}.}
  \item{n}{This is the number of samples to be drawn from the posterior
    distribution.}
  \item{CPUs}{This argument accepts an integer that specifies the number
    of central processing units (CPUs) of the multicore computer or
    computer cluster. This argument defaults to \code{CPUs=1}, in which
    parallel processing does not occur.}
  \item{Type}{This argument specifies the type of parallel processing to
    perform, accepting either \code{Type="PSOCK"} or
    \code{Type="MPI"}.}
}
\details{
  Sampling Importance Resampling (SIR) was introduced in Gordon, et
  al. (1993), and is the original particle filtering algorithm (and this
  family of algorithms is also known as Sequential Monte Carlo). A
  distribution is approximated with importance weights, which are
  approximations to the relative posterior densities of the particles,
  and the sum of the weights is one. In this terminology, each sample in
  the distribution is a ``particle''. SIR is a sequential or recursive
  form of importance sampling. As in importance sampling, the
  expectation of a function can be approximated as a weighted
  average. The optimal proposal distribution is the target distribution.

  In the \code{LaplacesDemon} package, the main use of the \code{SIR}
  function is to produce posterior samples for iterative quadrature,
  Laplace Approximation, or Variational Bayes, and \code{SIR} is called
  behind-the-scenes by the \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}}, or \code{\link{VariationalBayes}}
  function.

  Iterative quadrature estimates the posterior mean and the associated
  covariance matrix. Assuming normality, this output characterizes the
  marginal posterior distributions. However, it is often useful to have
  posterior samples, in which case the \code{SIR} function is used to
  draw samples. The number of samples, \code{n}, should increase with
  the number and intercorrelations of the parameters. Otherwise,
  multimodal posterior distributions may occur.

  Laplace Approximation estimates the posterior mode and the associated
  covariance matrix. Assuming normality, this output characterizes the
  marginal posterior distributions. However, it is often useful to have
  posterior samples, in which case the \code{SIR} function is used to
  draw samples. The number of samples, \code{n}, should increase with
  the number and intercorrelations of the parameters. Otherwise,
  multimodal posterior distributions may occur.

  Variational Bayes estimates both the posterior mean and
  variance. Assuming normality, this output characterizes the marginal
  posterior distributions. However, it is often useful to have posterior
  samples, in which case the \code{SIR} function is used to draw
  samples. The number of samples, \code{n}, should increase with the
  number of intercorrelations of the parameters. Otherwise, multimodal
  posterior distributions may occur.

  SIR is also commonly used when considering a mild change in a prior
  distribution. For example, suppose a model was updated in
  \code{\link{LaplacesDemon}}, and it had a least-informative prior
  distribution, but the statistician would like to estimate the impact
  of changing to a weakly-informative prior distribution. The change is
  made in the model specification function, and the posterior means and
  covariance are supplied to the \code{SIR} function. The returned
  samples are estimates of the posterior, given the different prior
  distribution. This is akin to sensitivity analysis (see the
  \code{SensitivityAnalysis} function).

  In other contexts (for which this function is not designed), SIR is
  used with dynamic linear models (DLMs) and state-space models (SSMs)
  for state filtering.

  Parallel processing may be performed when the user specifies
  \code{CPUs} to be greater than one, implying that the specified number
  of CPUs exists and is available. Parallelization may be performed on a
  multicore computer or a computer cluster. Either a Simple Network of
  Workstations (SNOW) or Message Passing Interface (MPI) is used. With
  small data sets and few samples, parallel processing may be slower,
  due to computer network communication. With larger data sets and more
  samples, the user should experience a faster run-time.
  
  This function was adapted from the \code{sir} function in the
  \code{LearnBayes} package.
}
\value{
  The \code{SIR} function returns a matrix of samples drawn from the
  posterior distribution.
}
\references{
  Gordon, N.J., Salmond, D.J., and Smith, A.F.M. (1993). "Novel Approach
  to Nonlinear/Non-Gaussian Bayesian State Estimation". \emph{IEEE
  Proceedings F on Radar and Signal Processing}, 140(2), p. 107--113.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{dmvn}},
  \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{SensitivityAnalysis}}, and
  \code{\link{VariationalBayes}}.
}
\keyword{High Performance Computing}
\keyword{Sampling Importance Resampling}
\keyword{Sequential Importance Resampling}
\keyword{Utility}
