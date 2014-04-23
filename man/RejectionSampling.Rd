\name{RejectionSampling}
\alias{RejectionSampling}
\title{Rejection Sampling}
\description{
  The \code{RejectionSampling} function implements rejection sampling
  of a target density given a proposal density.
}
\usage{
RejectionSampling(Model, Data, mu, S, df=Inf, logc, n=1000, CPUs=1, Type="PSOCK")
}
\arguments{
  \item{Model}{This is a model specification function. For more
    information, see \code{\link{LaplaceApproximation}}.}
  \item{Data}{This is a list of data. For more information, see
    \code{\link{LaplaceApproximation}}.}
  \item{mu}{This is a mean vector \eqn{\mu}{mu} for the multivariate
    normal or multivariate t proposal density.}
  \item{S}{This is a convariance matrix \eqn{\Sigma}{Sigma} for the
    multivariate normal or multivariate t proposal density.}
  \item{df}{This is a scalar degrees of freedom parameter
    \eqn{\nu}{nu}. It defaults to infinity, in which case the
    multivariate normal density is used.}
  \item{logc}{This is the logarithm of the rejection sampling constant.}
  \item{n}{This is the number of independent draws to be simulated from
    the proposal density.}
  \item{CPUs}{This argument accepts an integer that specifies the number
    of central processing units (CPUs) of the multicore computer or
    computer cluster. This argument defaults to \code{CPUs=1}, in which
    parallel processing does not occur.}
  \item{Type}{This argument specifies the type of parallel processing to
    perform, accepting either \code{Type="PSOCK"} or
    \code{Type="MPI"}.}
}
\details{
  Rejection sampling (von Neumann, 1951) is a Monte Carlo method for
  drawing independent samples from a distribution that is proportional
  to the target distribution, \eqn{f(x)}, given a sampling distribution,
  \eqn{g(x)}, from which samples can readily be drawn, and for which
  there is a finite constant \eqn{c}.

  Here, the target distribution, \eqn{f(x)}, is the result of the
  \code{Model} function. The sampling distribution, \eqn{g(x)}, is
  either a multivariate normal or multivariate t-distribution. The
  parameters of \eqn{g(x)} (\code{mu}, \code{S}, and \code{df}) are used
  to create random draws, \eqn{\theta}{theta}, of the sampling
  distribution, \eqn{g(x)}. These draws, \eqn{\theta}{theta}, are used
  to evaluate the target distribution, \eqn{f(x)}, via the \code{Model}
  specification function. The evaluations of the target distribution,
  sampling distribution, and the constant are used to create a
  probability of acceptance for each draw, by comparing to a vector of
  \eqn{n} uniform draws, \eqn{u}. Each draw, \eqn{\theta}{theta} is
  accepted if \deqn{u \le \frac{f(\theta|\textbf{y})}{cg(\theta)}}{u <=
    f(theta|y)/[cg(theta)]}

  Before beginning rejection sampling, a goal of the user is to find the
  bounding constant, \eqn{c}, such that \eqn{f(\theta|\textbf{y}) \le
  cg(\theta)}{f(theta|y) <= cg(theta)} for all
  \eqn{\theta}{theta}. These are all expressed in logarithms, so the
  goal is to find \eqn{\log f(\theta|\textbf{y}) - \log g(\theta) \le
  \log(c)}{log[f(theta|y)] - log[g(theta)] <= log(c)} for all
  \eqn{\theta}{theta}. This is done by maximizing \eqn{\log
  f(\theta|\textbf{y}) - \log g(\theta)}{log[f(theta|y)] -
  log[g(theta)]} over all \eqn{\theta}{theta}. By using, say,
  \code{\link{LaplaceApproximation}} to find the modes of the
  parameters of interest, and using the resultant \code{LP}, the mode
  of the logarithm of the joint posterior distribution, as
  \eqn{\log(c)}{log(c)}.

  The \code{RejectionSampling} function performs one iteration of
  rejection sampling. Rejection sampling is often iterated, then called
  the rejection sampling algorithm, until a sufficient number or
  proportion of \eqn{\theta}{theta} is accepted. An efficient
  rejection sampling algorithm has a high acceptance rate. However,
  rejection sampling becomes less efficient as the model dimension (the
  number of parameters) increases.
  
  Extensions of rejection sampling include Adaptive Rejection
  Sampling (ARS) (either derivative-based or derivative-free) and
  Adaptive Rejection Metropolis Sampling (ARMS), as in Gilks et
  al. (1995). The random-walk Metropolis algorithm (Metropolis et al.,
  1953) combined the rejection sampling (a method of Monte Carlo
  simulation) of von Neumann (1951) with Markov chains.

  Parallel processing may be performed when the user specifies
  \code{CPUs} to be greater than one, implying that the specified number
  of CPUs exists and is available. Parallelization may be performed on a
  multicore computer or a computer cluster. Either a Simple Network of
  Workstations (SNOW) or Message Passing Interface (MPI) is used. With
  small data sets and few samples, parallel processing may be slower,
  due to computer network communication. With larger data sets and more
  samples, the user should experience a faster run-time.

  This function is similar to the \code{rejectsampling} function in the
  \code{LearnBayes} package.
}
\value{
  The \code{RejectionSampling} function returns an object of class
  \code{rejection}, which is a matrix of accepted, independent,
  simulated draws from the target distribution.
}
\references{
  Gilks, W.R., Best, N.G., Tan, K.K.C. (1995). "Adaptive Rejection
  Metropolis Sampling within Gibbs Sampling". Journal of the Royal
  Statistical Society. Series C (Applied Statistics), Vol. 44, No. 4,
  p. 455--472.

  Metropolis, N., Rosenbluth, A.W., Rosenbluth, M.N., and Teller,
  E. (1953). "Equation of State Calculations by Fast Computing
  Machines". Journal of Chemical Physics, 21, p. 1087-1092.

  von Neumann, J. (1951). "Various Techniques Used in Connection with
  Random Digits. Monte Carlo Methods". National Bureau Standards, 12,
  p. 36--38.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{dmvn}},
  \code{\link{dmvt}},
  \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}}, and
  \code{\link{VariationalBayes}}.
}
\examples{
library(LaplacesDemon)
### Suppose an output object of class laplace is called Fit:
#rs <- RejectionSampling(Model, MyData, mu=Fit$Summary1[,1],
#     S=Fit$Covar, df=Inf, logc=Fit$LP.Final, n=1000)
#rs
}
\keyword{High Performance Computing}
\keyword{Monte Carlo}
