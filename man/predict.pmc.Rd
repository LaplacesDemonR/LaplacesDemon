\name{predict.pmc}
\alias{predict.pmc}
\title{Posterior Predictive Checks}
\description{
  This may be used to predict either new, unobserved instances of
  \eqn{\textbf{y}}{y} (called \eqn{\textbf{y}^{new}}{y[new]}) or
  replicates of \eqn{\textbf{y}}{y} (called
  \eqn{\textbf{y}^{rep}}{y[rep]}), and then perform posterior
  predictive checks. Either \eqn{\textbf{y}^{new}}{y[new]} or
  \eqn{\textbf{y}^{rep}}{y[rep]} is predicted given an object of
  class \code{demonoid}, the model specification, and data.
}
\usage{\method{predict}{pmc}(object, Model, Data, CPUs=1, Type="PSOCK", \dots)}
\arguments{
  \item{object}{An object of class \code{pmc} is required.}
  \item{Model}{The model specification function is required.}
  \item{Data}{A data set in a list is required. The dependent variable
    is required to be named either \code{y} or \code{Y}.}
  \item{CPUs}{This argument accepts an integer that specifies the number
    of central processing units (CPUs) of the multicore computer or
    computer cluster. This argument defaults to \code{CPUs=1}, in which
    parallel processing does not occur.}
  \item{Type}{This argument specifies the type of parallel processing to
    perform, accepting either \code{Type="PSOCK"} or
    \code{Type="MPI"}.}
  \item{\dots}{Additional arguments are unused.}
}
\details{
  This function passes each iteration of marginal posterior samples
  along with data to \code{Model}, where the fourth component in the
  return list is labeled \code{yhat}, and is a vector of expectations of
  \eqn{\textbf{y}}{y}, given the samples, model specification, and
  data. Stationary samples are used if detected, otherwise
  non-stationary samples will be used. To predict
  \eqn{\textbf{y}^{rep}}{y[rep]}, simply supply the data set used to
  estimate the model. To predict \eqn{\textbf{y}^{new}}{y[new]}, supply
  a new data set instead (though for some model specifications, this
  cannot be done, and \eqn{\textbf{y}_{new}}{y[new]} must be specified
  in the \code{Model} function). If the new data set does not have
  \eqn{\textbf{y}}{y}, then create \code{y} in the list and set it
  equal to something sensible, such as \code{mean(y)} from the
  original data set.

  The variable \code{y} must be a vector. If instead it is matrix
  \code{Y}, then it will be converted to vector \code{y}. The vectorized
  length of \code{y} or \code{Y} must be equal to the vectorized length
  of \code{yhat}, the fourth component of the return list of the
  \code{Model} function.

  Parallel processing may be performed when the user specifies
  \code{CPUs} to be greater than one, implying that the specified number
  of CPUs exists and is available. Parallelization may be performed on a
  multicore computer or a computer cluster. Either a Simple Network of
  Workstations (SNOW) or Message Passing Interface is used (MPI). With
  small data sets and few samples, parallel processing may be slower,
  due to computer network communication. With larger data sets and more
  samples, the user should experience a faster run-time.

  For more information on posterior predictive checks, see
  \url{https://web.archive.org/web/20150215050702/http://www.bayesian-inference.com/posteriorpredictivechecks}.
}
\value{
  This function returns an object of class \code{pmc.ppc} (where
  ppc stands for posterior predictive checks). The returned object is
  a list with the following components:

  \item{y}{This stores the vectorized form of \eqn{\textbf{y}}{y}, the
    dependent variable.}
  \item{yhat}{This is a \eqn{N \times S}{N x S} matrix, where \eqn{N} is
    the number of records of \eqn{\textbf{y}}{y} and \eqn{S} is the
    number of posterior samples.}
  \item{Deviance}{This is a vector of predictive deviance.}
}
\author{Statisticat, LLC.}
\seealso{
  \code{\link{PMC}}
}
\keyword{High Performance Computing}
\keyword{Posterior Predictive Checks}
\keyword{Predict}
