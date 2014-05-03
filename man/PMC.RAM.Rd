\name{PMC.RAM}
\alias{PMC.RAM}
\title{PMC RAM Estimate}
\description{
  This function estimates the random-access memory (RAM) required to
  update a given model and data with \code{\link{PMC}}.

  \emph{Warning:} Unwise use of this function may crash a computer, so
  please read the details below.
}
\usage{
PMC.RAM(Model, Data, Iterations, Thinning, M, N)
}
\arguments{
  \item{Model}{This is a model specification function. For more
    information, see \code{\link{PMC}}.}
  \item{Data}{This is a list of Data. For more information, see
    \code{\link{PMC}}.}
  \item{Iterations}{This is the number of iterations for which
    \code{\link{PMC}} would update. For more information,
    see \code{\link{PMC}}.}
  \item{Thinning}{This is the amount of thinning applied to the samples
    in \code{\link{PMC}}.For more information, see \code{\link{PMC}}.}
  \item{M}{This is the number of mixture components in
    \code{\link{PMC}}.}
  \item{N}{This is the number of samples in \code{\link{PMC}}.}
}
\details{
  The \code{PMC.RAM} function uses the
  \code{\link{object.size}} function to estimate the size in MB of RAM
  required to update in \code{\link{PMC}} for a given model and data,
  and for a number of iterations and specified thinning. When RAM is
  exceeded, the computer will crash. This function can be useful when
  trying to estimate how many samples and iterations to update a
  model without crashing the computer. However, when estimating the
  required RAM, \code{PMC.RAM} actually creates several large
  objects, such as \code{post} (see below). If too many iterations are
  given as an argument to \code{PMC.RAM}, for example, then it
  will crash the computer while trying to estimate the required RAM.

  The best way to use this function is as follows. First, prepare the
  model specification and list of data. Second, observe how much RAM the
  computer is using at the moment, as well as the maximum available
  RAM. The majority of the difference of these two is the amount of RAM
  the computer may dedicate to updating the model. Next, use this
  function with a small number of iterations. Note the estimated
  RAM. Increase the number of iterations, and again note the
  RAM. Continue to increase the number of iterations until, say,
  arbitrarily within 90\% of the above-mentioned difference in RAM.

  The computer operating system uses RAM, as does any other software
  running at the moment. R is currently using RAM, and other functions
  in the \code{LaplacesDemon} package, and any other package that is
  currently activated, are using RAM. There are numerous small objects
  that are not included in the returned list, that use RAM. For example,
  perplexity is a small vector, etc.

  A potentially large objects that is not included is a matrix used for
  estimating \code{\link{LML}}.
}
\value{
  \code{PMC.RAM} returns a list with several components. Each component
  is an estimate in MB for an object. The list has the following
  components:
  \item{alpha}{This is the estimated size in MB of RAM required for the
  matrix of mixture probabilities by iteration.}
  \item{Covar}{This is the estimated size in MB of RAM required for the
    covariance matrix or matrices.}
  \item{Data}{This is the estimated size in MB of RAM required for the
    list of data.}
  \item{Deviance}{This is the estimated size in MB of RAM required for
    the deviance vector before thinning.}
  \item{Initial.Values}{This is the estimated size in MB of RAM required
    for the matrix or vector of initial values.}
  \item{LH}{This is the estimated size in MB of RAM required for the
    \eqn{N \times T \times M}{N x T x M} array \code{LH}, where \eqn{N}
    is the number of samples, \eqn{T} is the number of iterations, and
    \eqn{M} is the number of mixture components. The \code{LH} array is
    not returned by \code{\link{PMC}}.}
  \item{LP}{This is the estimated size in MB of RAM required for the
    \eqn{N \times T \times M}{N x T x M} array \code{LP}, where \eqn{N}
    is the number of samples, \eqn{T} is the number of iterations, and
    \eqn{M} is the number of mixture components. The \code{LP} array is
    not returned by \code{\link{PMC}}.}
  \item{Model}{This is the estimated size in MB of RAM required for the
    model specification function.}
  \item{Monitor}{This is the estimated size in MB of RAM required for
    the \eqn{N \times J}{N x J} matrix \code{Monitor}, where \eqn{N} is
    the number of unthinned samples and J is the number of monitored
    variables. Although it is thinned later in the algorithm, the full
    matrix is created.}
  \item{Posterior1}{This is the estimated size in MB of RAM required for
    the \eqn{N \times J \times T \times M}{N x J x T x M} array
    \code{Posterior1}, where \eqn{N} is the number of samples, \eqn{J}
    is the number of parameters, \eqn{T} is the number of iterations,
    and \eqn{M} is the number of mixture components.}
  \item{Posterior2}{This is the estimated size in MB of RAM required for
    the \eqn{N \times J}{N x J} matrix \code{Posterior2}, where \eqn{N}
    is the number of samples and \eqn{J} is the number of initial values
    or parameters. Although this is thinned later, at one point it is
    un-thinned.}
  \item{Summary}{This is the estimated size in MB of RAM required for
    the summary table.}
  \item{W}{This is the estimated size in MB of RAM required for the
    matrix of importance weights.}
  \item{Total}{This is the estimated size in MB of RAM required in total
    to update with \code{\link{PMC}} for a given model and data, and for
    a number of iterations, specified thinning, mixture components, and
    number of samples.}
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{BigData}},
  \code{\link{LML}},
  \code{\link{object.size}}, and
  \code{\link{PMC}}.
}
\keyword{Memory}

