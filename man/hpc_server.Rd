\name{hpc_server}
\alias{server_Listening}
\title{Server Listening}
\description{
  This function is not intended to be called directly by the user. It is
  an internal-only function to prevent cluster problems while using the
  \code{INCA} algorithm in the \code{LaplacesDemon.hpc} function.
}
\usage{
server_Listening(n=2, port=19009)
}
\arguments{
  \item{n}{This is the number of CPUs. For more information, see
    \code{\link{LaplacesDemon.hpc}}.}
  \item{port}{This is a port for server listening, and defaults to
    port \code{19009}.}
}
\details{
  For the \code{INCA} algorithm, a server has been built into the
  \code{LaplacesDemon.hpc} function. The server exchanges information
  between processes, and has been designed to be portable. The
  \code{server_Listening} function is run as a separate process via the
  \code{system} function, when \code{INCA} is selected in
  \code{LaplacesDemon.hpc}.

  Socket connections and the \code{serialize} function are used as per
  the \pkg{Snow} package to update a single proposal covariance matrix
  given all parallel chains. The sockets are opened/closed in each
  process with a small random sleep time to avoid collisions during
  connections to the internal server of
  \code{LaplacesDemon.hpc}. Blocking sockets are used to synchronize
  processes.
}
\author{Silvere Vialet-Chabrand \email{silvere@vialet-chabrand.com}}
\seealso{
  \code{\link{LaplacesDemon}} and
  \code{\link{LaplacesDemon.hpc}}.
}
\keyword{High Performance Computing}
\keyword{Parallel Chains}