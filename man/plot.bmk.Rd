\name{plot.bmk}
\alias{plot.bmk}
\title{Plot Hellinger Distances}
\description{
  This function plots Hellinger distances in an object of class \code{bmk}.
}
\usage{
\method{plot}{bmk}(x, col=colorRampPalette(c("black","red"))(100),
     title="", PDF=FALSE, Parms=NULL, \dots)
}
\arguments{
  \item{x}{This required argument is an object of class \code{bmk}. See
    the \code{\link{BMK.Diagnostic}} function for more information.}
  \item{col}{This argument specifies the colors of the cells. By
    default, the \code{colorRampPalette} function colors large Hellinger
    distances as \code{red}, small as \code{black}, and provides 100
    color gradations.}
  \item{title}{This argument specifies the title of the plot, and the
    default does not include a title.}
  \item{PDF}{Logical. When \code{TRUE}, the plot is saved as a .pdf
    file.}
  \item{Parms}{
    This argument accepts a vector of quoted strings to be matched for
    selecting parameters for plotting. This argument defaults to
    \code{NULL} and selects every parameter for plotting. Each quoted
    string is matched to one or more parameter names with the
    \code{grep} function. For example, if the user specifies
    \code{Parms=c("eta", "tau")}, and if the parameter names
    are beta[1], beta[2], eta[1], eta[2], and tau, then all parameters
    will be selected, because the string \code{eta} is within
    \code{beta}. Since \code{grep} is used, string matching uses
    regular expressions, so beware of meta-characters, though these are
    acceptable: ".", "[", and "]".}
  \item{\dots}{Additional arguments are unused.}
}
\details{
  The \code{plot.bmk} function plots the Hellinger distances in an
  object of class \code{bmk}. This is useful for quickly finding
  portions of chains with large Hellinger distances, which indicates
  non-stationarity and non-convergence.
}
\seealso{\code{\link{BMK.Diagnostic}}}
\examples{
library(LaplacesDemon)
N <- 1000 #Number of posterior samples
J <- 10 #Number of parameters
Theta <- matrix(runif(N*J),N,J)
colnames(Theta) <- paste("beta[", 1:J, "]", sep="")
for (i in 2:N) {Theta[i,1] <- Theta[i-1,1] + rnorm(1)}
HD <- BMK.Diagnostic(Theta, batches=10)
plot(HD, title="Hellinger distance between batches")
}
\keyword{Plot}