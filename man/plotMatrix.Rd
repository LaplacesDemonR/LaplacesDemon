\name{plotMatrix}
\alias{plotMatrix}
\title{Plot a Numerical Matrix}
\description{
  This function plots a numerical matrix, and is often used to plot the
  following matrices: correlation, covariance, distance, and
  precision.
}
\usage{
plotMatrix(x, col=colorRampPalette(c("red","black","green"))(100),
     cex=1, circle=TRUE, order=FALSE, zlim=NULL, title="", PDF=FALSE, \dots)
}
\arguments{
  \item{x}{This required argument is a numerical matrix, or an
    object of class \code{bayesfactor}, \code{demonoid}, \code{iterquad},
    \code{laplace}, \code{pmc}, \code{posteriorchecks}, or \code{vb}.
    See more information below regarding these classes. One component of
    a blocked proposal covariance matrix must be pointed to explicitly,
    rather than to the object of class \code{demonoid}.}
  \item{col}{This argument specifies the colors of the circles. By
    default, the \code{colorRampPalette} function colors strong positive
    correlation as \code{green}, zero correlation as \code{black}, and
    strong negative correlation as \code{red}, and provides 100 color
    gradations.}
  \item{cex}{When \code{circle=TRUE}, this argument specifies the size
    of the marginal text, the names of the parameters or variables, and
    defaults to 1.}
  \item{circle}{Logical. When \code{TRUE}, each element in the numeric
    matrix is represented with a circle, and a larger circle is assigned
    to elements that are farther from zero. Also, when \code{TRUE}, the
    gradation scale does not appear to the right of the plot.}
  \item{order}{Logical. This argument defaults to \code{FALSE}, and
    presents the parameters or variables in the same order as in the
    numeric matrix. When \code{TRUE}, the parameters or variables
    are ordered using principal components analysis.}
  \item{zlim}{When \code{circle=FALSE}, the gradation scale may be
    constrained to an interval by \code{zlim}, such as
    \code{zlim=c(-1,1)}, and only values within the interval are
    plotted.}
  \item{title}{This argument specifies the title of the plot, and the
    default does not include a title. When \code{x} is of class
    \code{posteriorchecks}, the title is changed to \code{Posterior
      Correlation}.}
  \item{PDF}{Logical. When \code{TRUE}, the plot is saved as a .pdf
    file.}
  \item{\dots}{Additional arguments are unused.}
}
\details{
  The \code{plotMatrix} function produces one of two styles of plots,
  depending on the \code{circle} argument. A \eqn{K \times K}{K x K}
  numeric matrix of \eqn{K} parameters or variables is plotted. The plot
  is a matrix of the same dimensions, in which each element is colored
  (and sized, when \code{circle=TRUE}) according to its value.

  Although \code{plotMatrix} does not provide the same detail as a
  numeric matrix, it is easier to discover elements of interest
  according to color (and size when \code{circle=TRUE}).
  
  The \code{plotMatrix} function is not inherently Bayesian, and does
  not include uncertainty in matrices. Nonetheless, it is included
  because it is a useful graphical presentation of a numeric matrices,
  and is recommended to be used with the posterior correlation matrix in
  an object of class \code{posteriorchecks}.

  When \code{x} is an object of class \code{bayesfactor}, matrix
  \code{B} is plotted. When \code{x} is an object of class
  \code{demonoid} (if it is a matrix), \code{iterquad}, \code{laplace},
  \code{pmc}, or \code{vb}, the covariance matrix \code{Covar} is
  plotted. When \code{x} is an object of class \code{posteriorchecks},
  the posterior correlation matrix is plotted.

  This is a modified version of the \code{circle.corr} function
  of Taiyun Wei.
}
\author{Taiyun Wei}
\seealso{\code{\link{PosteriorChecks}}}
\examples{
library(LaplacesDemon)
### Although it is most commonly used with an object of class
### posteriorchecks, it is applied here to a different correlation matrix.
data(mtcars)
plotMatrix(cor(mtcars), col=colorRampPalette(c("green","gray10","red"))(100),
     cex=1, circle=FALSE, order=TRUE)
plotMatrix(cor(mtcars), col=colorRampPalette(c("green","gray10","red"))(100),
     cex=1, circle=TRUE, order=TRUE)}
\keyword{Plot}