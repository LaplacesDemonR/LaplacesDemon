\name{dist.Skew.Discrete.Laplace}
\alias{dsdlaplace}
\alias{psdlaplace}
\alias{qsdlaplace}
\alias{rsdlaplace}
\title{Skew Discrete Laplace Distribution: Univariate}
\description{
  These functions provide the density, distribution function, quantile
  function, and random generation for the univariate, skew discrete
  Laplace distribution with parameters \eqn{p}{p} and \eqn{q}{q}.
}
\usage{
dsdlaplace(x, p, q, log=FALSE)
psdlaplace(x, p, q)
qsdlaplace(prob, p, q)
rsdlaplace(n, p, q)
}
\arguments{
  \item{x}{This is a vector of data.}
  \item{p}{This is a scalar or vector of parameter \eqn{p \in [0,1]}{p
      in [0,1]}.}
  \item{q}{This is a scalar or vector of parameter \eqn{q \in [0,1]}{q
      in [0,1]}.}
  \item{prob}{This is a probability scalar or vector.}
  \item{n}{This is the number of observations, which must be a positive
    integer that has length 1.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
    \item Application: Discrete Univariate
    \item Density 1: \eqn{p(\theta) = \frac{(1-p)(1-q)}{1-pq}p^\theta;
      \theta=0,1,2,3,\dots}{p(theta) = (1-p)(1-q) / (1-pq)p^theta,
      theta=0,1,2,3,...}
    \item Density 2: \eqn{p(\theta) = \frac{(1-p)(1-q)}{1-pq}q^{|\theta|};
      x=0,-1,-2,-3,\dots}{p(theta) = (1-p)(1-q) / (1-pq)q^(|theta|),;
      x=0,-1,-2,-3,...}
    \item Inventor: Kozubowski, T.J. and Inusah, S. (2006)
    \item Notation 1: \eqn{\theta \sim \mathcal{DL}(p, q)}{theta ~ DL(p,
      q)}
    \item Notation 2: \eqn{p(\theta) = \mathcal{DL}(\theta | p,
      q)}{p(theta) = DL(theta | p, q)}
    \item Parameter 1: \eqn{p \in [0,1]}{p in [0,1]}
    \item Parameter 2: \eqn{q \in [0,1]}{q in [0,1]}
    \item Mean 1: \eqn{E(\theta) =
      \frac{1}{1-p}-\frac{1}{1-q}=\frac{p}{1-p}-\frac{q}{1-q}}{E(theta)
      = (1 / (1-p)) - (1 / (1-q)) = (p / (1-p)) - (q / (1-q))}
    \item Mean 2: \eqn{E(|\theta|) =
      \frac{q(1-p)^2+p(1-q)^2}{(1-qp)(1-q)(1-p)}}{E(|theta|) =
      (q(1-p)^2+p(1-q)^2) / ((1-qp)(1-q)(1-p))}
    \item Variance: \eqn{var(\theta) =
      \frac{1}{(1-p)^2(1-q)^2}[\frac{q(1-p)^3(1+q)+p(1-q)^3(1+p)}{1-pq}-(p-q)^2]}{var(theta)
      = (1 / ((1-p)^2(1-q)^2))[(q(1-p)^3(1+q)+p(1-q)^3(1+p)) / (1-pq) -
      (p-q)^2]}
    \item Mode: 
}

This is a discrete form of the skew-Laplace distribution. The symmetric
discrete Laplace distribution occurs when \eqn{p=q}{p=q}. DL(p,0) is a
geometric distribution, and DL(0,q) is a geometric distribution of
non-positive integers. The distribution is degenerate when
DL(0,0). Since the geometric distribution is a discrete analog of the
exponential distribution, the distribution of the difference of two
geometric variables is a discrete Laplace distribution.


These functions are similar to those in the \code{DiscreteLaplace} package.
}
\value{
  \code{dslaplace} gives the density,
  \code{pslaplace} gives the distribution function,
  \code{qslaplace} gives the quantile function, and
  \code{rslaplace} generates random deviates.
}
\references{
  Kozubowski, T.J. and Inusah, S. (2006). "A Skew Laplace Distribution
  on Integers". \emph{AISM}, 58, p. 555--571.
}
\seealso{
  \code{\link{dalaplace}},
  \code{\link{dexp}},
  \code{\link{dlaplace}},
  \code{\link{dlaplacep}}, and
  \code{\link{dslaplace}}.
}
\examples{
library(LaplacesDemon)
x <- dsdlaplace(1,0.5,0.5)
x <- psdlaplace(1,0.5,0.5)
x <- qsdlaplace(0.5,0.5,0.5)
x <- rsdlaplace(5,0.5,0.5)

#Plot Probability Functions
     x <- c(-3:3)
     plot(x, dsdlaplace(x,0.5,0.5), ylim=c(0,0.6), type="l", main="Probability Function",
          ylab="density", col="red")
     lines(x, dsdlaplace(x,0.3,0.6), type="l", col="green")
     lines(x, dsdlaplace(x,0.9,0.1), type="l", col="blue")
     legend(-2.5, 0.5, expression(paste(p==0.5, ", ", q==0.5),
          paste(p==0.3, ", ", q==0.6),
          paste(p==0.9, ", ", q==0.1)),
          lty=c(1,1,1), col=c("red","green","blue"))
}
\keyword{Distribution}