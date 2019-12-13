\name{dist.Zellner}
\alias{dhyperg}
\alias{dzellner}
\alias{rzellner}
\title{Hyperprior-g Prior and Zellner's g-Prior}
\description{
  These functions provide the density of the hyper-g prior (Liang et
  al., 2008), and both the density and random generation of Zellner's
  g-prior (Zellner, 1986).
}
\usage{
dhyperg(g, alpha=3, log=FALSE)
dzellner(beta, g, sigma, X, log=FALSE)
rzellner(n, g, sigma, X)
}
\arguments{
  \item{alpha}{This is a positive scale hyperhyperparameter that is
    proper when \eqn{\alpha > 2}{alpha > 2}. The default is
    \code{alpha=3}.}
  \item{beta}{This is regression effects \eqn{\beta}{beta}, a vector of
    length \eqn{J}.}
  \item{g}{This is hyperparameter \eqn{g}, a positive scalar.}
  \item{n}{This is the number of random deviates to generate.}
  \item{sigma}{This is the residual standard deviation
    \eqn{\sigma}{sigma}, a positive scalar.}
  \item{X}{This is a full-rank \eqn{N \times J}{N x J} design matrix
    \eqn{\textbf{X}}{X} for \eqn{N} records and \eqn{J} predictors,
    where \eqn{J+1 < N}. Zellner's g-prior has been extended (elsewhere)
    via singular value decomposition (SVD) to the case where
    \eqn{J > N}.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
}
\details{
  \itemize{
  \item Application: Continuous Multivariate
  \item Density: \eqn{p(\theta) = \frac{1}{(2\pi)^{J/2}|(g
      \sigma^2(\textbf{X}^T \textbf{X})^{-1})^{-1}|^{1/2}}
    \exp(-\frac{1}{2}(\theta - \mu)'(g \sigma^2(\textbf{X}^T
    \textbf{X})^{-1})^{-1}(\theta - \mu))}{p(theta) =
    (1/((2*pi)^(J/2)*|Sigma|^(1/2))) * 
    exp(-(1/2)*(theta-mu)'*Sigma^(-1)*(theta-mu))}
  \item Inventor: Zellner, A. (1986)
  \item Notation 1: \eqn{\theta \sim
    \mathrm{N}_J(0, g \sigma^2(\textbf{X}^T \textbf{X})^{-1})}{theta ~
    N[J](0, g sigma^2 (X^T X)^(-1))}
  \item Notation 2: \eqn{p(\theta) = \mathrm{N}_J(\theta | g, \sigma^2,
    \textbf{X})}{p(theta) = N[J](theta, g, sigma^2, X)}
  \item Parameter 1: location parameter \eqn{\beta}{beta}
  \item Parameter 2: scale parameter \eqn{g > 0}
  \item Parameter 3: scale parameter \eqn{\sigma^2 > 0}{sigma^2 > 0}
  \item Mean:
  \item Variance:
  \item Mode:
}
Zellner's g-prior is a popular, data-dependent, elliptical, improper,
least-informative prior distribution on regression effects
\eqn{\beta}{beta} in a Gaussian regression model. It is a particular
form in the conjugate Normal-Gamma family. Zellner's g-prior is also
used for estimating Bayes factors (for hypothesis testing) with a
simpler form, as well as in model selection and variable selection. The
marginal posterior distribution of regression effects \eqn{\beta}{beta}
is multivariate t.

One of many nice properties of Zellner's g-prior is that it adapts
automatically to near-collinearity between different
predictors. Zellner's g-prior puts most of its prior mass in the
direction that causes the regression coefficients of correlated
predictors to be smoothed away from each other. When coupled with model
selection, Zellner's g-prior discourages highly collinear predictors
from entering the models simultaneously by inducing a negative
correlation between the coefficients. However, when it is desirable for
collinear predictors to enter simultaneously, a modification has been
proposed (though not included here) in which
\eqn{(\textbf{X}^T \textbf{X})^{-1}}{(X^T X)^(-1)} is replaced with
\eqn{(\textbf{X}^T \textbf{X})^\lambda}{(X^T X)^lambda}. For more
information, see Krishna et al. (2009).

For variable selection, large values of \eqn{g}, with a prior mean of
zero for \eqn{\beta}{beta}, encourage models with few, large
coefficients. Conversely, small values of \eqn{g} encourage saturated
models with many, small coefficients.

The design matrix \eqn{\textbf{X}}{X} is converted to Fisher's
information matrix, which is used as a covariance matrix for
\eqn{\beta}{beta}. This is computationally efficient, because each
element of the covariance matrix does not need to be estimated as a
parameter. When \eqn{\textbf{X}}{X} is nearly singular, regression
effects \eqn{\beta}{beta} may be poorly estimated.

Hyperparameter \eqn{g} acts as an inverse relative prior sample size, or
as a dimensionality penalty. Zellner (1986) recommended that a
hyperprior distribution is assigned to \eqn{g} so that it is estimated
from the data, although in practice \eqn{g} has often been fixed, usually
to \eqn{N} when no information is available, since it has the
interpretation of adding prior information equivalent to one
observation. A variety of hyperpriors have been suggested for \eqn{g},
such as in Bove and Held (2011), Liang et al. (2008), and Maruyama and
George (2011). \eqn{g} becomes diffuse as it approaches infinity, and
the Bayes factor approaches zero. The hyper-g prior of Liang et al.
(2008) is proper when \eqn{\alpha > 2}{alpha > 2}, and any value in
the interval \eqn{(2,4]} may be reasonable.
}
\value{
  \code{dhyperg} gives the density of the hyper-g prior of Liang et
  al. (2008), \code{dzellner} gives the density of Zellner's g-prior,
  and \code{rzellner} generates random deviates.
}
\references{
  Bove, D.S. and Held, L. (2011). "Hyper-g Priors for Generalized
  Linear Models". \emph{Bayesian Analysis}, 6(3), p. 387--410.

  Krishna, A., Bondell, H.D., and Ghosh, S.K. (2009). "Bayesian Variable
  Selection Using an Adaptive Powered Correlation Prior". \emph{Journal
  of Statistical Planning Inference}, 139(8), p. 2665-2674..

  Liang, F., Paulo, R., Molina, G., Clyde, M.A., and Berger,
  J.O. (2008). "Mixtures of g Priors for Bayesian Variable
  Selection". \emph{Journal of the American Statistical Association},
  103, p. 410--423.

  Maruyama, Y. and George, E.I. (2011). "Fully Bayes Factors with a
  Generalised g-Prior". \emph{Annals of Statistics}, 39, p. 2740--2765.

  Zellner, A. (1986). "On Assessing Prior Distributions and Bayesian
  Regression Analysis with g-Prior Distributions". In \emph{Bayesian
  Inference and Decision Techniques: Essays in Honor of Bruno de
  Finetti}, p. 233--243. Elsevier: Amsterdam, North Holland.
}
\seealso{
  \code{\link{BayesFactor}} and
  \code{\link{dmvt}}
}
\examples{
library(LaplacesDemon)
set.seed(667)
beta <- rnorm(10)
g <- 100
sigma <- 2
X <- cbind(1,matrix(rnorm(100*9),100,9))
dhyperg(g, alpha=3)
dzellner(beta, g, sigma, X)
rzellner(1, g, sigma, X)
}
\keyword{Distribution}
