\name{Elicitation}
\alias{delicit}
\alias{elicit}
\title{Prior Elicitation}
\description{
  Prior elicitation is the act of inducing personal opinion to be
  expressed by the probabilities the person associates with an event
  (Savage, 1971). The \code{elicit} function elicits personal opinion
  and the \code{delicit} function estimates probability density to be
  used with model specification in the
  \code{\link{IterativeQuadrature}}, \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}}, \code{\link{LaplacesDemon.hpc}},
  \code{\link{PMC}}, or \code{\link{VariationalBayes}} functions.
}
\usage{
delicit(theta, x, a=-Inf, b=Inf, log=FALSE)
elicit(n, cats, cat.names, show.plot=FALSE)
}
\arguments{
  \item{theta}{This is a scalar or vector of parameters for which the
    density is estimated with respect to the kernel density estimate of
    \code{x}.}
  \item{x}{This is the elicited vector.}
  \item{a}{This is an optional lower bound for support.}
  \item{b}{This is an optional upper bound for support.}
  \item{log}{Logical. If \code{log=TRUE}, then the logarithm of the
    density is returned.}
  \item{n}{This is the number of chips.}
  \item{cats}{This is a vector of \eqn{k} categories, bins, or
    intervals. When the variable is continuous, the mid-point of each
    category is used. For example, if the continuous interval [0,1] has
    5 equal-sized categories, then \code{cats=c(0.1,0.3,0.5,0.7,0.9)}.}
  \item{cat.names}{This is a vector of category names. For example, if
    the continuous interval [0,1] has 5 equal-sized categories, then one
    way or naming the categories may be \code{cat.names=c("0:<.2",
      ".2:<.4", ".4:<.6", ".6:<.8", ".8:1")}.}
  \item{show.plot}{Logical. If \code{show.plot=TRUE}, then a barplot is
    shown after each allocation of chips.}
}
\details{
  The \code{elicit} function elicits a univariate, discrete,
  non-conjugate, informative, prior probability distribution by
  offering a number of chips (specified as \code{n} by the statistician)
  for the user to allocate into categories specified by the
  statistician. The results of multiple elicitations (meaning, with
  multiple people), each the output of \code{elicit}, may be combined
  with the \code{c} function in base R.

  This discrete distribution is included with the data for
  a model and supplied to a model specification function, where in turn
  it is supplied to the \code{delicit} function, which estimates the
  density at the current value of the prior distribution,
  \eqn{p(\theta)}{p(theta)}. The prior distribution may be either
  continuous or discrete, will be proper, and may have bounded support
  (constrained to an interval).

  For a minimal example, a statistician elicits the prior probability
  distribution for a regression effect, \eqn{\beta}. Non-statisticians
  would not be asked about expected parameters, but could be asked about
  how much \eqn{\textbf{y}}{y} would be expected to change given a
  one-unit change in \eqn{\textbf{x}}{x}. After consulting with others
  who have prior knowledge, the support does not need to be bounded,
  and their guesses at the range result in the statistician creating
  5 catgories from the interval [-1,4], where each interval has a width
  of one. The statistician schedules time with 3 people, and each person
  participates when the statistician runs the following R code:

  \code{x <- elicit(n=10, cats=c(-0.5, 0.5, 1.5, 2.5, 3.5),
    cat.names=c("-1:<0", "0:<1", "1:<2", "2:<3", "3:4"), show.plot=TRUE)}

  Each of the 3 participants receives 10 chips to allocate among the
  5 categories according to personal beliefs in the probability of the
  regression effect. When the statistician and each participant accept
  their elicited distribution, all 3 vectors are combined into one
  vector. In the model form, the prior is expressed as

  \deqn{p(\beta) \sim \mathcal{EL}}{p(beta) ~ EL}

  and the code for the model specification is

  \code{elicit.prior <- delicit(beta, x, log=TRUE)}

  This method is easily extended to priors that are multivariate,
  correlated, or conditional.

  As an alternative, Hahn (2006) also used a categorical approach,
  eliciting judgements about the relative likelihood of each category,
  and then minimizes the KLD (for more information on KLD, see the
  \code{\link{KLD}} function).
}
\references{
  Hahn, E.D. (2006). "Re-examining Informative Prior Elicitation Through
  the Lens of Markov chain Monte Carlo Methods". \emph{Journal of the
  Royal Statistical Society}, A 169 (1), p. 37--48.
  
  Savage, L.J. (1971). "Elicitation of Personal Probabilities and
  Expectations". \emph{Journal of the American Statistical Association},
  66(336), p. 783--801.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{de.Finetti.Game}},
  \code{\link{KLD}},
  \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{LaplacesDemon.hpc}},
  \code{\link{PMC}}, and
  \code{\link{VariationalBayes}}.
}
\examples{
library(LaplacesDemon)
x <- c(1,2,2,3,3,3,4,7,8,8,9,10) #Elicited with elicit function
theta <- seq(from=-5,to=15,by=.1)
plot(theta, delicit(theta,x), type="l", xlab=expression(theta),
     ylab=expression("p(" * theta * ")"))
}
\keyword{Elicitation}
\keyword{Personal Probability}
\keyword{Prior}
\keyword{Subjective Probability}
\keyword{Utility}