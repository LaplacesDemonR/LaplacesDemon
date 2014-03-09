\name{LossMatrix}
\alias{LossMatrix}
\title{Loss Matrix}
\description{
  A loss matrix is useful in Bayesian decision theory for selecting the
  Bayes action, the optimal Bayesian decision, when there are a discrete
  set of possible choices (actions) and a discrete set of possible
  outcomes (states of the world). The Bayes action is the action that
  minimizes expected loss, which is equivalent to maximizing expected
  utility.
}
\usage{
LossMatrix(L, p.theta)
}
\arguments{
  \item{L}{This required argument accepts a \eqn{S \times A}{S x A}
    matrix or \eqn{S \times A \times N}{S x A x N} array of losses,
    where \eqn{S} is the number of states of the world, \eqn{A} is the
    number of actions, and \eqn{N} is the number of samples. These
    losses have already been estimated, given a personal loss
    function. One or more personal losses has already been estimated
    for each combination of possible actions
    \eqn{a=1,\dots,A}{a=1,...,A} and possible states
    \eqn{s=1,\dots,S}{s=1,...,S}.}
  \item{p.theta}{This required argument accepts a
    \eqn{S \times A}{S x A} matrix or
    \eqn{S \times A \times N}{S x A x N} array of state prior
    probabilities, where \eqn{S} is the number of states of the world,
    \eqn{A} is the number of actions, and \eqn{N} is the number of
    samples. The sum of each column must equal one.}
}
\details{
  Bayesian inference is often tied to decision theory (Bernardo and
  Smith, 2000), and decision theory has long been considered the
  foundations of statistics (Savage, 1954).

  Before using the \code{LossMatrix} function, the user should have
  already considered all possible actions (choices), states of the world
  (outcomes unknown at the time of decision-making), chosen a loss
  function \eqn{L(\theta, \alpha)}{L(theta, alpha)}, estimated loss, and
  elicited prior probabilities \eqn{p(\theta | x)}{p(theta | x)}.
  
  Although possible actions (choices) for the decision-maker and
  possible states (outcomes) may be continuous or discrete, the loss
  matrix is used for discrete actions and states. An example of a
  continuous action may be that a decision-maker has already decided to
  invest, and the remaining, current decision is how much to invest.
  Likewise, an example of continuous states of the world (outcomes) may
  be how much profit or loss may occur after a given continuous unit of
  time.

  The coded example provided below is taken from Berger (1985, p. 6-7)
  and described here. The set of possible actions for a decision-maker
  is to invest in bond ZZZ or alternatively in bond XXX, as it is called
  here. A real-world decision should include a mutually exhaustive list
  of actions, such as investing in neither, but perhaps the
  decision-maker has already decided to invest and narrowed the options
  down to these two bonds.

  The possible states of the world (outcomes unknown at the time of
  decision-making) are considered to be two states: either the chosen
  bond will not default or it will default. Here, the loss function is
  a negative linear identity of money, and hence a loss in element
  \code{L[1,1]} of -500 is a profit of 500, while a loss in
  \code{L[2,1]} of 1,000 is a loss of 1,000.

  The decision-maker's dilemma is that bond ZZZ may return a higher
  profit than bond XXX, however there is an estimated 10\% chance, the
  prior probability, that bond ZZZ will default and return a substantial
  loss. In contrast, bond XXX is considered to be a sure-thing and 
  return a steady but smaller profit. The Bayes action is to choose the
  first action and invest in bond ZZZ, because it minimizes expected
  loss, even though there is a chance of default.

  A more realistic application of a loss matrix may be to replace the
  point-estimates of loss with samples given uncertainty around the
  estimated loss, and replace the point-estimates of the prior
  probability of each state with samples given the uncertainty of the
  probability of each state. The loss function used in the example is
  intuitive, but a more popular monetary loss function may be
  \eqn{-\log(E(W | R))}{-log(E(W | R))}, the negative log of the
  expectation of wealth, given the return. There are many alternative
  loss functions.

  Although isolated decision-theoretic problems exist such as the
  provided example, decision theory may also be applied to the results
  of a probability model (such as from
  \code{\link{IterativeQuadrature}}, \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}}, \code{\link{PMC}}), or
  \code{\link{VariationalBayes}}, contingent on how
  a decision-maker is considering to use the information from the
  model. The statistician may pass the results of a model to a client,
  who then considers choosing possible actions, given this
  information. The statistician should further assist the client with
  considering actions, states of the world, then loss functions, and
  finally eliciting the client's prior probabilities (such as with the
  \code{\link{elicit}} function).

  When the outcome is finally observed, the information from this
  outcome may be used to refine the priors of the next such decision. In
  this way, Bayesian learning occurs.
}
\value{
  The \code{LossMatrix} function returns a list with two components:
  \item{BayesAction}{This is a numeric scalar that indicates the action
    that minimizes expected loss.}
  \item{E.Loss}{This is a vector of expected losses, one for each
    action.}
}
\references{
  Berger, J.O. (1985). "Statistical Decision Theory and Bayesian
  Analysis", Second Edition. Springer: New York, NY.

  Bernardo, J.M. and Smith, A.F.M. (2000). "Bayesian Theory". John
  Wiley \& Sons: West Sussex, England.

  Savage, L.J. (1954). "The Foundations of Statistics". John Wiley \&
  Sons: West Sussex, England.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{elicit}},
  \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{PMC}}, and
  \code{\link{VariationalBayes}}.
}
\examples{
library(LaplacesDemon)
### Point-estimated loss and state probabilities
L <- matrix(c(-500,1000,-300,-300), 2, 2)
rownames(L) <- c("s[1]: !Defaults","s[2]: Defaults")
colnames(L) <- c("a[1]: Buy ZZZ", "a[2]: Buy XXX")
L
p.theta <- matrix(c(0.9, 0.1, 1, 0), 2, 2)
Fit <- LossMatrix(L, p.theta)

### Point-estimated loss and samples of state probabilities
L <- matrix(c(-500,1000,-300,-300), 2, 2)
rownames(L) <- c("s[1]: Defaults","s[2]: !Defaults")
colnames(L) <- c("a[1]: Buy ZZZ", "a[2]: Buy XXX")
L
p.theta <- array(runif(4000), dim=c(2,2,1000)) #Random probabilities,
#just for a quick example. And, since they must sum to one:
for (i in 1:1000) {
     p.theta[,,i] <- p.theta[,,i] / matrix(colSums(p.theta[,,i]),
          dim(p.theta)[1], dim(p.theta)[2], byrow=TRUE)}
Fit <- LossMatrix(L, p.theta)
Fit

### Point-estimates of loss may be replaced with samples as well.
}
\keyword{Decision Theory}
\keyword{Utility}
