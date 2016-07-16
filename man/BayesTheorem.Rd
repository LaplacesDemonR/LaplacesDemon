\name{BayesTheorem}
\alias{BayesTheorem}
\title{Bayes' Theorem}
\description{
  Bayes' theorem shows the relation between two conditional
  probabilities that are the reverse of each other. This theorem is
  named after Reverend Thomas Bayes (1702-1761), and is also referred to
  as Bayes' law or Bayes' rule (Bayes and Price, 1763). Bayes' theorem
  expresses the conditional probability, or `posterior probability', of
  an event \eqn{A} after \eqn{B} is observed in terms of the `prior
  probability' of \eqn{A}, prior probability of \eqn{B}, and the
  conditional probability of \eqn{B} given \eqn{A}. Bayes' theorem is
  valid in all common interpretations of probability. This function
  provides one of several forms of calculations that are possible with
  Bayes' theorem.
}
\usage{
BayesTheorem(PrA, PrBA)
}
\arguments{
  \item{PrA}{This required argument is the prior probability of \eqn{A},
    or \eqn{\Pr(A)}{Pr(A)}.}
  \item{PrBA}{This required argument is the conditional probability of
    \eqn{B} given \eqn{A} or \eqn{\Pr(B | A)}{Pr(B | A)}, and is known
    as the data, evidence, or likelihood.}
}
\value{
  The \code{BayesTheorem} function returns the conditional probability
  of \eqn{A} given \eqn{B}, known in Bayesian inference as the
  posterior. The returned object is of class \code{bayestheorem}.
}
\details{
  Bayes' theorem provides an expression for the conditional probability
  of \eqn{A} given \eqn{B}, which is equal to

  \deqn{\Pr(A | B) = \frac{\Pr(B | A)\Pr(A)}{\Pr(B)}}{Pr(A | B) = (Pr(B
    | A)Pr(A)) / Pr(B)}

  For example, suppose one asks the question: what is the probability of
  going to Hell, conditional on consorting (or given that a person
  consorts) with Laplace's Demon. By replacing \eqn{A} with \eqn{Hell}
  and \eqn{B} with \eqn{Consort}, the question becomes

  \deqn{\Pr(\mathrm{Hell} | \mathrm{Consort}) =
  \frac{\Pr(\mathrm{Consort} |
  \mathrm{Hell})\Pr(\mathrm{Hell})}{\Pr(\mathrm{Consort})}}{Pr(Hell |
  Consort) = (Pr(Consort | Hell)Pr(Hell)) / Pr(Consort)}

  Note that a common fallacy is to assume that \eqn{\Pr(A | B) = \Pr(B |
    A)}{Pr(A | B) = Pr(B | A)}, which is called the conditional
  probability fallacy.

  Another way to state Bayes' theorem (and this is the form in the
  provided function) is

  \deqn{\Pr(A_i | B) = \frac{\Pr(B | A_i)\Pr(A_i)}{\Pr(B | A_i)\Pr(A_i)
  +\dots+ \Pr(B | A_n)\Pr(A_n)}}{Pr(A[i] | B) = (Pr(B | A[i])Pr(A[i])) /
  (Pr(B | A[i])Pr(A[i]) +...+ Pr(B | A[n])Pr(A[n]))}

  Let's examine our \emph{burning} question, by replacing
  \eqn{A_i}{A[i]} with Hell or Heaven, and replacing \eqn{B} with
  Consort

  \itemize{
    \item \eqn{\Pr(A_1) = \Pr(\mathrm{Hell})}{Pr(A[1] = Pr(Hell)}
    \item \eqn{\Pr(A_2) = \Pr(\mathrm{Heaven})}{Pr(A[2] = Pr(Heaven)}
    \item \eqn{\Pr(B) = \Pr(\mathrm{Consort})}{Pr(B) = Pr(Consort)}
    \item \eqn{\Pr(A_1 | B) = \Pr(\mathrm{Hell} |
    \mathrm{Consort})}{Pr(A[1] | B) = Pr(Hell | Consort)}
    \item \eqn{\Pr(A_2 | B) = \Pr(\mathrm{Heaven} |
    \mathrm{Consort})}{Pr(A[2] | B) = Pr(Heaven | Consort)}
    \item \eqn{\Pr(B | A_1) = \Pr(\mathrm{Consort} |
    \mathrm{Hell})}{Pr(B | A[1]) = Pr(Consort | Heaven)}
    \item \eqn{\Pr(B | A_2) = \Pr(\mathrm{Consort} |
    \mathrm{Heaven})}{Pr(B | A[2]) = Pr(Consort | Heaven)}
}

  Laplace's Demon was conjured and asked for some data. He was glad to
  oblige.

  \itemize{
    \item 6 people consorted out of 9 who went to Hell.
    \item 5 people consorted out of 7 who went to Heaven.
    \item 75\% of the population goes to Hell.
    \item 25\% of the population goes to Heaven.
  }

  Now, Bayes' theorem is applied to the data. Four pieces are worked out
  as follows

  \itemize{
    \item \eqn{\Pr(\mathrm{Consort} | \mathrm{Hell}) = 6/9 =
    0.666}{Pr(Consort | Hell) = 6/9 = 0.666}
    \item \eqn{\Pr(\mathrm{Consort} | \mathrm{Heaven}) = 5/7 =
    0.714}{Pr(Consort | Heaven) = 5/7 = 0.714}
    \item \eqn{\Pr(\mathrm{Hell}) = 0.75}{Pr(Hell) = 0.75}
    \item \eqn{\Pr(\mathrm{Heaven}) = 0.25}{Pr(Heaven) = 0.25}
  }

  Finally, the desired conditional probability \eqn{\Pr(\mathrm{Hell} |
    \mathrm{Consort})}{Pr(Hell | Consort)} is calculated using Bayes'
  theorem

  \itemize{
    \item \eqn{\Pr(\mathrm{Hell} | \mathrm{Consort}) =
    \frac{0.666(0.75)}{0.666(0.75) + 0.714(0.25)}}{Pr(Hell | Consort) =
    0.666(0.75) / (0.666(0.75) + 0.714(0.25))}
    \item \eqn{\Pr(\mathrm{Hell} | \mathrm{Consort}) = 0.737}{Pr(Hell |
    Consort) = 0.737}
}
    The probability of someone consorting with Laplace's Demon and going
    to Hell is 73.7\%, which is less than the prevalence of 75\% in the
    population. According to these findings, consorting with Laplace's
    Demon does not increase the probability of going to Hell.

    For an introduction to model-based Bayesian inference, see the
    accompanying vignette entitled ``Bayesian Inference'' or
    \url{https://web.archive.org/web/20150206004608/http://www.bayesian-inference.com/bayesian}.
}
\references{
  Bayes, T. and Price, R. (1763). "An Essay Towards Solving a Problem in
  the Doctrine of Chances". By the late Rev. Mr. Bayes, communicated by
  Mr. Price, in a letter to John Canton, M.A. and F.R.S.
  \emph{Philosophical Transactions of the Royal Statistical Society of
  London}, 53, p. 370--418.
}
\author{Statisticat, LLC.}
\seealso{
  \code{\link{IterativeQuadrature}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{PMC}}, and
  \code{\link{VariationalBayes}}.
}
\examples{
# Pr(Hell|Consort) =
PrA <- c(0.75,0.25)
PrBA <- c(6/9, 5/7)
BayesTheorem(PrA, PrBA)
}
\keyword{Bayesian}
\keyword{Probability}

