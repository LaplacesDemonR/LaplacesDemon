\name{de.Finetti.Game}
\alias{de.Finetti.Game}
\title{de Finetti's Game}
\description{
  The \code{de.Finetti.Game} function estimates the interval of a
  subjective probability regarding a possible event in the near future.
}
\usage{
de.Finetti.Game(width)
}
\arguments{
  \item{width}{This is the maximum acceptable width of the interval for
    the returned subjective probability. The user must specify a width
    between 0 and 1.}
}
\details{
  This function is a variation on the game introduced by de Finetti,
  who is one of the main developers of subjective probability,
  along with Ramsey and Savage. In the original context, de Finetti
  proposed a gamble regarding life on Mars one billion years ago.

  The frequentist interpretation of probability defines the probability
  of an event as the limit of its relative frequency in a large number
  of trials. Frequentist inference is undefined, for example, when there
  are no trials from which to calculate a probability. By defining
  probability relative to frequencies of physical events, frequentists
  attempt to objectify probability. However, de Finetti asserts that the
  frequentist (or objective) interpretation always reduces to a subjective
  interpretation of probability, because probability is a human
  construct and does not exist independently of humans in
  nature. Therefore, probability is a degree of belief, and is called
  subjective or personal probability.
}
\value{
  The \code{de.Finetti.Game} function returns a vector of length
  two. The respective elements are the lower and upper bounds of the
  subjective probability of the participant regarding the possible event
  in the near future.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{elicit}}
}
\keyword{Personal Probability}
\keyword{Subjective Probability}
\keyword{Utility}
