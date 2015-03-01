\name{data.demonchoice}
\alias{demonchoice}
\title{Demon Choice Data Set}
\usage{data(demonchoice)}
\description{
  This data set is for discrete choice models and consists of the choice
  of commuting route to school: arterial, two-lane, or freeway. There
  were 151 Pennsylvania commuters who started from a residential complex
  in State College, PA, and commute to downtown State College.
}
\format{
  This data frame contains 151 rows of individual choices and 9
  columns. The following data dictionary describes each variable or
  column.
  \describe{
    \item{\code{Choice}}{This is the route choice: four-lane arterial
      (35 MPH speed limit), two-lane highway (35 MPH speed limit, with
      one lane in each direction), or a limited-access four-lane freeway
      (55 MPH speed liimit.)}
    \item{\code{HH.Income}}{This is an ordinal variable of annual
      household income of the commuter in USD. There are four
      categories: 1 is less than 20,000 USD, 2 is 20,000-29,999 USD, 3
      is 30,000-39,999 USD, and 4 is 40,000 USD or greater.}
    \item{\code{Vehicle.Age}}{This is the age in years of the vehicle of
      the commuter.}
    \item{\code{Stop.Signs.Arterial}}{This is the number of stop signs
      along the arterial route.}
    \item{\code{Stop.Signs.Two.Lane}}{This is the number of stop signs
      along the two-lane route.}
    \item{\code{Stop.Signs.Freeway}}{This is the number of stop signs
      along the freeway route.}
    \item{\code{Distance.Arterial}}{This is distance in miles of the
      arterial route.}
    \item{\code{Distance.Two.Lane}}{This is the distance in miles of
      the two-lane route.}
    \item{\code{Distance.Freeway}}{This is the distance in miles of
      the freeway route.}
    }
}
\source{Washington, S., Congdon, P., Karlaftis, M., and Mannering,
  F. (2009). "Bayesian Multinomial Logit: Theory and Route Choice
  Example". Transportation Research Record, 2136, p. 28--36.}
\keyword{datasets}
