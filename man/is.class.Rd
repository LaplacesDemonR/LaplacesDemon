\name{is.class}
\alias{is.bayesfactor}
\alias{is.blocks}
\alias{is.bmk}
\alias{is.demonoid}
\alias{is.demonoid.hpc}
\alias{is.demonoid.ppc}
\alias{is.demonoid.val}
\alias{is.hangartner}
\alias{is.heidelberger}
\alias{is.importance}
\alias{is.iterquad}
\alias{is.iterquad.ppc}
\alias{is.juxtapose}
\alias{is.laplace}
\alias{is.laplace.ppc}
\alias{is.miss}
\alias{is.pmc}
\alias{is.pmc.ppc}
\alias{is.pmc.val}
\alias{is.posteriorchecks}
\alias{is.raftery}
\alias{is.rejection}
\alias{is.sensitivity}
\alias{is.vb}
\alias{is.vb.ppc}
\title{Logical Check of Classes}
\description{
  These functions each provide a logical test of the class of an object.
}
\usage{
is.bayesfactor(x)
is.blocks(x)
is.bmk(x)
is.demonoid(x)
is.demonoid.hpc(x)
is.demonoid.ppc(x)
is.demonoid.val(x)
is.hangartner(x)
is.heidelberger(x)
is.importance(x)
is.iterquad(x)
is.iterquad.ppc(x)
is.juxtapose(x)
is.laplace(x)
is.laplace.ppc(x)
is.miss(x)
is.pmc(x)
is.pmc.ppc(x)
is.pmc.val(x)
is.posteriorchecks(x)
is.raftery(x)
is.rejection(x)
is.sensitivity(x)
is.vb(x)
is.vb.ppc(x)
}
\arguments{
  \item{x}{This is an object that will be subjected to a logical test of
    its class.}
}
\details{
  Functions in Laplace's Demon often assigns a class to an output
  object. For example, after updating a model with the
  \code{\link{LaplacesDemon}} or \code{\link{LaplacesDemon.hpc}}
  function, an output object is created. The output object is of class
  \code{demonoid} or \code{demonoid.hpc}, respectively. Likewise, after
  passing a model to the \code{\link{LaplaceApproximation}} function, an
  output object is created, and it is of class \code{laplace}. The class
  of these and other objects may be logically tested.

  By assigning a class to an output object, the package is able to
  discern which other functions are appropriate for it. For example,
  after updating a model with \code{\link{LaplacesDemon}}, which creates
  an object of class \code{demonoid}, the user may desire to plot its
  output. Since it is assigned a class, the user may use the generic
  \code{plot} function, which internally selects the
  \code{\link{plot.demonoid}} function, which differs from
  \code{\link{plot.laplace}} for objects of class \code{laplace}.

  For more information on object classes, see the \code{\link{class}}
  function.
}
\value{
  The \code{is.bayesfactor} function returns a logical value indicating
  whether or not the supplied object is of class \code{bayesfactor}.
  
  The \code{is.blocks} function returns a logical value indicating
  whether or not the supplied object is of class \code{blocks}.

  The \code{is.bmk} function returns a logical value indicating
  whether or not the supplied object is of class \code{bmk}.
  
  The \code{is.demonoid} function returns a logical value indicating
  whether or not the supplied object is of class \code{demonoid}.

  The \code{is.demonoid.hpc} function returns a logical value indicating
  whether or not the supplied object is of class \code{demonoid.hpc}.

  The \code{is.demonoid.ppc} function returns a logical value indicating
  whether or not the supplied object is of class \code{demonoid.ppc}.

  The \code{is.demonoid.val} function returns a logical value indicating
  whether or not the supplied object is of class \code{demonoid.val}.

  The \code{is.hangartner} function returns a logical value indicating
  whether or not the supplied object is of class \code{hangartner}.

  The \code{is.heidelberger} function returns a logical value indicating
  whether or not the supplied object is of class \code{heidelberger}.

  The \code{is.importance} function returns a logical value indicating
  whether or not the supplied object is of class \code{importance}.

  The \code{is.iterquad} function returns a logical value indicating
  whether or not the supplied object is of class \code{iterquad}.

  The \code{is.iterquad.ppc} function returns a logical value indicating
  whether or not the supplied object is of class \code{iterquad.ppc}.

  The \code{is.juxtapose} function returns a logical value indicating
  whether or not the supplied object is of class \code{juxtapose}.

  The \code{is.laplace} function returns a logical value indicating
  whether or not the supplied object is of class \code{laplace}.

  The \code{is.laplace.ppc} function returns a logical value indicating
  whether or not the supplied object is of class \code{laplace.ppc}.

  The \code{is.miss} function returns a logical value indicating
  whether or not the supplied object is of class \code{miss}.

  The \code{is.pmc} function returns a logical value indicating
  whether or not the supplied object is of class \code{pmc}.

  The \code{is.pmc.ppc} function returns a logical value indicating
  whether or not the supplied object is of class \code{pmc.ppc}.

  The \code{is.pmc.val} function returns a logical value indicating
  whether or not the supplied object is of class \code{pmc.val}.

  The \code{is.posteriorchecks} function returns a logical value
  indicating whether or not the supplied object is of class
  \code{posteriorchecks}.

  The \code{is.raftery} function returns a logical value indicating
  whether or not the supplied object is of class \code{raftery}.

  The \code{is.rejection} function returns a logical value indicating
  whether or not the supplied object is of class \code{rejection}.

  The \code{is.sensitivity} function returns a logical value indicating
  whether or not the supplied object is of class \code{sensitivity}.

  The \code{is.vb} function returns a logical value indicating
  whether or not the supplied object is of class \code{vb}.

  The \code{is.vb.ppc} function returns a logical value indicating
  whether or not the supplied object is of class \code{vb.ppc}.
}
\author{Statisticat, LLC. \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{BayesFactor}},
  \code{\link{Blocks}},
  \code{\link{BMK.Diagnostic}},
  \code{\link{class}},
  \code{\link{Hangartner.Diagnostic}},
  \code{\link{Heidelberger.Diagnostic}},
  \code{\link{Importance}},
  \code{\link{IterativeQuadrature}},
  \code{\link{Juxtapose}},
  \code{\link{LaplaceApproximation}},
  \code{\link{LaplacesDemon}},
  \code{\link{LaplacesDemon.hpc}},
  \code{\link{MISS}},
  \code{\link{PMC}},
  \code{\link{PosteriorChecks}},
  \code{\link{predict.demonoid}},
  \code{\link{predict.laplace}},
  \code{\link{predict.pmc}},
  \code{\link{predict.vb}},
  \code{\link{Raftery.Diagnostic}},
  \code{\link{RejectionSampling}},
  \code{\link{SensitivityAnalysis}},
  \code{\link{Validate}}, and
  \code{\link{VariationalBayes}}.
}
\keyword{Utility}