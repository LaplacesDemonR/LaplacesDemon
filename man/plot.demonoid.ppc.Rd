\name{plot.demonoid.ppc}
\alias{plot.demonoid.ppc}
\title{Plots of Posterior Predictive Checks}
\description{
  This may be used to plot, or save plots of, samples in an object of
  class \code{demonoid.ppc}. A variety of plots is provided.
}
\usage{\method{plot}{demonoid.ppc}(x, Style=NULL, Data=NULL, Rows=NULL,
     PDF=FALSE, \dots)}
\arguments{
  \item{x}{
    This required argument is an object of class \code{demonoid.ppc}.}
  \item{Style}{
    This optional argument specifies one of several styles of plots,
    and defaults to \code{NULL} (which is the same as
    \code{"Density"}). Styles of plots are indicated in quotes. Optional
    styles include \code{"Covariates"}, \code{"Covariates, Categorical DV"},
    \code{"Density"}, \code{"DW"}, \code{"DW, Multivariate, C"},
    \code{"ECDF"}, \code{"Fitted"}, \code{"Fitted, Multivariate, C"},
    \code{"Fitted, Multivariate, R"}, \code{"Jarque-Bera"},
    \code{"Jarque-Bera, Multivariate, C"},
    \code{"Mardia"}, \code{"Predictive Quantiles"},
    \code{"Residual Density"}, \code{"Residual Density, Multivariate, C"},
    \code{"Residual Density, Multivariate, R"},
    \code{"Residuals"}, \code{"Residuals, Multivariate, C"},
    \code{"Residuals, Multivariate, R"}, \code{"Space-Time by Space"},
    \code{"Space-Time by Time"}, \code{"Spatial"},
    \code{"Spatial Uncertainty"}, \code{"Time-Series"},
    \code{"Time-Series, Multivariate, C"}, and
    \code{"Time-Series, Multivariate, R"}. Details are given below.}
  \item{Data}{
    This optional argument accepts the data set used when updating the
    model. Data is required only with certain plot styles, including
    \code{"Covariates"}, \code{"Covariates, Categorical DV"},
    \code{"DW, Multivariate, C"}, \code{"Fitted, Multivariate, C"},
    \code{"Fitted, Multivariate, R"},
    \code{"Jarque-Bera, Multivariate, C"}, \code{"Mardia"},
    \code{"Residual Density, Multivariate, C"},
    \code{"Residual Density, Multivariate, R"},
    \code{"Residuals, Multivariate, C"},
    \code{"Residuals, Multivariate, R"},
    \code{"Space-Time by Space"}, \code{"Space-Time by Time"},
    \code{"Spatial"}, \code{"Spatial Uncertainty"},
    \code{"Time-Series, Multivariate, C"}, and
    \code{"Time-Series, Multivariate, R"}.}
  \item{Rows}{
    This optional argument is for a vector of row numbers that 
    specify the records associated by row in the object of class
    \code{demonoid.ppc}. Only these rows are plotted. The default is to
    plot all rows. Some plots do not allow rows to be specified.}
  \item{PDF}{
    This logical argument indicates whether or not the user wants
    Laplace's Demon to save the plots as a .pdf file.}
  \item{\dots}{Additional arguments are unused.}
}
\details{
  This function can be used to produce a variety of posterior predictive
  plots, and the style of plot is selected with the \code{Style}
  argument. Below are some notes on the styles of plots.

  \code{Covariates} requires \code{Data} to be specified, and also
  requires that the covariates are named \code{X} or \code{x}. A plot
  is produced for each covariate column vector against yhat, and is
  appropriate when y is not categorical.

  \code{Covariates, Categorical DV} requires \code{Data} to be
  specified, and also requires that the covariates are named \code{X} or
  \code{x}. A plot is produced for each covariate column vector against
  yhat, and is appropriate when y is categorical.
  
  \code{Density} plots show the kernel density of the posterior
  predictive distribution for each selected row of y (all are selected
  by default). A vertical red line indicates the position of the
  observed y along the x-axis. When the vertical red line is close to
  the middle of a normal posterior predictive distribution, then there
  is little discrepancy between y and the posterior predictive
  distribution. When the vertical red line is in the tail of the
  distribution, or outside of the kernel density altogether, then
  there is a large discrepancy between y and the posterior predictive
  distribution. Large discrepancies may be considered outliers, and
  moreover suggest that an improvement in model fit should be
  considered.

  \code{DW} plots the distributions of the Durbin-Watson (DW) test
  statistics (Durbin and Watson, 1950), both observed
  (\eqn{d^{obs}}{d.obs} as a transparent, black density) and replicated
  (\eqn{d^{rep}}{d.rep} as a transparent, red density). The distribution
  of \eqn{d^{obs}}{d.obs} is estimated from the model, and
  \eqn{d^{rep}}{d.rep} is simulated from normal residuals without
  autocorrelation, where the number of simulations are the same as the
  observed number. This DW test may be applied to the residuals of
  univariate time-series models (or otherwise ordered residuals) to
  detect first-order autocorrelation. Autocorrelated residuals are not
  independent. The DW test is applicable only when the residuals are
  normally-distributed, higher-order autocorrelation is not present, and
  y is not used also as a lagged predictor. The DW test statistic,
  \eqn{d^{obs}}{d[obs]}, occurs in the interval (0,4), where 0 is
  perfect positive autocorrelation, 2 is no autocorrelation, and 4 is
  perfect negative autocorrelation. The following summary is reported on
  the plot: the mean of \eqn{d^{obs}}{d[obs]} (and its 95\% probability
  interval), the probability that \eqn{d^{obs} > d^{rep}}{d[obs] >
  d[rep]}, and whether or not autocorrelation is found. Positive
  autocorrelation is reported when the observed process is greater than
  the replicated process in 2.5\% of the samples, and negative
  autocorrelation is reported when the observed process is greater than
  the replicated process in 97.5\% of the samples.

  \code{DW, Multivariate, C} requires \code{Data} to be
  specified, and also requires that variable \code{Y} exist in the data
  set with exactly that name. These plots compare each column-wise
  vector of residuals with a univariate Durbin-Watson test, as in
  \code{DW} above. This plot is appropriate when Y is multivariate, not
  categorical, and residuals are desired to be tested column-wise for
  first-order autocorrelation.
  
  \code{ECDF} (Empirical Cumulative Distribution Function) plots compare
  the ECDF of y with three ECDFs of yhat based on the 2.5\%, 50\%
  (median), and 97.5\% of its distribution. The ECDF(y) is defined as
  the proportion of values less than or equal to y. This plot is
  appropriate when y is univariate and at least ordinal.

  \code{Fitted} plots compare y with the probability interval of its
  replicate, and provide loess smoothing. This plot is appropriate when
  y is univariate and not categorical.

  \code{Fitted, Multivariate, C} requires \code{Data} to be
  specified, and also requires that variable \code{Y} exists in the data
  set with exactly that name. These plots compare each column-wise
  vector of y in Y with its replicates and provide loess smoothing.
  This plot is appropriate when Y is multivariate, not categorical, and
  desired to be seen column-wise.
  
  \code{Fitted, Multivariate, R} requires \code{Data} to be
  specified, and also requires that variable \code{Y} exists in the data
  set with exactly that name. These plots compare each row-wise
  vector of y in Y with its replicates and provide loess smoothing.
  This plot is appropriate when Y is multivariate, not categorical, and
  desired to be seen row-wise.

  \code{Jarque-Bera} plots the distributions of the Jarque-Bera (JB)
  test statistics (Jarque and Bera, 1980), both observed
  (\eqn{JB^{obs}}{JB.obs} as a transparent black density) and replicated
  (\eqn{JB^{rep}}{JB.rep} as a transparent red density). The
  distribution of \eqn{JB^{obs}}{JB.obs} is estimated from the model,
  and \eqn{JB^{rep}}{JB.rep} is simulated from normal residuals, where
  the number of simulations are the same as the observed number. This
  Jarque-Bera test may be applied to the residuals of 
  univariate models to test for normality. The Jarque-Bera test does not
  test normality per se, but whether or not the distribution has
  kurtosis and skewness that match a normal distribution, and is
  therefore a test of the moments of a normal distribution. The
  following summary is reported on the plot: the mean of
  \eqn{JB^{obs}}{JB[obs]} (and its 95\% probability interval), the
  probability that \eqn{JB^{obs} > JB^{rep}}{JB[obs] > JB[rep]}, and
  whether or not normality is indicated. Non-normality is reported when
  the observed process is greater than the replicated process in either
  2.5\% or 97.5\% of the samples.

  \code{Jarque-Bera, Multivariate, C} requires \code{Data} to be
  specified, and also requires that variable \code{Y} exist in the data
  set with exactly that name. These plots compare each column-wise
  vector of residuals with a univariate Jarque-Bera test, as in
  \code{Jarque-Bera} above. This plot is appropriate when Y is
  multivariate, not categorical, and residuals are desired to be
  tested column-wise for normality.

  \code{Mardia} plots the distributions of the skewness (K3) and
  kurtosis (K4) test statistics (Mardia, 1970), both observed
  (\eqn{K3^{obs}}{K3.obs} and \eqn{K4^{obs}}{K4.obs} as transparent
  black density) and replicated (\eqn{K3^{rep}}{K3.rep} and
  \eqn{K4^{rep}}{K4.rep} as transparent red density). The distributions
  of \eqn{K3^{obs}}{K3.obs} and \eqn{K4^{obs}}{K4.obs} are estimated
  from the model, and both \eqn{K3^{rep}}{K3.rep} \eqn{K4^{rep}}{K4.rep}
  are simulated from multivariate normal residuals, where the number of
  simulations are the same as the observed number. This Mardia's test
  may be applied to the residuals of multivariate models to test for
  multivariate normality. Mardia's test does not test for multivariate
  normality per se, but whether or not the distribution has kurtosis and
  skewness that match a multivariate normal distribution, and is 
  therefore a test of the moments of a multivariate normal
  distribution. The following summary is reported on the plots: the
  means of \eqn{K3^{obs}}{K3[obs]} and \eqn{K4^{obs}}{K4[obs]} (and
  the associated 95\% probability intervals), the probabilities that
  \eqn{K3^{obs} > K3^{rep}}{K3[obs] > K3[rep]} and
  \eqn{K4^{obs} > K4^{rep}}{K4[obs] > K4[rep]}, and whether or not
  multivariate normality is indicated. Non-normality is reported when 
  the observed process is greater than the replicated process in either
  2.5\% or 97.5\% of the samples. \code{Mardia} requires \code{Data} to
  be specified, and also requires that variable \code{Y} exist in the
  data set with exactly that name. \code{Y} must be a \eqn{N \times P}{N
  x P} matrix of \eqn{N}{N} records and \eqn{P}{P} variables. Source
  code was modified from the deprecated package QRMlib.
  
  \code{Predictive Quantiles} plots compare y with the predictive
  quantile (PQ) of its replicate. This may be useful in looking for
  patterns with outliers. Instances outside of the gray lines are
  considered outliers.

  \code{Residual Density} plots the residual density of the median of
  the samples. A vertical red line occurs at zero. This plot may be
  useful for inspecting a distributional assumption of residual
  variance. This plot is appropriate when y is univariate and
  continuous.

  \code{Residual Density, Multivariate C} requires \code{Data} to be
  specified, and also requires that variable \code{Y} exist in the data
  set with exactly that name. These are column-wise plots of residual
  density, given the median of the samples. These plots may be useful
  for inspecting a distributional assumption of residual variance.
  This plot is appropriate when Y is multivariate, continuous, and
  densities are desired to be seen column-wise.

  \code{Residual Density, Multivariate R} requires \code{Data} to be
  specified, and also requires that variable \code{Y} exist in the data
  set with exactly that name. These are row-wise plots of residual
  density, given the median of the samples. These plots may be useful
  for inspecting a distributional assumption of residual variance.
  This plot is appropriate when Y is multivariate, continuous, and
  densities are desired to be seen row-wise.

  \code{Residuals} plots compare y with its residuals. The probability
  interval is plotted as a line. This plot is appropriate when y
  is univariate.

  \code{Residuals, Multivariate, C} requires \code{Data} to be
  specified, and also requires that variable \code{Y} exist in the data
  set with exactly that name. These are plots of each column-wise
  vector of residuals. The probability interval is plotted as a
  line. This plot is appropriate when Y is multivariate, not
  categorical, and the residuals are desired to be seen column-wise.

  \code{Residuals, Multivariate, R} requires \code{Data} to be
  specified, and also requires that variable \code{Y} exist in the data
  set with exactly that name. These are plots of each row-wise
  vector of residuals. The probability interval is plotted as a
  line. This plot is appropriate when Y is multivariate, not
  categorical, and the residuals are desired to be seen row-wise.
  
  \code{Space-Time by Space} requires \code{Data} to be specified, and
  also requires that the following variables exist in the data set with
  exactly these names: \code{latitude}, \code{longitude}, \code{S}, and
  \code{T}. These space-time plots compare the S x T matrix Y with the S
  x T matrix Yrep, producing one time-series plot per point s in space,
  for a total of S plots. Therefore, these are time-series plots for
  each point s in space across T time-periods. See \code{Time-Series}
  plots below.

  \code{Space-Time by Time} requires \code{Data} to be specified, and
  also requires that the following variables exist in the data set with
  exactly these names: \code{latitude}, \code{longitude}, \code{S}, and
  \code{T}. These space-time plots compare the S x T matrix Y with the S
  x T matrix Yrep, producing one spatial plot per time-period, and T
  plots will be produced. See \code{Spatial} plots below.

  \code{Spatial} requires \code{Data} to be specified, and also requires
  that the following variables exist in the data set with exactly these
  names: \code{latitude} and \code{longitude}. This spatial plot shows
  yrep plotted according to its coordinates, and is color-coded so that
  higher values of yrep become more red, and lower values become more
  yellow.

  \code{Spatial Uncertainty} requires \code{Data} to be specified, and
  also requires that the following variables exist in the data set with
  exactly these names: \code{latitude} and \code{longitude}. This
  spatial plot shows the probability interval of yrep plotted according
  to its coordinates, and is color-coded so that wider probability
  intervals become more red, and lower values become more yellow.

  \code{Time-Series} plots compare y with its replicate, including the
  median and probability interval quantiles. This plot is appropriate
  when y is univariate and ordered by time.

  \code{Time-Series, Multivariate, C} requires \code{Data} to be
  specified, and also requires that variable \code{Y} exist in the data
  set with exactly that name. These plots compare each column-wise
  time-series in Y with its replicate, including the median and
  probability interval quantiles. This plot is appropriate when y is
  multivariate and each time-series is indexed by column in Y.

  \code{Time-Series, Multivariate, R} requires \code{Data} to be
  specified, and also requires that variable \code{Y} exist in the data
  set with exactly that name. These plots compare each row-wise
  time-series in Y with its replicate, including the median and
  probability interval quantiles. This plot is appropriate when y is
  multivariate and each time-series is indexed by row in Y, such as is
  typically true in panel models.
}
\references{
  Durbin, J., and Watson, G.S. (1950). "Testing for Serial Correlation
  in Least Squares Regression, I." \emph{Biometrika}, 37, p. 409--428.
  
  Jarque, C.M. and Bera, A.K. (1980). "Efficient Tests for Normality,
  Homoscedasticity and Serial Independence of Regression Residuals".
  \emph{Economics Letters}, 6(3), p. 255--259.

  Mardia, K.V. (1970). "Measures of Multivariate Skewness and Kurtosis
  with Applications". \emph{Biometrika}, 57(3), p. 519--530.
}
\author{Statisticat, LLC \email{software@bayesian-inference.com}}
\seealso{
  \code{\link{LaplacesDemon}} and
  \code{\link{predict.demonoid}}.
}
\examples{### See the LaplacesDemon function for an example.}
\keyword{Plot}