\name{plot.ensemble}
\alias{plot.ensemble}
\alias{plot.ensemble.rsf}
\title{Plot Ensemble Estimates}
\description{
    Plot ensemble survival curves and ensemble estimates of mortality.
}
\usage{
    \method{plot}{ensemble}(x, plots.one.page = TRUE, ...)
}
\arguments{
    \item{x}{An object of class \code{(rsf, grow)} or \code{(rsf,
	predict)}.}
    \item{plots.one.page}{Should plots be placed on one page?}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{  
  Four plots are produced.  Going from top to bottom, left to right: (1)
  Ensemble survival function for each individual.  Thick red line is
  overall ensemble survival, thick green line is Nelson-Aalen
  estimator. (2) Comparison of the population ensemble survival function
  to the Nelson-Aalen estimator. (3) Brier score (0=perfect, 1=poor, and
  0.25=guessing) stratified by ensemble mortality.  Based on the method
  described in Gerds et al. (2008) in which the censoring distribution
  is estimated using the Kaplan-Meier estimator.  Stratification is into
  4 groups corresponding to the 0-25, 25-50, 50-75 and 75-100 percentile
  values of mortality.  Red line is the overall (non-stratified) Brier
  score.  (4) Plot of mortality versus observed time.  Points in blue
  correspond to events, black points are censored observations.

  Note that when \code{x} is of class \code{(rsf, predict)} not all
  plots will be produced.

  Applies to competing risk analyses but plots are non-event specific.
  Use \code{competing.risk} for event-specific curves and for a more
  comprehensive analysis in such cases.

  Whenever possible, out-of-bag (OOB) values are used.
}
\author{
  Hemant Ishwaran
  Udaya B. Kogalur
}
\references{
    Gerds T.A., Cai T. and Schumacher M. (2008).  The performance of
    risk prediction models, \emph{Biometrical J.}, 4:457-479.
  
    Graf E., Schmoor C., Sauerbrei W. and Schumacher M. (1999).
    Assessment and comparison of prognostic classification
    schemes for survival data, \emph{Statist. in Med.},
    18:2529-2545.
}
\seealso{
  \code{rsf},
  \code{predict.rsf}.
}
\examples{
data(veteran, package = "randomSurvivalForest") 
v.out <- rsf(Surv(time, status) ~ ., veteran, ntree = 1000)
plot.ensemble(v.out)

# plot of ensemble survival for a single individual
surv.ensb <- t(exp(-v.out$oob.ensemble))
plot(v.out$timeInterest, surv.ensb[, 1])
}
\keyword{file}
