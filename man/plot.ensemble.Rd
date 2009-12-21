\name{plot.ensemble}
\alias{plot.ensemble}
\title{Plot of Ensemble Estimates}
\description{
    Plot ensemble survival curves and ensemble estimates of mortality.
}
\usage{
    plot.ensemble(x, plots.one.page = TRUE, ...)
}
\arguments{
    \item{x}{An object of class \code{(rsf, grow)} or \code{(rsf,
	predict)}.}
    \item{plots.one.page}{Logical.  Should plots be placed on one page?
        Default is TRUE.}
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

  For competing risk analyses plots correspond to unconditional values
  (i.e., they are non-event specific).  Use \code{competing.risk} for
  event-specific curves and for a more comprehensive analysis in such
  cases.

  Whenever possible, out-of-bag (OOB) values are used.
}
\author{
    Hemant Ishwaran \email{hemant.ishwaran@gmail.com} and
    Udaya B. Kogalur \email{ubk2101@columbia.edu}
}
\references{

    T.A. Gerds, T. Cai and M. Schumacher (2008).  The performance of
    risk prediction models, \emph{Biometrical Jour.}, 4:457-479.
  
    E. Graf, C. Schmoor, W. Sauerbrei and M. Schumacher M (1999).
    Assessment and comparison of prognostic classification
    schemes for survival data, \emph{Statistics in Medicine},
    18:2529-2545.

    H. Ishwaran, U.B. Kogalur (2007).  Random survival forests for R,
    \emph{Rnews}, 7/2:25-31.
}
\seealso{
  \code{competing.risk},
  \code{rsf},
  \code{predict.rsf}.
}
\examples{
data(veteran, package = "randomSurvivalForest") 
v.out <- rsf(Survrsf(time, status)~., veteran, ntree = 1000)
plot.ensemble(v.out)

# plot of ensemble survival for a single individual
surv.ensb <- t(exp(-v.out$oob.ensemble))
plot(v.out$timeInterest, surv.ensb[, 1])

}
\keyword{file}
