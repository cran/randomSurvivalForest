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
  This shows the ensemble survival function for each individual in the
  data.  Thick red line is overall ensemble survival, thick green line
  is Nelson-Aalen estimator. (2) This is a comparison of the population
  ensemble survival function to the Nelson-Aalen estimator. (3) The
  Brier score, a value between 0 and 1 (where 0=perfect, 1=poor, and
  0.25=guessing) is plotted at each of the unique event times with plot
  stratified by ensemble mortality value (see Graf et al. for more
  background on Brier scores).  Stratification is into 4 groups
  corresponding to the 0-25, 25-50, 50-75 and 75-100 percentile values
  of mortality.  Green line is non-stratified score.  (4) Plot of
  estimated mortality versus observed time.  Points in blue correspond
  to events, black points are censored observations.

  Note that when \code{x} is of class \code{(rsf, predict)} not all
  plots will be produced.
}
\author{
    Hemant Ishwaran \email{hemant.ishwaran@gmail.com} and
    Udaya B. Kogalur \email{ubk2101@columbia.edu}
}
\references{
    H. Ishwaran and Udaya B. Kogalur (2006).  Random Survival
    Forests.  \emph{Cleveland Clinic Technical Report}.   

    E. Graf, C. Schmoor, W. Sauerbrei and M. Schumacher M (1999).
    Assessment and comparison of prognostic classification
    schemes for survival data, \emph{Statistics in Medicine},
    18:2529-2545.
}
\seealso{
  \code{rsf},
  \code{predict.rsf}.
}
\examples{
  data(veteran, package = "randomSurvivalForest") 
  v.out <- rsf(Survrsf(time, status)~., veteran, ntree = 1000)
  plot.ensemble(v.out)
}
\keyword{file}