\name{plot.ensemble}
\alias{plot.ensemble}
\title{Plot of Ensemble Estimates}
\description{
    Plot ensemble survival curves and ensemble estimates of mortality.
}
\usage{
    plot.ensemble(x, main = deparse(substitute(x)), ...)
}
\arguments{
    \item{x}{An object of class \code{randomSurvivalForest},
    as that created by the function \code{rsf}.}
    \item{main}{Plot title.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
  Four plots are produced.  Going from top to bottom, left to right: (1)
  This shows the ensemble survival function for each individual in the
  data.  Thick red line is overall ensemble survival, thick green line
  is Nelson-Aalen estimator. (2) This is a comparison of the population
  ensemble survival function to the Nelson-Aalen estimator. (3) The
  Brier score, a value between 0 and 1 (where 0=perfect, 1=poor, and
  0.25=guessing; see reference below) is plotted at each of the unique
  event times with plot stratified by ensemble mortality value.
  Stratification is into 4 groups corresponding to the 0-25, 25-50,
  50-75 and 75-100 percentile values of mortality.  Green line is
  non-stratified score.  (4) Plot of estimated mortality versus observed
  time.  Points in blue correspond to events, black points are censored
  observations.
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
  \code{print.rsf},
  \code{plot.variable},
  \code{plot.error},
  \code{plot.proximity}.
}
\examples{
  data(veteran, package = "randomSurvivalForest") 
  v.out <- rsf(Survrsf(time, status)~., veteran, ntree = 1000)
  plot.ensemble(v.out)
}
\keyword{file}
