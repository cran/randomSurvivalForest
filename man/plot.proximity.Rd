\name{plot.proximity}
\alias{plot.proximity}
\alias{plot.proximity.rsf}
\title{Extract and Plot the Proximity Matrix}
\description{
    Multidimensional scaling plot of the proximity matrix.
}
\usage{
    \method{plot}{proximity}(x, plot = TRUE, ...)
}
\arguments{
    \item{x}{An object of class \code{(rsf, grow)} or \code{(rsf,
	predict)}. }
    \item{plot}{Should proximity be plotted.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
  Extracts the proximity information from the object \code{x} and
  transforms this to a symmetric proximity matrix.  Dissimilarities
  between points are converted into distances using the multidimensional
  scaling function \code{cmdscale} and then plotted.  Overlayed on the
  plot are mortality values, rescaled from 1-100, with 1 indicating low
  mortality, and 100 indicating high mortality.  Mortality values will
  be well separated in successful analyses.  Note that points in blue
  correspond to events, whereas black points are censored observations.
}
\value{
   Invisibly, the proximity matrix with entries transformed to relative
   frequencies.
}
\author{
  Hemant Ishwaran
  Udaya B. Kogalur
}
\seealso{
  \code{rsf},
  \code{predict.rsf}.
}
\examples{
  data(pbc, package = "randomSurvivalForest") 
  pbc.prox.out <- rsf(Surv(days,status) ~ ., pbc, ntree = 100, proximity = TRUE)
  plot.proximity(pbc.prox.out)
}
\keyword{file}
