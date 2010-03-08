\name{plot.proximity}
\alias{plot.proximity}
\title{Plot of Proximity}
\description{
    Multidimensional scaling plot of proximity matrix.
}
\usage{
    plot.proximity(x, plot = TRUE, ...)
}
\arguments{
    \item{x}{An object of class \code{(rsf, grow)} or \code{(rsf,
	predict)}. }
    \item{plot}{Logical.  If TRUE, proximity is plotted.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
  Extracts proximity information from \code{x} and transforms this to a
  symmetric proximity matrix.  Dissimilarities between points are then
  converted into distances using the R multidimensional scaling function
  \code{cmdscale} and then plotted.  Overlayed on the plot are mortality
  values, rescaled from 1-100, with 1 indicating low mortality, and 100
  indicating high mortality.  Mortality values will be well separated in
  successful analyses.  Note that points in blue correspond to events,
  whereas black points are censored observations.
}
\value{
   Invisibly, the proximity matrix with entries transformed to relative
   frequencies.
}
\author{
    Hemant Ishwaran \email{hemant.ishwaran@gmail.com} and
    Udaya B. Kogalur \email{kogalurshear@gmail.com}
}
\seealso{
  \code{rsf},
  \code{predict.rsf}.
}
\examples{
  data(pbc, package = "randomSurvivalForest") 
  pbc.prox.out <- rsf(Survrsf(days,status)~., pbc, ntree = 100, proximity = TRUE)
  plot.proximity(pbc.prox.out)
}
\keyword{file}
