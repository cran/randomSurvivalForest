\name{plot.proximity}
\alias{plot.proximity}
\title{Plot of Proximity}
\description{
    Multidimensional scaling plot of proximity matrix.
}
\usage{
    plot.proximity(x, ...)
}
\arguments{
    \item{x}{An object of class \code{(rsf, grow)} or \code{(rsf,
	predict)}. Note that \code{proximity}=TRUE must be used in
        the original \code{rsf} call.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
  Extracts proximity information from \code{x} and converts this to a
  symmetric proximity matrix.  Dissimilarities between points are then
  converted into distances using the R multidimensional scaling function
  \code{cmdscale} and then plotted.  Overlayed on the plot are mortality
  values, rescaled from 1-100, with 1 indicating low mortality, and 100
  indicating high mortality.  Mortality values will be well separated in
  successful analyses.
}
\value{
   Invisibly, the proximity matrix.
}
\author{
    Hemant Ishwaran \email{hemant.ishwaran@gmail.com} and
    Udaya B. Kogalur \email{ubk2101@columbia.edu}
}
\references{
    H. Ishwaran and Udaya B. Kogalur (2006).  Random Survival
    Forests.  \emph{Cleveland Clinic Technical Report}.
}
\seealso{
  \code{rsf},
  \code{predict.rsf}.
}
\examples{
  data(pbc, package = "randomSurvivalForest") 
  pbc.prox.out <- rsf(Survrsf(days,status)~., pbc, ntree = 1000, proximity = TRUE)
  plot.proximity(pbc.prox.out)
}
\keyword{file}
