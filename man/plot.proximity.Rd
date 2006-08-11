\name{plot.proximity}
\alias{plot.proximity}
\title{Plot of Proximity}
\description{
    Multidimensional scaling plot of proximity matrix.
}
\usage{
    plot.proximity(x, main = deparse(substitute(x)), ...)
}
\arguments{
    \item{x}{An object of class \code{randomSurvivalForest},
    as that created by the function \code{rsf}.}
    \item{main}{Plot title.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
  Extracts proximity information from \code{x} and converts it to a
  symmetric proximity matrix.  Dissimilarities between points are
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
  \code{print.rsf},
  \code{plot.ensemble},
  \code{plot.variable},
  \code{plot.error},
  \code{cmdscale}.
}
\examples{
  data(pbc, package = "randomSurvivalForest") 
  pbc.prox.out <- rsf(Survrsf(days,status)~., pbc, ntree = 1000,
                      proximity = TRUE)
  plot.proximity(pbc.prox.out)
}
\keyword{file}
