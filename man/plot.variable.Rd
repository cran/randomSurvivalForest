\name{plot.variable}
\alias{plot.variable}
\title{Plot of Ensemble Survival Effect of Predictors}
\description{
  Plot of ensemble mortality for each predictor.
}
\usage{
    plot.variable(x,
                  plots.per.page = 4,
                  granule = 25,
                  prednames = NULL,
                  main = deparse(substitute(x)),
                  ...)
}
\arguments{
    \item{x}{An object of class \code{randomSurvivalForest},
        as that created by the function \code{rsf}.}
    \item{plots.per.page}{Integer value controlling page layout.}
    \item{granule}{Integer value controlling whether a plot for a
        specific predictor should be given as a boxplot or scatter 
        plot.  Larger values coerce boxplots.}
    \item{prednames}{Character vector of predictor names.  Only these
        predictors will be plotted.}
    \item{main}{Plot title.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
    The ensemble mortality is plotted against each predictor used in
    growing the forest.  Ensemble mortality values (vertical axis)
    should be interpreted in terms of total number of deaths.  For
    example, if individual \code{i} has a mortality value of 100, then
    if all individuals had the same predictor as \code{i}, there would
    be 100 deaths in the dataset (on average).
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
  \code{plot.error},
  \code{plot.proximity}.
}
\examples{
  data(veteran,package = "randomSurvivalForest") 
  v.out <- rsf(Survrsf(time,status)~., veteran, ntree = 1000)
  plot.variable(v.out)
  plot.variable(v.out, plots.per.page = 2, granule = 5,
                prednames = c("trt", "karno", "age"))
}
\keyword{file}
