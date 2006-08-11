\name{print.rsf}
\alias{print.rsf}
\title{Print Summary Output of Analysis}
\description{
    Print summary output from Random Survival Forests analysis.
    
}
\usage{
    print.rsf(x, ...)
}
\arguments{
    \item{x}{An object of class \code{randomSurvivalForest},
    as that created by the function \code{rsf}.}
    \item{...}{Further arguments passed to or from other methods.}
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
  \code{plot.ensemble},
  \code{plot.variable},
  \code{plot.error},
  \code{plot.proximity}.
}
\examples{
  data(veteran, package = "randomSurvivalForest") 
  v.out=rsf(Survrsf(time, status)~.,veteran,ntree = 1000)
  print.rsf(v.out)
}
\keyword{file}
