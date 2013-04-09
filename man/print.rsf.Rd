\name{print.rsf}
\alias{print.rsf}
\title{Print Summary Output of Analysis}
\description{
    Print summary output from a Random Survival Forests analysis.
    Note that this is the default print method for the package.
}
\usage{
    \method{print}{rsf}(x, ...)
}
\arguments{
    \item{x}{An object of class \code{(rsf, grow)} or \code{(rsf,
	predict)}.}
    \item{...}{Further arguments passed to or from other methods.}
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
  data(veteran, package = "randomSurvivalForest") 
  v.out=rsf(Surv(time, status)~.,veteran, ntree = 1000)
  print.rsf(v.out)
}
\keyword{file}
