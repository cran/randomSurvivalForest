\name{plot.error}
\alias{plot.error}
\alias{plot.rsf}
\title{Plot of Error Rate and Variable Importance}
\description{
  Plot out-of-bag (OOB) error rate for the ensemble as a
  function of number of trees in the forest.  Also plots importance
  values for predictors.  Note this is the default plot method
  for the package.
}
\usage{
    plot.error(x, ...)
    plot.rsf(x, ...)
}
\arguments{
    \item{x}{An object of class \code{(rsf, grow)} or \code{(rsf,
	predict)}.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
  Plot of OOB error rate, with the b-th value being the error rate for
  the ensemble computed using the first b trees.  Error rate is 1-C,
  where C is Harrell's concordance index.  Rates given are between 0 and
  1, with 0.5 representing the benchmark value of a procedure based on
  random guessing.  A value of 0 is perfect.

  In the orginal call, if \code{importance}=TRUE (the default setting),
  then variable importance (VIMP) values for predictors are plotted.  A
  matrix with up to 3 columns is also printed.  First column are VIMP
  values, second column are standardized VIMP values (divided by the
  maximum importance value), third column is the vector
  \code{predictorWt} (only printed if values are distinct).
}
\author{
    Hemant Ishwaran \email{hemant.ishwaran@gmail.com} and
    Udaya B. Kogalur \email{ubk2101@columbia.edu}
}
\references{
    H. Ishwaran, Udaya B. Kogalur, Eugene H. Blackstone and Michael S. Lauer
    (2007).  Random Survival Forests. \emph{Cleveland Clinic Technical
    Report}.

    L. Breiman (2001). Random forests, \emph{Machine Learning}, 45:5-32.

    F.E. Harrell et al. (1982).  Evaluating the yield of medical tests,
    \emph{J. Amer. Med. Assoc.}, 247, 2543-2546.
}
\seealso{
  \code{rsf},
  \code{predict.rsf}.
}
\examples{
  data(veteran, package = "randomSurvivalForest") 
  v.out <- rsf(Survrsf(time, status)~., veteran, ntree = 1000)
  plot.error(v.out, veteran)
}
\keyword{file}
