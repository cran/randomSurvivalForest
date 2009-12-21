\name{plot.error}
\alias{plot.error}
\alias{plot.rsf}
\title{Plot of Error Rate and Variable Importance}
\description{
  Plot the out-of-bag (OOB) error rates for the ensemble CHF.  Also plot
  variable importance (VIMP) for predictors.  This is the default plot
  method for the package.
}
\usage{
    plot.error(x, sorted = TRUE, ...)
    plot.rsf(x, sorted = TRUE, ...)
}
\arguments{
    \item{x}{An object of class \code{(rsf, grow)} or \code{(rsf,
	predict)}.}
    \item{sorted}{Should variables be sorted by importance values?}
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

  For competing risks, error rates and VIMPs are given for the ensemble
  CHF and for the ensemble conditional CHF for each event type.
}
\author{
    Hemant Ishwaran \email{hemant.ishwaran@gmail.com} and
    Udaya B. Kogalur \email{ubk2101@columbia.edu}
}
\references{
    H. Ishwaran, U.B. Kogalur, E.H. Blackstone and M.S. Lauer
    (2008).  Random survival forests, \emph{Ann. App.
    Statist.}, 2:841-860.

    L. Breiman (2001). Random forests, \emph{Machine Learning}, 45:5-32.

    F.E. Harrell et al. (1982).  Evaluating the yield of medical tests,
    \emph{J. Amer. Med. Assoc.}, 247, 2543-2546.
}
\seealso{
  \code{rsf},
  \code{predict.rsf}.
}
\examples{
\dontrun{
  data(veteran, package = "randomSurvivalForest") 
  v.out <- rsf(Survrsf(time, status)~., veteran, ntree = 1000, nsplit = 3)
  plot.error(v.out)
}
}
\keyword{file}
