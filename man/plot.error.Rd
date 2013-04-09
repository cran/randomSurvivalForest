\name{plot.error}
\alias{plot.error}
\alias{plot.error.rsf}
\alias{plot.rsf}
\title{Plot Error Rate and Variable Importance}
\description{
  Plot the out-of-bag (OOB) error rates for the ensemble CHF.  Also plot
  variable importance (VIMP) for predictors.  This is the default plot
  method for the package.
}
\usage{
    \method{plot}{error}(x, sorted = TRUE, ...)
    \method{plot}{rsf}(x, sorted = TRUE, ...)
}
\arguments{
    \item{x}{An object of class \code{(rsf, grow)} or \code{(rsf,
	predict)}.}
    \item{sorted}{Should variables be sorted by importance values?}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
  Plot of OOB error rate, with the bth value being the error rate for
  the ensemble computed using the first b trees.  Error rate is 1-C,
  where C is Harrell's concordance index.  Rates given are between 0 and
  1, with 0.5 representing the benchmark value of a procedure based on
  random guessing.  A value of 0 is perfect.

  Plots VIMP for the x-predictors and prints these values using a matrix
  with up to 3 columns.  The first column is VIMP, the second column is
  standardized VIMP (divided by the maximum importance value), the third
  column is \code{predictorWt} (this is only printed if its values are
  distinct).

  For competing risks, error rates and VIMPs are given for the ensemble
  CHF and for the ensemble conditional CHF for each event type.  For
  more details see Ishwaran et al. (2010).
}
\author{
  Hemant Ishwaran
  Udaya B. Kogalur
}
\references{
    Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S.
    (2008).  Random survival forests, \emph{Ann. App. Statist.}, 2:841-860.

    Ishwaran H., Kogalur U.B., Moore R.D., Gange S.J. and Lau
    B.M. (2010).  Random survival forests for competing risks.

    Breiman L. (2001). Random forests, \emph{Machine Learning}, 45:5-32.

    Harrell F.E. et al. (1982).  Evaluating the yield of medical tests,
    \emph{J. Amer. Med. Assoc.}, 247, 2543-2546.
}
\seealso{
  \code{rsf},
  \code{predict.rsf}.
}
\examples{
\dontrun{
  data(veteran, package = "randomSurvivalForest") 
  v.out <- rsf(Surv(time, status) ~ ., veteran, ntree = 1000, nsplit = 3)
  plot.error(v.out)
}
}
\keyword{file}
