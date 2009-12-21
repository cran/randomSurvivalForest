\name{breast}
\docType{data}
\alias{breast}
\title{German Breast Cancer Study Group 2}
\description{
  Breast cancer survival data.
}
\format{A data frame containing:
  \tabular{ll}{
    hormone \tab hormonal therapy: no, yes\cr
    age     \tab age (in years)\cr
    menostat\tab menopausal status: pre=premenopausal, post=postmenopausal\cr
    size    \tab tumor size (in mm)\cr
    grade   \tab tumor grade: I < II < III\cr
    pnodes  \tab number of positive nodes\cr
    progrec \tab progesterone receptor (in fmol)\cr
    estrec  \tab estrogen receptor (in fmol)\cr
    time    \tab recurrence free survival time (in days)\cr
    cens    \tab censoring status: 0=censored, 1=recurrence\cr
  }
}
\references{
  M. Schumacher, G. Basert, H. Bojar,  K. Huebner, M. Olschewski, 
  W. Sauerbrei, C. Schmoor, C. Beyerle, R.L.A. Neumann and H.F. Rauschecker
  for the German Breast Cancer Study Group (1994),
  Randomized \eqn{2\times2} trial evaluating hormonal treatment
  and the duration of chemotherapy in node-positive breast cancer patients.
  \emph{Journal of Clinical Oncology}, \bold{12}, 2086--2093.

  W. Sauerbrei and P. Royston (1999). Building multivariable prognostic
  and diagnostic models: transformation of the predictors by using 
  fractional polynomials. \emph{Journal of the Royal Statistics Society
  Series A}, Volume \bold{162}(1), 71--94.
}
\examples{data(breast, package = "randomSurvivalForest")}
\keyword{datasets}
