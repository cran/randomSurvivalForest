\name{recid}
\docType{data}
\alias{recid}
\title{Recidivism Data}
\description{
  This data comprises 432 male inmates released from Maryland state
  prisons in the early 1970's.  Dates of arrest, which constitutes an
  event, were recorded for each inmate up to 1 year following release.
  A major goal of the analysis was to determine if those who received
  financial aid had lower arrest rates than those who did not.
  Moreover, the effect of having financial aid is likely to be more
  important immediately following arrest and less important as time
  passes.
}
\format{A data frame containing:
    \tabular{ll}{
        week   \tab survival time \cr
        arrest \tab censoring variable \cr
        fin    \tab financial aid, 0=no 1=yes \cr
        age    \tab age in years \cr
        race   \tab race (1=black, 0=otherwise) \cr
        wexp   \tab full time work experience (0=no, 1=yes) \cr
        mar    \tab married (0=no, 1=yes) \cr
        paro   \tab released on parole (0=no, 1=yes) \cr
        prio   \tab number of prior convictions
    }
}
\source{Rossi, Berk and Lenihan, 1980.}
\references{P.H. Rossi and R.A. Berk and K.J. Lenihan (1980).
    \emph{Money, Work and Crime: Some Experimental Results},
     Academic Press:New York.}
}
\examples{data(recid, package = "randomSurvivalForest")}
\keyword{datasets}
