\name{follic}
\docType{data}
\alias{follic}
\title{Follicular Cell Lymphoma}
\description{
    Competing risk data set involving follicular cell lymphoma.
}
\format{A data frame containing:
    \tabular{ll}{
        age      \tab age\cr
	hgb      \tab hemoglobin (g/l)\cr
        clinstg	 \tab clinical stage: 1=stage I, 2=stage II\cr
	ch       \tab chemotherapy: Y=Yes, N=No\cr
	rt       \tab radiotherapy: Y=Yes, N=No\cr
        time     \tab first failure time\cr
        status   \tab censoring status: 0=censored, 1=relapse, 2=death
    }
  }

\source{Follicular cell lymphoma data set. Table 1.4b, \emph{Competing
    Risks: A Practical Perspective.}  541 patients identified as having
    follicular type lymphoma, with early stage disease (I or II), and
    treated with radiation alone (RT) or with radiation and chemotherapy
    (CMT).  }
\references{
    Pintilie M., (2006) \emph{Competing Risks: A
      Practical Perspective.} West Sussex: John Wiley and Sons.
}
\examples{data(follic, package = "randomSurvivalForest")}
\keyword{datasets}
