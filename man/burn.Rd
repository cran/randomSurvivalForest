\name{burn}
\docType{data}
\alias{burn}
\title{Burn Data}
\description{
    Data on 154 burn patients.
}
\format{A data frame containing:
    \tabular{ll}{
        Z1 \tab treatment: 0=routine bathing, 1=Body cleansing \cr
        Z2 \tab gender: 0=male, 1=female \cr
        Z3 \tab race: 0=nonwhite 1=white \cr
        Z4 \tab percentage of total surface area burned \cr
        Z5 \tab burn site indicator head: 1=yes, 0=no \cr
        Z6 \tab burn site indicator buttock: 1=yes, 0=no \cr
        Z7 \tab burn site indicator trunk: 1=yes, 0=no \cr
        Z8 \tab burn site indicator upper leg: 1=yes, 0=no \cr
        Z9 \tab burn site indicator lower leg: 1=yes, 0=no \cr
        Z10 \tab respiratory tract: 1=yes, 0=no \cr
        Z11 \tab type of burn: 1=chemical, 2=scald, 3=electric, 4=flame \cr
        T1  \tab time to excision or on study time \cr
        D1  \tab excision indicator: 1=yes, 0=no \cr
        T2  \tab time to prophylactic antibiotic treatment or on study time \cr
        D2  \tab prophylactic antibiotic treatment: 1=yes 0=no \cr
        T3  \tab time to straphylocous aureaus infection or on study time \cr
        D3  \tab straphylocous aureaus infection: 1=yes 0=no 
      }
}
\source{
  Klein and Moeschberger (1997).
  \emph{Survival Analysis: Techniques for Censored and Truncated Data},
  Springer: New York.
}
\references{
   Ichida et al. (1993) \emph{Stat. Med.}, 12:301-310.
}
\examples{data(burn, package = "randomSurvivalForest")}
\keyword{datasets}
