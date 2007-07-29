\name{Survrsf}
\alias{Survrsf}
\title{Definition of Time and Censoring Information}
\description{
    \code{Survrsf} defines the time and censoring information
    that will be extracted from the data.
}
\usage{
    Survrsf(t, d)
}
\arguments{
    \item{t}{survival time}
    \item{d}{censoring indicator (0=censored, 1=death)}
}
\value{
     An object of class \code{rsf.formula}.
}
\author{
    Hemant Ishwaran \email{hemant.ishwaran@gmail.com} and
    Udaya B. Kogalur \email{ubk2101@columbia.edu}
}
\references{
    H. Ishwaran, Udaya B. Kogalur, Eugene H. Blackstone and Michael S. Lauer
    (2007).  Random Survival Forests. \emph{Cleveland Clinic Technical
    Report}.
}
\seealso{
  \code{rsf}.
}
\examples{
  data(veteran, package = "randomSurvivalForest") 
  v.out <- rsf(Survrsf(time, status)~., veteran, ntree = 1000)
}
\keyword{file}
