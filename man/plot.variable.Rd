\name{plot.variable}
\alias{plot.variable}
\title{Plot of Ensemble Survival Effect of Variables}
\description{
  Plot of ensemble mortality, predicted median survival, or predicted
  survival time for each variable.  Users can select between marginal
  and partial plots.
}
\usage{
    plot.variable(x,
                  plots.per.page = 4,
                  granule = 5,
                  sorted = TRUE,
                  type = c("mort", "rel.freq", "surv", "time")[1],
                  partial = FALSE,
                  predictorNames = NULL,
                  npred = NULL,
                  npts = 25,
                  subset = NULL,
                  ...)
}
\arguments{
    \item{x}{An object of class \code{(rsf, grow)} or \code{(rsf,
	predict)}.}
    \item{plots.per.page}{Integer value controlling page layout.}
    \item{granule}{Integer value controlling whether a plot for a
        specific variable should be given as a boxplot or scatter 
        plot. Larger values coerce boxplots.}
    \item{sorted}{Should variables be sorted by importance values (only
        applies if importance values are available)?  Default is TRUE.}
      \item{type}{Select type of value to be plotted on the vertical
	axis. See details.}
    \item{partial}{Logical. Should partial plots be created?  Default
        is FALSE.}
    \item{predictorNames}{Character vector of variable names.  Only these
        variables will be plotted.  Default is all.}
    \item{npred}{Number of variables to be plotted (only applies when
	\code{predictorNames}=NULL).  Default is all.}
    \item{npts}{Maximum number of points used when generating partial
      plots for continuous variables.}
    \item{subset}{Indices indicating which rows of the predictor matrix
      to be used (note: this applies to the \emph{processed} predictor
      matrix, \code{predictors} of the object).  Default is to use all
      rows.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
    Either mortality, relative frequency of mortality, predicted median
    survival, or predicted survival times are plotted on the vertical
    axis (y-value) against x-variables on the horizontal axis.  The
    choice of x-variables are specified by \code{predictorNames} or
    using \code{npred}.  The choice of y-value is controlled by
    \code{type}.  There are 4 different choices: \code{mort} is ensemble
    mortality, \code{rel.freq} is standardized mortality, \code{surv} is
    predicted median survival (predicted survival at the median survival
    time), \code{time} is the predicted survival time (this last option
    only applies to partial plots, however).  For continuous variables,
    points are colored so that blue corresponds to events, whereas black
    points represent censored observations.
  
    Ensemble mortality values should be interpreted in terms of total
    number of deaths.  For example, if individual \code{i} has a
    mortality value of 100, then if all individuals were the same as
    \code{i}, we would expect to find 100 deaths on average in the data.
    If \code{type} is set to \code{rel.freq}, then mortality values are
    divided by an adjusted sample size, defined as the maximum of the
    sample size and the maximum mortality value.  The standardized
    mortality values no longer indicate total deaths, but instead
    reflect relative mortality.

    Partial plots are created when \code{partial}=TRUE.  Interpretation
    for these are different than marginal plots.  The y-value for a
    variable \eqn{X}, evaluated at \eqn{X=x}, is
    \deqn{
    \tilde{f}(x) = \frac{1}{n} \sum_{i=1}^n \hat{f}(x, x_{i,O}),
    }
    where \eqn{x_{i,O}} represents the value for all other variables
    other than \eqn{X} for individual \eqn{i} and \eqn{\hat{f}} is the
    predicted value. Generating partial plots can be very slow.
    Choosing a small value for \code{npts} can speed up computational
    times as this restricts the number of distinct \eqn{x} values used
    in computing \eqn{\tilde{f}}.

    For continuous variables, red points are used to indicate partial
    values and dashed red lines represent a lowess smoothed error bar of
    +/- two standard errors.  Black dashed line is the lowess estimate
    of the partial values.  For discrete variables, partial values are
    indicated using boxplots with whiskers extending out approximately
    two standard errors from the mean.  Standard errors are meant only
    to be a guide and should be interpreted with caution.

    Partial plots can be slow.  Setting type to \code{time} can greatly
    speed things up.  Setting \code{npts} to a smaller number should
    also be tried.
}
\author{
    Hemant Ishwaran \email{hemant.ishwaran@gmail.com} and
    Udaya B. Kogalur \email{ubk2101@columbia.edu}
}
\references{
    H. Ishwaran, U.B. Kogalur (2007).  Random survival forests for R,
    \emph{Rnews}, 7/2:25-31.

    J.H. Friedman (2001). Greedy function approximation: a gradient
    boosting machine, \emph{Ann. of Stat.}, 5:1189-1232.

    A. Liaw and M. Wiener (2002). Classification and regression by
    randomForest, \emph{R News}, 2:18-22.    
}
\seealso{
  \code{rsf},
  \code{predict.rsf}.
}
\examples{
#------------------------------------------------------------------------
# Some examples: veteran data.

data(veteran, package = "randomSurvivalForest") 
v.out <- rsf(Survrsf(time,status)~., veteran, forest = TRUE, ntree = 1000)
plot.variable(v.out, plots.per.page = 3)
plot.variable(v.out, plots.per.page = 2, predictorNames = c("trt", "karno", "age"))
plot.variable(v.out, type = "rel.freq", partial = TRUE, plots.per.page = 2, npred=3)

\dontrun{
#------------------------------------------------------------------------
# Fast partial plots using 'time' type.
# Top 8 predictors from PBC data.

data(pbc, package = "randomSurvivalForest") 
pbc.out <- rsf(Survrsf(days,status)~., pbc, ntree = 1000, forest = TRUE)
plot.variable(pbc.out, type = "time", partial = TRUE, npred=8)
}

}
\keyword{file}
