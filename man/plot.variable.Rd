\name{plot.variable}
\alias{plot.variable}
\title{Plot of Ensemble Survival Effect of Variables}
\description{
  Plot of ensemble mortality for each variable.  Users can select
  between marginal and partial plots.
}
\usage{
    plot.variable(x,
                  plots.per.page = 4,
                  granule = 5,
                  sort = TRUE,
                  rel.freq = FALSE,
                  partial = FALSE,
                  predictorNames = NULL,
                  n.pred = NULL,
                  n.pts = 25,
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
    \item{sort}{Should variables be sorted by importance values (only
        applies if importance values are available)?  Default is TRUE.}
    \item{rel.freq}{Should y-axis be in terms of number of deaths
	(the default), or should relative frequencies be used?} 
    \item{partial}{Logical. Should partial plots be created?  Default
        is FALSE.}
    \item{predictorNames}{Character vector of variable names.  Only these
        variables will be plotted.  Default is all.}
    \item{n.pred}{Number of variables to be plotted (only applies when
	\code{predictorNames}=NULL).  Default is all.}
    \item{n.pts}{Maximum number of points used when generating partial
        plots for continuous variables.}
    \item{subset}{Logical vector indicating which records to use.
        If NULL, then all records used.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
    Ensemble mortality is plotted against the value of a variable.
    Ensemble mortality values (vertical axis) should be interpreted in
    terms of total number of deaths.  For example, if individual
    \code{i} has a mortality value of 100, then if all individuals were
    the same as \code{i}, we would expect to find 100 deaths on average
    in the data.  If \code{rel.freq}=TRUE, then mortality values are
    divided by an adjusted sample size, defined as the maximum of the
    sample size and the maximum mortality value.  The standardized
    mortality values no longer indicate total deaths, but instead
    reflect relative mortality.

    The default is to create marginal plots.  Thus, each point
    represents the estimated mortality of an individual \code{i} against
    the value for \code{i}'s variable.  For continuous variables, points
    are colored so that blue corresponds to events, whereas black points
    represent censored observations.

    Partial plots are created when \code{partial}=TRUE.  In this
    case, the mortality value being plotted for a variable \eqn{X}
    evaluated at \eqn{X=x} is 
    \deqn{
    \tilde{f}(x) = \frac{1}{n} \sum_{i=1}^n \hat{f}(x, x_{i,O}),
    }
    where \eqn{x_{i,O}} represents the value for all other variables
    other than \eqn{X} for individual \eqn{i} and \eqn{\hat{f}} is the
    ensemble mortality predictor. Generating partial plots can be very
    slow.  Choosing a small value for \code{n.pts} can speed up
    computational times as this restricts the number of distinct \eqn{x}
    values used in computing \eqn{\tilde{f}}.

    For continuous variables, red points are used to indicate partial
    values and dashed red lines represent a lowess smoothed error bar of
    +/- two standard errors.  Black dashed line is the lowess estimate
    of the partial values.  For discrete variables, partial values are
    indicated using boxplots with whiskers extending out approximately
    two standard errors from the mean.  Standard errors are meant only
    to be a guide and should be interpreted with caution.
}
\author{
    Hemant Ishwaran \email{hemant.ishwaran@gmail.com} and
    Udaya B. Kogalur \email{ubk2101@columbia.edu}
}
\references{
    H. Ishwaran, Udaya B. Kogalur, Eugene H. Blackstone and Michael S. Lauer
    (2007).  Random Survival Forests. \emph{Cleveland Clinic Technical
    Report}.

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
  data(veteran, package = "randomSurvivalForest") 
  v.out <- rsf(Survrsf(time,status)~., veteran, forest = TRUE, ntree = 1000)
  plot.variable(v.out, plots.per.page = 3)
  plot.variable(v.out, plots.per.page = 2, predictorNames = c("trt", "karno", "age"))
  plot.variable(v.out, rel.freq = TRUE, partial = TRUE, plots.per.page = 2, n.pred=3)
}
\keyword{file}
