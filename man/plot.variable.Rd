\name{plot.variable}
\alias{plot.variable}
\alias{plot.variable.rsf}
\title{Plot Survival Effect of Variables}
\description{
  Plot of ensemble mortality, predicted survival, or predicted survival
  time against a given x-variable.  Users can select between marginal
  and partial plots.
}
\usage{
    \method{plot}{variable}(x, plots.per.page = 4, granule = 5, sorted = TRUE,
                  type = c("mort", "rel.freq", "surv", "time")[1],
                  partial = FALSE, predictorNames = NULL, npred = NULL,
                  npts = 25, subset = NULL, percentile = 50, ...)
}
\arguments{
    \item{x}{An object of class \code{(rsf, grow)} or \code{(rsf,
	predict)}.}
    \item{plots.per.page}{Integer value controlling page layout.}
    \item{granule}{Integer value controlling whether a plot for a
        specific variable should be given as a boxplot or scatter 
        plot. Larger values coerce boxplots.}
    \item{sorted}{Should variables be sorted by importance values (only
        applies if importance values are available)?}
      \item{type}{Select type of value to be plotted on the vertical
	axis. See details.}
    \item{partial}{Should partial plots be created?}
    \item{predictorNames}{Character vector of x-variables to be plotted.
        Default is all.}
    \item{npred}{Number of variables to be plotted.  Default is all.}
    \item{npts}{Maximum number of points used when generating partial
      plots for continuous variables.}
    \item{subset}{Indices indicating which rows of the predictor matrix
      to be used (note: this applies to the \emph{processed} predictor
      matrix, \code{predictors} of the object).  Default is to use all
      rows.}
    \item{percentile}{Percentile of follow up time used for plotting
      predicted survival.  See details below.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
    Either mortality, relative frequency of mortality, predicted
    survival, or predicted survival times are plotted on the vertical
    axis (y-value) against x-variables on the horizontal axis.  The
    choice of x-variables can be specified using \code{predictorNames}.
    The choice of y-value is controlled by \code{type}.  There are 4
    different choices: (1) "\code{mort}" is ensemble mortality; (2)
    "\code{rel.freq}" is standardized mortality; (3) "\code{surv}" is
    predicted survival at a given time point (the default is the median
    follow up time, but this can be set using the option
    \code{percentile}); (4) "\code{time}" is the predicted survival time
    (this last option only applies to partial plots).  For continuous
    variables, points are colored with blue corresponding to events, and
    black to censored observations.
  
    Ensemble mortality should be interpreted in terms of total number of
    deaths.  For example, if \code{i} has a mortality value of 100, then
    if all individuals were the same as \code{i}, the expected number of
    deaths would be 100.  If \code{type}="\code{rel.freq}", then
    mortality values are divided by an adjusted sample size, defined as
    the maximum of the sample size and the maximum mortality value.
    Standardized mortality values do not indicate total deaths, but
    rather relative mortality.

    Partial plots are created when \code{partial}=TRUE.  Interpretation
    for these are different than marginal plots.  The partial value for a
    variable \eqn{X}, evaluated at \eqn{X=x}, is
    \deqn{f(x) = \frac{1}{n} \sum_{i=1}^n \hat{f}(x, x_{i,o}),}
    where \eqn{\hat{f}} is the predicted value and where for each
    individual \eqn{i}, \eqn{x_{i,o}} represents the value for all other
    variables other than \eqn{X}.  For continuous variables, red points
    are used to indicate partial values and dashed red lines represent
    an error bar of +/- two standard errors.  A black dashed line
    indicates the lowess estimate of the partial values.  For discrete
    variables, partial values are indicated using boxplots with whiskers
    extending out approximately two standard errors from the mean.
    Standard errors are provided only as a guide and should be interpreted
    with caution.

    Partial plots can be slow.  Setting type="\code{time}" can improve
    matters.  Setting \code{npts} to a smaller number should also be
    tried.
    

    For competing risk analyses, plots correspond to unconditional values
    (i.e., they are non-event specific).  Use \code{competing.risk} for
    event-specific curves and for a more comprehensive analysis in such
    cases.
}
\author{
  Hemant Ishwaran
  Udaya B. Kogalur
}
\references{
    Ishwaran H., Kogalur U.B. (2007).  Random survival forests for R,
    \emph{Rnews}, 7(2):25-31.

    Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S.
    (2008).  Random survival forests, \emph{Ann. App. Statist.}, 2:841-860.

    Friedman J.H. (2001). Greedy function approximation: a gradient
    boosting machine, \emph{Ann. of Statist.}, 5:1189-1232.

    Liaw A. and Wiener M. (2002). Classification and regression by
    randomForest, \emph{R News}, 2:18-22.    
}
\seealso{
  \code{rsf},
  \code{predict.rsf}.
}
\examples{
# Some examples applied to veteran data.
data(veteran, package = "randomSurvivalForest") 
v.out <- rsf(Surv(time,status) ~ ., veteran, nsplit = 10, ntree = 1000)
plot.variable(v.out, plots.per.page = 3)
plot.variable(v.out, plots.per.page = 2,
        predictorNames = c("trt", "karno", "age"))
plot.variable(v.out, type = "surv", npred = 1, percentile = 50)
plot.variable(v.out, type = "rel.freq", partial = TRUE,
        plots.per.page = 2, npred=3)

\dontrun{
# Fast partial plots using 'time' type.
# Top 8 predictors from PBC data.
data(pbc, package = "randomSurvivalForest") 
pbc.out <- rsf(Surv(days,status) ~ ., pbc, ntree = 1000, nsplit = 3)
plot.variable(pbc.out, type = "time", partial = TRUE, npred=8)
}
}
\keyword{file}
