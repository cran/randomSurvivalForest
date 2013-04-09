\name{vimp}
\alias{vimp}
\title{VIMP for Single or Grouped Variables}
\description{
  Calculate variable importance (VIMP) for a single variable or group of
  variables.
}
\usage{
  vimp(object, predictorNames = NULL, subset = NULL,
       joint = TRUE, rough = FALSE,
       importance = c("randomsplit", "permute", "none")[1], 
       seed = NULL, do.trace = FALSE, ...)
}
\arguments{
    \item{object}{An object of class \code{(rsf, grow)} or \code{(rsf,
      forest)}.}
    \item{predictorNames}{Character vector of x-variable names to be
      considered.  If NULL (the default) all variables are used.  Only
      x-variables listed in the object predictor matrix will be used.}
    \item{subset}{Indices indicating which rows of the predictor
      matrix to be used (note: this applies to the \emph{object}
      predictor matrix, \code{predictors}).  Default is to use all rows.}
    \item{joint}{Should joint-VIMP or individual VIMP be calculated?
      See details below.}
    \item{rough}{Should fast approximation be used?}
    \item{importance}{Type of VIMP.}
    \item{seed}{Seed (negative integer) for random number generator.}
    \item{do.trace}{Should trace output be enabled? 
      Integer values can also be passed.  A positive value
      causes output to be printed each \code{do.trace} iteration.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
    Using a previously grown forest, and restricting the data to that
    indicated by \code{subset}, calculate the VIMP for variables listed
    in \code{predictorNames}.  If \code{joint}=TRUE, then joint-VIMP for
    the group of variables is calculated.  This equals the amount that
    prediction error changes when the group of variables are
    simultaneously "noised-up" (perturbed).  If \code{joint}=FALSE, the VIMP for
    each variable is calculated separately.
    
    Depending upon the setting for \code{importance}, VIMP is determined
    by using either random daughter assignment ("\code{randomsplit}") or
    random permutation ("\code{permute}") noising-up of the variable(s).
    A NULL value is returned if \code{importance}="\code{none}".

    For competing risk data, VIMP and error rates are given for the
    ensemble CHF and the conditional CHF (CCHF) for each event type.
    For more details see Ishwaran et al. (2010).
}
\value{
    A list with the following components:
    \item{err.rate}{OOB error rate for the (unperturbed) ensemble
      restricted to the subsetted data.}
    \item{err.perturb.rate}{OOB error rate for the perturbed ensemble
      restricted to the subsetted data.   Dimension depends upon the
      option \code{joint}.}
    \item{importance}{Variable importance (VIMP). Dimension depends
      upon the option \code{joint}.} 
}
\author{
  Hemant Ishwaran
  Udaya B. Kogalur
}
\references{
    Ishwaran H. (2007).  Variable importance in binary regression
    trees and forests,  \emph{Electronic J. Statist.}, 1:519-537.

    Ishwaran H., Kogalur U.B., Moore R.D., Gange S.J. and Lau
    B.M. (2010).  Random survival forests for competing risks.
}
\seealso{
  \code{find.interaction}.
}
\examples{
#------------------------------------------------------------------------
# Example of paired-VIMP. 
# Veteran data.

data(veteran, package = "randomSurvivalForest") 
v.out <- rsf(Survrsf(time,status) ~ ., veteran, ntree = 1000)
vimp(v.out, c("karno","celltype"))$importance

\dontrun{
#------------------------------------------------------------------------
# Individual VIMP for data restricted to events only.
# PBC data.

data(pbc, package = "randomSurvivalForest") 
rsf.out <- rsf(Surv(days,status) ~ ., pbc, ntree = 1000, nsplit = 3)
o.r <- rev(order(rsf.out$importance))
imp <- rsf.out$importance[o.r]
imp.events <- rep(0, length(imp))
events <- which(rsf.out$cens == 1)
imp.events <-
 vimp(rsf.out, names(imp), events, joint = FALSE)$importance
imp.all <- as.data.frame(cbind(imp.events = imp.events, imp = imp))
print(round(imp.all, 3))

#------------------------------------------------------------------------
# Estimate variability of VIMP in two ways (PBC data):
# (i)  Monte Carlo:  Estimates variability of the procedure
# (ii) Bootstrap:    Estimates statistical variability

data(pbc, package = "randomSurvivalForest") 
rsf.out <- rsf(Surv(days,status) ~ ., pbc, ntree = 1000, nsplit = 3)
o.r <- rev(order(rsf.out$importance))
imp.names <- names(rsf.out$importance[o.r])
subset.index <- 1:nrow(rsf.out$predictors)
imp.mc <- imp.boot <- NULL
for (k in 1:100) {
  cat("iteration:", k , "\n")
  imp.mc <-
    cbind(imp.mc, vimp(rsf.out, imp.names, joint = FALSE)$importance)
  imp.boot <-
    cbind(imp.boot, vimp(rsf.out, imp.names,
    subset = sample(subset.index, replace = TRUE), joint = FALSE)$importance)
}
imp.mc <- as.data.frame(cbind(imp.mean = apply(imp.mc, 1, mean),
                    imp.sd = apply(imp.mc, 1, sd)))
imp.boot <- as.data.frame(cbind(imp.mean = apply(imp.boot, 1, mean),
                    imp.sd = apply(imp.boot, 1, sd)))
print(round(imp.mc, 3))
print(round(imp.boot, 3))
}
}
\keyword{file}

