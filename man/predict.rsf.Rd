\name{predict.rsf}
\alias{predict.rsf}
\title{Random Survival Forest Prediction}
\description{
    Prediction on test data using Random Survival Forests.  
}
\usage{
  \method{predict}{rsf}(object = NULL, test = NULL,
              importance = c("randomsplit", "permute", "none")[1],
              na.action = c("na.omit", "na.impute")[1],
              outcome = c("train", "test")[1],
              proximity = FALSE, split.depth = FALSE, seed = NULL,
              do.trace = FALSE, ...)
}
\arguments{
    \item{object}{An object of class \code{(rsf, grow)} or \code{(rsf,
	forest)}.}
    \item{test}{Data frame containing test data. Missing values allowed.}
    \item{importance}{Method used to compute variable importance (VIMP).  Only
        applies when test data contains survival outcomes.}
    \item{na.action}{Action to be taken if the data contains NA's.  Possible
        values are "\code{na.omit}", which removes the entire record if
        even one of its entries is NA, and "\code{na.impute}", which
        imputes the test data.  See details below.}
    \item{outcome}{Data frame used in calculating the ensemble.  By
      default this is always the training data, but see details below.}
    \item{proximity}{Should proximity measure between test observations
      be calculated? Can be large.}
    \item{split.depth}{Return minimal depth for each variable for each
        test set individual?}
    \item{seed}{Seed (negative integer) for random number generator.}
    \item{do.trace}{Logical. Should trace output be enabled? Integer
      values can also be passed.  A positive value causes output to
      be printed each \code{do.trace} iteration.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
   The test data is dropped down the grow-forest (i.e., the forest grown
   from the training data) yielding the ensemble cumulative hazard
   function (CHF) for each individual in the test data evaluated at each
   unique death time point from the grow data.  If survival outcome
   information is present in the test data, the overall error rate and
   VIMP for each variable is also returned.  Setting
   \code{na.action}="\code{na.impute}" imputes missing test data
   (x-variables or outcomes).  Imputation uses the grow-forest such that
   only training data is used when imputing test data to avoid biasing
   error rates and VIMP (Ishwaran et al. 2008).

  For competing risks, the ensemble conditional CHF (CCHF) is computed
  for each event type in addition to the ensemble CHF.

  If \code{outcome}="\code{test}", the ensemble is calculated by
  specifically using survival information from the test data (survival
  information must be present).  In this case, the terminal nodes from
  the grow-forest are recalculated using survival data from the test
  set.  This yields a modified predictor in which the topology of the
  forest is based solely on the training data, but where the predicted
  value is based on the test data.  Error rates and VIMP are calculated
  by bootstrapping the test data and using out-of-bagging to ensure
  unbiased estimates.  See Examples 2 and 3 below for illustration.
}    
\value{
     An object of class \code{(rsf, predict)}, which is a list with the
     following components:
     \item{call}{The original grow call to \code{rsf}.}
     \item{forest}{The grow forest.}
     \item{ntree}{Number of trees in the grow forest.}
     \item{leaf.count}{Number of terminal nodes for each tree in the
         grow forest.  Vector of length \code{ntree}.}
     \item{timeInterest}{Sorted unique event times from grow (training) data.
         Ensemble values given for these time points only.}
     \item{n}{Sample size of test data (depends upon NA's, see \code{na.action}).}
     \item{ndead}{Number of deaths in test data (can be NULL).}
     \item{time}{Vector of survival times from test data (can be NULL).}
     \item{cens}{Vector of censoring indicators from test data (can be NULL).}
     \item{predictorNames}{Character vector of variable names.}
     \item{predictors}{Data frame comprising x-variables used for
       prediction.}
     \item{ensemble}{Matrix containing the ensemble CHF for the test
         data. Each row corresponds to the CHF for an individual in the
         test set evaluated at each of the time points in
         \code{timeInterest}.  For competing risks, a 3-D array where
         the 3rd dimension is for the ensemble CHF and each of the
         CCHFs, respectively.}
     \item{poe}{Matrix containing the ensemble probability
         of an event (POE) for each test set individual: used to estimate the
         CIF.  Rows correspond to each of the event types.  Applies only
	 to competing risk data.  NULL otherwise.}
     \item{mortality}{Vector containing ensemble mortality for each
         individual in the test data.  Ensemble mortality should
         be interpreted in terms of total number of training deaths if
         \code{outcome}="\code{train}".}
     \item{err.rate}{Vector of length \code{ntree} of the test-set error
         rate.  For competing risks, a matrix of test-set errors with
         rows corresponding to the ensemble CHF and each of the CCHFs,
         respectively. Can be NULL.  If \code{outcome}="\code{test}" 
         the test-set error is non-cumulative (i.e., it is for the full
	 forest).}
     \item{importance}{VIMP of each variable in the test data.  For
         competing risks, a matrix with rows corresponding to the
         ensemble CHF and each of the CCHFs, respectively.  Can be
         NULL.}
     \item{proximity}{If \code{proximity}=TRUE, a matrix recording
	 proximity of the inputs from test data is computed.  Value
	 returned is a vector of the lower diagonal of the matrix.  Use
	 \code{plot.proximity} to extract this information.}
     \item{imputedIndv}{Vector of indices of records in test data with
         missing values.  Can be NULL.}
     \item{imputedData}{Data frame containing the imputed test data.  First
         two columns are censoring and survival time, respectively. The
         remaining columns are the x-variables.  Row i contains imputed
         outcomes and x-variables for row \code{imputedIndv}[i] of
        \code{predictors}.  Can be NULL.}
      \item{splitDepth}{Matrix where entry [i][j] is the mean minimal
	  depth for variable [j] for case [i] in the test data. Used for
	 variable selection (see \code{max.subtree}).  Can be NULL.}
}
\author{
  Hemant Ishwaran
  Udaya B. Kogalur
}
\references{
    Breiman L. (2001). Random forests, \emph{Machine Learning}, 45:5-32.

    Ishwaran H., Kogalur U.B. (2007).  Random survival forests for R,
    \emph{Rnews}, 7(2):25-31.

    Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S.
    (2008).  Random survival forests, \emph{Ann. App. Statist.}, 2:841-860.
    
    Ishwaran H., Kogalur U.B., Moore R.D., Gange S.J. and Lau
    B.M. (2010).  Random survival forests for competing risks.
}	
\seealso{
  \code{rsf}.
}
\examples{
#------------------------------------------------------------
# Example 1:  Typical call (veteran data)

data(veteran, package = "randomSurvivalForest")
pt.train <- sample(1:nrow(veteran), round(nrow(veteran)*0.80))
veteran.out <- rsf(Surv(time, status) ~ ., data = veteran[pt.train , ])
veteran.pred <- predict(veteran.out, veteran[-pt.train , ])

\dontrun{
#------------------------------------------------------------
# Example 2:  Get out-of-bag error rate using the training
# data as test data (pbc example)

data(pbc, package = "randomSurvivalForest")
pbc.grow <- rsf(Surv(days, status) ~ ., pbc, nsplit = 3)
pbc.pred <- predict(pbc.grow, pbc, outcome = "test")
cat("GROW error rate  :", round(pbc.grow$err.rate[1000], 3))
cat("PRED error rate  :", round(pbc.pred$err.rate, 3))

#------------------------------------------------------------
# Example 3:  Verify reproducibility of forest (pbc data)

#primary call
data(pbc, package = "randomSurvivalForest")
pt.train <- sample(1:nrow(pbc), round(nrow(pbc)*0.50))
pbc.out <- rsf(Surv(days, status) ~ ., nsplit = 3, 
               data = pbc[pt.train, ])

#make separate predict calls using the outcome option
pbc.train <- predict(pbc.out, pbc[-pt.train, ], outcome = "train")
pbc.test <- predict(pbc.out, pbc[-pt.train, ], outcome = "test")

#check forest reproducibilility by comparing predicted survival curves
timeInterest <- pbc.out$timeInterest
surv.train <- exp(-pbc.train$ensemble)
surv.test <- exp(-pbc.test$ensemble)
matplot(timeInterest, t(surv.train - surv.test), type = "l")

#test reproducibility by repeating B times
#compute l1-difference in predicted survival
B <- 25
l1.valid <- rep(NA, B)
for (b in 1:B) {
 cat("Replication:", b, "\n")
 pt.train <- sample(1:nrow(pbc), round(nrow(pbc)*0.50))
 pbc.out <- rsf(Surv(days, status) ~ ., nsplit = 3, 
                 data = pbc[pt.train, ])
 surv.train <- exp(-predict(pbc.out, pbc[-pt.train, ],
                 outcome = "train")$ensemble)
 surv.test <- exp(-predict(pbc.out, pbc[-pt.train, ],
                 outcome = "test")$ensemble)
 l1.valid <-
   mean(apply(abs(surv.train - surv.test), 1, mean, na.rm = TRUE), na.rm = TRUE)
}
cat("l1-reproducibility:", round(mean(l1.valid, na.rm = TRUE), 3), "\n")
}
}
\keyword{survival}
\keyword{tree}
