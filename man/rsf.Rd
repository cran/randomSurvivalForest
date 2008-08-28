\name{rsf.default}
\alias{rsf.default}
\alias{rsf}
\alias{randomSurvivalForest}
\title{Random Survival Forest Entry Point}
\description{
  Random Survival Forests (RSF) (Ishwaran, Kogalur, Blackstone and
  Lauer, 2008) is an extension of Breiman's Random Forests (Breiman,
  2001) to right-censored survival analysis settings.  A forest of
  survival trees is grown and used to estimate an ensemble cumulative
  hazard function (CHF).  Trees can be grown using different survival
  tree splitting rules.  An \dQuote{out-of-bag} estimate of Harrell's
  concordance index (Harrell, 1982) is provided for assessing prediction
  accuracy of the CHF.  Variable importance (VIMP) can be computed for
  single, as well as grouped variables, as a means to filter variables
  and to assess variable predictiveness.  RSF can be used to predict on
  test data. Missing data (x-variables, survival times, censoring
  indicators) can be imputed on both training and test data. Note this
  is the default generic method for the package.
}
\usage{
\method{rsf}{default}(formula,
    data = NULL,
    ntree = 1000,
    mtry = NULL,
    nodesize = NULL,
    splitrule = c("logrank", "conserve", "logrankscore", "random")[1],
    nsplit = 0,
    importance = c("randomsplit", "permute", "none")[1],
    big.data = FALSE,
    na.action = c("na.omit", "na.impute")[1],
    nimpute = 1,
    predictorWt = NULL,
    forest = FALSE,
    proximity = FALSE,
    varUsed = NULL,  
    seed = NULL,
    do.trace = FALSE,
    ...)
}
\arguments{
    \item{formula}{A symbolic description of the model to be fit. 
        Details for model specification are given below.}
    \item{data}{Data frame containing the data used in the formula.
        Missing values allowed.  See \code{na.action} for details.}
    \item{ntree}{Number of trees to grow.  This should not be set
        to a number too small, in order to ensure that every input
        row gets predicted at least a few times.}
    \item{mtry}{Number of variables randomly sampled at each split.
        The default is sqrt(\code{p}), where \code{p} equals the number
        of variables.}
    \item{nodesize}{Minimum number of deaths with unique survival
        times required for a terminal node.  Default is roughly
        min(3,round(0.632*\code{ndead})). Larger values cause smaller
	trees to be grown.}
    \item{splitrule}{Splitting rule used to grow trees. See details
      below.}
    \item{nsplit}{Non-negative integer value.  If non-zero, a random
      version of the specified tree splitting rule is implemented.  This can
      significantly increase speed.  See details below.}
    \item{importance}{Method used to compute variable importance.
	See details below.}    
    \item{big.data}{Logical.  Set this value to true when the number of
      variables \code{p} is \emph{very} large, or the data is very
      large.  See details below.}
    \item{na.action}{Action to be taken if the data contain NA's.
        Possible values are \code{na.omit} and \code{na.impute}.
        Default is \code{na.omit}, which removes the entire record if
        even one of its entries is NA (for x-variables this applies only
        to those specifically listed in 'formula').  The action
        \code{na.impute} implements a sophisticated tree imputation
        technique.  See details below.}
    \item{nimpute}{Number of iterations of missing data algorithm.}
    \item{predictorWt}{Vector of non-negative weights where entry
        \code{k}, after normalizing, is the probability of selecting
        variable \code{k} as a candidate for splitting.  Default is to
        use uniform weights.  Vector must be of dimension \code{p},
        where \code{p} equals the number of variables.}
    \item{forest}{Logical. Should the forest object be returned?  Used
        for prediction on new data.  Default is FALSE.}
    \item{proximity}{Logical. Should proximity measure between
        observations be calculated?  Creates an \code{n}x\code{n}
        matrix (which can be huge).  Default is FALSE.}
    \item{varUsed}{Analyzes which variables are used (split upon) in the topology
        of the forest.  Default is NULL.  Possible values are \code{all.trees},
	\code{by.tree}. See details below.}
    \item{seed}{Seed for random number generator.  Must be a negative
        integer (the R wrapper handles incorrectly set seed values).}
    \item{do.trace}{Logical. Should trace output be enabled? Default is
        FALSE. Integer values can also be passed.  A positive value
	causes output to be printed each \code{do.trace} iteration.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{  
  Four primary splitting rules are available for growing a survival
  forest.  The default rule, \code{logrank}, splits tree nodes by
  maximization of the log-rank test statistic (Segal, 1988; Leblanc and
  Crowley, 1993).  A second rule, \code{conserve}, splits nodes by
  finding daughters closest to the conservation of events principle (see
  Naftel, Blackstone and Turner, 1985).  A third rule,
  \code{logrankscore}, uses a standardized log-rank statistic (Hothorn
  and Lausen, 2003).  A fourth rule, \code{random}, implements pure
  random splitting.  For each node, a variable is randomly selected from
  a random set of \code{mtry} variables and the node is split using a
  random split point (Lin and Jeon, 2006).

  A random version of the \code{logrank}, \code{conserve} and
  \code{logrankscore} splitting rules can be invoked using
  \code{nsplit}.  If \code{nsplit} is set to a non-zero positive
  integer, then a maximum of \code{nsplit} split points are chosen
  randomly for each of the \code{mtry} variables within a node (this is
  in contrast to deterministic splitting where all possible split points
  for each of the \code{mtry} variables are considered).  The splitting
  rule is applied to these random split points and the node is split on
  that variable and random split point maximizing survival difference
  (as measured by the splitting rule).  Note that \code{nsplit} has no
  effect if the splitting rule is \code{random}.

  In terms of performance, a detailed study carried out by Ishwaran et
  al. (2008) found \code{logrank} and \code{logrankscore} to be the most
  accurate in terms of prediction error, followed by \code{conserve}.
  Setting \code{nsplit}=1 and using \code{logrank} splitting gave
  performance close to \code{logrank}, but with significantly shorter
  computational times.  Accuracy can be good further improved without
  overly compromising speed by using larger values of \code{nsplit}.

  In addition to random splitting, computation times for \emph{very}
  large data sets can be improved by discretizing continuous variables
  and/or the observed survival times.  Discretization does not have to
  be overly granular for substantial gains to be seen.  Users may also
  consider setting \code{big.data}=TRUE for data with a large number of
  variables.  This bypasses the large overhead R needs to create design
  matrices and parse formula.  Be aware, however, that variables are not
  processed and are interpreted \emph{as is} when this option is turned
  on.  Think of the data frame as containing time and censoring
  information and the rest of the data as the pre-processed design
  matrix when this option is on.  In particular, transformations used in
  the formula (such as logs etc.) are ignored.

  A typical RSF formula has the form \code{Survrsf(time, censoring) ~
  terms}, where \code{time} is survival time and \code{censoring} is a
  binary censoring indicator.  Note that censoring must be coded as
  0=censored and 1=death (event) and \code{time} must be strictly
  positive.

  Variables encoded as factors are treated as such.  If the factor is
  ordered, then splits are similar to real valued variables.  If the
  factor is unordered, a split will move a subset of the levels in the
  parent node to the left daughter, and the complementary subset to the
  right daughter.  All possible complementary pairs are considered and
  apply to factors with an unlimited number of levels.  However, there
  is an optimization check to ensure that the number of splits attempted
  is not greater than the number of cases in a node (this internal check
  will override the nsplit value in random splitting mode if nsplit is
  large enough).  Note that when predicting on test data involving
  factors, the factor labels in the test data must be the same as in the
  grow (training) data.  Consider setting labels that are unique in the
  test data to missing.
 
  Other than factors, all other x-variables are coerced and treated as
  being real valued.
  
  Variable importance (VIMP) is computed similar to Breiman (2001),
  although there are two ways to perturb a variable to determine its
  VIMP: \code{randomsplit}, \code{permute}.  The default method is
  \code{randomsplit} which works as follows.  To calculate VIMP for a
  variable \code{x}, out-of-bag (OOB) cases are dropped down the
  bootstrap (in-bag) survival tree.  A case is assigned a daughter node
  randomly whenever an \code{x}-split is encountered.  An OOB ensemble
  cumulative hazard function (CHF) is computed from the forest of such
  trees and its OOB error rate calculated.  The VIMP for \code{x} is the
  difference between this and the OOB error rate for the original forest
  (without random node assignment using \code{x}).  If \code{permute} is
  used, then \code{x} is randomly permuted in OOB data and dropped down
  the in-bag tree.  See Ishwaran et al. (2008) for further details.
  
  Prediction error is measured by 1-C, where C is Harrell's concordance
  index.  Prediction error is between 0 and 1, and measures how well the
  ensemble correctly ranks (classifies) any two individuals in terms of
  survival.  A value of 0.5 is no better than random guessing. A value
  of 0 is perfect.  Because VIMP is based on the concordance index, VIMP
  indicates how much misclassification increases, or decreases, for a
  new test case if a given variable were not available for that case
  (given that the forest was grown using that variable).
  
  Setting \code{na.action} to \code{na.impute} implements a tree
  imputation method whereby missing data (x-variables or outcomes) are
  imputed dynamically as a tree is grown by randomly sampling from the
  distribution within the current node (Ishwaran et al. 2008).  OOB data
  is not used in imputation to avoid biasing prediction error and VIMP
  estimates.  Final imputation for integer valued variables and
  censoring indicators use a maximal class rule, whereas continuous
  variables and survival time use a mean rule.  Records in which all
  outcome and x-variable information are missing are removed.  Variables
  having all missing values are removed.  The algorithm can be iterated
  by setting \code{nimpute} to a positive integer greater than 1.  A few
  iterations should be used in heavy missing data settings to improve
  accuracy of imputed values (see Ishwaran et al., 2008).  Note if the
  algorithm is iterated, a side effect is that missing values in
  returned objects \code{predictors}, \code{time} and \code{cens} are
  replaced by imputed values.  Further, imputed objects such as
  \code{imputedData} are set to NULL.

  If \code{varUsed}=\code{all.trees}, a vector of size \code{p} is
  returned.  Each element contains a count of the number of times a
  split has occurred on this variable.  If
  \code{varUsed}=\code{by.tree}, a matrix of size \code{ntree}x\code{p}
  is returned.  Each element [i][j] contains a count of the number of
  times a split has occurred on variable [j] in tree [i].
}   
\value{
     An object of class \code{(rsf, grow)}, which is a list with the
     following components:
     \item{call}{The original call to \code{rsf}.}
     \item{formula}{The formula used in the call.}
     \item{n}{Sample size of the data (depends upon NA's, see \code{na.action}).}
     \item{ndead}{Number of deaths.}
     \item{ntree}{Number of trees grown.}
     \item{mtry}{Number of variables randomly selected for splitting at each node.}
     \item{nodesize}{Minimum size of terminal nodes.}
     \item{splitrule}{Splitting rule used.}
     \item{nsplit}{Number of randomly selected split points.}
     \item{time}{Vector of length \code{n} of survival times.}
     \item{cens}{Vector of length \code{n} of censoring information (0=censored, 1=death).}
     \item{timeInterest}{Sorted unique event times.  Ensemble values are given
         for these time points only.}
     \item{predictorNames}{A character vector of the variable names used in
         growing the forest.}
     \item{predictorWt}{Vector of non-negative weights used for randomly
         sampling variables for splitting.}
     \item{predictors}{Data frame comprising x-variables used to grow the forest.}
     \item{ensemble}{A matrix of the bootstrap ensemble CHF with each row
         corresponding to an individual's CHF evaluated at each of the
         time points in \code{timeInterest}.}
     \item{oob.ensemble}{Same as \code{ensemble}, but based on the OOB CHF.}
     \item{mortality}{A vector of length \code{n}, with each value
         containing the bootstrap ensemble mortality for an
         individual in the data.  Ensemble mortality values should
         be interpreted in terms of total number of deaths.}
     \item{oob.mortality}{Same as \code{mortality}, but based on \code{oob.ensemble}.}
     \item{err.rate}{Vector of length \code{ntree} containing OOB error
         rates for the ensemble, with the b-th element being the error
         rate for the ensemble formed using the first b trees.
         Error rates are measured using 1-C, where C is Harrell's
	 concordance index.}
     \item{leaf.count}{Number of terminal nodes for each tree in the
         forest. Vector of length \code{ntree}.  A value of zero indicates
         a rejected tree (sometimes occurs when imputing missing data).
         Values of one indicate tree stumps.}
     \item{importance}{VIMP for each variable.}
     \item{forest}{If \code{forest}=TRUE, the forest object is returned.
         This object can then be used for prediction with new test data
	 sets.}
     \item{proximity}{If \code{proximity}=TRUE, a matrix of dimension
         \code{n}x\code{n} recording the frequency pairs of data points
         occur within the same terminal node. Value returned is a
         vector of the lower diagonal of the matrix.  Use
         \code{plot.proximity()} to extract this information.}
     \item{varUsed}{Count of the number of times a variable is used in
         growing the forest.  Can be a vector, matrix, or NULL.}
     \item{imputedIndv}{Vector of indices for cases with missing
         values. Can be NULL.}
     \item{imputedData}{Data frame comprising imputed data.  First two
        columns are censoring and survival time, respectively.
        Remaining columns are the x-variables.  Row i contains imputed
        outcomes and x-variables for row \code{imputedIndv}[i] of
        \code{predictors}.  Can be NULL.}
}
\note{
    The key deliverable is the matrix \code{ensemble} containing the
    bootstrap ensemble CHF function for each individual evaluated at a
    set of distinct time points (an OOB ensemble, \code{oob.ensemble},
    is also returned).  The vector \code{mortality} (likewise
    \code{oob.mortality}) is a weighted sum over the columns of
    \code{ensemble}, weighted by the number of individuals at risk at
    the different time points.  Entry \code{i} of the vector represents
    the estimated total mortality of individual \code{i} in terms of
    total number of deaths.  In other words, if \code{i} has a mortality
    value of 100, then if all individuals had the same x-values as
    \code{i}, there would be on average 100 deaths in the dataset.
    
    Different R wrappers are provided with the package to aid in
    interpreting the ensemble.
}
\author{
    Hemant Ishwaran \email{hemant.ishwaran@gmail.com} and
    Udaya B. Kogalur \email{ubk2101@columbia.edu}
}
\references{
    L. Breiman (2001). Random forests, \emph{Machine Learning}, 45:5-32.

    F.E. Harrell et al. (1982).  Evaluating the yield of medical tests,
    \emph{J. Amer. Med. Assoc.}, 247:2543-2546.

    T. Hothorn and B. Lausen (2003). On the exact distribution of maximally selected
    rank statistics, \emph{Comp. Stat. Data Anal.}, 43:121-137.    

    H. Ishwaran, U.B. Kogalur, E.H. Blackstone and M.S. Lauer
    (2008).  Random survival forests, \emph{To appear in Ann. App.
    Statist.}.

    H. Ishwaran, U.B. Kogalur (2007).  Random survival forests for R,
    \emph{Rnews}, 7/2:25-31.

    H. Ishwaran (2007).  Variable importance in binary regression
    trees and forests,  \emph{Electronic J. Statist.}, 1:519-537.

    M. LeBlanc and J. Crowley (1993).  Survival trees by goodness of split,
    \emph{J. Amer. Stat. Assoc.}, 88:457-467.
    
    A. Liaw and M. Wiener (2002). Classification and regression by randomForest,
    \emph{R News}, 2:18-22.

    Y. Lin and Y. Jeon (2006).  Random forests and adaptive nearest
    neighbors, \emph{J. Amer. Stat. Assoc.}, 101:578-590.

    D.C. Naftel, E.H. Blackstone and M.E. Turner (1985). Conservation of
    events, unpublished notes.

    M. R. Segal. (1988).  Regression trees for censored data, 
    \emph{Biometrics},  44:35-47.
}
\seealso{
  \code{plot.ensemble},
  \code{plot.variable},
  \code{plot.error},
  \code{plot.proximity},
  \code{predict.rsf},
  \code{print.rsf},
  \code{find.interaction},
  \code{pmml2rsf},
  \code{rsf2pmml},
  \code{Survrsf}.
}
\examples{
#------------------------------------------------------------
# Example 1:  Veteran's Administration lung cancer trial from
# Kalbfleisch & Prentice.  Randomized trial of two treatment
# regimens for lung cancer.  Minimal argument call.  Print
# results, then plot error rate and importance values.

data(veteran, package = "randomSurvivalForest")
veteran.out <- rsf(Survrsf(time, status)~., data = veteran)
print(veteran.out)
plot(veteran.out)

#------------------------------------------------------------
# Example 2:  Richer argument call (veteran data).
# Forest is saved by setting 'forest' option to true
# (see 'rsf.predict' for more details about prediction).
# Coerce variable 'celltype' as a factor, and karnofsky score
# as an ordered factor to illustrate factor useage in RSF.
# Use random splitting with 'nsplit'.
# Use 'varUsed' option.

data(veteran, package = "randomSurvivalForest")
veteran.f <- as.formula(Survrsf(time, status)~.)
veteran$celltype <- factor(veteran$celltype,
    labels=c("squamous", "smallcell",  "adeno",  "large"))
veteran$karno <- factor(veteran$karno, ordered = TRUE)
ntree <- 200
mtry <- 2
nodesize <- 3
splitrule <- "logrank"
nsplit <- 10
varUsed <- "by.tree"
forest <- TRUE
proximity <- TRUE
do.trace <- 1
veteran2.out <- rsf(veteran.f, veteran, ntree,
       mtry, nodesize, splitrule, nsplit,
       varUsed = varUsed, forest = forest, 
       proximity = proximity, do.trace = do.trace)
print(veteran2.out)
plot.proximity(veteran2.out)

# Take a peek at the forest ...
head(veteran2.out$forest$nativeArray)

# Average number of times a variable was split on.
apply(veteran2.out$varUsed,2,mean)

# Partial plot of top variable.
plot.variable(veteran2.out, partial = TRUE, npred=1)


\dontrun{
#------------------------------------------------------------
# Example 3:  Veteran data (again).
# Consider Karnofsky performance score. Compare to Kaplan-Meier.
# Assumes "survival" library is loaded.

if (library("survival", logical.return = TRUE))
{
        data(veteran, package = "randomSurvivalForest")
        veteran3.out <- rsf(Survrsf(time, status)~karno,
                       veteran,
                       ntree = 1000)
        plot.ensemble(veteran3.out)
        par(mfrow = c(1,1))
        plot(survfit(Surv(time, status)~karno, data = veteran))
}

#------------------------------------------------------------
# Example 4:  Primary biliary cirrhosis (PBC) of the liver.
# Data found in Appendix D.1 of Fleming and Harrington, Counting
# Processes and Survival Analysis, Wiley, 1991 (modified so
# that age is in days and sex and stage variables are not
# missing for observations 313-418).  

data(pbc, package = "randomSurvivalForest") 
pbc.out <- rsf(Survrsf(days,status)~., pbc, ntree = 1000)
print(pbc.out)

#------------------------------------------------------------
# Example 5:  Same as Example 4, but with imputation for
# missing values.

data(pbc, package = "randomSurvivalForest") 
pbc2.out <- rsf(Survrsf(days,status)~., pbc, ntree = 1000,
                na.action="na.impute")
# summary of analysis
print(pbc2.out)
# Combine original data + imputed data.
pbc.imputed.data <- cbind(status=pbc2.out$cens, days=pbc2.out$time,
                          pbc2.out$predictors)
pbc.imputed.data[pbc2.out$imputedIndv,] <- pbc2.out$imputedData
tail(pbc)
tail(pbc.imputed.data)
# Iterate the missing data algorithm.
# Use logrank random splitting (with nsplit=5) to increase speed.
# Use trace to track algorithm in detail.
# Note that a side effect of iterating is that the original data
# are replaced by imputed values.
pbc3.out <- rsf(Survrsf(days,status)~., pbc, ntree = 1000, nsplit=5, 
         na.action="na.impute", nimpute=3, do.trace = TRUE)
pbc.iterate.imputed.data <- cbind(status=pbc3.out$cens,
         days=pbc3.out$time, pbc3.out$predictors)

#------------------------------------------------------------
# Example 6:  Compare Cox regression to RSF (PBC data).
# Compute OOB estimate of Harrell's concordance 
# index for Cox regression using B = 100 bootstrap draws.
# Assumes "Hmisc" and "survival" libraries are loaded. 

if (library("survival", logical.return = TRUE) 
    & library("Hmisc", logical.return = TRUE))
{
  data(pbc, package = "randomSurvivalForest")
  pbc3.out <- rsf(Survrsf(days,status)~., pbc, mtry = 2, ntree = 1000)
  B <- 100 
  cox.err <- rep(NA, B) 
  cox.f <- as.formula(Surv(days,status)~.)  
  pbc.data <- pbc[apply(is.na(pbc), 1, sum) == 0,] ##remove NA's 
  cat("Out-of-bag Cox Analysis ...", "\n")
  for (b in 1:B) {
    cat("Cox bootstrap:", b, "\n") 
    bag.sample <- sample(1:nrow(pbc.data),
                         nrow(pbc.data),
                         replace = TRUE) 
    oob.sample <- setdiff(1:nrow(pbc.data), bag.sample)
    train <- pbc.data[bag.sample,]
    test <- pbc.data[oob.sample,]
    cox.out <- coxph(cox.f, train)
    cox.out <- tryCatch({coxph(cox.f, train)}, error=function(ex){NULL})
    if (is.list(cox.out)) {
      cox.predict <- predict(cox.out, test)
      cox.err[b] <- rcorr.cens(cox.predict, 
              Surv(pbc.data$days[oob.sample],
              pbc.data$status[oob.sample]))[1]
    }
  }
  cat("Error rates:", "\n")
  cat("Random Survival Forests:", pbc3.out$err.rate[pbc3.out$ntree], "\n")
  cat("         Cox Regression:", mean(cox.err, na.rm = TRUE), "\n")
}
}
}
\keyword{survival}
\keyword{tree}
