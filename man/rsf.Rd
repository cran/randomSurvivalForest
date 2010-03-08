\name{rsf.default}
\alias{rsf.default}
\alias{rsf}
\alias{randomSurvivalForest}
\title{Random Survival Forest Entry Point}
\description{  
  Prediction and variable selection for right censored survival and
  competing risk data using Random Survival Forests (RSF) (Ishwaran,
  Kogalur, Blackstone and Lauer, 2008).  A random forest (Breiman, 2001)
  of survival trees is grown and used to estimate an ensemble cumulative
  hazard function (CHF) and ensemble conditional cumulative hazard function
  (CCHF) (the latter in the case of competing risks).  Trees can be
  grown using different survival tree splitting rules.  An
  \dQuote{out-of-bag} estimate of Harrell's concordance index (Harrell,
  1982) is provided for assessing prediction accuracy.  Variable
  importance (VIMP) for single, as well as grouped variables, can be
  used to filter variables and to assess variable predictiveness.
  Missing data (x-variables, survival times, censoring indicators) can
  be imputed on both training and test data. This is the default
  generic method for the package.
}
\usage{
\method{rsf}{default}(formula, data = NULL, ntree = 1000, mtry = NULL, nodesize = NULL,
    splitrule = NULL, nsplit = 0, importance = c("randomsplit", "permute", "none")[1], 
    big.data = FALSE, na.action = c("na.omit", "na.impute")[1], nimpute = 1, predictorWt = NULL,
    forest = FALSE, proximity = FALSE, varUsed = NULL, split.depth = FALSE, seed = NULL, 
    do.trace = FALSE, ...)
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
    \item{nodesize}{Minimum number of deaths with unique survival times
        required for a terminal node.  Default is approximately 3 for
        right-censoring and 6 for competing risk data. Larger values create
        smaller trees.}
    \item{splitrule}{Splitting rule used to grow trees. See details
      below.}
    \item{nsplit}{Non-negative integer value.  If non-zero, the
      specified tree splitting rule is randomized which can
      significantly increase speed.  See details below.}
    \item{importance}{Method used to compute variable importance.
	See details below.}    
    \item{big.data}{Set this value to TRUE when the number of
        variables \code{p} is \emph{very} large, or the sample size is
        very large.  See details below.}
    \item{na.action}{Action to be taken if the data contains NA's.
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
    \item{forest}{Should the forest object be returned?  Used for
        prediction on new data.  Default is FALSE.}
    \item{proximity}{Should the proximity between observations be
        calculated?  Creates an \code{n}x\code{n} matrix (which can be
        huge).  Default is FALSE.}
    \item{varUsed}{Analyzes which variables are used (split upon) in the
        forest. Default is NULL.  Possible values are \code{all.trees} 
	and \code{by.tree}. See details below.}
    \item{split.depth}{Return minimal depth for each variable for each case.
        Default is FALSE. Used for variable selection: see details below.}
    \item{seed}{Seed for random number generator.  Must be a negative
        integer (the R wrapper handles incorrectly set seed values).}
    \item{do.trace}{Should trace output be enabled? Default is
        FALSE. Integer values can also be passed.  A positive value
	causes output to be printed each \code{do.trace} iteration.}
    \item{...}{Further arguments passed to or from other methods.}
}
\details{
  ---> Splitting Rules:
  
  Four primary splitting rules are available for growing a survival
  forest for right-censored data: \code{logrank}, \code{conserve},
  \code{logrankscore}, and \code{random}.
  
  The default rule, \code{logrank}, splits tree nodes by maximization of
  the log-rank test statistic (Segal, 1988; Leblanc and Crowley, 1993).
  \code{conserve}, splits nodes by finding daughters closest to the
  conservation of events principle (see Naftel, Blackstone and Turner,
  1985).  \code{logrankscore}, uses a standardized log-rank statistic
  (Hothorn and Lausen, 2003).  Finally, \code{random}, implements pure
  random splitting.  For each node, a variable is randomly selected from
  a random set of \code{mtry} variables and the node is split using a
  random split point (Cutler and Zhao, 2001; Lin and Jeon, 2006).  Note,
  however, that because random splitting promotes splits near the edges,
  node splitting can terminate early resulting in extremely unbalanced
  trees.  To correct this, the definition of \code{nodesize} is taken in
  this setting (and this setting alone) to equal the minimum number of
  unique deaths within a node required to split the node.

  A random version of the \code{logrank}, \code{conserve} and
  \code{logrankscore} splitting rules can be invoked using
  \code{nsplit}.  If \code{nsplit} is set to a non-zero positive
  integer, then a maximum of \code{nsplit} split points are chosen
  randomly for each of the \code{mtry} variables within a node (this is
  in contrast to deterministic splitting, i.e. \code{nsplit}=0, where
  all possible split points for each of the \code{mtry} variables are
  considered).  The splitting rule is applied to these random split
  points and the node is split on that variable and random split point
  maximizing survival difference (as measured by the splitting rule).
  
  In terms of performance, a detailed study carried out by Ishwaran et
  al. (2008) found \code{logrank} and \code{logrankscore} to be the most
  accurate in terms of prediction error, followed by \code{conserve}.
  Setting \code{nsplit}=1 and using \code{logrank} splitting gave
  performance close to \code{logrank}, but with significantly shorter
  computational times.  Accuracy can be further improved without overly
  compromising speed by using larger values of \code{nsplit}.  

  Trees tend to favor splits on continuous variables (Loh and Shih,
  1997), so it is good practice to use the \code{nsplit} option when the
  data contains a mix of continuous and discrete variables.  Using a
  reasonably small value, such as \code{nsplit}=3, mitigates bias.
  
  ---> Large Data Sets:
  
  Computation times for \emph{very} large data sets can be improved by
  discretizing continuous variables and/or the observed survival times;
  in addition to using random splitting.  Discretization does not have
  to be overly granular for substantial gains to be seen.  Users may
  also consider setting \code{big.data}=TRUE for data with a large
  number of variables.  This bypasses the large overhead R needs to
  create design matrices and parse formula.  Be aware, however, that
  variables are not processed and are interpreted \emph{as is} under
  this option.  Think of the data frame as containing time and censoring
  information and the rest of the data as the pre-processed design
  matrix.  In particular, transformations used in the formula (such as
  logs etc.) are ignored.

  
  ---> Formula:
  
  A typical RSF formula has the form \code{Survrsf(time, censoring) ~
  terms}, where \code{time} is survival time and \code{censoring} is a
  binary censoring indicator.  Censoring must be coded as a non-negative
  integer with 0 reserved for censoring and (usually) 1=death (event).
  Also, \code{time} must be strictly positive.  Other permissible
  formula are \code{Surv(time, censoring) ~ terms}.


  ---> Factors and Variable Types:
  
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
  test data to missing to avoid issues.
 
  Other than factors, all other x-variables are coerced and treated as
  being real valued.


  ---> Variable Importance:
  
  Variable importance (VIMP) is computed similar to Breiman (2001),
  although there are two ways to perturb a variable to determine its
  VIMP: \code{randomsplit}, \code{permute}.  The default method is
  \code{randomsplit} which works as follows.  To calculate VIMP for a
  variable \code{x}, out-of-bag (OOB) cases are dropped down the
  in-bag (bootstrap) survival tree.  A case is assigned a daughter node
  randomly whenever an \code{x}-split is encountered.  An OOB ensemble
  cumulative hazard function (CHF) is computed from the forest of such
  trees and its OOB error rate calculated.  The VIMP for \code{x} is the
  difference between this and the OOB error rate for the original forest
  (without random node assignment using \code{x}).  If \code{permute} is
  used, then \code{x} is randomly permuted in OOB data and dropped down
  the in-bag tree.  See Ishwaran et al. (2008) for further details.


  ---> Predition Error:
  
  Prediction error is measured by 1-C, where C is Harrell's concordance
  index.  Prediction error is between 0 and 1, and measures how well the
  ensemble correctly ranks (classifies) two random individuals in terms
  of survival.  A value of 0.5 is no better than random guessing. A
  value of 0 is perfect.  Because VIMP is based on the concordance
  index, VIMP indicates how much misclassification increases, or
  decreases, for a new test case if a given variable were not available
  for that case (given that the forest was grown using that variable).

  ---> Competing Risks:
  
  The implementation is similar to right-censoring but with the
  following caveats:

  (1) Censoring must be coded as a non-negative integer where 0
  indicates right-censoring and non-zero values indicate different event
  types.  While \code{0,1,2,..,J} is standard, events can be coded
  non-sequentially, although 0 must always be used for censoring.
    
  (2) The default splitting rule is \code{logrankCR}, a modified log-rank
  splitting rule tailored for competing risks.  Over-riding this by
  manually selecting any split rule other than \code{logrankCR} or
  \code{random} will result in a right-censored
  analysis in which all (non-censored) events are treated as if they are
  the same type (indeed, they will coerced as such).  Note that
  \code{nsplit} works as in right-censoring.

  (3) The ensemble (see below) is a 3-D array in which the 3rd dimension
  is reserved for the ensemble CHF and each of the J ensemble
  conditional CHFs (CCHFs).  The wrapper \code{competing.risk} can be
  used to process the ensemble and to generate event-specific cumulative
  incidence functions (CIF) and subsurvival functions (see Gray (1988)
  for background and definitions).

  (4) The cases within a terminal node are used to estimate both the
  unconditional survival function and the event-specific subsurvival
  functions and for this reason \code{nodesize} should generally be set
  larger than in right-censored data settings.
  
  
  ---> Missing Data and Imputation:

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
  algorithm is iterated, a side effect is that missing values in the
  returned objects \code{predictors}, \code{time} and \code{cens} are
  replaced by imputed values.  Further, imputed objects such as
  \code{imputedData} are set to NULL.  See the examples below.
  Also see the function \code{impute.rsf} for a fast impute interface.


  ---> Miscellanea:
  
  Setting \code{varUsed}=\code{all.trees} returns a vector of size
  \code{p} where each element is a count of the number of times a split
  occurred on a variable.  If \code{varUsed}=\code{by.tree}, a matrix of
  size \code{ntree}x\code{p} is returned.  Each element [i][j] is the
  count of the number of times a split occurred on variable [j] in tree
  [i].

  Setting \code{split.depth}=TRUE returns a matrix of size
  \code{n}x\code{p} where entry [i][j] is the mean minimal depth for
  variable [j] for case [i].  Used to select variables at the
  case-level.  See \code{max.subtree} for more details regarding minimal
  depth.
}   
\value{
     An object of class \code{(rsf, grow)} with the following components:
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
     \item{cens}{Vector of length \code{n} of censoring information (0=censored).  }
     \item{timeInterest}{Sorted unique event times.  Ensemble values are given
         for these time points only.}
     \item{predictorNames}{A character vector of the variable names used in
         growing the forest.}
     \item{predictorWt}{Vector of non-negative weights used for randomly
         sampling variables for splitting.}
     \item{predictors}{Data frame comprising x-variables used to grow
       the forest.}
     \item{ensemble}{Matrix for the in-bag ensemble CHF with each
         row corresponding to an individual's CHF evaluated at each of
         the time points in \code{timeInterest}.  For competing risks, a
         3-D array where the 3rd dimension is for the ensemble CHF and
         each of the CCHFs, respectively.} 
     \item{oob.ensemble}{Same as \code{ensemble}, but based on OOB
         data.}
     \item{poe}{Matrix for the in-bag ensemble probability
         of an event (poe) for each individual: used to estimate the
         CIF.  Rows correspond to each of the event types.  Applies only
	 to competing risk data.  NULL otherwise.}
     \item{oob.poe}{Same as \code{poe}, but based on OOB data.}
     \item{mortality}{A vector of length \code{n} for the in-bag
         ensemble mortality for an individual in the data.  Ensemble
         mortality values should be interpreted in terms of total number
         of deaths.}
     \item{oob.mortality}{Same as \code{mortality}, but based on \code{oob.ensemble}.}
     \item{err.rate}{Vector of length \code{ntree} containing OOB error
         rates for the ensemble, with the b-th element being the error
         rate for the ensemble formed using the first b trees.  Error
         rates are measured using 1-C, where C is Harrell's concordance
         index. For competing risks, a matrix with rows corresponding to
         the ensemble CHF and each of the CCHFs, respectively.}
     \item{leaf.count}{Number of terminal nodes for each tree in the
         forest. Vector of length \code{ntree}.  A value of zero indicates
         a rejected tree (sometimes occurs when imputing missing data).
         Values of one indicate tree stumps.}
       \item{importance}{Vector recording VIMP for each variable.  For
         competing risks, a matrix with rows corresponding to the ensemble CHF
         and each of the CCHFs, respectively.}
     \item{forest}{If \code{forest}=TRUE, the forest object is returned.
         This object can then be used for prediction with new test data
	 sets and is required for other R-wrappers.}
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
         outcomes and x-variables for row j of \code{predictors}, where
	 j=\code{imputedIndv}[i].  See the examples below. Can be NULL.}
     \item{splitDepth}{Matrix of size \code{n}x\code{p} where entry
         [i][j] is the mean minimal depth for variable [j] for case [i].
         Can be NULL.}
}
\note{
  The key deliverable is the matrix \code{ensemble} (and its OOB
  counterpart, \code{oob.ensemble}) containing the ensemble CHF function
  for each individual evaluated at a set of distinct time points.  The
  vector \code{mortality} (likewise \code{oob.mortality}) is a weighted
  sum over the columns of \code{ensemble}, weighted by the number of
  individuals at risk at the different time points.  Entry [i] of the
  vector represents the estimated total mortality of individual \code{i}
  in terms of total number of deaths.  In other words, if \code{i} has a
  mortality value of 100, then if all individuals had the same x-values
  as \code{i}, there would be on average 100 deaths in the dataset.
    
  Different R-wrappers are provided to aid in parsing the ensemble.
}
\author{
    Hemant Ishwaran \email{hemant.ishwaran@gmail.com} and
    Udaya B. Kogalur \email{kogalurshear@gmail.com}
}
\references{
    Breiman L. (2001). Random forests, \emph{Machine Learning}, 45:5-32.

    Cutler A. and Zhao G. (2001). Pert-Perfect random tree ensembles.
    \emph{Comp. Sci. Statist.}, 33: 490-497.

    Gray R.J. (1988).  A class of k-sample tests for comparing the
    cumulative incidence of a competing risk, \emph{Ann. Statist.},
    16: 1141-1154.

    Harrell F.E. et al. (1982).  Evaluating the yield of medical tests,
    \emph{J. Amer. Med. Assoc.}, 247:2543-2546.

    Hothorn T. and Lausen B. (2003). On the exact distribution of maximally selected
    rank statistics, \emph{Comp. Statist. Data Anal.}, 43:121-137.    

    Ishwaran H. (2007).  Variable importance in binary regression
    trees and forests,  \emph{Electronic J. Statist.}, 1:519-537.

    Ishwaran H., Kogalur U.B. (2007).  Random survival forests for R,
    \emph{Rnews}, 7/2:25-31.

    Ishwaran H., Kogalur U.B., Blackstone E.H. and Lauer M.S.
    (2008).  Random survival forests, \emph{Ann. App. Statist.}, 2:841-860.
    
    Ishwaran H., Kogalur U.B., Gorodeski E.Z, Minn A.J. and 
    Lauer M.S. (2010).  High-dimensional variable selection for survival
    data.  \emph{J. Amer. Statist. Assoc. (in press)}.

    Ishwaran H., Kogalur U.B., Chen X. and Minn A.J. (2010).  Random
    survival forests for high-dimensional data.  

    Ishwaran H., Kogalur U.B. and Lau B.M. (2010).  Random survival
    forests for competing risks.
    
    LeBlanc M. and Crowley J. (1993).  Survival trees by goodness of split,
    \emph{J. Amer. Statist. Assoc.}, 88:457-467.
    
    Liaw A. and Wiener M. (2002). Classification and regression by randomForest,
    \emph{R News}, 2:18-22.

    Lin Y. and Jeon Y. (2006).  Random forests and adaptive nearest
    neighbors, \emph{J. Amer. Statist. Assoc.}, 101:578-590.

    Loh W.-Y and Shih Y.-S (1997).  Split selection methods for
    classification trees, \emph{Statist. Sinica}, 7:815-840, 1997.

    Naftel N.C., Blackstone E.H. and Turner M.E. (1985). Conservation of
    events, unpublished notes.
  
    Segal M.R. (1988).  Regression trees for censored data, 
    \emph{Biometrics},  44:35-47.
}
\seealso{
  \code{competing.risk},
  \code{find.interaction},
  \code{impute.rsf},
  \code{max.subtree},
  \code{plot.ensemble},
  \code{plot.variable},
  \code{plot.error},
  \code{plot.proximity},
  \code{pmml2rsf},
  \code{predict.rsf},
  \code{print.rsf},
  \code{rsf2rfz},
  \code{rsf2pmml},
  \code{Survrsf},
  \code{varSel},
  \code{vimp}.
}
\examples{
#------------------------------------------------------------
# Example 1:  Veteran's Administration lung cancer data
# Randomized trial of two treatment regimens for lung cancer
# See Kalbfleisch & Prentice

data(veteran, package = "randomSurvivalForest")
veteran.out <- rsf(Survrsf(time, status) ~ ., data = veteran)
print(veteran.out)
plot(veteran.out)

#------------------------------------------------------------
# Example 2:  More detailed call (veteran data)
# Use random splitting by invoking 'nsplit' option
# Use 'varUsed' option
# Use forest=TRUE option

#read the data, set various options
data(veteran, package = "randomSurvivalForest")
veteran.f <- as.formula(Survrsf(time, status) ~ .)
ntree <- 200
nsplit <- 3
varUsed <- "by.tree"
forest <- TRUE

# coerce 'celltype' as a factor and 'karnofsky score'
# as an ordered factor to illustrate factor useage
veteran$celltype <- factor(veteran$celltype,
    labels=c("squamous", "smallcell",  "adeno",  "large"))
veteran$karno <- factor(veteran$karno, ordered = TRUE)

# grow call
veteran2.out <- rsf(veteran.f, veteran, 
       ntree = ntree, nsplit = nsplit,
       varUsed = varUsed, forest = forest) 

# plot of ensemble survival for a single individual
surv.ensb <- t(exp(-veteran2.out$oob.ensemble))
plot(veteran2.out$timeInterest, surv.ensb[, 1])

# take a peek at the forest
head(veteran2.out$forest$nativeArray)

# average number of times a variable was split 
apply(veteran2.out$varUsed, 2, mean)

# partial plot of top variable
plot.variable(veteran2.out, partial = TRUE, npred = 1)


\dontrun{
#------------------------------------------------------------
# Example 3:  Competing risks

# Follicular Cell Lymphoma
data(follic, package = "randomSurvivalForest")
follic.out <- rsf(Surv(time, status) ~ ., follic, nsplit = 3, ntree = 400)
print(follic.out)
plot(follic.out, sorted = FALSE)
competing.risk(follic.out)

# Hodgkin's disease
data(hd, package = "randomSurvivalForest")
hd.out <- rsf(Surv(time, status) ~ ., hd, nsplit = 3, ntree = 400)
print(hd.out)
plot(hd.out, sorted = FALSE)
competing.risk(hd.out)

#------------------------------------------------------------
# Example 4:  Primary biliary cirrhosis (PBC) of the liver
# See Appendix D.1 of Fleming and Harrington

data(pbc, package = "randomSurvivalForest") 
pbc.out <- rsf(Surv(days, status) ~ ., pbc, nsplit = 3)
print(pbc.out)

#------------------------------------------------------------
# Example 5:  Same as Example 4, but with data imputation
# Also see the R-wrapper "impute.rsf"

# rsf call with imputation
data(pbc, package = "randomSurvivalForest") 
pbc2.out <- rsf(Surv(days, status)~., pbc, 
                nsplit = 3, na.action="na.impute")
print(pbc2.out)

# here's a nice wrapper to combine original data + imputed data
combine.impute <- function(object) {
  imputed.data <- cbind(cens = object$cens,
                        time = object$time,
                        object$predictors)
  if (!is.null(object$imputedIndv)) {
    imputed.data[object$imputedIndv, ] <- object$imputedData
  }
  colnames(imputed.data)[c(2,1)] <- all.vars(object$formula)[1:2]
  imputed.data
}

# combine original data + imputed data
pbc.imputed.data <- combine.impute(pbc2.out)

# iterate the missing data algorithm
# compare to non-iterated algorithm
pbc3.out <- rsf(Surv(days, status)~., pbc, nsplit=5,  
         na.action="na.impute", nimpute = 3)
pbc.iterate.imputed.data <- combine.impute(pbc3.out)
tail(pbc.imputed.data)
tail(pbc.iterate.imputed.data)


#------------------------------------------------------------
# Example 6:  German breast cancer data
# Variable selection using minimal depth

data(breast, package = "randomSurvivalForest")
breast.out <- rsf(Surv(time, cens) ~ . , breast,
            forest = TRUE, nsplit = 3)

# use varSel to select variables
# see the help file of varSel for details/examples

breast.vs <- varSel(object=breast.out)

#------------------------------------------------------------
# Example 7:  Compare Cox regression to RSF using PBC data
# OOB estimate of C-index for Cox based on 100 bootstraps
# Assumes "Hmisc" and "survival" libraries are loaded

if (library("survival", logical.return = TRUE) 
    & library("Hmisc", logical.return = TRUE))
{
  data(pbc, package = "randomSurvivalForest")
  pbc3.out <- rsf(Surv(days, status)~., pbc, nsplit = 10, mtry = 2)
  B <- 100 
  cox.err <- rep(NA, B) 
  cox.f <- as.formula(Surv(days, status)~.)  
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
