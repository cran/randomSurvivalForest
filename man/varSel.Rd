\name{varSel}
\alias{varSel}
\title{Variable selection using RSF}
\description{
    Variable selection using minimal depth theory for trees.
}
\usage{
    varSel(formula = NULL,          
           data = NULL,
           object = NULL,      
           method = c("vh", "vhVIMP", "md")[3],
           ntree = (if (method == "md") 1000 else 500),
           mvars = (if (!is.null(data) & method != "md") min(1000, round(ncol(data)/5)) else NULL),
           mtry = (if (!is.null(data) & method == "md") max(sqrt(ncol(data)), ncol(data)/3) else NULL),
           nodesize = (if (method == "vh" | method == "vhVIMP") 2 else NULL),
           nsplit = 10,
           predictorWt = NULL,
           big.data = FALSE,
           na.action = c("na.omit", "na.impute")[1],
           do.trace = 0,       
           always.use = NULL,  
           nrep = 50,        
           K = 5,             
           nstep = 1,         
           verbose = TRUE,
           ...
           )

}
\arguments{
   \item{formula}{A symbolic description of the model to be fit. 
       Must be specified unless \code{object} is given.}
   \item{data}{Data frame containing the data used in the formula.
       Missing values allowed. Must be specified unless \code{object}
       is given.}
   \item{object}{An object of class \code{(rsf, grow)}.  Requires
       \code{forest}=TRUE in the original \code{rsf} call.
       Can be NULL.}
   \item{method}{Variable selection method: vh=variable hunting;
       vhVIMP=variable hunting with VIMP (variable importance);
       md=minimal depth.  See details below.}
   \item{ntree}{Number of trees to grow.}
   \item{mvars}{Number of randomly selected variables used in the
        variable hunting algorithm.}  
   \item{mtry}{Number of variables randomly sampled at each split.
        Should be large when the goal is variable selection.}
   \item{nodesize}{Minimum number of deaths with unique survival
       times required for a terminal node.  Should be small if number
       of variables is large.}
   \item{nsplit}{Non-negative integer value.  If non-zero, the
       specified tree splitting rule is randomized which significantly
       increases speed.}
  \item{predictorWt}{Vector of non-negative weights specifying
        the probability of selecting a variable for splitting.
        Must be of dimension equal to the number of variables.  Default
	(NULL) invokes a data-adaptive method.}
  \item{big.data}{Only set this value to TRUE when the sample size is
      very large.}
  \item{na.action}{Action to be taken if the data contains NA's.}
  \item{do.trace}{Should trace output be enabled? Default is
      FALSE. A positive integer value causes output to be printed
      each \code{do.trace} iteration.}
  \item{always.use}{Character vector of variable names to be always
      included in the model selection procedure and in the final
      selected model.}
  \item{nrep}{Number of Monte Carlo iterations of the variable hunting algorithm.}
  \item{K}{Integer value specifying the K-fold size used in the variable hunting
      algorithm.}
  \item{nstep}{Integer value controlling the step size used in the
    forward selection process of the variable hunting algorithm.
    Increasing this will encourage more variables to be selected.}
  \item{verbose}{Set to TRUE to get verbose output.}
  \item{...}{Further arguments passed to or from other methods.}
}
\details{
    Variable selection using minimal depth theory for trees (Ishwaran et
    al., 2009).  The option \code{method} allows for two different
    approaches: (1) minimal depth: uses all data and all variables
    simultaenously; and (2) variable hunting: uses K-fold Monte Carlo
    validation, random selection of variables, and regularized forward
    selection.

    ---> Minimal Depth variable selection (\code{method}="md")

    The maximal subtree (Ishwaran et al., 2009) for a variable \code{x}
    is the largest subtree whose root node splits on \code{x} (all
    parent nodes of \code{x}'s maximal subtree have nodes that split on
    variables other than \code{x}).  The minimal depth of a maximal
    subtree equals the shortest distance (the depth) from the root node
    to the parent node of the maximal subtree (zero is the smallest
    value possible).  The smaller the minimal depth, the more impact
    \code{x} has on prediction.

    Variables are selected using an adaptive threshold based on minimal
    depth (Ishwaran et al., 2009) coupled with minor supervision using
    VIMP.

    Set \code{mtry} to larger values when the number of variables is
    high.
    

    ---> Variable Hunting (\code{method}="vh" or \code{method}="vhVIMP")

    Variable hunting is used for problems where the number of variables
    is magnitudes larger than the sample size and the sample size is
    reasonably small.  Microarray data is a good example.

    Using training data from random K-fold subsampling, a forest is fit
    to a randomly selected set of variables of size \code{mvars}
    (variables are chosen with probability proportional to weights
    determined using an initial forest fit on the training data).  The
    subset of variables are ordered by increasing minimal depth and
    added sequentially (starting from a minimal model) until joint VIMP
    no longer increases (signifying the final model; Ishwaran et al.,
    2009).  A forest is refit with these variables and applied to test
    data to estimate prediction error and VIMP.  The process is repeated
    \code{nrep} times.  Final selected variables are the top P ranked
    variables, where P is the average model size and variables are
    ranked by average minimal depth.

    A rough rule for choosing \code{mvars} is to set it equal to some
    fraction of the number of variables.

    The same algorithm is used when \code{method}="vhVIMP", but
    variables are ordered using VIMP (including the final model).  This
    is faster, but not as accurate.

    If \code{method}="vh", and the number of variables is large, set
    \code{nsplit} to a fairly large number, such as 10, to ensure that
    tree splitting is not overly influenced by noisy variables.  

    ---> Miscellanea

    If \code{big.data}=TRUE, and variable hunting is used, the training
    data is chosen to be of size n/K, where n=sample size (i.e., the
    size of the training data is swapped with the test data).  This
    speeds up the algorithm.  Increasing K also helps.
    
    For efficiency, transformations used in the formula (such as logs
    etc.) are ignored.  Variables are interpreted \emph{as is}.

    Can be used for competing risk data.  Variable selection is based on
    the ensemble CHF.
}
\value{
  A list with the following components:
  \item{err.rate}{Prediction error for the forest (a vector of
      length \code{nrep} if variable hunting used).}
  \item{modelSize}{Number of variables selected.}
  \item{topvars}{Character vector of names of the final selected variables.}  
  \item{varselect}{Matrix of values used in determining the set of
      selected variables.}
  \item{rsf.out}{Refitted forest using the final set of selected variables.  
      NULL if \code{big.data}=TRUE.}
}
\author{
    Hemant Ishwaran  \email{hemant.ishwaran@gmail.com} and
    Udaya B. Kogalur \email{kogalurshear@gmail.com}
}
\references{
    H. Ishwaran, U.B. Kogalur, E.Z. Gorodeski, A.J. Minn and 
    M.S. Lauer (2009).  High-dimensional variable selection for survival
    data.  J. Amer. Stat. Assoc. (in press).

    H. Ishwaran, U.B. Kogalur, X. Chen and  A.J. Minn (2009).  Random
    survival forests for high-dimensional data.  
}
\seealso{
  \code{max.subtree},
  \code{rsf}.
}
\examples{
\dontrun{
#------------------------------------------------------------
# Minimal depth variable selection: pbc data with noise

data(pbc, package = "randomSurvivalForest") 
vs <- varSel(Surv(days, status) ~ ., pbc)

# As dimension increases, mtry should increase
pbc.noise <- cbind(pbc, noise = matrix(rnorm(nrow(pbc) * 1000), nrow(pbc)))
vs.bigp <- varSel(Surv(days, status) ~ ., pbc.noise, mtry = 100)


#------------------------------------------------------------
# Variable hunting: van de Vijver microarray breast cancer
# Note: nrep is small for illustration; typical values are nrep = 100

data(vdv, package = "randomSurvivalForest")
vh <- varSel(Surv(Time, Censoring) ~ ., vdv, method = "vh", nrep = 10, nstep = 5)

# Same analysis, but using predefined weights for selecting a gene 
# for node splitting.  We illustrate this using univarate cox p-values.

if (library("survival", logical.return = TRUE) 
    & library("Hmisc", logical.return = TRUE))
{
  cox.weights <- function(rsf.f, rsf.data) {
    event.names <- all.vars(rsf.f)[1:2]
    p <- ncol(rsf.data) - 2
    event.pt <- match(event.names, names(rsf.data))
    predictor.pt <- setdiff(1:ncol(rsf.data), event.pt)
    sapply(1:p, function(j) {
      cox.out <- coxph(rsf.f, rsf.data[, c(event.pt, predictor.pt[j])])
      pvalue <- summary(cox.out)$coef[5]
      if (is.na(pvalue)) 1.0 else 1/(pvalue + 1e-100)
    })
  }       

  data(vdv, package = "randomSurvivalForest")
  rsf.f <- as.formula(Surv(Time, Censoring) ~ .)
  cox.wts <- cox.weights(rsf.f, vdv)
  vh.cox <- varSel(rsf.f, vdv, method = "vh", nstep = 5, predictorWt = cox.wts)

}
}}
\keyword{survival}
\keyword{tree}
