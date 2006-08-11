##**********************************************************************
##**********************************************************************
##
##  RANDOM SURVIVAL FOREST 1.0.0
##
##  Copyright 2006, Cleveland Clinic
##
##  This program is free software; you can redistribute it and/or
##  modify it under the terms of the GNU General Public License
##  as published by the Free Software Foundation; either version 2
##  of the License, or (at your option) any later version.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public
##  License along with this program; if not, write to the Free
##  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
##  Boston, MA  02110-1301, USA.
##
##  Project funded by:
##    National Institutes of Health, HL072771-01
##
##    Michael S Lauer, MD, FACC, FAHA
##    Cleveland Clinic Lerner College of Medicine of CWRU
##    9500 Euclid Avenue
##    Cleveland, OH 44195
##
##    email:  lauerm@ccf.org
##    phone:   216-444-6798
##
##  Written by:
##    Hemant Ishwaran, Ph.D.
##    Dept of Quantitative Health Sciences/Wb4
##    Cleveland Clinic Foundation
##    9500 Euclid Avenue
##    Cleveland, OH 44195
##
##    email:  hemant.ishwaran@gmail.com
##    phone:  216-444-9932
##    URL:    www.bio.ri.ccf.org/Resume/Pages/Ishwaran/ishwaran.html
##    --------------------------------------------------------------
##    Udaya B. Kogalur, Ph.D.
##    Kogalur Shear Corporation
##    5425 Nestleway Drive, Suite L1
##    Clemmons, NC 27012
##
##    email:  ubk2101@columbia.edu
##    phone:  919-824-9825
##    URL:    www.kogalur-shear.com
##
##**********************************************************************
##**********************************************************************

#################################################################
# Primary R function for Random Surival Forests
# ---------------------------------------------------------------
# Description:
#    'rsf' implements Ishwaran and Kogalur's random survival
#    forests algorithm for right censored survival data.  Uses
#    a recursive tree growing procedure with different splitting
#    rules for growing an ensemble cumulative hazard function.
#    An out-of-bag (OOB) estimate of Harrell's concordance index
#    is provided for assessing prediction.
#################################################################

rsf <- function(
    formula = formula(data),
    data = sys.parent(),
    ntree = 100,
    mtry = NULL,
    nodesize = NULL,
    ntime = NULL,
    splitrule = c("conserve", "logrank"),
    proximity = FALSE,
    seed = NULL,
    ntime.upper = 1000,
    trace = FALSE,
    ...)
{
    ### model details
    splitrule <- match.arg(splitrule)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula","data"),names(mf),0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- T
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf,parent.frame())
    mt <- attr(mf,"terms")
    survNames <- all.vars(formula)

    ### remove na's and issue warning
    if (any(is.na(data))){
      cat("\n","*** Warning *** found NA's ... removing all records with NA's","\n")
      if (survNames[3] == ".") {
        data <- na.omit(data)
      }
      else {
        data <- na.omit(data[,is.element(names(data), survNames)])
      }  
    }

    ### process data for native code
    if (survNames[3] == "."){
        cov.names <- names(data)[
            !is.element(names(data),
            survNames[1:2])
        ]
        predictors <- model.matrix(
           as.formula(paste("~", paste(cov.names, collapse = "+"), sep = "")),
           data[is.element(names(data), cov.names)])
    }
    else {
        cov.names <- names(data)[is.element(names(data), c(survNames[-c(1,2)]))]
        predictors <- model.matrix(
            as.formula(paste("~", paste(cov.names, collapse = "+"), sep = "")),
            data[is.element(names(data), cov.names)])
    }
    predictors <- as.matrix(predictors[,!is.element(colnames(predictors), "(Intercept)")])
    cov.names <- colnames(predictors)
    Died <- data[,is.element(names(data), survNames[2])]
    Time <- data[,is.element(names(data), survNames[1])]

    ### internal checks
    if (min(Died) < 0 || max(Died) > 1) 
        stop("censoring variable must be coded as 0 [censored] and 1 [death]")
    if (sum(Died == 0) == length(Died)) 
        stop("no deaths in data: convert 0's to 1's and try again")
 
    ### prepare variables for native code
    ### check that option parameters are correctly specified
    ### cap the number of unique points at a maximum of ntime.upper
    ### work out individuals at risk (used later for mortality calculation)
    rownames(predictors) <- colnames(predictors) <- NULL
    ncov <- length(cov.names)
    if (!is.null(mtry)){if (mtry > ncov) mtry <- ncov}
    if (is.null(mtry)) mtry <- max(floor(sqrt(ncov)), 1)
    if (is.null(nodesize)) nodesize <- min(sum(Died == 1), 3)
    if (is.null(seed) || abs(seed) < 1) seed <- rnorm(1)*1000
    seed <- as.integer(-1 * abs(seed))
    tunique <- sort(unique(Time[Died == 1]))
    if (ntime.upper > 5 & ntime.upper < length(tunique))
        tunique <- tunique[sort(unique(as.integer(seq(1, length(tunique),
                                       length = ntime.upper))))]
    N <- length(tunique)
    if (is.null(ntime) || ntime <= 1) ntime <- N
    if (N > ntime) tunique <- tunique[seq(1, N, length = ntime)]
    Risk <- apply(cbind(1:length(tunique)),
                  1,
                  function(i, tau, tunq){sum(tau >= tunq[i])},
                  tau = Time, tunq = tunique)
    Risk <- Risk - c(Risk[-1],0)
    
    #############################################################
    # Parameters passed to native C code are as follows:
    #############################################################
    # 00 - C code library name
    # 01 - trace (0 = do not send trace, !0 = send trace)
    # 02 - memory useage (0 = no proximity, !0 = proximity)
    # 03 - random seed for repeatability (int < 0) 
    # 04 - split rule to be used 
    #      (1 = Log-Rank, 2 = Conservation of Events)
    # 05 - number of covariates to be randomly selected for
    #      growing tree (int > 0)
    # 06 - number of bootstrap samples to be used (int B > 0)
    # 07 - minimum number of deaths in a node (int > 0)
    # 08 - number of observations (int n)
    # 09 - time (vector double)
    # 10 - status (0 = censored, 1 = death) (vector int)
    # 11 - number of time points of interest (int)
    # 12 - time points of interest (vector double)
    # 13 - number of predictors (int p > 0 )
    # 14 - [n x p] matrix of predictor observations
    #############################################################

    RSF <- .C("rsf",
        as.integer(if (trace == FALSE) 0 else 1),
        as.integer(if (proximity == FALSE) 0 else 1),
        as.integer(seed),
        as.integer(if (splitrule == "logrank") 1 else 2), 
        as.integer(mtry),
        as.integer(ntree),
        as.integer(nodesize),
        as.integer(dim(predictors)[1]),
        as.double(Time),
        as.integer(Died),
        as.integer(length(tunique)),
        as.double(tunique),
        as.integer(dim(predictors)[2]),
        as.numeric(predictors),
        ensemble = as.double(rep(0, dim(predictors)[1]*length(tunique))),
        err.rate = as.double(rep(0, ntree)),
        leaf.count = as.integer(rep(0,ntree)),
        proximity = as.integer(if (proximity == FALSE) 0 else
            rep(0,(dim(predictors)[1]+1)*dim(predictors)[1]/2)))
    mortality <- apply(matrix(RSF$ensemble, 
                       nrow = length(Died),
                       byrow = FALSE),
                       1,
                       function(x, wt) {sum(x*wt)},
                       wt = Risk)
    out <- list(call = match.call(),
        formula = formula,
        terms = mt,
        n = length(Died),
        ndead = sum(Died == 1),
        ntree = ntree,
        mtry = mtry,
        nodesize = nodesize,
        splitrule = splitrule,
        Time = Time,
        Died = Died,
        Time.unique = tunique,
        prednames = cov.names,
        predictors = predictors,
        ensemble = matrix(RSF$ensemble, nrow = length(Died), byrow = FALSE),
        mortality = if (max(mortality) <= length(Died)) mortality else
                        round(mortality*(length(Died))/(1*(max(mortality) == 0)+max(mortality))),
        err.rate = RSF$err.rate,
        leaf.count = RSF$leaf.count,
        proximity = (if (proximity == FALSE) NULL else RSF$proximity))
    class(out) <- "randomSurvivalForest"
    return(out)
}
