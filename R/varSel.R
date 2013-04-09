####**********************************************************************
####**********************************************************************
####
####  RANDOM SURVIVAL FOREST 3.6.4
####
####  Copyright 2013, Cleveland Clinic Foundation
####
####  This program is free software; you can redistribute it and/or
####  modify it under the terms of the GNU General Public License
####  as published by the Free Software Foundation; either version 2
####  of the License, or (at your option) any later version.
####
####  This program is distributed in the hope that it will be useful,
####  but WITHOUT ANY WARRANTY; without even the implied warranty of
####  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
####  GNU General Public License for more details.
####
####  You should have received a copy of the GNU General Public
####  License along with this program; if not, write to the Free
####  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
####  Boston, MA  02110-1301, USA.
####
####  Written by:
####    Hemant Ishwaran, Ph.D.
####    Director of Statistical Methodology
####    Professor, Division of Biostatistics
####    Clinical Research Building, Room 1058
####    1120 NW 14th Street
####    University of Miami, Miami FL 33136
####
####    email:  hemant.ishwaran@gmail.com
####    URL:    http://web.ccs.miami.edu/~hishwaran
####    --------------------------------------------------------------
####    Udaya B. Kogalur, Ph.D.
####    Adjunct Staff
####    Dept of Quantitative Health Sciences
####    Cleveland Clinic Foundation
####    
####    Kogalur & Company, Inc.
####    5425 Nestleway Drive, Suite L1
####    Clemmons, NC 27012
####
####    email:  commerce@kogalur.com
####    URL:    http://www.kogalur.com
####    --------------------------------------------------------------
####
####**********************************************************************
####**********************************************************************

varSel <-
  function(
           formula = NULL,#          
           data = NULL,#
           object = NULL,#      
           method = c("vh", "vhVIMP", "md")[3],#
           ntree = (if (method == "md") 1000 else 500),#
           mvars = (if (!is.null(data) & method != "md") min(1000, round(ncol(data)/5)) else NULL),#
           mtry = (if (!is.null(data) & method == "md") max(sqrt(ncol(data)), ncol(data)/3) else NULL),#
           nodesize = (if (method == "vh" | method == "vhVIMP") 2 else NULL),#
           nsplit = 10,#
           predictorWt = NULL,#
           big.data = FALSE,#
           na.action = c("na.omit", "na.impute")[1],#
           do.trace = 0,#
           always.use = NULL,#  
           nrep = 50,#        
           K = 5,#             
           nstep = 1,#         
           verbose = TRUE,#
           ...)
{
  n.vimp     <- 3      #add top positive n.vimp variables for MD 
  very.big.P <- 5000   #max no. genes in prefiltering for VH
  get.imp <- function(f.o) {
    if (!is.null(dim(f.o$importance))) {
      c(cbind(f.o$importance)[1, ])
    }
    else {
      c(f.o$importance)
    }
  }
  get.imp.all <- function(f.o, pretty = TRUE) {
    if (cr.flag) {
      if (is.null(dim(f.o$importance))) {
        warning("subsampling is coercing a right-censored analysis for CR data\n")
        NA
      }
      else {
        imp.all <- t(f.o$importance)
        if (pretty) colnames(imp.all) <- paste("vimp.", colnames(imp.all), sep = "")
        imp.all
      }
    }
    else {
      imp.all <- cbind(c(f.o$importance))
      if (pretty) colnames(imp.all) <- "vimp"
      imp.all
    }
  }
  get.err <- function(f.o) {
    ntree <- f.o$ntree
    if (!is.null(dim(f.o$err.rate))) {
        cbind(f.o$err.rate)[1, ntree]
    }
    else {
      f.o$err.rate[ntree]
    }
  }
  SD <- function(x) {
    if (all(is.na(x))) {
      NA
    }
    else {
      sd(x, na.rm = TRUE)
    }
  }
  LENGTH <- function(x, y) {
    (length(x) > 0 & length(y) > 0)
  }
  Mtry <- function(x, y) {
    mty <- round((length(x) - length(y))/3)
    if (mty == 0) round(length(x)/3) else mty
  }
  resample <- function(x, size, ...) {
    if(length(x) <= 1) {
      if(!missing(size) && size == 0) x[FALSE] else x
    }
    else {
      sample(x, size, ...)
    }
  }
 cv.folds <- function (n, folds = 10) {
   split(resample(1:n), rep(1:folds, length = n))
 }
 permute.rows <-function(x) {
   n <- nrow(x)
   p <- ncol(x)
   mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
   matrix(t(x)[order(mm)], n, p, byrow = TRUE)
 }
 balanced.folds <- function(y, nfolds = min(min(table(y)), 10)) {
   y[is.na(y)] <- resample(y[!is.na(y)], size = sum(is.na(y)), replace = TRUE)
   totals <- table(y)
   if (length(totals) < 2) {
     return(cv.folds(length(y), nfolds))
   }
   else {
     fmax <- max(totals)
     nfolds <- min(nfolds, fmax)     
     nfolds <- max(nfolds, 2)
     folds <- as.list(seq(nfolds))
     yids <- split(seq(y), y) 
     bigmat <- matrix(NA, ceiling(fmax/nfolds) * nfolds, length(totals))
     for(i in seq(totals)) {
       if(length(yids[[i]])>1){bigmat[seq(totals[i]), i] <- sample(yids[[i]])}
       if(length(yids[[i]])==1){bigmat[seq(totals[i]), i] <- yids[[i]]}
     }
     smallmat <- matrix(bigmat, nrow = nfolds)
     smallmat <- permute.rows(t(smallmat)) 
     res <- vector("list", nfolds)
     for(j in 1:nfolds) {
       jj <- !is.na(smallmat[, j])
       res[[j]] <- smallmat[jj, j]
     }
     return(res)
   }
 }
 rsf.var.hunting <- function(train.id, var.pt, nstep) {
  if (verbose) cat("\t", paste("selecting variables using", mName), "...\n")
  drop.var.pt <- setdiff(var.columns, var.pt)
  if (sum(data[train.id, names(data) == all.vars(rsf.all.f)[2]], na.rm = TRUE) < 2) {
    stop("training data has insufficient deaths: K is probably set too high\n")
  }
  rsf.filter.out  <- rsf(rsf.all.f,
              data=(if (LENGTH(var.pt, drop.var.pt)) data[train.id, -drop.var.pt]
                    else data[train.id, ]),
              ntree = ntree,
              nodesize = nodesize,
              mtry = Mtry(var.columns, drop.var.pt),
              nsplit = nsplit,
              na.action = na.action,
              big.data = TRUE,
              forest = TRUE,
              do.trace = do.trace)
   imp <- get.imp(rsf.filter.out)
   names(imp) <- rsf.filter.out$predictorNames
   if (method == "vhVIMP") {
     VarStrength <- sort(imp, decreasing = TRUE)
     lower.VarStrength <- min(VarStrength) - 1 
     n.lower <- min(2, length(VarStrength))    
     avg.depth <- m.depth <- NA
     sig.vars.old <- names(VarStrength)[1]
   }
   else {
     v <- max.subtree(rsf.filter.out)
     if (is.null(v$order)) { 
       VarStrength <- lower.VarStrength <- 0
       avg.depth <- m.depth <- NA
       sig.vars.old <- names(sort(imp, decreasing = TRUE))[1]
     }
     else {       
       m.depth <- VarStrength <- v$order[, 1]
       avg.depth <- floor(mean(apply(v$nodesAtDepth, 2, function(x){sum(!is.na(x))}), na.rm=TRUE))
       exact.threshold <- v$threshold
       top.vimp <- imp[order(imp, decreasing = TRUE)[min(length(imp), n.vimp)]]
       n.lower <- max(min(2, length(VarStrength)), 
                      sum(VarStrength <= exact.threshold | (imp >= top.vimp & imp > 0)))     
       VarStrength <- max(VarStrength) - VarStrength
       names(m.depth) <- names(VarStrength) <- rsf.filter.out$predictorNames
       VarStrength <- sort(VarStrength, decreasing = TRUE)
       lower.VarStrength <- -1  
       sig.vars.old <- names(VarStrength)[1]
     }
   }
   nstep <- max(round(length(rsf.filter.out$predictorNames)/nstep), 1)
   imp.old <- 0
  for (b in 1:nstep) {
    if (b == 1) {
      if (sum(VarStrength > lower.VarStrength) == 0) {
        sig.vars <- sig.vars.old
        break
      }
      n.upper <- max(which(VarStrength > lower.VarStrength), n.lower)
      threshold <- unique(round(seq(n.lower, n.upper, length = nstep)))
      if (length(threshold) < nstep) {
        threshold <- c(threshold, rep(max(threshold), nstep - length(threshold)))
      }
    }
    sig.vars <- names(VarStrength)[1:threshold[b]]
    if (!is.null(always.use)) {
      sig.vars <- unique(c(sig.vars, always.use))
    }
    if (length(sig.vars) <= 1) { 
      sig.vars <- sig.vars.old
      break
    }
    imp <- get.imp(vimp(rsf.filter.out, sig.vars, joint=TRUE))
    if (verbose) cat("\t iteration: ", b,
                     "  # vars:",     length(sig.vars),
                     "  joint-vimp:",  round(imp, 3),
                     "\r")
    #break when joint vimp no longer increases (strict inequality is a safety feature)
    if (imp  <= imp.old) {
      sig.vars <- sig.vars.old
      break
    }
    else {
      var.pt <- var.columns[match(sig.vars, predictorNames)]
      sig.vars.old <- sig.vars
      imp.old <- imp
    }
  }
  var.pt <- var.columns[match(sig.vars, predictorNames)]
  drop.var.pt <- setdiff(var.columns, var.pt)
  rsf.out  <- rsf(rsf.all.f,
              data=(if (LENGTH(var.pt, drop.var.pt)) data[train.id, -drop.var.pt]
                    else data[train.id, ]),
              ntree = ntree,
              nsplit = nsplit,
              na.action = na.action,                  
              big.data = TRUE,
              forest = TRUE,
              do.trace = do.trace)
  return(list(rsf.out=rsf.out, sig.vars=rsf.out$predictorNames, avg.depth=avg.depth, m.depth=m.depth))
 }
 if (!is.null(object)) {
    if (sum(inherits(object, c("rsf", "grow"), TRUE) == c(1, 2)) != 2) {
       stop("This function only works for objects of class `(rsf, grow)'")
    }
    if (is.null(object$forest)) 
      stop("Forest is empty!  Re-run grow call with forest set to 'TRUE'")
    formula <- object$formula
 }
 else {
    if (is.null(formula) | is.null(data))
      stop("Need to specify 'formula' and 'data'")
 }
 method.names <- c("vh", "vhVIMP", "md")
 method.idx <- which(method.names == method)
 if (length(method.idx) != 1) {
      stop("Invalid method:  ", method)
 }
 fNames <- all.vars(formula, max.names=1e7)
 rsf.all.f <- as.formula(paste("Surv(", fNames[1], ",", fNames[2], ") ~ .", sep=""))
 if (!is.null(object)) {
   data <- as.data.frame(cbind(time = object$time, cens = object$cens, object$predictors))
   colnames(data)[c(1,2)] <- c(fNames[1], fNames[2])
   if (is.null(mvars)) mvars <- min(1000, round(ncol(data)/5))
 }
 if (fNames[3] != ".") {
   data <- data[, is.element(names(data), fNames)]
 }
 if (method == "vh") {
   mName <- "Variable Hunting"
 }
  else if (method == "vhVIMP") {
    mName <- "Variable Hunting (VIMP)"
 }
 else {
   mName <- "Minimal Depth"
 }
 n <- nrow(data)
 P <- ncol(data) - 2
 survival.time <- data[ , names(data) == fNames[1]]
 censoring.status <- data[ , names(data) == fNames[2]]
 eventType <- unique(na.omit(censoring.status))
 eventType <- eventType[eventType > 0]
 cr.flag <- (length(eventType) > 1)
 var.columns <- (1:ncol(data))[-c(which(names(data) == fNames[1]),
             which(names(data) == fNames[2]))]
 predictorNames <- names(data)[var.columns]
 if (!is.null(always.use)) {
   always.use.pt <- var.columns[match(always.use, predictorNames)]
 }
 else {
   always.use.pt <- NULL
 }
 if (!is.null(predictorWt)) {
   if (any(predictorWt < 0) | length(predictorWt) != P | all(predictorWt == 0)) {
     predictorWt <- rep(1/P, P)
   }
   else {
     predictorWt <-predictorWt/sum(predictorWt)
   }
 }
  if (!is.null(mtry)) {
    mtry <- round(mtry)
    if (mtry < 1 | mtry > P) mtry <- max(1, min(mtry, P))
  }
 if (method == "md") {
   if (verbose) cat("minimal depth variable selection ...\n")
   if (!is.null(object)) {
     v <- max.subtree(object)
     pe <- get.err(object)
     ntree <- object$ntree
     mtry <- object$mtry
     nsplit <- object$nsplit
     nodesize <- object$nodesize
     imp <- get.imp(object)
     imp.all <- get.imp.all(object)
     if (is.null(imp)) imp <- rep(0, P)
   }
   else {
     rsf.out <- rsf(rsf.all.f,
                    data,
                    ntree = ntree,
                    mtry = mtry,
                    nsplit = nsplit,
                    nodesize = nodesize,
                    na.action = na.action,
                    big.data = TRUE,
                    forest = TRUE,
                    do.trace = do.trace,
                    predictorWt = predictorWt)
      v <- max.subtree(rsf.out)
      pe <- get.err(rsf.out)
      imp <- get.imp(rsf.out)
      imp.all <- get.imp.all(rsf.out)
      nodesize <- rsf.out$nodesize
      n <- nrow(rsf.out$predictors)
      if (verbose) cat("forest analysis complete ...\n")
   }
   depth <- v$order[, 1]
   top.vimp <- imp[order(imp, decreasing = TRUE)[min(length(imp), n.vimp)]]
   top.var.pt <- (depth <= v$threshold | (imp >= top.vimp & imp > 0))
   modelSize <- sum(top.var.pt)
   o.r.m <- order(depth, decreasing = FALSE)
   top.var.pt <- top.var.pt[o.r.m]
   varselect <- as.data.frame(cbind(depth = depth, vimp = imp.all))[o.r.m, ]
   topvars <- unique(c(always.use, rownames(varselect)[top.var.pt]))
   if (!big.data) {
     if (verbose) cat("finalizing forests ...\n")
     var.pt <- var.columns[match(topvars, predictorNames)]
     var.pt <- unique(c(var.pt, always.use.pt))
     drop.var.pt <- setdiff(var.columns, var.pt)
     rsf.out  <- rsf(rsf.all.f,
               data=(if (LENGTH(var.pt, drop.var.pt)) data[, -drop.var.pt]
                    else data),
               ntree = ntree,
               nsplit = nsplit,
               na.action = na.action,
               big.data = TRUE,
               forest = TRUE,
               do.trace = do.trace)
    }
   else {
     rsf.out <- NULL
   }
   if (verbose) {
     cat("\n\n")
     cat("-----------------------------------------------------------\n")
     cat("var. selection :", mName, "\n")
     cat("no. variables  :", P, "\n")
     cat("no. samples    :", n, "\n")
     cat("no. deaths     :", sum(censoring.status != 0), "\n")
     cat("ntree          :", ntree, "\n")
     cat("nsplit         :", nsplit, "\n")
     if (!is.null(mtry))     cat("mtry           :", mtry, "\n")
     if (!is.null(nodesize)) cat("nodesize       :", nodesize, "\n")
     cat("big data       :", big.data, "\n")
     cat("model size     :", modelSize, "\n")
     cat("PE (full model):", round(mean(100*pe), 4), "\n")
     cat("\n\n")
     cat("Top variables:\n")
     print(round(varselect[top.var.pt, ], 3))
     cat("-----------------------------------------------------------\n")
   }
   return(list(err.rate=pe,
             modelSize=modelSize,
             topvars=topvars,
             varselect=varselect,
             rsf.out=rsf.out
       ))
 }  
 pred.results <- dim.results <- forest.depth <- rep(0, nrep)
 var.signature <- NULL
 var.depth <- matrix(NA, nrep, P)
 var.vimp <- array(NA, c(nrep, P, 1 + (cr.flag)*length(eventType)))
 very.big  <- (P > very.big.P)
 if (mvars < P & very.big & is.null(predictorWt)) {
     if (verbose) cat("Determining selection weights for variables...\n")
     rsf.wts.out  <- rsf(rsf.all.f,
              ntree = max(100, ntree),
              nodesize = max(3, nodesize),
              nsplit = 1,
              data=data,
              na.action = na.action,
              big.data = TRUE)
 }
 for (m in 1:nrep) {
   if (verbose & nrep>1) cat("---------------------  Iteration:", m, "  ---------------------\n")
   all.folds <- balanced.folds(censoring.status, K)
   if (big.data) {
     train.id <- all.folds[[1]]
     test.id <- all.folds[[2]]
   }
   else {
     test.id <- all.folds[[1]]
     train.id <- setdiff(1:n, test.id)
   }
   if (mvars < P & is.null(predictorWt)) {
     if (!very.big) {
        rsf.wts.out  <- rsf(rsf.all.f,
              data=data[train.id, ],
              ntree = min(100, ntree),
              mtry = mtry,
              nsplit = nsplit,
              na.action = na.action,
              big.data = TRUE)
     }
     wts <- pmax(get.imp(rsf.wts.out), 0)
     if (any(wts > 0)) {
       var.pt <- unique(resample(var.columns, mvars, replace = TRUE, prob = wts))
     }
       else {
         var.pt <- var.columns[1:P]
     }
   }
   else {
     var.pt <- var.columns[1:P]
   }
   if (!is.null(predictorWt)) {
     var.pt <- unique(resample(var.columns, mvars, replace = TRUE, prob = predictorWt))
   }
   if (!is.null(always.use)) {
     var.pt <- unique(c(var.pt, always.use.pt))
   }
   object <- rsf.var.hunting(train.id = train.id, var.pt = var.pt, nstep = nstep)
   rsf.out <- object$rsf.out
   sig.vars <- object$sig.vars
   if (method == "vh") {
     forest.depth[m] <- object$avg.depth
     var.depth[m, match(names(object$m.depth), predictorNames)] <- object$m.depth
   }
   #RSF-predictor
   pred.out <- predict(rsf.out,  data[test.id, ])
   pred.results[m] <- get.err(pred.out)
   dim.results[m] <- length(sig.vars)
   var.signature <- c(var.signature, sig.vars)
   #test set vimp
   test.var.pt <- match(pred.out$predictorNames, predictorNames)
   if (m < nrep) {
     var.vimp[m, test.var.pt, ] <- get.imp.all(pred.out, FALSE)
   }
   else {
     imp.all.temp <- get.imp.all(pred.out, TRUE)
     if (!all(is.na(imp.all.temp))) { 
       dimnames(var.vimp) <- list(NULL, NULL, colnames(imp.all.temp))
     }
     var.vimp[m, test.var.pt, ] <- imp.all.temp
   }
   if (verbose) {
     cat("\t                                                                \r")
     cat("\t PE:", round(pred.results[m], 4), "     dim:", dim.results[m], "\n")
   }
 }
 if (is.null(nodesize)) nodesize <- min(3, round(0.632*sum(censoring.status == 1)))
 pred.results <- c(na.omit(pred.results))
 var.freq.all.temp <- 100*tapply(var.signature, var.signature, length)/nrep
 freq.pt <- match(names(var.freq.all.temp), predictorNames)
 var.freq.all <- rep(0, P)
 var.freq.all[freq.pt] <- var.freq.all.temp
 if (method == "vh") {
   var.depth.all <- apply(var.depth, 2, mean, na.rm = T)
 }
 var.vimp.all <- apply(var.vimp, c(2, 3), mean, na.rm = T)
 if (method == "vh") {
   varselect <- cbind(depth = var.depth.all, vimp = var.vimp.all, rel.freq = var.freq.all)
 }
 else {
   varselect <- cbind(vimp = var.vimp.all, rel.freq = var.freq.all)
 }
 o.r.f <- order(var.freq.all, var.vimp.all[, 1], decreasing = TRUE)
 rownames(varselect) <- predictorNames
 varselect <- varselect[o.r.f,, drop = FALSE]
 modelSize <- round(mean(dim.results))  
 topvars <- unique(c(always.use, rownames(varselect)[1:modelSize]))
 if (!big.data) {
   if (verbose) cat("finalizing forests ...\n")
   var.pt <- var.columns[match(rownames(varselect)[1:modelSize], predictorNames)]
   drop.var.pt <- setdiff(var.columns, var.pt)
   rsf.out  <- rsf(rsf.all.f,
               data=(if (LENGTH(var.pt, drop.var.pt)) data[, -drop.var.pt]
                    else data),
               na.action = na.action,
               big.data = TRUE,
               forest = TRUE,
               ntree = ntree,
               nodesize = nodesize,
               nsplit = nsplit,
               do.trace = do.trace)
 }
 else {
   rsf.out <- NULL
 }
 if (verbose) {
   cat("\n\n")
   cat("-----------------------------------------------------------\n")
   cat("var. selection  :", mName, "\n")
   cat("no. variables   :", P, "\n")
   cat("no. samples     :", n, "\n")
   cat("no. deaths      :", sum(censoring.status != 0), "\n")
   cat("K-fold          :", K, "\n")
   cat("no. reps        :", nrep, "\n")
   cat("nstep           :", nstep, "\n")
   cat("ntree           :", ntree, "\n")
   cat("mvars           :", mvars, "\n")
   if (!is.null(mtry))     cat("mtry            :", mtry, "\n")
   if (!is.null(nodesize)) cat("nodesize        :", nodesize, "\n")
   cat("nsplit          :", nsplit, "\n")
   cat("big data        :", big.data, "\n")
   if (method == "vh") {
     cat("depth ratio     :", round(mean(mvars/(2^forest.depth)), 4), "\n")
   }
   cat("model size      :", round(mean(dim.results), 4), "+/-", round(SD(dim.results), 4), "\n")
   cat("PE              :", round(mean(100*pred.results), 4), "+/-", round(SD(100*pred.results), 4), "\n")
   cat("\n\n")
   cat("Top variables:\n")
   print(round(varselect[1:modelSize,, drop = FALSE], 3))
   cat("-----------------------------------------------------------\n")
 }
 return(list(err.rate=pred.results,
             modelSize=modelSize,
             topvars=topvars,
             varselect=varselect,
             rsf.out=rsf.out
       ))
}
