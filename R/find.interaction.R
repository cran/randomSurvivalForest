####**********************************************************************
####**********************************************************************
####
####  RANDOM SURVIVAL FOREST 3.6.3
####
####  Copyright 2009, Cleveland Clinic Foundation
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
####  ----------------------------------------------------------------
####  Project Partially Funded By:
####    --------------------------------------------------------------
####    National Institutes of Health,  Grant HHSN268200800026C/0001
####
####    Michael S. Lauer, M.D., FACC, FAHA 
####    National Heart, Lung, and Blood Institute
####    6701 Rockledge Dr, Room 10122
####    Bethesda, MD 20892
####
####    email:  lauerm@nhlbi.nih.gov
####
####    --------------------------------------------------------------
####    Case Western Reserve University/Cleveland Clinic  
####    CTSA Grant:  UL1 RR024989, National Center for
####    Research Resources (NCRR), NIH
####
####    --------------------------------------------------------------
####    Dept of Defense Era of Hope Scholar Award, Grant W81XWH0910339
####    Andy Minn, M.D., Ph.D.
####    Department of Radiation and Cellular Oncology, and
####    Ludwig Center for Metastasis Research
####    The University of Chicago, Jules F. Knapp Center, 
####    924 East 57th Street, Room R318
####    Chicago, IL 60637
#### 
####    email:  aminn@radonc.uchicago.edu
####
####    --------------------------------------------------------------
####    Bryan Lau, Ph.D.
####    Department of Medicine, Johns Hopkins School of Medicine,
####    Baltimore, Maryland 21287
####
####    email:  blau1@jhmi.edu
####
####  ----------------------------------------------------------------
####  Written by:
####    --------------------------------------------------------------
####    Hemant Ishwaran, Ph.D.
####    Dept of Quantitative Health Sciences/Wb4
####    Cleveland Clinic Foundation
####    9500 Euclid Avenue
####    Cleveland, OH 44195
####
####    email:  hemant.ishwaran@gmail.com
####    phone:  216-444-9932
####    URL:    www.bio.ri.ccf.org/Resume/Pages/Ishwaran/ishwaran.html
####
####    --------------------------------------------------------------
####    Udaya B. Kogalur, Ph.D.
####    Dept of Quantitative Health Sciences/Wb4
####    Cleveland Clinic Foundation
####    
####    Kogalur Shear Corporation
####    5425 Nestleway Drive, Suite L1
####    Clemmons, NC 27012
####
####    email:  ubk2101@columbia.edu
####    phone:  919-824-9825
####    URL:    www.kogalur-shear.com
####    --------------------------------------------------------------
####
####**********************************************************************
####**********************************************************************

find.interaction <- function (
    object,# 
    predictorNames = NULL,#
    method = c("maxsubtree", "vimp")[1],#
    sorted = TRUE,#
    npred = NULL,#
    subset = NULL,#                              
    nrep  = 1,#
    rough = FALSE,#
    importance = c("randomsplit", "permute")[1],#
    seed = NULL,#
    do.trace = FALSE,#
    ...) {
 
    ## Check that 'object' is of the appropriate type.
    if (is.null(object)) stop("Object is empty!")
    if (sum(inherits(object, c("rsf", "grow"), TRUE) == c(1, 2)) != 2    &
        sum(inherits(object, c("rsf", "forest"), TRUE) == c(1, 2)) != 2)
       stop("This function only works for objects of class `(rsf, grow)' or '(rsf, forest)'.")
    if (sum(inherits(object, c("rsf", "grow"), TRUE) == c(1, 2)) == 2) {
      if (is.null(object$forest)) 
        stop("Forest is empty!  Re-run grow call with forest set to 'TRUE'.")
    }

    ### check method option
    if (method != "maxsubtree" &  method != "vimp")
       stop("Invalid choice for 'method':  " + importance)

    ### check importance option
    if (importance != "randomsplit" &  importance != "permute")
       stop("Invalid choice for 'importance':  " + importance)

    ### number of event types
    n.event  <- length(unique(na.omit(object$cens)[na.omit(object$cens) > 0]))
    if (n.event > 1) n.event <- n.event + 1
    
    ### variable name details
    ### special treatment needed if user passes predictor names
    ### should predictors be sorted by importance?
    ### determine number of variables to be paired-up
    cov.names <- object$predictorNames
    n.interact <- n.cov <- length(cov.names)
    if (!is.null(predictorNames)) {
      if (sum(is.element(cov.names, predictorNames)) == 0) {
           cat("Variables do not match original analysis:","\n")
           print(cov.names)
           stop()
      }
      predictorNames <- unique(predictorNames[is.element(predictorNames, cov.names)])
      if (length(predictorNames) == 1)
        stop("Pairwise comparisons require more than one candidate variable.")
      cov.names <- predictorNames
      n.interact <- length(predictorNames)
    }
    if (sorted) {
      if (!is.null(object$importance)) {
        if (n.event == 1) {
          o.r <- order(object$importance, decreasing = TRUE)
        }
        else {
          o.r <- order(object$importance[1, ], decreasing = TRUE)
        }
        cov.names <- cov.names[o.r]
      }
    }
    if (!is.null(npred)) {
      n.interact <- min(n.cov, max(round(npred), 1))
    }
    if (n.interact == 1) stop("Pairwise comparisons require more than one candidate variable.")


    if (method == "vimp") {
      ### VIMP approach
      ### interaction loop: call vimp
      rownames.interact.imp <- NULL
      interact.imp <- vector("list", n.event)
      for (k in 1:(n.interact-1)) {
        n.joint.cov <- n.interact-k
        imp <- matrix(0,1+n.joint.cov, n.event)
        imp.joint <- matrix(0,n.joint.cov, n.event)
        for (l in (k+1):n.interact) {
          cat("Pairing",cov.names[k],"with",cov.names[l],"\n")
          for (m in 1:nrep) {
            rsfOutput  <-  vimp(object, cov.names[c(k,l)], importance=importance,
                     subset=subset, joint=FALSE, rough=rough, seed=seed, do.trace=do.trace)
            rsfOutput.joint  <-  vimp(object, cov.names[c(k,l)], importance=importance,
                     subset=subset, rough=rough, seed=seed, do.trace=do.trace)
            imp[1,] <- imp[1,]+rbind(rsfOutput$importance)[,1]
            imp[l-k+1,] <- imp[l-k+1,]+rbind(rsfOutput$importance)[,2]
            imp.joint[l-k,] <- imp.joint[l-k,]+rsfOutput.joint$importance
          }
        }
        imp[1,] <- imp[1,]/n.joint.cov
        imp <- imp/nrep
        imp.joint <- imp.joint/nrep
        for (N in 1:n.event) {
          interact.imp[[N]] <- rbind(interact.imp[[N]],
                 cbind(imp.joint[,N],(imp[1,N]+imp[,N])[-1],imp.joint[,N]-(imp[1,N]+imp[,N])[-1]))
        }
        rownames.interact.imp <- c(rownames.interact.imp,
                                 paste(cov.names[k],":",cov.names[(k+1):n.interact],
                                 sep=""))
      }
      if (n.event == 1) {
        interact.imp <- interact.imp[[1]]
        colnames(interact.imp) <- c("Paired","Additive","Difference")
        rownames(interact.imp) <- rownames.interact.imp
      }
      else {
        for (N in 1:n.event) {
          colnames(interact.imp[[N]]) <- c("Paired","Additive","Difference")
          rownames(interact.imp[[N]]) <- rownames.interact.imp
        }
        names(interact.imp) <- c("CHF", paste("condCHF.", 1:(n.event - 1), sep = ""))
      }
      n.pairs <- length(rownames.interact.imp)

      ### output table
      cat("\n")
      cat("                      Technique used: ", method,              "\n", sep="")
      cat("                    No. of variables: ", n.cov,               "\n", sep="")
      cat("           Variables sorted by VIMP?: ", sorted,              "\n", sep="")
      cat("   No. of variables used for pairing: ", n.interact,          "\n", sep="")
      cat("    Total no. of paired interactions: ", n.pairs,             "\n", sep="")
      cat("            Monte Carlo replications: ", nrep,                "\n", sep="")
      cat("    Type of noising up used for VIMP: ", importance,          "\n", sep="")
      cat("                  Fast approximation: ", rough,               "\n", sep="")
      cat("\n")
      if (n.event == 1) print(round(interact.imp,4)) else print(interact.imp)

      ### return the goodies
      invisible(interact.imp)
    }
    else {
      ### maximal subtree approach
      v <- max.subtree(object, sub.order = TRUE, max.order = 1)
      subOrder <- v$subOrder
      if (sorted) {
        o.r <- order(diag(subOrder), decreasing = FALSE)
        subOrder <- subOrder[o.r, o.r]
      }
      cov.pt <- is.element(colnames(subOrder), cov.names[1:n.interact])
      subOrder <- subOrder[cov.pt, cov.pt]
      
      ### output table
      cat("\n")
      cat("                      Technique used: ", method,              "\n", sep="")
      cat("                    No. of variables: ", n.cov,               "\n", sep="")
      cat("  Variables sorted by minimal depth?: ", sorted,              "\n", sep="")
      cat("\n")
      print(round(subOrder, 2))

      ### return the goodies
      invisible(subOrder)
    }
     
}
