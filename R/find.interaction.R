##**********************************************************************
##**********************************************************************
##
##  RANDOM SURVIVAL FOREST 3.2.3
##
##  Copyright 2008, Cleveland Clinic Foundation
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

find.interaction <- function (
    object, 
    predictorNames = NULL,
    sorted = TRUE,
    npred = NULL, 
    subset = NULL,                              
    nrep  = 1,
    rough = FALSE,
    importance = c("randomsplit", "permute")[1],
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

    ### get ntree
    ntree <- length(unique(object$nativeArray[,1]))

    ### check importance option
    if (importance != "randomsplit" &  importance != "permute")
       stop("Invalid choice for 'importance' option:  " + importance)
        
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
      sorted <- F
      npred <- NULL
      n.interact <- length(predictorNames)
    }
    if (sorted) {
      if (!is.null(object$importance)) {
        o.r <- rev(order(object$importance))
        cov.names <- cov.names[o.r]
      }
    }
    if (!is.null(npred)) {
      n.interact <- min(n.cov, max(round(npred), 1))
    }
    if (n.interact == 1) stop("Pairwise comparisons require more than one candidate variable.")

    
    ### interaction loop: call interaction.rsf
    interact.imp <- rownames.interact.imp <- NULL
    for (k in 1:(n.interact-1)) {
      n.joint.cov <- n.interact-k
      imp <- rep(0,1+n.joint.cov)
      imp.joint <- rep(0,n.joint.cov)
      for (l in (k+1):n.interact) {
        cat("Pairing",cov.names[k],"with",cov.names[l],"\n")
        for (m in 1:nrep) {
          rsfOutput  <-  interaction.rsf(object, cov.names[c(k,l)], importance=importance,
                                         subset=subset, joint=FALSE, rough=rough)
          rsfOutput.joint  <-  interaction.rsf(object, cov.names[c(k,l)], importance=importance,
                                           subset=subset, rough=rough)
          imp[1] <- imp[1]+rsfOutput$importance[1]
          imp[l-k+1] <- imp[l-k+1]+rsfOutput$importance[2]
          imp.joint[l-k] <- imp.joint[l-k]+rsfOutput.joint$importance
        }
      }
      imp[1] <- imp[1]/n.joint.cov
      imp <- imp/nrep
      imp.joint <- imp.joint/nrep
      interact.imp <- rbind(interact.imp,
                 cbind(imp.joint,(imp[1]+imp)[-1],imp.joint-(imp[1]+imp)[-1]))
      rownames.interact.imp <- c(rownames.interact.imp,
                                 paste(cov.names[k],":",cov.names[(k+1):n.interact],
                                 sep=""))
    }
    colnames(interact.imp) <- c("Paired","Additive","Difference")
    rownames(interact.imp) <- rownames.interact.imp

    ### output table
    cat("\n")
    cat("                    No. of variables: ", n.cov,               "\n", sep="")
    cat("                   Variables sorted?: ", sorted,              "\n", sep="")
    cat("   No. of variables used for pairing: ", n.interact,          "\n", sep="")
    cat("    Total no. of paired interactions: ", dim(interact.imp)[1],"\n", sep="")
    cat("            Monte Carlo replications: ", nrep,                "\n", sep="")
    cat("                                VIMP: ", importance,          "\n", sep="")
    cat("                  Fast approximation: ", rough,               "\n", sep="")
    cat("\n")
    print(round(interact.imp,4))

    ### return the goodies
    invisible(interact.imp)
}
