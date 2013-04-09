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

competing.risk <- function (x, plot = TRUE, ...) {
  if (sum(inherits(x, c("rsf", "grow"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(x, c("rsf", "predict"), TRUE) == c(1, 2)) != 2)
    stop("This function only works for objects of class `(rsf, grow)' or '(rsf, predict)'.")
  if (sum(inherits(x, c("rsf", "predict"), TRUE) == c(1, 2)) == 2) {
      rsfPred <- TRUE
    }
    else {
      rsfPred <- FALSE
  }
  if (length(dim(x$ensemble)) == 2) {
    if (rsfPred) {
      ensb <- array(x$ensemble, dim = c(nrow(x$ensemble), ncol(x$ensemble), 1), 1)
    }
    else {
      ensb <- array(x$oob.ensemble, dim = c(nrow(x$oob.ensemble), ncol(x$oob.ensemble), 1), 1)
    }
    n.event <- 1
  }
  else {
    if (rsfPred) ensb <- x$ensemble else ensb <- x$oob.ensemble 
    n.event <- dim(ensb)[3]
  }
  surv.ensb.avg <- apply(exp(-ensb[,, 1]), 2, mean, na.rm = TRUE)
  if (n.event > 1) {
    cif.ensb <- sub.ensb <- array(0, dim = c(dim(ensb)[c(1,2)], n.event - 1))
    cond.mort <- matrix(0, n.event - 1, dim(ensb)[1])
    cif.ensb.avg <- sub.ensb.avg <- csurv.ensb.avg <- NULL
    for (j in 1:(n.event-1)) {
      poe.j <- if (rsfPred) c(x$poe[j, ]) else c(x$oob.poe[j, ]) 
      h.j   <- ensb[,, 1 + j]
      S.j   <- poe.j*exp(-h.j)
      cif.ensb[,, j]   <- poe.j - S.j
      sub.ensb[,, j]   <- S.j
      cond.mort[j, ] <- apply(h.j, 1, sum, na.rm = TRUE)
      cif.ensb.avg    <- cbind(cif.ensb.avg, apply(cif.ensb[,,j], 2, mean, na.rm = TRUE))
      sub.ensb.avg    <- cbind(sub.ensb.avg, apply(S.j, 2, mean, na.rm = TRUE))
      csurv.ensb.avg  <- cbind(csurv.ensb.avg, apply(exp(-h.j), 2, mean, na.rm = TRUE))
    }
    dimnames(cif.ensb)[[3]] <- paste("CIF.", 1:(n.event-1), sep = "")
    dimnames(sub.ensb)[[3]] <- paste("subS.", 1:(n.event-1), sep = "")
    rownames(cond.mort) <- paste("event.", 1:(n.event-1), sep = "")
  }
  matPlot <- function(M, ylab = "Survival (%)", legend = "ensemble survival", pos = 2) {
     m <- dim(cbind(M))[2]
     if (!rsfPred) legend <- paste("oob", legend)     
     if (m > 1) legend <- paste(legend, 1:m, "  ")
     matplot(x$timeInterest, 100*M, xlab = "Time", ylab = ylab, type = "l",
             col = (1:m), lty = 1,
             lwd = c(3, 1.5)[1 + 1 * (m > 10)])
    legend(c("topleft", "topright")[pos], legend = legend, col = (1:m), lty = 1,
           xjust = c(0, 1)[pos],
           lwd = c(3, 1.5)[1 + 1 * (m > 10)],
           cex = c(0.75, 0.65)[1 + 1 * (m > 10)])
  }
  if (plot) {
    if (n.event == 1) {
      old.par <- par(no.readonly = TRUE) 
      matPlot(surv.ensb.avg)
      par(old.par) 
      return(NULL)
    }
    else {
      old.par <- par(no.readonly = TRUE) 
      par(mfrow = c(2,2))
      matPlot(cif.ensb.avg,   "Probability (%)", "ensemble CIF", 1)
      matPlot(sub.ensb.avg,   "Probability (%)", "ensemble subsurvival")
      matPlot(csurv.ensb.avg, "Survival (%)", "ensemble cond. surv.")
      matPlot(surv.ensb.avg)
      par(old.par) 
    }
  }
  if (n.event > 1) {
    invisible(list(cif.ensb = cif.ensb, sub.ensb = sub.ensb, cond.mortality = cond.mort))
  } 
}
