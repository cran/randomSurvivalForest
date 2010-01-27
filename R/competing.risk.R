####**********************************************************************
####**********************************************************************
####
####  RANDOM SURVIVAL FOREST 3.6.1
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

competing.risk <- function (x, plot = TRUE, ...) {

  if (sum(inherits(x, c("rsf", "grow"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(x, c("rsf", "predict"), TRUE) == c(1, 2)) != 2)
    stop("This function only works for objects of class `(rsf, grow)' or '(rsf, predict)'.")

  # number of event types
  if (length(dim(x$ensemble)) == 2) {
    ensb <- array(x$ensemble, dim = c(nrow(x$ensemble), ncol(x$ensemble), 1), 1)
    n.event <- 1
    surv.names <- c("surv")
    legend.names <- c("ensemble DF  ")
  }
  else {
    ensb <- x$ensemble
    n.event <- dim(ensb)[3]
    surv.names <- c("dist", paste("subdist", 1:(n.event-1)))
    legend.names <- c("ensemble Surv  ", paste("ensemble subSurv ", 1:(n.event-1), "  "))
  }

  # ensemble CIF 
  # ensemble subsurvival functions
  # ensemble conditional survival functions
  # ensemble survival function
  surv.ensb.avg <- apply(exp(-ensb[,, 1]), 2, mean, na.rm = TRUE)
  if (n.event > 1) {
    cif.ensb <- sub.ensb <- array(0, dim = c(dim(ensb)[c(1,2)], n.event - 1))
    cif.ensb.avg <- sub.ensb.avg <- csurv.ensb.avg <- NULL
    for (j in 1:(n.event-1)) {
      poe.j <- c(x$poe[j, ])
      h.j   <- ensb[,, 1+j]
      S.j   <- poe.j*exp(-h.j)
      cif.ensb[,, j]   <- poe.j - S.j
      sub.ensb[,, j]   <- S.j
      cif.ensb.avg    <- cbind(cif.ensb.avg, apply(cif.ensb[,,j], 2, mean, na.rm = TRUE))
      sub.ensb.avg    <- cbind(sub.ensb.avg, apply(S.j, 2, mean, na.rm = TRUE))
      csurv.ensb.avg  <- cbind(csurv.ensb.avg, apply(exp(-h.j), 2, mean, na.rm = TRUE))
    }
  }
    
  # plots
  matPlot <- function(M, ylab = "Survival (%)", legend = "ensemble survival", pos = 1) {
     m <- dim(cbind(M))[2]
     if (m > 1) legend <- paste(legend, 1:m, "  ")
     matplot(x$timeInterest, 100*M, xlab = "Time", ylab = ylab, type = "l",
             col = (1:m), lty = 1, lwd = 3)
    legend(c("topright", "bottomright")[pos], legend = legend, col = (1:m), lty = 1, lwd = 3)
  }
  if (plot) {
    if (n.event == 1) {
      matPlot(surv.ensb.avg)
      return(NULL)
    }
    else {
      old.par <- par(no.readonly = TRUE)#save par settings
      par(mfrow = c(2,2))
      matPlot(cif.ensb.avg,   "Probability (%)", "ensemble CIF", 2)
      matPlot(sub.ensb.avg,   "Probability (%)", "ensemble subsurvival")
      matPlot(csurv.ensb.avg, "Survival (%)", "ensemble cond. surv.")
      matPlot(surv.ensb.avg)
      par(old.par)#restore par 
    }
  }
    
  # return the goodies
  invisible(list(cif.ensb = cif.ensb, sub.ensb = sub.ensb))
  
}
