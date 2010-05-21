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

plot.proximity <- function (x, plot = TRUE, ...) {

   ### check that object is interpretable
   if (sum(inherits(x, c("rsf", "grow"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(x, c("rsf", "predict"), TRUE) == c(1, 2)) != 2)
    stop("This function only works for objects of class `(rsf, grow)' or '(rsf, predict)'")

  ### proximity checks
  if (is.null(x$proximity))
    stop("Proximity is NULL.  Re-run analysis with proximity set to 'TRUE'.")
  if (length(x$proximity) == 1) return()

  ### generate plot
  proximity <- matrix(0, x$n, x$n)
  count <- 0
  for (k in 1:x$n){
    proximity[k,1:k] <- x$proximity[(count+1):(count+k)]
    proximity[1:k,k] <- proximity[k,1:k]
    count <- count+k
  }
  proximity <- proximity/diag(proximity)
  if (plot) {
    loc <- cmdscale(sqrt(1-proximity))
    old.par <- par(no.readonly = TRUE)
    par(mfrow = c(1,1))
    plot(loc[,1], -loc[,2], type="n", xlab="", ylab="", main="cmdscale(Ensemble Mortality)")
    if (sum(x$cens == 1) > 1 & sum(x$cens == 0) > 1) {
      text(loc[,1], -loc[,2],
           pmax(1, ceiling(100*(x$mortality -
                      min(x$mortality, na.rm = TRUE))/diff(range(x$mortality, na.rm = TRUE)))),
           col = 3*x$cens+1,
           cex = 0.8)
    }
    else {
      text(loc[,1], -loc[,2],
           pmax(1, ceiling(100*(x$mortality -
                      min(x$mortality, na.rm = TRUE))/diff(range(x$mortality, na.rm = TRUE)))),
           cex = 0.8)
    }
    par(old.par)
  }
  invisible(proximity)
}
