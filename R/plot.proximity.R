##**********************************************************************
##**********************************************************************
##
##  RANDOM SURVIVAL FOREST 2.0.0
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

plot.proximity <- function (x, ...) {

   ### check that object is interpretable
   if (sum(inherits(x, c("rsf", "grow"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(x, c("rsf", "predict"), TRUE) == c(1, 2)) != 2)
    stop("This function only works for objects of class `(rsf, grow)' or '(rsf, predict)'.")

  ### proximity checks
  if (is.null(x$proximity))
    stop("Proximity is NULL.  Re-run rsf (grow) analysis with proximity set to 'TRUE'.")
  if (length(x$proximity) == 1) return()

  ### generate plot
  proxm <- matrix(0, x$n, x$n)
  count <- 0
  for (k in 1:x$n){
    proxm[k,1:k] <- x$proximity[(count+1):(count+k)]
    proxm[1:k,k] <- proxm[k,1:k]
    count <- count+k
  }
  proxm <- proxm/diag(proxm)
  loc <- cmdscale(sqrt(1-proxm))
  old.par <- par(no.readonly = TRUE)
  par(mfrow = c(1,1))
  plot(loc[,1], -loc[,2], type="n", xlab="", ylab="", main="cmdscale(Ensemble Mortality)")
  text(loc[,1], -loc[,2],
       pmax(1, ceiling(100*(x$mortality-min(x$mortality))/diff(range(x$mortality)))),
       cex=0.8)
  par(old.par)
  invisible(proximity=proxm)
}