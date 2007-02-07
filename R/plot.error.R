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

plot.error <- function (x, ...) {

  if (sum(inherits(x, c("rsf", "grow"), TRUE) == c(1, 2)) != 2 &
      sum(inherits(x, c("rsf", "predict"), TRUE) == c(1, 2)) != 2)
    stop("This function only works for objects of class `(rsf, grow)' or '(rsf, predict)'.")

    old.par <- par(no.readonly = TRUE)
    if (is.null(x$importance)) {
      err <- x$err.rate
      par(mfrow = c(1,1))
      plot(1:length(err), 
           err,
           xlab = "Number of Trees",
           ylab = "Error Rate",
           type = "l")
    }
    else {
      err <- x$err.rate
      imp <- x$importance 
      pred.order <- order(imp)
      n.pred <- length(imp)
      if (n.pred > 100) {
        dotchart.labels <- rep("",n.pred)
        pretty.pt <- pretty(1:n.pred, n=100)
        dotchart.labels[pretty.pt] <- x$predictorNames[pred.order][pretty.pt]
      }
      else {
        dotchart.labels <- x$predictorNames[pred.order]
      }
      par(mfrow = c(1,2))
      plot(1:length(err), 
           err,
           xlab = "Number of Trees",
           ylab = "Error Rate",
           type = "l")
      dotchart(imp[pred.order], dotchart.labels,
               xlab="Importance",
               bg="blue")
      imp.out=as.data.frame(cbind(imp,x$predictorWt),row.names=x$predictorNames)[rev(pred.order),]
      colnames(imp.out) <- c("Importance","predictorWt")
      cat("\n")
      print(round(imp.out,4), justify="right", print.gap=3)
    } 
    par(old.par)      
}