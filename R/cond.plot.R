###########################################################################
# cond.plot                                                               #
#                                                                         #
# The purpose of the cond.plot function is to provide several styles of   #
# conditional plots with base graphics.                                   #
###########################################################################

cond.plot <- function(x, y, z, Style="smoothscatter")
     {
     ### Initial Checks
     if(missing(x)) stop("x is a required argument.")
     if(missing(z)) z <- rep(1, length(x))
     if(Style %in% c("boxplot","scatter","smoothscatter")) {
          if(missing(y)) stop("y is a required argument.")
          if(length(x) != length(y)) stop("length differs for x and y.")}
     if(Style == "boxplot") x <- round(x)
     ### Object Names
     xname <- deparse(substitute(x))
     yname <- deparse(substitute(y))
     zname <- deparse(substitute(z))
     ### Set Plot Window
     znum <- length(unique(z))
     if(znum == 1) par(mfrow=c(1,1))
     else if(znum == 2) par(mfrow=c(2,1))
     else if(znum > 2 & znum < 5) par(mfrow=c(2,2))
     else if(znum > 4 & znum < 7) par(mfrow=c(2,3))
     else par(mfrow=c(3,3))
     if(Style == "boxplot") { ### Conditional Boxplot
          for (i in 1:znum) {
               boxplot(y[which(z == i)] ~ x[which(z == i)], col="red",
                    varwidth=TRUE, main=paste(zname, " = ", i, sep=""),
                    ylab=yname)} #xlab=xname doesn't work here?
          }
     if(Style == "densover") { ### Conditional Density Overlay
          par(mfrow=c(1,1))
          dens <- list()
          for (i in 1:znum) dens[[i]] <- density(x[which(z == i)])
          xmin <- min(sapply(dens, function(x) min(x$x)))
          xmax <- max(sapply(dens, function(x) max(x$x)))
          ymax <- max(sapply(dens, function(x) max(x$y)))
          plot(dens[[1]], xlim=c(xmin, xmax), ylim=c(0,ymax), col=1,
               main="Density Overlay",
               xlab=paste("f(", xname, " | ", zname, ")", sep=""))
          for (i in 2:znum) lines(dens[[i]], col=i)
          }
     else if(Style == "hist") { ### Conditional Histogram
          for (i in 1:znum) {
               hist(x[which(z == i)], col="red", prob=TRUE,
                    main=paste(zname, " = ", i, sep=""),
                    xlab=xname)
               lines(density(x[which(z == i)]))}
          }
     else if(Style == "scatter") { ### Conditional Scatterplot
          mycol <- rgb(100, 100, 100, 50, maxColorValue=255)
          for (i in 1:znum) {
               plot(x[which(z == i)], y[which(z == i)], xlim=range(x),
                    ylim=range(y), col=mycol, pch=16, xlab=xname,
                    ylab=yname, main=paste(zname, " = ", i, sep=""))
               panel.smooth(x[which(z == i)], y[which(z == i)],
                    col=mycol, pch=16)}
          }
     else if(Style == "smoothscatter") { #Conditional smoothScatter
          colramp=colorRampPalette(c("black","red"))
          lowess.na <- function(x, y=NULL, f=2/3, ...) {
               x1 <- subset(x, is.finite(x) & is.finite(y))
               y1 <- subset(y, is.finite(x) & is.finite(y))
               lowess.na <- lowess(x1, y1, f, ...)}
          for (i in 1:znum) {
          smoothScatter(x[which(z == i)], y[which(z == i)],
               colramp=colramp, nrpoints=0, xlab=xname, ylab=yname,
               main=paste(zname, " = ", i, sep=""))
               lines(lowess.na(x[which(z == i)], y[which(z == i)]),
                    col="white")}
          }
     return(invisible())
     }

#End
