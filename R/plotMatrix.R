###########################################################################
# plotMatrix                                                              #
#                                                                         #
# The purpose of the plotMatrix function is to plot a numerical matrix.   #
###########################################################################

plotMatrix <- function(x, col=colorRampPalette(c("red","black","green"))(100),
     cex=1, circle=TRUE, order=FALSE, zlim=NULL, title="", PDF=FALSE, ...)
     {
     ### Initial Checks
     if(missing(x)) stop("x is a required argument.")
     if(identical(class(x), "bayesfactor")) {
          x <- x$B
          title <- "Bayes Factors"}
     else if(identical(class(x), "demonoid")) {
          if(is.null(x$Covar) | is.list(x$Covar) | is.vector(x$Covar))
               stop("Covar=NULL.")
          x <- x$Covar
          title <- "Covariance"}
     else if(identical(class(x), "iterquad")) {
          if(is.null(x$Covar)) stop("Covar=NULL.")
          x <- x$Covar
          title <- "Covariance"}
     else if(identical(class(x), "laplace")) {
          if(is.null(x$Covar)) stop("Covar=NULL.")
          x <- x$Covar
          title <- "Covariance"}
     else if(identical(class(x), "pmc")) {
          if(is.null(x$Covar)) stop("Covar=NULL.")
          x <- x$Covar
          title <- "Covariance"}
     else if(identical(class(x), "posteriorchecks")) {
          x <- x$Posterior.Correlation
          title <- "Posterior Correlation"}
     else if(identical(class(x), "vb")) {
          if(is.null(x$Covar)) stop("Covar=NULL.")
          x <- x$Covar
          title <- "Covariance"}
     else if(!is.matrix(x)) x <- as.matrix(x)
     min <- min(x)
     max <- max(x)
     yLabels <- rownames(x)
     xLabels <- colnames(x)
     if(is.null(rownames(x))) xLabels <- 1:nrow(x)
     if(is.null(colnames(x))) yLabels <- 1:ncol(x)
     if(!is.null(zlim)) {
          if(length(zlim) != 2) stop("zlim must have length 2.")
          if(zlim[1] >= zlim[2])
               stop("zlim[1] must be lower than zlim[2].")
          min <- zlim[1]
          max <- zlim[2]}
     if(circle == TRUE & order == TRUE) {
          if(!nrow(x) == ncol(x)) 
               stop("The matrix must be square if order is TRUE.")
          x.eigen <- eigen(x)$vectors[, 1:2]
          e1 <- x.eigen[, 1]
          e2 <- x.eigen[, 2]
          alpha <- ifelse(e1 > 0, atan(e2/e1), atan(e2/e1) + pi)
          x <- x[order(alpha), order(alpha)]
          yLabels <- rownames(x)
          xLabels <- colnames(x)}
     ### Plot Matrix
     if(PDF == TRUE) pdf("plotMatrix.pdf")
     if(circle == FALSE) {
          ### Layout and Colors
          layout(matrix(data=c(1,2), nrow=1, ncol=2), widths=c(4,1),
               heights=c(1,1))
          ColorRamp <- col
          ColorLevels <- seq(min, max, length=length(ColorRamp))
          ### Reverse y-axis
          reverse <- nrow(x):1
          yLabels <- yLabels[reverse]
          x <- x[reverse,]
          ### Data Map
          par(mar = c(3,5,2.5,2))
          image(1:length(xLabels), 1:length(yLabels), t(x), col=ColorRamp,
               xlab="", ylab="", axes=FALSE, zlim=c(min,max))
          if(!is.null(title)) title(main=title)
          axis(BELOW<-1, at=1:length(xLabels), labels=xLabels, cex.axis=0.7)
          axis(LEFT <-2, at=1:length(yLabels), labels=yLabels,
               las=HORIZONTAL<-1, cex.axis=0.7)
          par(mar=c(3,2.5,2.5,2))
          image(1, ColorLevels,
               matrix(data=ColorLevels, ncol=length(ColorLevels), nrow=1),
               col=ColorRamp, xlab="", ylab="", xaxt="n")
          layout(1)}
     else {
          col <- col[length(col):1]
          ### Scale Covariance/Precision Matrices
          maxel <- max(abs(x))
          if(maxel > 1) x <- x * (1/maxel)
          ### Plot Setup
          par(mar=c(0, 0, 2, 0), bg="white")
          plot.new()
          plot.window(c(0, ncol(x)), c(0, nrow(x)), asp=1)
          cex.x <- 1 / {log(length(xLabels))/5 + 1}
          cex.y <- 1 / {log(length(yLabels))/5 + 1}
          xlabwidth <- max(strwidth(yLabels, cex=cex))
          ylabwidth <- max(strwidth(xLabels, cex=cex))
          plot.window(c(-xlabwidth + 0.5, ncol(x) + 0.5),
               c(0, nrow(x) + 1 + ylabwidth),
               asp=1, xlab="", ylab="")
          bg <- "gray10"
          rect(0.5, 0.5, ncol(x) + 0.5, nrow(x) + 0.5, col=bg)
          text(rep(-xlabwidth/2, nrow(x)), nrow(x):1, xLabels,
               col="black", cex=cex.x)
          text(1:ncol(x), rep(nrow(x) + 1 + ylabwidth/2, ncol(x)),
               yLabels, srt=90, col="black", cex=cex.y)
          ### Add grid
          lines <- "gray30"
          segments(rep(0.5, nrow(x) + 1), 0.5 + 0:nrow(x),
               rep(ncol(x) + 0.5, nrow(x) + 1),
               0.5 + 0:nrow(x), col=lines)
          segments(0.5 + 0:ncol(x), rep(0.5, ncol(x) + 1),
               0.5 + 0:ncol(x), rep(nrow(x) + 0.5, ncol(x)), col=lines)
          ### Assign circles' fill color
          nc <- length(col)
          if(nc==1) bg <- rep(col, prod(dim(x)))
          else {
               ff <- seq(-1,1, length=nc+1)
               bg2 <- rep(0, prod(dim(x)))
               for (i in 1:prod(dim((x)))) {
                    bg2[i] <- rank(c(ff[2:nc], as.vector(x)[i]),
                         ties.method="random")[nc]}
               bg <- (col[nc:1])[bg2]}
          ### Plot n*m circles using vector language, suggested by Yihui Xie
          ### the area of circles denotes the absolute value of coefficient
          symbols(rep(1:ncol(x), each=nrow(x)), rep(nrow(x):1, ncol(x)),
             add=TRUE, inches=F,
             circles=as.vector(sqrt(abs(x))/2), bg=as.vector(bg))
          }
     if(PDF == TRUE) dev.off()
     }

#End
