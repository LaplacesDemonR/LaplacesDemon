###########################################################################
# plot.bmk                                                                #
#                                                                         #
# The purpose of the plot.bmk function is to plot an object of class bmk. #
###########################################################################

plot.bmk <- function(x, col=colorRampPalette(c("black","red"))(100),
     title="", PDF=FALSE, Parms=NULL, ...)
     {
     ### Initial Checks
     if(missing(x)) stop("x is a required argument.")
     if(!identical(class(x), "bmk"))
          stop("x must be an object of class bmk.")
     ### Selecting Parms
     if(!is.null(Parms)) {
          Parms <- sub("\\[","\\\\[",Parms)
          Parms <- sub("\\]","\\\\]",Parms)
          Parms <- sub("\\.","\\\\.",Parms)
          if(length(grep(Parms[1], rownames(x))) == 0)
               stop("Parameter in Parms does not exist.")
          keeprows <- grep(Parms[1], rownames(x))
          if(length(Parms) > 1) {
               for (i in 2:length(Parms)) {
                    if(length(grep(Parms[i], rownames(x))) == 0)
                         stop("Parameter in Parms does not exist.")
                    keeprows <- c(keeprows,
                         grep(Parms[i], rownames(x)))}}
          x.temp <- as.matrix(x[keeprows,])
          rownames(x.temp) <- rownames(x)[keeprows]
          x <- x.temp
          rm(x.temp)}
     ### Initial Settings
     min <- 0
     max <- 1
     if(is.null(rownames(x))) xLabels <- 1:nrow(x)
     else xLabels <- colnames(x)
     if(is.null(colnames(x))) yLabels <- 1:ncol(x)
     else yLabels <- rownames(x)
     ### plot.bmk
     if(PDF == TRUE) pdf("plot.bmk.pdf")
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
     layout(1)
     if(PDF == TRUE) dev.off()
     }

#End
