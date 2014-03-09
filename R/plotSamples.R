###########################################################################
# plotSamples                                                             #
#                                                                         #
# The purpose of this function is to provide basic plots that are         #
# extended to include samples. This takes a N x S matrix X with N records #
# and S simulations.                                                      #
###########################################################################

plotSamples <- function(X, Style="KDE", LB=0.025, UB=0.975, Title=NULL)
     {
     ### Initial Checks
     if(missing(X)) stop("The X argument is required.")
     if(!is.matrix(X)) stop("X is required to be a matrix.")
     Xname <- deparse(substitute(X))
     N <- nrow(X)
     S <- ncol(X)
     if({LB < 0} | {LB >= 0.5})
          stop("LB must be in the interval [0, 0.5).")
     if({UB <= 0.5} | {UB > 1})
          stop("UB must be in the interval (0.5, 1].")
     if(Style == "barplot") {
          ### Qauntiles
          qLB <- apply(X, 1, quantile, prob=LB, na.rm=TRUE)
          qMed <- apply(X, 1, quantile, prob=0.5, na.rm=TRUE)
          qUB <- apply(X, 1, quantile, prob=UB, na.rm=TRUE)
          ### Plot
          barplot(qUB, names.arg=rownames(X),
               ylim=c(min(qLB), max(qUB)),
               col=rgb(255, 0, 0, 50, maxColorValue=255),
               main=Title, xlab=Xname, ylab="Value")
          barplot(qMed, add=TRUE,
               col=rgb(255, 0, 0, 75, maxColorValue=255))
          barplot(qLB, add=TRUE,
               col=rgb(255, 0, 0, 100, maxColorValue=255))
          }
     else if(Style == "dotchart") {
          ### Qauntiles
          qLB <- apply(X, 1, quantile, prob=LB, na.rm=TRUE)
          qMed <- apply(X, 1, quantile, prob=0.5, na.rm=TRUE)
          qUB <- apply(X, 1, quantile, prob=UB, na.rm=TRUE)
          ### Plot
          dotchart(qMed, xlim=range(c(qLB,qUB)), main=Title)
          for (i in 1:nrow(X)) lines(x=c(qLB[i],qUB[i]), y=c(i,i))
          }
     else if(Style == "hist") {
          ### Histogram counts by Sample
          br <- hist(X, plot=FALSE)
          H <- matrix(0, length(br$counts), S)
          for (s in 1:S)
               H[,s] <- hist(X[,s], breaks=br$breaks, plot=FALSE)$counts
          ### Qauntiles
          qLB <- apply(H, 1, quantile, prob=LB, na.rm=TRUE)
          qMed <- apply(H, 1, quantile, prob=0.5, na.rm=TRUE)
          qUB <- apply(H, 1, quantile, prob=UB, na.rm=TRUE)
          ### Plot
          barplot(qUB, names.arg=br$mids, ylim=c(0, max(qUB)),
               col=rgb(255, 0, 0, 50, maxColorValue=255),
               main=Title, xlab=Xname, ylab="Frequency")
          barplot(qMed, add=TRUE,
               col=rgb(255, 0, 0, 75, maxColorValue=255))
          barplot(qLB, add=TRUE,
               col=rgb(255, 0, 0, 100, maxColorValue=255))
          }
     else if(Style == "KDE") {
          ### KDE by Sample
          D.y <- density(X[,1], na.rm=TRUE)
          D.x <- D.y$x
          D.y <- D.y$y
          D.y <- matrix(D.y, length(D.y), S)
          for (s in 2:S) D.y[,s] <- density(X[,s])$y
          ### Quantiles
          qLB <- apply(D.y, 1, quantile, prob=LB, na.rm=TRUE)
          qMed <- apply(D.y, 1, quantile, prob=0.5, na.rm=TRUE)
          qUB <- apply(D.y, 1, quantile, prob=UB, na.rm=TRUE)
          ### Plot
          plot(D.x, qUB, col=rgb(255, 0, 0, 50, maxColorValue=255),
               type="l", ylim=c(0,max(qUB)), main=Title, xlab=Xname,
               ylab="Density")
          polygon(x=D.x, y=qUB, col=rgb(255, 0, 0, 50, maxColorValue=255),
               border=NULL)
          polygon(x=D.x, y=qMed, col=rgb(255, 0, 0, 75, maxColorValue=255),
               border=NULL)
          polygon(x=D.x, y=qLB, col=rgb(255, 0, 0, 100, maxColorValue=255),
               border=NULL)
          }
     else if(Style == "Time-Series") {
          ### Qauntiles
          qLB <- apply(X, 1, quantile, prob=LB, na.rm=TRUE)
          qMed <- apply(X, 1, quantile, prob=0.5, na.rm=TRUE)
          qUB <- apply(X, 1, quantile, prob=UB, na.rm=TRUE)
          ### Plot
          plot(1:nrow(X), qMed, ylim=range(c(qLB,qUB)), col="white",
                main=Title, xlab="Time", ylab="Value")
          polygon(c(1:nrow(X),rev(1:nrow(X))),c(qLB,rev(qUB)),
                col=rgb(255, 0, 0, 50, maxColorValue=255),
                border=FALSE)
          lines(1:nrow(X), qMed, col="red")
          }
     else stop("Style is unknown.")
     return(invisible())
     }

#End

