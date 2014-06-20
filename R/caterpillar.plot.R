###########################################################################
# caterpillar.plot                                                        #
#                                                                         #
# The purpose of the caterpillar.plot function is to provide a            #
# caterpillar plot of the posterior summaries in an object of class       #
# demonoid, laplace, or pmc, or also of S x J matrix of S samples and J   #
# variables.                                                              #
###########################################################################

caterpillar.plot <- function(x, Parms=NULL, Title=NULL)
     {
     ### Initial Checks
     if(missing(x)) stop("The x argument is required.")
     par(mfrow=c(1,1))
     if(identical(class(x), "demonoid")) {
          if(any(is.na(x$Summary2))) {
               x <- x$Summary1
               x.lab <- "All Samples"}
          else {
               x <- x$Summary2
               x.lab <- "Stationary Samples"}
          if(!is.null(Parms)) {
               if(is.character(Parms)) {
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
                                   grep(Parms[i], rownames(x)))}
                         }
                    }
               if(is.numeric(Parms)) keeprows <- Parms
               temp <- x
               x <- matrix(x[keeprows,], length(keeprows), ncol(temp))
               rownames(x) <- rownames(temp)[keeprows]
               colnames(x) <- colnames(temp)
               }
          ### Setup
          x.rows <- nrow(x)
          x.lim <- c(min(x[,5]), max(x[,7]))
          y.lim <- c(0, x.rows+1)
          ### Basic Plot
          plot(0, 0, ylim=y.lim, xlim=x.lim, main=Title, sub="",
               xlab=x.lab, ylab="", type="n", ann=TRUE, yaxt="n")
          abline(v=0, col="gray")
          ### Add Medians
          points(x[,6], x.rows:1, pch=20)
          ### Add Horizontal Lines for 2.5%-97.5% Quantiles
          for (i in 1:x.rows) {
               lines(x[i,c(5,7)], c(x.rows-i+1, x.rows-i+1))}
          ### Add y-axis labels
          yy <- x.rows:1
          cex.labels <- 1 / {log(x.rows)/5 + 1}
          axis(2, labels=rownames(x), tick=FALSE, las=1, at=yy,
               cex.axis=cex.labels)
          }
     else if(identical(class(x), "demonoid.hpc")) {
          Chains <- length(x)
          x.temp <- list()
          for (i in 1:Chains) {x.temp[[i]] <- x[[i]][["Summary1"]]}
          x <- x.temp; remove(x.temp)
          x.lab <- "All Samples"
          if(!is.null(Parms)) {
               if(is.character(Parms)) {
                    Parms <- sub("\\[","\\\\[",Parms)
                    Parms <- sub("\\]","\\\\]",Parms)
                    Parms <- sub("\\.","\\\\.",Parms)
                    if(length(grep(Parms[1], rownames(x[[1]]))) == 0)
                         stop("Parameter in Parms does not exist.")
                    keeprows <- grep(Parms[1], rownames(x[[1]]))
                    if(length(Parms) > 1) {
                         for (i in 2:length(Parms)) {
                              if(length(grep(Parms[i], rownames(x[[1]]))) == 0)
                                   stop("Parameter in Parms does not exist.")
                              keeprows <- c(keeprows,
                                   grep(Parms[i], rownames(x[[1]])))}
                         }
                    }
               if(is.numeric(Parms)) keeprows <- Parms
               temp <- x
               for (i in 1:Chains) {
                    x[[i]] <- matrix(x[[i]][keeprows,], length(keeprows),
                         ncol(temp[[1]]))}
               rownames(x[[1]]) <- rownames(temp[[1]])[keeprows]
               colnames(x[[1]]) <- colnames(temp[[1]])
               }
          ### Setup
          x.rows <- nrow(x[[1]])
          x.lim <- c(min(x[[1]][,5]), max(x[[1]][,7]))
          y.lim <- c(0, x.rows+1)
          ### Basic Plot
          plot(0, 0, ylim=y.lim, xlim=x.lim, main=Title, sub="",
               xlab=x.lab, ylab="", type="n", ann=TRUE, yaxt="n")
          abline(v=0, col="gray")
          ### Add Medians
          points(x[[1]][,6], x.rows:1, pch=20)
          for (i in 2:Chains) {points(x[[i]][,6], (x.rows:1)-(i/10), col=i, pch=20)}
          ### Add Horizontal Lines for 2.5%-97.5% Quantiles
          for (i in 1:x.rows) {
               lines(x[[1]][i,c(5,7)], c(x.rows-i+1, x.rows-i+1))}
          for (j in 2:Chains) {for (i in 1:x.rows) {
               lines(x[[j]][i,c(5,7)], c(x.rows-i+1-(j/10), x.rows-i+1-(j/10)), col=j)}}
          ### Add y-axis labels
          yy <- x.rows:1
          cex.labels <- 1 / {log(x.rows)/5 + 1}
          axis(2, labels=rownames(x[[1]]), tick=FALSE, las=1, at=yy,
               cex.axis=cex.labels)
          }
     else if(identical(class(x), "iterquad")) {
          if(any(is.na(x$Posterior))) {
               x <- x$Summary1
               x.lab <- "Point-Estimates"}
          else {
               x <- x$Summary2[1:length(x$Initial.Values),]
               x.lab <- "SIR Samples"}
          if(is.null(Parms)) {
               keeprows <- Parms <- 1:length(x$Initial.Values)}
          else {
               if(is.numeric(Parms)) keeprows <- Parms
               if(is.character(Parms)) {
                    Parms <- sub("\\[","\\\\[",Parms)
                    Parms <- sub("\\]","\\\\]",Parms)
                    Parms <- sub("\\.","\\\\.",Parms)
                    if(length(grep(Parms[1], rownames(x))) == 0)
                         stop("Parameter in Parms does not exist.")
                    keeprows <- grep(Parms[1], rownames(x))
                    if(length(Parms) > 1) {
                         for (i in 2:length(Parms)) {
                              if(length(grep(Parms[i],
                                   rownames(x))) == 0)
                                   stop("Parameter in Parms does not exist.")
                              keeprows <- c(keeprows,
                                   grep(Parms[i], rownames(x)))}
                         }
                    }
               }
          temp <- x
          x <- matrix(x[keeprows,], length(keeprows),
               ncol(temp))
          rownames(x) <- rownames(temp)[keeprows]
          colnames(x) <- colnames(temp)
          if(x.lab != "SIR Samples") Modes <- x[,1]
          else Modes <- x[,6]
          if(x.lab != "SIR Samples") LB <- x[,3]
          else LB <- x[,5]
          if(x.lab != "SIR Samples") UB <- x[,4]
          else UB <- x[,7]
          ### Setup
          x.rows <- length(Modes)
          x.lim <- c(min(LB), max(UB))
          y.lim <- c(0, x.rows+1)
          ### Basic Plot
          plot(0, 0, ylim=y.lim, xlim=x.lim, main=Title, sub="",
               xlab=x.lab, ylab="", type="n", ann=TRUE, yaxt="n")
          abline(v=0, col="gray")
          ### Add Modes
          points(Modes, x.rows:1, pch=20)
          ### Add Horizontal Lines for 2.5%-97.5% Quantiles
          for (i in 1:x.rows) {
               lines(c(LB[i], UB[i]), c(x.rows-i+1, x.rows-i+1))}
          ### Add y-axis labels
          yy <- x.rows:1
          cex.labels <- 1/{log(x.rows)/5 + 1}
          axis(2, labels=rownames(x), tick=FALSE, las=1,
               at=yy, cex.axis=cex.labels)
          }
     else if(identical(class(x), "laplace")) {
          if(any(is.na(x$Posterior))) {
               x <- x$Summary1
               x.lab <- "Point-Estimates"}
          else {
               x <- x$Summary2[1:length(x$Initial.Values),]
               x.lab <- "SIR Samples"}
          if(is.null(Parms)) {
               keeprows <- Parms <- 1:length(x$Initial.Values)}
          else {
               if(is.numeric(Parms)) keeprows <- Parms
               if(is.character(Parms)) {
                    Parms <- sub("\\[","\\\\[",Parms)
                    Parms <- sub("\\]","\\\\]",Parms)
                    Parms <- sub("\\.","\\\\.",Parms)
                    if(length(grep(Parms[1], rownames(x))) == 0)
                         stop("Parameter in Parms does not exist.")
                    keeprows <- grep(Parms[1], rownames(x))
                    if(length(Parms) > 1) {
                         for (i in 2:length(Parms)) {
                              if(length(grep(Parms[i],
                                   rownames(x))) == 0)
                                   stop("Parameter in Parms does not exist.")
                              keeprows <- c(keeprows,
                                   grep(Parms[i], rownames(x)))}
                         }
                    }
               }
          temp <- x
          x <- matrix(x[keeprows,], length(keeprows),
               ncol(temp))
          rownames(x) <- rownames(temp)[keeprows]
          colnames(x) <- colnames(temp)
          if(x.lab != "SIR Samples") Modes <- x[,1]
          else Modes <- x[,6]
          if(x.lab != "SIR Samples") LB <- x[,3]
          else LB <- x[,5]
          if(x.lab != "SIR Samples") UB <- x[,4]
          else UB <- x[,7]
          ### Setup
          x.rows <- length(Modes)
          x.lim <- c(min(LB), max(UB))
          y.lim <- c(0, x.rows+1)
          ### Basic Plot
          plot(0, 0, ylim=y.lim, xlim=x.lim, main=Title, sub="",
               xlab=x.lab, ylab="", type="n", ann=TRUE, yaxt="n")
          abline(v=0, col="gray")
          ### Add Modes
          points(Modes, x.rows:1, pch=20)
          ### Add Horizontal Lines for 2.5%-97.5% Quantiles
          for (i in 1:x.rows) {
               lines(c(LB[i], UB[i]), c(x.rows-i+1, x.rows-i+1))}
          ### Add y-axis labels
          yy <- x.rows:1
          cex.labels <- 1/{log(x.rows)/5 + 1}
          axis(2, labels=rownames(x), tick=FALSE, las=1,
               at=yy, cex.axis=cex.labels)
          }
     else if(identical(class(x), "pmc")) {
          x <- x$Summary
          x.lab <- "All Samples"
          if(!is.null(Parms)) {
               if(is.character(Parms)) {
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
                                   grep(Parms[i], rownames(x)))}
                         }
                    }
               if(is.numeric(Parms)) keeprows <- Parms
               temp <- x
               x <- matrix(x[keeprows,], length(keeprows), ncol(temp))
               rownames(x) <- rownames(temp)[keeprows]
               colnames(x) <- colnames(temp)
               }
          ### Setup
          x.rows <- nrow(x)
          x.lim <- c(min(x[,5]), max(x[,7]))
          y.lim <- c(0, x.rows+1)
          ### Basic Plot
          plot(0, 0, ylim=y.lim, xlim=x.lim, main=Title, sub="",
               xlab=x.lab, ylab="", type="n", ann=TRUE, yaxt="n")
          abline(v=0, col="gray")
          ### Add Medians
          points(x[,6], x.rows:1, pch=20)
          ### Add Horizontal Lines for 2.5%-97.5% Quantiles
          for (i in 1:x.rows) {
               lines(x[i,c(5,7)], c(x.rows-i+1, x.rows-i+1))}
          ### Add y-axis labels
          yy <- x.rows:1
          cex.labels <- 1 / {log(x.rows)/5 + 1}
          axis(2, labels=rownames(x), tick=FALSE, las=1, at=yy,
               cex.axis=cex.labels)
          }
      else if(identical(class(x), "vb")) {
          if(any(is.na(x$Posterior))) {
               x <- x$Summary1
               x.lab <- "Point-Estimates"}
          else {
               x <- x$Summary2[1:length(x$Initial.Values),]
               x.lab <- "SIR Samples"}
          if(is.null(Parms)) {
               keeprows <- Parms <- 1:length(x$Initial.Values)}
          else {
               if(is.numeric(Parms)) keeprows <- Parms
               if(is.character(Parms)) {
                    Parms <- sub("\\[","\\\\[",Parms)
                    Parms <- sub("\\]","\\\\]",Parms)
                    Parms <- sub("\\.","\\\\.",Parms)
                    if(length(grep(Parms[1], rownames(x))) == 0)
                         stop("Parameter in Parms does not exist.")
                    keeprows <- grep(Parms[1], rownames(x))
                    if(length(Parms) > 1) {
                         for (i in 2:length(Parms)) {
                              if(length(grep(Parms[i],
                                   rownames(x))) == 0)
                                   stop("Parameter in Parms does not exist.")
                              keeprows <- c(keeprows,
                                   grep(Parms[i], rownames(x)))}
                         }
                    }
               }
          temp <- x
          x <- matrix(x[keeprows,], length(keeprows),
               ncol(temp))
          rownames(x) <- rownames(temp)[keeprows]
          colnames(x) <- colnames(temp)
          if(x.lab != "SIR Samples") Modes <- x[,1]
          else Modes <- x[,6]
          if(x.lab != "SIR Samples") LB <- x[,3]
          else LB <- x[,5]
          if(x.lab != "SIR Samples") UB <- x[,4]
          else UB <- x[,7]
          ### Setup
          x.rows <- length(Modes)
          x.lim <- c(min(LB), max(UB))
          y.lim <- c(0, x.rows+1)
          ### Basic Plot
          plot(0, 0, ylim=y.lim, xlim=x.lim, main=Title, sub="",
               xlab=x.lab, ylab="", type="n", ann=TRUE, yaxt="n")
          abline(v=0, col="gray")
          ### Add Modes
          points(Modes, x.rows:1, pch=20)
          ### Add Horizontal Lines for 2.5%-97.5% Quantiles
          for (i in 1:x.rows) {
               lines(c(LB[i], UB[i]), c(x.rows-i+1, x.rows-i+1))}
          ### Add y-axis labels
          yy <- x.rows:1
          cex.labels <- 1/{log(x.rows)/5 + 1}
          axis(2, labels=rownames(x), tick=FALSE, las=1,
               at=yy, cex.axis=cex.labels)
          }
     else {
          x <- as.matrix(x)
          x.hpd <- p.interval(x, HPD=TRUE, MM=FALSE, prob=0.95)
          x.median <- apply(x, 2, median)
          x <- cbind(colMeans(x), sqrt(.colVars(x)),
               apply(x, 2, MCSE), ESS(x), x.hpd[,1],
               apply(x, 2, median), x.hpd[,2])
          rownames(x) <- rownames(x.hpd)
          colnames(x) <- c("Mean","SD","MCSE","ESS","LB","Median","UB")
          x.lab <- "All Samples"
          ### Setup
          x.rows <- nrow(x)
          x.lim <- c(min(x[,5]), max(x[,7]))
          y.lim <- c(0, x.rows+1)
          ### Basic Plot
          plot(0, 0, ylim=y.lim, xlim=x.lim, main=Title, sub="",
               xlab=x.lab, ylab="", type="n", ann=TRUE, yaxt="n")
          abline(v=0, col="gray")
          ### Add Medians
          points(x[,6], x.rows:1, pch=20)
          ### Add Horizontal Lines for 2.5%-97.5% Quantiles
          for (i in 1:x.rows) {
               lines(x[i,c(5,7)], c(x.rows-i+1, x.rows-i+1))}
          ### Add y-axis labels
          yy <- x.rows:1
          cex.labels <- 1 / {log(x.rows)/5 + 1}
          axis(2, labels=rownames(x), tick=FALSE, las=1, at=yy,
               cex.axis=cex.labels)
          }
     }

#End
