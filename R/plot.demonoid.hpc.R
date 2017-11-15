###########################################################################
# plot.demonoid.hpc                                                       #
#                                                                         #
# The purpose of the plot.demonoid.hpc function is to plot an object of   #
# class demonoid.hpc.                                                     #
###########################################################################

plot.demonoid.hpc <- function(x, BurnIn=0, Data=NULL, PDF=FALSE, Parms=NULL,
                              FileName = paste0("laplacesDemon-plot_", format(Sys.time(), "%Y-%m-%d_%T"), ".pdf"), ...)
     {
     ### Initial Checks
     if(missing(x)) stop("The x argument is required.")
     if(class(x) != "demonoid.hpc")
          stop("x must be of class demonoid.hpc.")
     Chains <- length(x)
     if(is.null(Data)) stop("The Data argument is NULL.")
     nn <- nrow(x[[1]]$Posterior1)
     if(BurnIn >= nn) BurnIn <- 0
     Stat.at <- BurnIn + 1
     ### Selecting Parms
     if(is.null(Parms)) {
          Posterior <- list()
          for (i in 1:Chains) {
               Posterior[[i]] <- x[[i]][["Posterior1"]]}
          }
     else {
          Parms <- sub("\\[","\\\\[",Parms)
          Parms <- sub("\\]","\\\\]",Parms)
          Parms <- sub("\\.","\\\\.",Parms)
          if(length(grep(Parms[1], colnames(x[[1]]$Posterior1))) == 0)
               stop("Parameter in Parms does not exist.")
          keepcols <- grep(Parms[1], colnames(x[[1]]$Posterior1))
          if(length(Parms) > 1) {
               for (i in 2:length(Parms)) {
                    if(length(grep(Parms[i],
                         colnames(x[[1]]$Posterior1))) == 0)
                         stop("Parameter in Parms does not exist.")
                    keepcols <- c(keepcols,
                         grep(Parms[i], colnames(x[[1]]$Posterior1)))}}
          Posterior <- list()
          for (i in 1:Chains) {
               Posterior[[i]] <- matrix(x[[i]][["Posterior1"]][,keepcols],
                    nn, length(keepcols))}
          }
     if(PDF == TRUE) {
          pdf(FileName)
          par(mfrow=c(3,3))
          }
     else {par(mfrow=c(3,3), ask=TRUE)}
     ### Plot Parameters
     for (j in 1:ncol(Posterior[[1]])) {
          plot(Stat.at:nn, Posterior[[1]][Stat.at:nn,j],
               ylim=c(min(matrix(sapply(Posterior, function(x) {
                         min = min(x[,j])}), nn, Chains)[Stat.at:nn,]),
                    max(matrix(sapply(Posterior, function(x) {
                         max = max(x[,j])}), nn, Chains)[Stat.at:nn,])),
               col=rgb(0,0,0,50,maxColorValue=255), type="l",
               xlab="Thinned Samples", ylab="Value",
               main=colnames(Posterior[[1]])[j])
          for (n in 2:Chains) {
               lines(Stat.at:nn, Posterior[[n]][Stat.at:nn,j],
                    col=rgb(col2rgb(n)[1],col2rgb(n)[2],col2rgb(n)[3],50,
                    maxColorValue=255))}
          plot(density(Posterior[[1]][Stat.at:nn,j]), col="white",
               xlab="Value", main=colnames(Posterior[[1]])[j])
          polygon(density(Posterior[[1]][Stat.at:nn,j]),
               col=rgb(0,0,0,50,maxColorValue=255), border=NA)
          for (n in 2:Chains) {
               polygon(density(Posterior[[n]][Stat.at:nn,j]),
               col=rgb(col2rgb(n)[1],col2rgb(n)[2],col2rgb(n)[3],50,
                    maxColorValue=255), border=NA)}
          abline(v=0, col="red", lty=2)
          ### Only plot an ACF if there's > 1 unique values
          if(!is.constant(Posterior[[1]][Stat.at:nn,j])) {
               z <- acf(Posterior[[1]][Stat.at:nn,j], plot=FALSE)
               se <- 1/sqrt(length(Posterior[[1]][Stat.at:nn,j]))
               plot(z$lag, z$acf, ylim=c(min(z$acf,-2*se),1),
                    col=rgb(0,0,0,50,maxColorValue=255), type="h",
                    main=colnames(Posterior[[1]])[j], xlab="Lag",
                    ylab="Correlation")
               abline(h=(2*se), col="red", lty=2)
               abline(h=(-2*se), col="red", lty=2)
               for (n in 2:Chains) {
                    z <- acf(Posterior[[n]][Stat.at:nn,j], plot=FALSE)
                    se <- 1/sqrt(length(Posterior[[n]][Stat.at:nn,j]))
                    lines(z$lag, z$acf, col=rgb(col2rgb(n)[1],
                         col2rgb(n)[2],col2rgb(n)[3],50,
                         maxColorValue=255))}
               }
          else {plot(0, 0, main=paste(colnames(Posterior[[1]])[j],
               "is a constant."))}
          }
     rm(Posterior)
     ### Plot Deviance
     Deviance <- list()
     for (i in 1:Chains) {Deviance[[i]] <- x[[i]][["Deviance"]]}
     plot(Stat.at:nn, Deviance[[1]][Stat.at:nn],
          ylim=c(min(sapply(Deviance, function(x) {min(x[Stat.at:nn])})),
               max(sapply(Deviance, function(x) {max(x[Stat.at:nn])}))),
          col=rgb(0,0,0,50,maxColorValue=255),
          type="l", xlab="Thinned Samples", ylab="Value", main="Deviance")
     for (n in 2:Chains) {
          lines(Stat.at:nn, Deviance[[n]][Stat.at:nn],
               col=rgb(col2rgb(n)[1], col2rgb(n)[2],col2rgb(n)[3],50,
                    maxColorValue=255))}
     plot(density(Deviance[[1]][Stat.at:nn]), col="white",
          xlab="Value", main="Deviance")
     polygon(density(Deviance[[1]][Stat.at:nn]),
          col=rgb(0,0,0,50,maxColorValue=255), border=NA)
     for (n in 2:Chains) {
          polygon(density(Deviance[[n]][Stat.at:nn]),
               col=rgb(col2rgb(n)[1], col2rgb(n)[2],col2rgb(n)[3],50,
                    maxColorValue=255), border=NA)}
     abline(v=0, col="red", lty=2)
     #### Only plot an ACF if there's > 1 unique values
     if(!is.constant(Deviance[[1]][Stat.at:nn])) {
          z <- acf(Deviance[[1]][Stat.at:nn], plot=FALSE)
          se <- 1/sqrt(length(Deviance[[1]][Stat.at:nn]))
          plot(z$lag, z$acf, ylim=c(min(z$acf,-2*se),1),
               col=rgb(0,0,0,50,maxColorValue=255), type="h",
               main="Deviance", xlab="Lag", ylab="Correlation")
          abline(h=(2*se), col="red", lty=2)
          abline(h=(-2*se), col="red", lty=2)
          for (n in 2:Chains) {
               z <- acf(Deviance[[n]][Stat.at:nn], plot=FALSE)
               se <- 1/sqrt(length(Deviance[[n]][Stat.at:nn]))
               lines(z$lag, z$acf, col=rgb(col2rgb(n)[1],
                      col2rgb(n)[2],col2rgb(n)[3],50,
                      maxColorValue=255))}
          }
     else {plot(0, 0, main="Deviance is a constant.")}
     rm(Deviance)
     #### Plot Monitored Variables
     J <- length(Data[["mon.names"]])
     Monitor <- list()
     for (i in 1:Chains) {
          Monitor[[i]] <- matrix(x[[i]][["Monitor"]], nn, J)}
     for (j in 1:J) {
          plot(Stat.at:nn, Monitor[[1]][Stat.at:nn,j],
               ylim=c(min(sapply(Monitor, function(x) {min(x[Stat.at:nn,j])})),
                    max(sapply(Monitor, function(x) {max(x[Stat.at:nn,j])}))),
               col=rgb(0,0,0,50,maxColorValue=255),
               type="l", xlab="Thinned Samples", ylab="Value",
               main=Data[["mon.names"]][j])
          for (n in 2:Chains) {
               lines(Stat.at:nn, Monitor[[n]][Stat.at:nn,j],
                    col=rgb(col2rgb(n)[1],col2rgb(n)[2],col2rgb(n)[3],50,
                    maxColorValue=255))}
          plot(density(Monitor[[1]][Stat.at:nn,j]), col="white",
               xlab="Value", main=Data[["mon.names"]][j])
          polygon(density(Monitor[[1]][Stat.at:nn,j]),
               col=rgb(0,0,0,50,maxColorValue=255), border=NA)
          for (n in 2:Chains) {
               polygon(density(Monitor[[n]][Stat.at:nn,j]),
                    col=rgb(col2rgb(n)[1],col2rgb(n)[2],col2rgb(n)[3],50,
                         maxColorValue=255), border=NA)}
          abline(v=0, col="red", lty=2)
          ### Only plot an ACF if there's > 1 unique values
          if(!is.constant(Monitor[[1]][Stat.at:nn,j])) {
               z <- acf(Monitor[[1]][Stat.at:nn,j], plot=FALSE)
               se <- 1/sqrt(length(Monitor[[1]][Stat.at:nn,j]))
               plot(z$lag, z$acf, ylim=c(min(z$acf,-2*se),1),
                    col=rgb(0,0,0,50,maxColorValue=255), type="h",
                    main=Data[["mon.names"]][j], xlab="Lag",
                    ylab="Correlation")
               abline(h=(2*se), col="red", lty=2)
               abline(h=(-2*se), col="red", lty=2)
               for (n in 2:Chains) {
                    z <- acf(Monitor[[n]][Stat.at:nn,j], plot=FALSE)
                    se <- 1/sqrt(length(Monitor[[n]][Stat.at:nn,j]))
                    lines(z$lag, z$acf, col=rgb(col2rgb(n)[1],
                         col2rgb(n)[2],col2rgb(n)[3],50,
                         maxColorValue=255))}
               }
          else {plot(0, 0, main=paste(Data[["mon.names"]][j],
                     "is a constant."))}
          }
     rm(Monitor)
     #### Diminishing Adaptation
     if(nrow(x[[1]]$CovarDHis) > 1) {
          Diff <- abs(diff(x[[1]]$CovarDHis))
          adaptchange <- matrix(NA, nrow(Diff), 3)
          for (i in 1:nrow(Diff)) {
               adaptchange[i,1:3] <- as.vector(quantile(Diff[i,],
                    probs=c(0.025, 0.500, 0.975)))}
          plot(adaptchange[,2], ylim=c(min(adaptchange), max(adaptchange)),
               type="l", col=rgb(0,0,0,50,maxColorValue=255),
               xlab="Adaptations", ylab="Absolute Difference",
               main="Proposal Variance", sub="Median=Red, 95% Bounds=Gray")
          for (n in 2:Chains) {
               Diff <- abs(diff(x[[n]]$CovarDHis))
               adaptchange <- matrix(NA, nrow(Diff), 3)
               for (i in 1:nrow(Diff)) {
                    adaptchange[i,1:3] <- as.vector(quantile(Diff[i,],
                         probs=c(0.025, 0.500, 0.975)))}
               lines(adaptchange[,2], col=rgb(col2rgb(n)[1],
                    col2rgb(n)[2],col2rgb(n)[3],50, maxColorValue=255))}
          }
     if(PDF == TRUE) dev.off()
     }

#End
