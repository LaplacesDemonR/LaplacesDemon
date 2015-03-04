###########################################################################
# plot.laplace.ppc                                                        #
#                                                                         #
# The purpose of the plot.laplace.ppc function is to plot an object of    #
# class laplace.ppc.                                                      #
###########################################################################

plot.laplace.ppc <- function(x, Style=NULL, Data=NULL, Rows=NULL,
     PDF=FALSE, ...)
     {
     ### Initial Checks
     if(missing(x)) stop("The x argument is required.")
     if(class(x) != "laplace.ppc") stop("x is not of class laplace.ppc.")
     if(is.null(Style)) Style <- "Density"
     if(is.null(Rows)) Rows <- 1:nrow(x[["yhat"]])
     ### Plots
     if(Style == "Covariates") {
          if(PDF == TRUE) {
               pdf("PPC.Plots.Covariates.pdf")
               par(mfrow=c(3,3))}
          else par(mfrow=c(3,3), ask=TRUE)
          if(is.null(Data))
               stop("Data is required for Style=Covariates.")
          if(is.null(Data[["X"]]) & is.null(Data[["x"]]))
               stop("X or x is required in Data.")
          if(is.null(Data[["X"]]))
               co <- matrix(Data[["x"]], length(Data[["x"]]), 1)
          else if(is.null(Data[["x"]])) co <- Data[["X"]]
          temp <- summary(x, Quiet=TRUE)$Summary
          mycol <- rgb(0, 100, 0, 50, maxColorValue=255)
          for (i in 1:ncol(co)) {
               plot(co[Rows,i], temp[Rows,5], col=mycol, pch=16, cex=0.75,
                    ylim=c(min(temp[Rows,c(1,4:6)]),max(temp[Rows,c(1,4:6)])),
                    xlab=paste("X[,",i,"]", sep=""),
                    ylab="yhat",
                    sub="Gray lines are yhat at 2.5% and 95%.")
               panel.smooth(co[Rows,i], temp[Rows,5], col=mycol, pch=16,
                    cex=0.75)}}
     if(Style == "Covariates, Categorical DV") {
          if(PDF == TRUE) {
               pdf("PPC.Plots.Covariates.Cat.pdf")
               par(mfrow=c(3,3))}
          else par(mfrow=c(3,3), ask=TRUE)
          if(is.null(Data))
               stop("Data is required for Style=Covariates.")
          if(is.null(Data[["X"]]) & is.null(Data[["x"]]))
               stop("X or x is required in Data.")
          if(is.null(Data[["X"]]))
               co <- matrix(Data[["x"]], length(Data[["x"]]), 1)
          else if(is.null(Data[["x"]])) co <- Data[["X"]]
          temp <- summary(x, Categorical=TRUE, Quiet=TRUE)$Summary
          ncat <- length(table(temp[,1]))
          mycol <- rgb(0, 100, 0, 50, maxColorValue=255)
          for (i in 1:ncol(co)) {for (j in 2:(ncat+1)) {
               plot(co[Rows,i], temp[Rows,j], col=mycol, pch=16, cex=0.75,
                    xlab=paste("X[,",i,"]", sep=""),
                    ylab=colnames(temp)[j])
               panel.smooth(co[Rows,i], temp[Rows,j], col=mycol, pch=16,
                    cex=0.75)}}}
     if(Style == "Density") {
          if(PDF == TRUE) {
               pdf("PPC.Plots.Density.pdf")
               par(mfrow=c(3,3))}
          else par(mfrow=c(3,3), ask=TRUE)
          for (j in 1:length(Rows)) {
               plot(density(x[["yhat"]][Rows[j],]),
                    main=paste("Post. Pred. Plot of yhat[", Rows[j],
                         ",]", sep=""), xlab="Value",
                    sub="Black=Density, Red=y")
               polygon(density(x[["yhat"]][Rows[j],]), col="black",
                    border="black")
               abline(v=x[["y"]][Rows[j]], col="red")}}
     if(Style == "DW") {
          if(PDF == TRUE) pdf("PPC.Plots.DW.pdf")
          par(mfrow=c(1,1))
          epsilon.obs <- x[["y"]] - x[["yhat"]]
          N <- nrow(epsilon.obs)
          S <- ncol(epsilon.obs)
          epsilon.rep <- matrix(rnorm(N*S), N, S)
          d.obs <- d.rep <-  rep(0, S)
          for (s in 1:S) {
               d.obs[s] <- sum(c(0,diff(epsilon.obs[,s]))^2, na.rm=TRUE) / 
                    sum(epsilon.obs[,s]^2, na.rm=TRUE)
               d.rep[s] <- sum(c(0,diff(epsilon.rep[,s]))^2, na.rm=TRUE) / 
                    sum(epsilon.rep[,s]^2, na.rm=TRUE)}
          result <- "no"
          if(mean(d.obs > d.rep, na.rm=TRUE) < 0.025) result <- "positive"
          if(mean(d.obs > d.rep, na.rm=TRUE) > 0.975) result <- "negative"
          d.d.obs <- density(d.obs, na.rm=TRUE)
          d.d.rep <- density(d.rep, na.rm=TRUE)
          plot(d.d.obs, xlim=c(0,4),
               ylim=c(0, max(d.d.obs$y, d.d.rep$y)), col="white",
               main="Durbin-Watson test", 
               xlab=paste("d.obs=", round(mean(d.obs, na.rm=TRUE),2), " (",
               round(as.vector(quantile(d.obs, probs=0.025, na.rm=TRUE)),2),
               ", ", round(as.vector(quantile(d.obs, probs=0.975, na.rm=TRUE)),
               2), "), p(d.obs > d.rep) = ", round(mean(d.obs > d.rep,
               na.rm=TRUE),3), " = ", result, " autocorrelation", sep=""))
          polygon(d.d.obs, col=rgb(0,0,0,50,maxColorValue=255), border=NA)
          polygon(d.d.rep, col=rgb(255,0,0,50,maxColorValue=255), border=NA)
          abline(v=2, col="red")}
     if(Style == "DW, Multivariate, C") {
          if(PDF == TRUE) {
               pdf("PPC.Plots.DW.M.pdf")
               par(mfrow=c(1,1))}
          else par(mfrow=c(1,1), ask=TRUE)
          if(is.null(Data))
               stop("Data is required for Style=Fitted, Multivariate, C.")
          if(is.null(Data[["Y"]])) stop("Y is required in Data.")
          M <- nrow(Data[["Y"]])
          J <- ncol(Data[["Y"]])
          epsilon.obs <- x[["y"]] - x[["yhat"]]
          N <- nrow(epsilon.obs)
          S <- ncol(epsilon.obs)
          epsilon.rep <- matrix(rnorm(N*S), N, S)
          d.obs <- d.rep <-  rep(0, S)
          for (j in 1:J) {
               for (s in 1:S) {
                    d.obs[s] <- sum(c(0,diff(epsilon.obs[((j-1)*M+1):(j*M),s]))^2, na.rm=TRUE) / 
                         sum(epsilon.obs[((j-1)*M+1):(j*M),s]^2, na.rm=TRUE)
                    d.rep[s] <- sum(c(0,diff(epsilon.rep[((j-1)*M+1):(j*M),s]))^2, na.rm=TRUE) / 
                         sum(epsilon.rep[((j-1)*M+1):(j*M),s]^2, na.rm=TRUE)}
               result <- "no"
               if(mean(d.obs > d.rep, na.rm=TRUE) < 0.025) result <- "positive"
               if(mean(d.obs > d.rep, na.rm=TRUE) > 0.975) result <- "negative"
               d.d.obs <- density(d.obs, na.rm=TRUE)
               d.d.rep <- density(d.rep, na.rm=TRUE)
               plot(d.d.obs, xlim=c(0,4),
                    ylim=c(0, max(d.d.obs$y, d.d.rep$y)), col="white",
                    main="Durbin-Watson test", 
                    xlab=paste("d.obs=", round(mean(d.obs, na.rm=TRUE),2), " (",
                    round(as.vector(quantile(d.obs, probs=0.025, na.rm=TRUE)),2),
                    ", ", round(as.vector(quantile(d.obs, probs=0.975, na.rm=TRUE)),
                    2), "), p(d.obs > d.rep) = ", round(mean(d.obs > d.rep,
                    na.rm=TRUE),3), " = ", result, " autocorrelation", sep=""),
                    sub=paste("Y[,",j,"]",sep=""))
               polygon(d.d.obs, col=rgb(0,0,0,50,maxColorValue=255),
                    border=NA)
               polygon(d.d.rep, col=rgb(255,0,0,50,maxColorValue=255),
                    border=NA)
               abline(v=2, col="red")}}
     if(Style == "ECDF") {
          if(PDF == TRUE) pdf("PPC.Plots.ECDF.pdf")
          par(mfrow=c(1,1))
          plot(ecdf(x[["y"]][Rows]), verticals=TRUE, do.points=FALSE,
               main="Cumulative Fit",
               xlab="y (black) and yhat (red; gray)",
               ylab="Cumulative Frequency")
          lines(ecdf(apply(x[["yhat"]][Rows,], 1, quantile, probs=0.975)),
               verticals=TRUE, do.points=FALSE, col="gray")
          lines(ecdf(apply(x[["yhat"]][Rows,], 1, quantile, probs=0.025)),
               verticals=TRUE, do.points=FALSE, col="gray")
          lines(ecdf(apply(x[["yhat"]][Rows,], 1, quantile, probs=0.500)),
               verticals=TRUE, do.points=FALSE, col="red")}
     if(Style == "Fitted") {
          if(PDF == TRUE) pdf("PPC.Plots.Fitted.pdf")
          par(mfrow=c(1,1))
          temp <- summary(x, Quiet=TRUE)$Summary
          plot(temp[Rows,1], temp[Rows,5], pch=16, cex=0.75,
               ylim=c(min(temp[Rows,4], na.rm=TRUE),
               max(temp[Rows,6], na.rm=TRUE)),
               xlab="y", ylab="yhat", main="Fitted")
          for (i in Rows) {
               lines(c(temp[Rows[i],1], temp[Rows[i],1]),
                    c(temp[Rows[i],4], temp[Rows[i],6]))}
          panel.smooth(temp[Rows,1], temp[Rows,5], pch=16, cex=0.75)}
     if(Style == "Fitted, Multivariate, C") {
          if(PDF == TRUE) {
               pdf("PPC.Plots.Fitted.M.pdf")
               par(mfrow=c(1,1))}
          else par(mfrow=c(1,1), ask=TRUE)
          if(is.null(Data))
               stop("Data is required for Style=Fitted, Multivariate, C.")
          if(is.null(Data[["Y"]])) stop("Y is required in Data.")
          temp <- summary(x, Quiet=TRUE)$Summary
          for (i in 1:ncol(Data[["Y"]])) {
               temp1 <- as.vector(matrix(temp[,1], nrow(Data[["Y"]]),
                    ncol(Data[["Y"]]))[,i])
               temp2 <- as.vector(matrix(temp[,4], nrow(Data[["Y"]]),
                    ncol(Data[["Y"]]))[,i])
               temp3 <- as.vector(matrix(temp[,5], nrow(Data[["Y"]]),
                    ncol(Data[["Y"]]))[,i])
               temp4 <- as.vector(matrix(temp[,6], nrow(Data[["Y"]]),
                    ncol(Data[["Y"]]))[,i])
               plot(temp1, temp3, pch=16, cex=0.75,
                    ylim=c(min(temp2, na.rm=TRUE),
                         max(temp4, na.rm=TRUE)),
                    xlab=paste("Y[,", i, "]", sep=""), ylab="yhat",
                    main="Fitted")
               for (j in 1:nrow(Data[["Y"]])) {
                    lines(c(temp1[j], temp1[j]),
                         c(temp2[j], temp4[j]))}
               panel.smooth(temp1, temp3, pch=16, cex=0.75)}}
     if(Style == "Fitted, Multivariate, R") {
          if(PDF == TRUE) {
               pdf("PPC.Plots.Fitted.M.pdf")
               par(mfrow=c(1,1))}
          else par(mfrow=c(1,1), ask=TRUE)
          if(is.null(Data))
               stop("Data is required for Style=Fitted, Multivariate, R.")
          if(is.null(Data[["Y"]])) stop("Y is required in Data.")
          temp <- summary(x, Quiet=TRUE)$Summary
          for (i in 1:nrow(Data[["Y"]])) {
               temp1 <- as.vector(matrix(temp[,1], nrow(Data[["Y"]]),
                    ncol(Data[["Y"]]))[i,])
               temp2 <- as.vector(matrix(temp[,4], nrow(Data[["Y"]]),
                    ncol(Data[["Y"]]))[i,])
               temp3 <- as.vector(matrix(temp[,5], nrow(Data[["Y"]]),
                    ncol(Data[["Y"]]))[i,])
               temp4 <- as.vector(matrix(temp[,6], nrow(Data[["Y"]]),
                    ncol(Data[["Y"]]))[i,])
               plot(temp1, temp3, pch=16, cex=0.75,
                    ylim=c(min(temp2, na.rm=TRUE),
                         max(temp4, na.rm=TRUE)),
                    xlab=paste("Y[,", i, "]", sep=""), ylab="yhat",
                    main="Fitted")
               for (j in 1:ncol(Data[["Y"]])) {
                    lines(c(temp1[j], temp1[j]),
                         c(temp2[j], temp4[j]))}
               panel.smooth(temp1, temp3, pch=16, cex=0.75)}}
     if(Style == "Jarque-Bera") {
          if(PDF == TRUE) pdf("PPC.Plots.Jarque.Bera.pdf")
          par(mfrow=c(1,1))
          epsilon.obs <- epsilon.rep <- x[["y"]][Rows] - x[["yhat"]][Rows,]
          kurtosis <- function(x) {  
               m4 <- mean((x-mean(x, na.rm=TRUE))^4, na.rm=TRUE) 
               kurt <- m4/(sd(x, na.rm=TRUE)^4)-3  
               return(kurt)}
          skewness <-  function(x) {
               m3 <- mean((x-mean(x, na.rm=TRUE))^3, na.rm=TRUE)
               skew <- m3/(sd(x, na.rm=TRUE)^3)
               return(skew)}
          JB.obs <- JB.rep <- rep(0, ncol(epsilon.obs))
          N <- nrow(epsilon.obs)
          for (s in 1:ncol(epsilon.obs)) {
               epsilon.rep[,s] <- rnorm(N, mean(epsilon.obs[,s],
                    na.rm=TRUE), sd(epsilon.obs[,s], na.rm=TRUE))
               K.obs <- kurtosis(epsilon.obs[,s])
               S.obs <- skewness(epsilon.obs[,s])
               K.rep <- kurtosis(epsilon.rep[,s])
               S.rep <- skewness(epsilon.rep[,s])
               JB.obs[s] <- (N/6)*(S.obs^2 + ((K.obs-3)^2)/4)
               JB.rep[s] <- (N/6)*(S.rep^2 + ((K.rep-3)^2)/4)}
          p <- round(mean(JB.obs > JB.rep, na.rm=TRUE), 3)
          result <- "Non-Normality"
          if((p >= 0.025) & (p <= 0.975)) result <- "Normality"
          d.obs <- density(JB.obs)
          d.rep <- density(JB.rep)
          plot(d.obs, xlim=c(min(d.obs$x,d.rep$x), max(d.obs$x,d.rep$x)),
               ylim=c(0, max(d.obs$y, d.rep$y)), col="white",
               main="Jarque-Bera Test",
               xlab="JB", ylab="Density",
               sub=paste("JB.obs=", round(mean(JB.obs, na.rm=TRUE),2),
               " (", round(as.vector(quantile(JB.obs, probs=0.025,
               na.rm=TRUE)),2), ",", round(as.vector(quantile(JB.obs,
               probs=0.975, na.rm=TRUE)),2), "), p(JB.obs > JB.rep) = ",
               p, " = ", result, sep=""))
          polygon(d.obs, col=rgb(0,0,0,50,maxColorValue=255), border=NA)
          polygon(d.rep, col=rgb(255,0,0,50,maxColorValue=255), border=NA)}
     if(Style == "Jarque-Bera, Multivariate, C") {
          if(PDF == TRUE) {
               pdf("PPC.Plots.Jarque.Bera.pdf")
               par(mfrow=c(1,1))}
          else par(mfrow=c(1,1), ask=TRUE)
          if(is.null(Data))
               stop("Data is required for Style=Jarque-Bera, Multivariate, C.")
          if(is.null(Data[["Y"]])) stop("Y is required in Data.")
          M <- nrow(Data[["Y"]])
          J <- ncol(Data[["Y"]])
          epsilon.obs <- epsilon.rep <- x[["y"]] - x[["yhat"]]
          kurtosis <- function(x) {  
               m4 <- mean((x-mean(x, na.rm=TRUE))^4, na.rm=TRUE) 
               kurt <- m4/(sd(x, na.rm=TRUE)^4)-3  
               return(kurt)}
          skewness <-  function(x) {
               m3 <- mean((x-mean(x, na.rm=TRUE))^3, na.rm=TRUE)
               skew <- m3/(sd(x, na.rm=TRUE)^3)
               return(skew)}
          JB.obs <- JB.rep <- rep(0, ncol(epsilon.obs))
          N <- nrow(epsilon.obs)
          for (j in 1:J) {
               for (s in 1:ncol(epsilon.obs)) {
                    e.obs <- matrix(epsilon.obs[,s], M, J)
                    e.rep <- rnorm(M, mean(e.obs[,j], na.rm=TRUE),
                         sd(e.obs[,j], na.rm=TRUE))
                    K.obs <- kurtosis(e.obs[,j])
                    S.obs <- skewness(e.obs[,j])
                    K.rep <- kurtosis(e.rep)
                    S.rep <- skewness(e.rep)
                    JB.obs[s] <- (N/6)*(S.obs^2 + ((K.obs-3)^2)/4)
                    JB.rep[s] <- (N/6)*(S.rep^2 + ((K.rep-3)^2)/4)}
               p <- round(mean(JB.obs > JB.rep, na.rm=TRUE), 3)
               result <- "Non-Normality"
               if((p >= 0.025) & (p <= 0.975)) result <- "Normality"
               d.obs <- density(JB.obs)
               d.rep <- density(JB.rep)
               plot(d.obs, xlim=c(min(d.obs$x,d.rep$x), max(d.obs$x,d.rep$x)),
                    ylim=c(0, max(d.obs$y, d.rep$y)), col="white",
                    main="Jarque-Bera Test",
                    xlab=paste("JB for Y[,",j,"]", sep=""), ylab="Density",
                    sub=paste("JB.obs=", round(mean(JB.obs, na.rm=TRUE),2),
                    " (", round(as.vector(quantile(JB.obs, probs=0.025,
                    na.rm=TRUE)),2), ",", round(as.vector(quantile(JB.obs,
                    probs=0.975, na.rm=TRUE)),2), "), p(JB.obs > JB.rep) = ",
                    p, " = ", result, sep=""))
               polygon(d.obs, col=rgb(0,0,0,50,maxColorValue=255),
                    border=NA)
               polygon(d.rep, col=rgb(255,0,0,50,maxColorValue=255),
                    border=NA)}}
     if(Style == "Mardia") {
          if(PDF == TRUE) pdf("PPC.Plots.Mardia.pdf")
          par(mfrow=c(2,1))
          if(is.null(Data))
               stop("Data is required for Style=Mardia, C.")
          if(is.null(Data[["Y"]]))
               stop("Variable Y is required for Style=Mardia, C.")
          epsilon.obs <- x[["y"]] - x[["yhat"]]
          M <- nrow(Data[["Y"]])
          J <- ncol(Data[["Y"]])
          K3.obs <- K3.rep <- K4.obs <- K4.rep <- rep(0, ncol(epsilon.obs))
          for (s in 1:ncol(epsilon.obs)) {
               e.obs <- matrix(epsilon.obs[,s], M, J)
               e.obs.mu <- colMeans(e.obs)
               e.obs.mu.mat <- matrix(e.obs.mu, M, J, byrow=TRUE)
               e.obs.stand <- e.obs - e.obs.mu.mat
               S.obs <- var(e.obs)
               A.obs <- t(chol(S.obs))
               A.inv.obs <- solve(A.obs)
               Z.obs <- t(A.inv.obs %*% t(e.obs.stand))
               Dij.obs <- Z.obs %*% t(Z.obs)
               D2.obs <- diag(Dij.obs)
               K3.obs[s] <- mean(as.vector(Dij.obs)^3)
               K4.obs[s] <- mean(D2.obs^2)
               e.rep <- rmvn(M, e.obs.mu.mat, S.obs)
               e.rep.mu <- colMeans(e.rep)
               e.rep.mu.mat <- matrix(e.rep.mu, M, J, byrow=TRUE)
               e.rep.stand <- e.rep - e.rep.mu.mat
               S.rep <- var(e.rep)
               A.rep <- t(chol(S.rep))
               A.inv.rep <- solve(A.rep)
               Z.rep <- t(A.inv.rep %*% t(e.rep.stand))
               Dij.rep <- Z.rep %*% t(Z.rep)
               D2.rep <- diag(Dij.rep)
               K3.rep[s] <- mean(as.vector(Dij.rep)^3)
               K4.rep[s] <- mean(D2.rep^2)}
          p.K3 <- round(mean(K3.obs > K3.rep), 3)
          p.K4 <- round(mean(K4.obs > K4.rep), 3)
          K3.result <- K4.result <- "Non-Normality"
          if((p.K3 >= 0.025) & (p.K3 <= 0.975)) K3.result <- "Normality"
          if((p.K4 >= 0.025) & (p.K4 <= 0.975)) K4.result <- "Normality"
          d.K3.obs <- density(K3.obs)
          d.K3.rep <- density(K3.rep)
          d.K4.obs <- density(K4.obs)
          d.K4.rep <- density(K4.rep)
          plot(d.K3.obs, xlim=c(min(d.K3.obs$x, d.K3.rep$x),
               max(d.K3.obs$x, d.K3.rep$x)),
               ylim=c(0, max(d.K3.obs$y, d.K3.rep$y)), col="white",
               main="Mardia's Test of MVN Skewness",
               xlab="Skewness Test Statistic (K3)", ylab="Density",
               sub=paste("K3.obs=", round(mean(K3.obs, na.rm=TRUE), 2),
                    " (", round(quantile(K3.obs, probs=0.025, na.rm=TRUE),
                    2), ", ", round(quantile(K3.obs, probs=0.975,
                    na.rm=TRUE), 2), "), p(K3.obs > K3.rep) = ",
                    p.K3, " = ", K3.result, sep=""))
          polygon(d.K3.obs, col=rgb(0,0,0,50,maxColorValue=255), border=NA)
          polygon(d.K3.rep, col=rgb(255,0,0,50,maxColorValue=255), border=NA)
          plot(d.K4.obs, xlim=c(min(d.K4.obs$x, d.K4.rep$x),
               max(d.K4.obs$x, d.K4.rep$x)),
               ylim=c(0, max(d.K4.obs$y, d.K4.rep$y)), col="white",
               main="Mardia's Test of MVN Kurtosis",
               xlab="Kurtosis Test Statistic (K4)", ylab="Density",
               sub=paste("K4.obs=", round(mean(K4.obs, na.rm=TRUE), 2),
                    " (", round(quantile(K4.obs, probs=0.025, na.rm=TRUE),
                    2), ", ", round(quantile(K4.obs, probs=0.975,
                    na.rm=TRUE), 2), "), p(K4.obs > K4.rep) = ",
                    p.K4, " = ", K4.result, sep=""))
          polygon(d.K4.obs, col=rgb(0,0,0,50,maxColorValue=255), border=NA)
          polygon(d.K4.rep, col=rgb(255,0,0,50,maxColorValue=255), border=NA)}
     if(Style == "Predictive Quantiles") {
          if(PDF == TRUE) pdf("PPC.Plots.PQ.pdf")
          par(mfrow=c(1,1))
          temp <- summary(x, Quiet=TRUE)$Summary
          mycol <- rgb(0, 100, 0, 50, maxColorValue=255)
          plot(temp[Rows,1], temp[Rows,7], ylim=c(0,1), col=mycol,
               pch=16, cex=0.75, xlab="y", ylab="PQ",
               main="Predictive Quantiles")
          panel.smooth(temp[Rows,1], temp[Rows,7], col=mycol, pch=16,
               cex=0.75)
          abline(h=0.025, col="gray")
          abline(h=0.975, col="gray")}
     if(Style == "Residual Density") {
          if(PDF == TRUE) pdf("PPC.Plots.Residual.Density.pdf")
          par(mfrow=c(1,1))
          epsilon <- x[["y"]] - x[["yhat"]]
          epsilon.summary <- apply(epsilon, 1, quantile,
               probs=c(0.025,0.500,0.975), na.rm=TRUE)
          dens <- density(epsilon.summary[2,Rows], na.rm=TRUE)
          plot(dens, col="black", main="Residual Density",
               xlab=expression(epsilon), ylab="Density")
          polygon(dens, col="black", border="black")
          abline(v=0, col="red")}
     if(Style == "Residual Density, Multivariate, C") {
          if(PDF == TRUE) {
               pdf("PPC.Plots.Residual.Density.pdf")
               par(mfrow=c(1,1))}
          else par(mfrow=c(1,1), ask=TRUE)
          if(is.null(Data))
               stop("Data is required for Style=Residual Density, Multivariate, C.")
          if(is.null(Data[["Y"]]))
               stop("Variable Y is required for Style=Residual Density, Multivariate, C.")
          epsilon <- x[["y"]] - x[["yhat"]]
          epsilon.summary <- apply(epsilon, 1, quantile,
               probs=c(0.025,0.500,0.975), na.rm=TRUE)
          epsilon.500 <- matrix(epsilon.summary[2,], nrow(Data[["Y"]]),
               ncol(Data[["Y"]]))
          for (i in 1:ncol(Data[["Y"]])) {
               dens <- density(epsilon.500[,i], na.rm=TRUE)
               plot(dens, col="black", main="Residual Density",
                    xlab=paste("epsilon[,", i, "]", sep=""),
                    ylab="Density")
               polygon(dens, col="black", border="black")
               abline(v=0, col="red")}}
     if(Style == "Residual Density, Multivariate, R") {
          if(PDF == TRUE) {
               pdf("PPC.Plots.Residual.Density.pdf")
               par(mfrow=c(1,1))}
          else par(mfrow=c(1,1), ask=TRUE)
          if(is.null(Data))
               stop("Data is required for Style=Residual Density, Multivariate, R.")
          if(is.null(Data[["Y"]]))
               stop("Variable Y is required for Style=Residual Density, Multivariate, R.")
          epsilon <- x[["y"]] - x[["yhat"]]
          epsilon.summary <- apply(epsilon, 1, quantile,
               probs=c(0.025,0.500,0.975), na.rm=TRUE)
          epsilon.500 <- matrix(epsilon.summary[2,], nrow(Data[["Y"]]),
               ncol(Data[["Y"]]))
          for (i in 1:nrow(Data[["Y"]])) {
               dens <- density(epsilon.500[i,], na.rm=TRUE)
               plot(dens, col="black", main="Residual Density",
                    xlab=paste("epsilon[", i, ",]", sep=""),
                    ylab="Density")
               polygon(dens, col="black", border="black")
               abline(v=0, col="red")}}
     if(Style == "Residuals") {
          if(PDF == TRUE) pdf("PPC.Plots.Residuals.pdf")
          par(mfrow=c(1,1))
          epsilon <- x[["y"]] - x[["yhat"]]
          epsilon.summary <- apply(epsilon, 1, quantile,
               probs=c(0.025,0.500,0.975), na.rm=TRUE)
          plot(epsilon.summary[2,Rows], pch=16, cex=0.75,
               ylim=c(min(epsilon.summary[,Rows], na.rm=TRUE),
                    max(epsilon.summary[,Rows], na.rm=TRUE)),
               xlab="y", ylab=expression(epsilon))
          lines(rep(0, ncol(epsilon.summary[,Rows])), col="red")
          for (i in Rows) {
               lines(c(i,i), c(epsilon.summary[1,Rows[i]],
                    epsilon.summary[3,Rows[i]]), col="black")}}
     if(Style == "Residuals, Multivariate, C") {
          if(PDF == TRUE) {
               pdf("PPC.Plots.Residuals.pdf")
               par(mfrow=c(1,1))}
          else par(mfrow=c(1,1), ask=TRUE)
          if(is.null(Data))
               stop("Data is required for Style=Residuals, Multivariate, C.")
          if(is.null(Data[["Y"]]))
               stop("Variable Y is required for Style=Residuals, Multivariate, C.")
          epsilon <- x[["y"]] - x[["yhat"]]
          epsilon.summary <- apply(epsilon, 1, quantile,
               probs=c(0.025,0.500,0.975), na.rm=TRUE)
          epsilon.025 <- matrix(epsilon.summary[1,], nrow(Data[["Y"]]),
               ncol(Data[["Y"]]))
          epsilon.500 <- matrix(epsilon.summary[2,], nrow(Data[["Y"]]),
               ncol(Data[["Y"]]))
          epsilon.975 <- matrix(epsilon.summary[3,], nrow(Data[["Y"]]),
               ncol(Data[["Y"]]))
          for (i in 1:ncol(Data[["Y"]])) {
               plot(epsilon.500[,i], pch=16, cex=0.75,
                    ylim=c(min(epsilon.025[,i], na.rm=TRUE),
                         max(epsilon.975[,i], na.rm=TRUE)),
                    xlab=paste("Y[,", i, "]", sep=""), ylab=expression(epsilon))
               lines(rep(0, nrow(epsilon.500)), col="red")
               for (j in 1:nrow(Data[["Y"]])) {
                    lines(c(j,j), c(epsilon.025[j,i],
                         epsilon.975[j,i]), col="black")}}}
     if(Style == "Residuals, Multivariate, R") {
          if(PDF == TRUE) {
               pdf("PPC.Plots.Residuals.pdf")
               par(mfrow=c(1,1))}
          else par(mfrow=c(1,1), ask=TRUE)
          if(is.null(Data))
               stop("Data is required for Style=Residuals, Multivariate, C.")
          if(is.null(Data[["Y"]]))
               stop("Variable Y is required for Style=Residuals, Multivariate, C.")
          epsilon <- x[["y"]] - x[["yhat"]]
          epsilon.summary <- apply(epsilon, 1, quantile,
               probs=c(0.025,0.500,0.975), na.rm=TRUE)
          epsilon.025 <- matrix(epsilon.summary[1,], nrow(Data[["Y"]]),
               ncol(Data[["Y"]]))
          epsilon.500 <- matrix(epsilon.summary[2,], nrow(Data[["Y"]]),
               ncol(Data[["Y"]]))
          epsilon.975 <- matrix(epsilon.summary[3,], nrow(Data[["Y"]]),
               ncol(Data[["Y"]]))
          for (i in 1:nrow(Data[["Y"]])) {
               plot(epsilon.500[i,], pch=16, cex=0.75,
                    ylim=c(min(epsilon.025[i,], na.rm=TRUE),
                         max(epsilon.975[i,], na.rm=TRUE)),
                    xlab=paste("Y[", i, ",]", sep=""), ylab=expression(epsilon))
               lines(rep(0, ncol(epsilon.500)), col="red")
               for (j in 1:ncol(Data[["Y"]])) {
                    lines(c(j,j), c(epsilon.025[i,j],
                         epsilon.975[i,j]), col="black")}}}
     if(Style == "Space-Time by Space") {
          if(PDF == TRUE) {
               pdf("PPC.Plots.SpaceTime.pdf")
               par(mfrow=c(1,1))}
          else par(mfrow=c(1,1), ask=TRUE)
          if(is.null(Data))
               stop("Data is required for Style=Space-Time by Space.")
          if(is.null(Data[["longitude"]]))
               stop("Variable longitude is required in Data.")
          if(is.null(Data[["latitude"]]))
               stop("Variable latitude is required in Data.")
          if(is.null(Data[["S"]])) stop("Variable S is required in Data.")
          if(is.null(Data[["T"]])) stop("Variable T is required in Data.")
          temp <- summary(x, Quiet=TRUE)$Summary
          for (s in 1:Data[["S"]]) {
               plot(matrix(temp[,1], Data[["S"]], Data[["T"]])[s,],
                    ylim=c(min(c(matrix(temp[,4], Data[["S"]], Data[["T"]])[s,],
                         matrix(temp[,1], Data[["S"]], Data[["T"]])[s,]), na.rm=TRUE),
                         max(c(matrix(temp[,6], Data[["S"]], Data[["T"]])[s,],
                         matrix(temp[,1], Data[["S"]], Data[["T"]])[s,]), na.rm=TRUE)),
                    type="l", xlab="Time", ylab="y",
                    main=paste("Space-Time at Space s=",s," of ",
                         Data[["S"]], sep=""),
                    sub="Actual=Black, Fit=Red, Interval=Transparent Red")
               polygon(c(1:Data[["T"]],rev(1:Data[["T"]])),
                    c(matrix(temp[,4], Data[["S"]], Data[["T"]])[s,],
                    rev(matrix(temp[,6], Data[["S"]], Data[["T"]])[s,])),
                    col=rgb(255, 0, 0, 50, maxColorValue=255), border=FALSE)
               lines(matrix(temp[,5], Data[["S"]], Data[["T"]])[s,], col="red")}}
     if(Style == "Space-Time by Time") {
          if(PDF == TRUE) {
               pdf("PPC.Plots.SpaceTime.pdf")
               par(mfrow=c(1,1))}
          else par(mfrow=c(1,1), ask=TRUE)
          if(is.null(Data))
               stop("Data is required for Style=Space-Time by Time.")
          if(is.null(Data[["longitude"]]))
               stop("Variable longitude is required in Data.")
          if(is.null(Data[["latitude"]]))
               stop("Variable latitude is required in Data.")
          if(is.null(Data[["S"]])) stop("Variable S is required in Data.")
          if(is.null(Data[["T"]])) stop("Variable T is required in Data.")
          Heat <- (1-(x[["y"]]-min(x[["y"]], na.rm=TRUE)) /
               max(x[["y"]]-min(x[["y"]], na.rm=TRUE), na.rm=TRUE)) * 99 + 1
          Heat <- matrix(Heat, Data[["S"]], Data[["T"]])
          for (t in 1:Data[["T"]]) {
               plot(Data[["longitude"]], Data[["latitude"]],
                    col=heat.colors(120)[Heat[,t]],
                    pch=16, cex=0.75, xlab="Longitude", ylab="Latitude",
                    main=paste("Space-Time at t=",t," of ", Data[["T"]],
                         sep=""), sub="Red=High, Yellow=Low")}}
     if(Style == "Spatial") {
          if(PDF == TRUE) pdf("PPC.Plots.Spatial.pdf")
          par(mfrow=c(1,1))
          if(is.null(Data)) stop("Data is required for Style=Spatial.")
          if(is.null(Data[["longitude"]]))
               stop("Variable longitude is required in Data.")
          if(is.null(Data[["latitude"]]))
               stop("Variable latitude is required in Data.")
          heat <- (1-(x[["y"]][Rows]-min(x[["y"]][Rows], na.rm=TRUE)) /
               max(x[["y"]][Rows]-min(x[["y"]][Rows], na.rm=TRUE), na.rm=TRUE)) * 99 + 1
          plot(Data[["longitude"]][Rows], Data[["latitude"]][Rows],
               col=heat.colors(120)[heat],
               pch=16, cex=0.75, xlab="Longitude", ylab="Latitude",
               main="Spatial Plot", sub="Red=High, Yellow=Low")}
     if(Style == "Spatial Uncertainty") {
          if(PDF == TRUE) pdf("PPC.Plots.Spatial.Unc.pdf")
          par(mfrow=c(1,1))
          if(is.null(Data))
               stop("Data is required for Style=Spatial Uncertainty.")
          if(is.null(Data[["longitude"]]))
               stop("Variable longitude is required in Data.")
          if(is.null(Data[["latitude"]]))
               stop("Variable latitude is required in Data.")
          heat <- apply(x[["yhat"]], 1, quantile, probs=c(0.025,0.975))
          heat <- heat[2,] - heat[1,]
          heat <- (1-(heat[Rows]-min(heat[Rows])) /
               max(heat[Rows]-min(heat[Rows]))) * 99 + 1
          plot(Data[["longitude"]][Rows], Data[["latitude"]][Rows],
               col=heat.colors(120)[heat],
               pch=16, cex=0.75, xlab="Longitude", ylab="Latitude",
               main="Spatial Uncertainty Plot",
               sub="Red=High, Yellow=Low")}
     if(Style == "Time-Series") {
          if(PDF == TRUE) pdf("PPC.Plots.TimeSeries.pdf")
          par(mfrow=c(1,1))
          temp <- summary(x, Quiet=TRUE)$Summary
          plot(Rows, temp[Rows,1],
               ylim=c(min(temp[Rows,c(1,4)], na.rm=TRUE),
                    max(temp[Rows,c(1,6)], na.rm=TRUE)),
               type="l", xlab="Time", ylab="y",
               main="Plot of Fitted Time-Series",
               sub="Actual=Black, Fit=Red, Interval=Transparent Red")
          polygon(c(Rows,rev(Rows)),c(temp[Rows,4],rev(temp[Rows,6])),
                col=rgb(255, 0, 0, 50, maxColorValue=255),
                border=FALSE)
          lines(Rows, temp[Rows,1])
          lines(Rows, temp[Rows,5], col="red")}
     if(Style == "Time-Series, Multivariate, C") {
          if(PDF == TRUE) {
               pdf("PPC.Plots.TimeSeries.pdf")
               par(mfrow=c(1,1))}
          else {par(mfrow=c(1,1), ask=TRUE)}
          if(is.null(Data))
               stop("Data is required for Style=Time-Series, Multivariate.")
          if(is.null(Data[["Y"]]))
               stop("Variable Y is required in Data.")
          temp <- summary(x, Quiet=TRUE)$Summary
          for (i in 1:ncol(Data[["Y"]])) {
          tempy <- matrix(temp[Rows,1], nrow(Data[["Y"]]),
               ncol(Data[["Y"]]))[,i]
          qLB <- matrix(temp[Rows,4], nrow(Data[["Y"]]),
               ncol(Data[["Y"]]))[,i]
          qMed <- matrix(temp[Rows,5], nrow(Data[["Y"]]),
               ncol(Data[["Y"]]))[,i]
          qUB <- matrix(temp[Rows,6], nrow(Data[["Y"]]),
               ncol(Data[["Y"]]))[,i]
          plot(1:length(tempy), tempy,
               ylim=c(min(Data[["Y"]][,i],
                    matrix(temp[Rows,4], nrow(Data[["Y"]]),
                         ncol(Data[["Y"]]))[,i], na.rm=TRUE),
                    max(Data[["Y"]][,i],
                    matrix(temp[Rows,6], nrow(Data[["Y"]]),
                         ncol(Data[["Y"]]))[,i], na.rm=TRUE)),
               type="l", xlab="Time", ylab="y",
               main=paste("Time-Series ", i, " of ", ncol(Data[["Y"]]), sep=""),
               sub="Actual=Black, Fit=Red, Interval=Transparent Red")
          polygon(c(1:length(tempy),rev(1:length(tempy))),c(qLB,rev(qUB)),
                col=rgb(255, 0, 0, 50, maxColorValue=255),
                border=FALSE)
          lines(1:length(tempy), tempy)
          lines(1:length(tempy), qMed, col="red")}}
     if(Style == "Time-Series, Multivariate, R") {
          if(PDF == TRUE) {
               pdf("PPC.Plots.TimeSeries.pdf")
               par(mfrow=c(1,1))}
          else par(mfrow=c(1,1), ask=TRUE)
          if(is.null(Data))
               stop("Data is required for Style=Time-Series, Multivariate.")
          if(is.null(Data[["Y"]]))
               stop("Variable Y is required in Data.")
          temp <- summary(x, Quiet=TRUE)$Summary
          for (i in 1:nrow(Data[["Y"]])) {
          tempy <- matrix(temp[Rows,1], nrow(Data[["Y"]]), ncol(Data[["Y"]]))[i,]
          qLB <- matrix(temp[Rows,4], nrow(Data[["Y"]]), ncol(Data[["Y"]]))[i,]
          qMed <- matrix(temp[Rows,5], nrow(Data[["Y"]]), ncol(Data[["Y"]]))[i,]
          qUB <- matrix(temp[Rows,6], nrow(Data[["Y"]]), ncol(Data[["Y"]]))[i,]
          plot(1:length(tempy), tempy,
               ylim=c(min(Data[["Y"]][i,],
                    matrix(temp[Rows,4], nrow(Data[["Y"]]),
                         ncol(Data[["Y"]]))[i,], na.rm=TRUE),
                    max(Data[["Y"]][i,],
                    matrix(temp[Rows,6], nrow(Data[["Y"]]),
                         ncol(Data[["Y"]]))[i,], na.rm=TRUE)),
               type="l", xlab="Time", ylab="y",
               main=paste("Time-Series ", i, " of ", nrow(Data[["Y"]]), sep=""),
               sub="Actual=Black, Fit=Red, Interval=Transparent Red")
          polygon(c(1:length(tempy),rev(1:length(tempy))),c(qLB,rev(qUB)),
                col=rgb(255, 0, 0, 50, maxColorValue=255),
                border=FALSE)
          lines(1:length(tempy), tempy)
          lines(1:length(tempy), qMed, col="red")}}
     if(PDF == TRUE) dev.off()
     }

#End
