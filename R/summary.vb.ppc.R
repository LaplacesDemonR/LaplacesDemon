###########################################################################
# summary.vb.ppc                                                          #
#                                                                         #
# The purpose of the summary.vb.ppc function is to summarize an object of #
# class vb.ppc (posterior predictive check).                              #
###########################################################################

summary.vb.ppc <- function(object=NULL, Categorical=FALSE, Rows=NULL,
     Discrep=NULL, d=0, Quiet=FALSE, ...)
     {
     if(is.null(object)) stop("The object argument is NULL.")
     y <- object$y
     yhat <- object$yhat
     Deviance <- object$Deviance
     monitor <- object$monitor
     if(is.null(Rows)) Rows <- 1:length(y)
     if(any(Rows > length(y)) || any(Rows <= 0)) {
          warning("Invalid Rows argument; All rows included.")
          Rows <- 1:length(y)}
     ### Create Continuous Summary Table for y and yhat
     if(Categorical == FALSE) {
          Summ <- matrix(NA, length(y), 8, dimnames=list(1:length(y),
               c("y","Mean","SD","LB","Median","UB","PQ","Discrep")))
          Summ[,1] <- y
          Summ[,2] <- round(rowMeans(yhat),3)
          Summ[,3] <- round(sqrt(.rowVars(yhat)),3)
          for(i in 1:length(y))
              {
              Summ[i,4] <- round(quantile(yhat[i,], probs=0.025,
                   na.rm=TRUE),3)
              Summ[i,5] <- round(quantile(yhat[i,], probs=0.500,
                   na.rm=TRUE),3)
              Summ[i,6] <- round(quantile(yhat[i,], probs=0.975,
                   na.rm=TRUE),3)
              Summ[i,7] <- round(mean(yhat[i,] >= y[i], na.rm=TRUE),3)
              }
          ### Discrepancy Statistics
          Concordance <- 1 - mean(({Summ[,7] < 0.025} | {Summ[,7] > 0.975}),
               na.rm=TRUE)
          Discrepancy.Statistic <- 0
          if(!is.null(Discrep) && {Discrep == "Chi-Square"}) {
               Summ[,8] <- round((y - rowMeans(yhat))^2 /
                    .rowVars(yhat),3)
               Discrepancy.Statistic <- round(sum(Summ[,8], na.rm=TRUE),3)}
          if(!is.null(Discrep) && {Discrep == "Chi-Square2"}) {
               chisq.obs <- chisq.rep <- yhat
               E.y <- E.yrep <- rowMeans(yhat, na.rm=TRUE)
               for (i in 1:nrow(yhat)) {
                    chisq.obs[i,] <- (y[i] - E.y[i])^2 / E.y[i]
                    chisq.rep[i,] <- (yhat[i,] - E.yrep[i])^2 / E.yrep[i]
                    }
               Summ[,8] <- round(rowMeans(chisq.rep > chisq.obs,
                    na.rm=TRUE),3)
               Discrepancy.Statistic <- round(mean((Summ[,8] < 0.025) |
                    (Summ[,8] > 0.975), na.rm=TRUE),3)}
          if(!is.null(Discrep) && {Discrep == "Kurtosis"}) {
               kurtosis <- function(x) {  
                    m4 <- mean((x-mean(x, na.rm=TRUE))^4, na.rm=TRUE) 
                    kurt <- m4/(sd(x, na.rm=TRUE)^4)-3  
                    return(kurt)}
               for (i in 1:length(y)) {Summ[i,8] <- round(kurtosis(yhat[i,]),3)}
               Discrepancy.Statistic <- round(mean(Summ[,8], na.rm=TRUE),3)}
          if(!is.null(Discrep) && {Discrep == "L.criterion"}) {
               Summ[,8] <- round(sqrt(.rowVars(yhat) +
                    (y - rowMeans(yhat))^2),3)
               Discrepancy.Statistic <- round(sum(Summ[,8], na.rm=TRUE),3)}
          if(!is.null(Discrep) && {Discrep == "MASE"}) {
               Summ[,8] <- round(abs(rowMeans(y - yhat, na.rm=TRUE) /
                    mean(abs(diff(y)), na.rm=TRUE)), 3)
               Discrepancy.Statistic <- round(mean(Summ[,8], na.rm=TRUE),3)}
          if(!is.null(Discrep) && {Discrep == "MSE"}) {
               Summ[,8] <- round(rowMeans((y - yhat)^2, na.rm=TRUE),3)
               Discrepancy.Statistic <- round(mean(Summ[,8], na.rm=TRUE),3)}
          if(!is.null(Discrep) && {Discrep == "PPL"}) {
               Summ[,8] <- round(.rowVars(yhat) + (d/(d+1)) *
                    (rowMeans(yhat) - y)^2,3)
               Discrepancy.Statistic <- round(sum(Summ[,8], na.rm=TRUE),3)}
          if(!is.null(Discrep) && {Discrep == "Quadratic Loss"}) {
               Summ[,8] <- round(rowMeans((y - yhat)^2, na.rm=TRUE),3)
               Discrepancy.Statistic <- round(mean(Summ[,8], na.rm=TRUE),3)}
          if(!is.null(Discrep) && {Discrep == "Quadratic Utility"}) {
               Summ[,8] <- round(rowMeans(-1*(y - yhat)^2, na.rm=TRUE),3)
               Discrepancy.Statistic <- round(mean(Summ[,8], na.rm=TRUE),3)}
          if(!is.null(Discrep) && {Discrep == "RMSE"}) {
               Summ[,8] <- round(sqrt(rowMeans((y - yhat)^2, na.rm=TRUE)),3)
               Discrepancy.Statistic <- round(mean(Summ[,8], na.rm=TRUE),3)}
          if(!is.null(Discrep) && {Discrep == "Skewness"}) {
               skewness <-  function(x) {
                    m3 <- mean((x-mean(x, na.rm=TRUE))^3, na.rm=TRUE)
                    skew <- m3/(sd(x, na.rm=TRUE)^3)
                    return(skew)}
               for (i in 1:length(y)) {Summ[i,8] <- round(skewness(yhat[i,]),3)}
               Discrepancy.Statistic <- round(mean(Summ[,8], na.rm=TRUE),3)}
          if(!is.null(Discrep) && {Discrep == "max(yhat[i,]) > max(y)"}) {
               for (i in 1:length(y)) {Summ[i,8] <- max(yhat[i,]) > max(y)}
               Discrepancy.Statistic <- round(mean(Summ[,8], na.rm=TRUE),3)}
          if(!is.null(Discrep) && {Discrep == "mean(yhat[i,]) > mean(y)"}) {
               for (i in 1:length(y)) {Summ[i,8] <- mean(yhat[i,]) > mean(y)}
               Discrepancy.Statistic <- round(mean(Summ[,8], na.rm=TRUE),3)}
          if(!is.null(Discrep) && {Discrep == "mean(yhat[i,] > d)"}) {
               for (i in 1:length(y)) {Summ[i,8] <- mean(yhat[i,] > d)}
               Discrepancy.Statistic <- round(mean(Summ[,8], na.rm=TRUE),3)}
          if(!is.null(Discrep) && {Discrep == "mean(yhat[i,] > mean(y))"}) {
               for (i in 1:length(y)) {Summ[i,8] <- mean(yhat[i,] > mean(y))}
               Discrepancy.Statistic <- round(mean(Summ[,8], na.rm=TRUE),3)}
          if(!is.null(Discrep) && {Discrep == "min(yhat[i,]) < min(y)"}) {
               for (i in 1:length(y)) {Summ[i,8] <- min(yhat[i,]) < min(y)}
               Discrepancy.Statistic <- round(mean(Summ[,8], na.rm=TRUE),3)}
          if(!is.null(Discrep) && {Discrep == "round(yhat[i,]) = d"}) {
               for (i in 1:length(y)) {
                    Summ[i,8] <- round(mean(round(yhat[i,]) == d,
                         na.rm=TRUE), 3)}
               Discrepancy.Statistic <- round(mean(Summ[,8], na.rm=TRUE),3)}
          if(!is.null(Discrep) && {Discrep == "sd(yhat[i,]) > sd(y)"}) {
               for (i in 1:length(y)) {Summ[i,8] <- sd(yhat[i,]) > sd(y)}
               Discrepancy.Statistic <- round(mean(Summ[,8], na.rm=TRUE),3)}
          L <- round(sqrt(.rowVars(yhat) + (y - rowMeans(yhat))^2), 3)
          S.L <- round(sd(L, na.rm=TRUE),3); L <- round(sum(L, na.rm=TRUE),3)
          ### Deviance
          Dbar <- round(mean(Deviance, na.rm=TRUE),3)
          pD <- round(var(Deviance, na.rm=TRUE) / 2,3)
          BPIC <- Dbar + 2*pD
          bpic <- matrix(c(Dbar, pD, BPIC), 1, 3)
          colnames(bpic) <- c("Dbar","pD","BPIC"); rownames(bpic) <- ""
          ### Create Summary Table for monitored variables
          Mon <- matrix(NA, nrow(monitor), 5,
               dimnames=list(c(rownames(monitor)),
               c("Mean","SD","LB","Median","UB")))
          for (i in 1:nrow(monitor)) {
               Mon[i,1] <- mean(monitor[i,])
               Mon[i,2] <- round(sd(monitor[i,]),3)
               Mon[i,3] <- round(quantile(monitor[i,], probs=0.025),3)
               Mon[i,4] <- round(quantile(monitor[i,], probs=0.500),3)
               Mon[i,5] <- round(quantile(monitor[i,], probs=0.975),3)
               }
          ### Create Output
          Summ.out <- list(BPIC=bpic,
               Concordance=Concordance,
               Discrepancy.Statistic=round(Discrepancy.Statistic,5),
               L.criterion=L,
               S.L=S.L,
               Summary=Summ[Rows,])
          if(Quiet == FALSE) {
               cat("Bayesian Predictive Information Criterion:\n")
               print(bpic)
               cat("Concordance: ", Concordance, "\n")
               cat("Discrepancy Statistic: ",
                    round(Discrepancy.Statistic,5), "\n")
               cat("L-criterion: ", L, ", S.L: ", S.L, sep="", "\n")
               cat("Monitors:\n")
               print(Mon)
               cat("\n\nRecords:\n")
               print(Summ[Rows,])}
          }
     ### Create Categorical Summary Table
     else {
          catcounts <- table(y)
          sumnames <- rep(NA, length(catcounts)+3)
          sumnames[1] <- "y"
          for (i in 1:length(catcounts)) {
               sumnames[i+1] <- paste("p(yhat=",names(catcounts)[i],")",sep="")}
          sumnames[length(sumnames)-1] <- "Lift"
          sumnames[length(sumnames)] <- "Discrep"
          Summ <- matrix(NA, length(y), length(sumnames),
               dimnames=list(1:length(y), sumnames))
          Summ[,1] <- y
          for (i in 1:length(catcounts)) {
               Summ[,i+1] <- rowSums(yhat == as.numeric(names(catcounts)[i])) /
                    ncol(yhat)}
          Summ[,{ncol(Summ)-1}] <- 1
          for (i in 1:length(y)) {
               Summ[i,{ncol(Summ)-1}] <- Summ[i,
                    grep(Summ[i,1],names(catcounts))+1] / 
                    {as.vector(catcounts[grep(Summ[i,1],names(catcounts))]) /
                    sum(catcounts)} - 1}
          ### Discrepancy Statistics
          Mean.Lift <- round(mean(Summ[,{ncol(Summ)-1}]),3)
          Discrepancy.Statistic <- 0
          if(!is.null(Discrep) && {Discrep == "p(yhat[i,] != y[i])"}) {
               for (i in 1:length(y)) { Summ[i,ncol(Summ)] <- 1 - 
                    Summ[i, grep(Summ[i,1],names(catcounts))+1]}
               Discrepancy.Statistic <- round(mean(Summ[,ncol(Summ)],
                    na.rm=TRUE),3)}
          ### Deviance
          Dbar <- round(mean(Deviance, na.rm=TRUE),3)
          pD <- round(var(Deviance, na.rm=TRUE) / 2,3)
          BPIC <- Dbar + 2*pD
          bpic <- matrix(c(Dbar, pD, BPIC), 1, 3)
          colnames(bpic) <- c("Dbar","pD","BPIC"); rownames(bpic) <- ""
          ### Create Summary Table for monitored variables
          Mon <- matrix(NA, nrow(monitor), 5,
               dimnames=list(c(rownames(monitor)),
               c("Mean","SD","LB","Median","UB")))
          for (i in 1:nrow(monitor)) {
               Mon[i,1] <- mean(monitor[i,])
               Mon[i,2] <- sd(monitor[i,])
               Mon[i,3] <- quantile(monitor[i,], probs=0.025)
               Mon[i,4] <- quantile(monitor[i,], probs=0.500)
               Mon[i,5] <- quantile(monitor[i,], probs=0.975)
               }
          ### Create Output
          Summ.out <- list(BPIC=bpic,
               Mean.Lift=Mean.Lift,
               Discrepancy.Statistic=round(Discrepancy.Statistic,5),
               Summary=Summ[Rows,])
          if(Quiet == FALSE) {
               cat("Bayesian Predictive Information Criterion:\n")
               print(bpic)
               cat("Mean Lift: ", Mean.Lift, "\n")
               cat("Discrepancy Statistic: ",
                    round(Discrepancy.Statistic,5), "\n")
               cat("Monitors:\n")
               print(Mon)
               cat("\n\nRecords: \n")
               print(Summ[Rows,])}
          }
     return(invisible(Summ.out))
     }

#End
