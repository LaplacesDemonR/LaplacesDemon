###########################################################################
# deburn                                                                  #
#                                                                         #
# The purpose of deburn() is to remove the user-specified burn-in from an #
# object of class demonoid.                                               #
###########################################################################

deburn <- function(x, BurnIn=0)
     {
     ### Initial Checks
     if(!identical(class(x), "demonoid"))
          stop("x is not an object of class demonoid.")
     S <- nrow(x$Posterior1)
     if(S < 22) stop("x has too few posterior samples.")
     BurnIn <- abs(round(BurnIn))
     if(BurnIn >= S) BurnIn <- S - 2
     LIV <- x$Parameters
     ### Remove Burn-in
     x$Posterior1 <- x$Posterior2 <- x$Posterior1[(BurnIn+1):S,]
     x$Deviance <- x$Deviance[(BurnIn+1):S]
     x$Monitor <- x$Monitor[(BurnIn+1):S,,drop=FALSE]
     x$Rec.BurnIn.Thinned <- 0
     x$Rec.BurnIn.UnThinned <- 0
     x$Thinned.Samples <- x$Thinned.Samples - BurnIn
     ### Summary1
     x$Summary1[1:LIV,1] <- colMeans(x$Posterior1)
     x$Summary1[1:LIV,2] <- sqrt(.colVars(x$Posterior1))
     x$Summary1[1:LIV,4] <- ESS(x$Posterior1)
     x$Summary1[1:LIV,5] <- apply(x$Posterior1, 2, quantile, c(0.025),
          na.rm=TRUE)
     x$Summary1[1:LIV,6] <- apply(x$Posterior1, 2, quantile, c(0.500),
          na.rm=TRUE)
     x$Summary1[1:LIV,7] <- apply(x$Posterior1, 2, quantile, c(0.975),
          na.rm=TRUE)
     for (i in 1:LIV) {
          temp <- try(MCSE(x$Posterior1[,i]), silent=TRUE)
          if(!inherits(temp, "try-error")) x$Summary1[i,3] <- temp
          else x$Summary1[i,3] <- MCSE(x$Posterior1[,i],
               method="sample.variance")}
     ### Deviance
     x$Summary1[LIV+1,1] <- mean(x$Deviance)
     x$Summary1[LIV+1,2] <- sd(x$Deviance)
     temp <- try(MCSE(x$Deviance), silent=TRUE)
     if(inherits(temp, "try-error"))
          temp <- MCSE(x$Deviance, method="sample.variance")
     x$Summary1[LIV+1,3] <- temp
     x$Summary1[LIV+1,4] <- ESS(x$Deviance)
     x$Summary1[LIV+1,5] <- as.numeric(quantile(x$Deviance, probs=0.025,
          na.rm=TRUE))
     x$Summary1[LIV+1,6] <- as.numeric(quantile(x$Deviance, probs=0.500,
          na.rm=TRUE))
     x$Summary1[LIV+1,7] <- as.numeric(quantile(x$Deviance, probs=0.975,
          na.rm=TRUE))
     ### Monitor
     Num.Mon <- ncol(x$Monitor)
     x$Summary1[LIV+1+1:Num.Mon,1] <- colMeans(x$Monitor)
     x$Summary1[LIV+1+1:Num.Mon,2] <- sqrt(.colVars(x$Monitor))
     x$Summary1[LIV+1+1:Num.Mon,4] <- ESS(x$Monitor)
     x$Summary1[LIV+1+1:Num.Mon,5] <- apply(x$Monitor, 2, quantile,
          c(0.025), na.rm=TRUE)
     x$Summary1[LIV+1+1:Num.Mon,6] <- apply(x$Monitor, 2, quantile,
          c(0.500), na.rm=TRUE)
     x$Summary1[LIV+1+1:Num.Mon,7] <- apply(x$Monitor, 2, quantile,
          c(0.975), na.rm=TRUE)
     for (i in 1:Num.Mon) {
          temp <- try(MCSE(x$Monitor[,i]), silent=TRUE)
          if(!inherits(temp, "try-error")) x$Summary1[LIV+1+i,3] <- temp
          else x$Summary1[LIV+1+i,3] <- MCSE(x$Monitor[,i],
               method="sample.variance")}
     ### Summary2
     x$Summary2 <- x$Summary1
     ### DIC
     x$DIC1 <- x$DIC2 <- c(mean(x$Deviance), var(x$Deviance)/2,
          mean(x$Deviance) + var(x$Deviance)/2)
     ### Output
     return(x)
     }

#End
