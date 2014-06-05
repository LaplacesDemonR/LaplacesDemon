###########################################################################
# Importance                                                              #
#                                                                         #
# The purpose of the Importance function is to compare the impact of      #
# design matrix X on replicates when each column vector (predictor) is    #
# sequentially removed.                                                   #
###########################################################################

Importance <- function(object, Model, Data, Categorical=FALSE, Discrep,
     d=0, CPUs=1, Type="PSOCK")
     {
     if(missing(object)) stop("The object argument is required.")
     if(missing(Model)) stop("The Model arguement is required.")
     if(missing(Data)) stop("The Data argument is required.")
     if(is.null(Data[["X"]])) stop("Data must have X.")
     if(missing(Discrep)) Discrep <- NULL
     X.orig <- Data[["X"]]
     cat("\nX has", ncol(X.orig), "variables")
     cat("\nEstimating the full model...")
     Pred <- predict(object, Model, Data)
     Summ <- summary(Pred, Categorical=Categorical, Discrep=Discrep, d=d,
          Quiet=TRUE)
     out <- matrix(0, ncol(X.orig) + 1, 4)
     out[1,1] <- Summ$BPIC[1,3]
     if(Categorical == FALSE) out[1,2] <- round(Summ$Concordance, 3)
     else out[1,2] <- round(Summ$Mean.Lift, 3)
     out[1,3] <- Summ$Discrepancy.Statistic
     if(Categorical == FALSE) {
          out[1,4] <- Summ$L.criterion
          S.L <- Summ$S.L}
     else S.L <- NA
     for (i in 1:ncol(X.orig)) {
          cat("\nEstimating without X[,", i, "]...", sep="")
          X.temp <- X.orig
          X.temp[,i] <- 0
          Data[["X"]] <- X.temp
          Pred <- predict(object, Model, Data, CPUs, Type)
          Summ <- summary(Pred, Categorical=Categorical,
               Discrep=Discrep, d=d, Quiet=TRUE)
          out[i+1,1] <- Summ$BPIC[1,3]
          if(Categorical == FALSE) out[i+1,2] <- round(Summ$Concordance, 3)
          else out[i+1,2] <- round(Summ$Mean.Lift, 3)
          out[i+1,3] <- Summ$Discrepancy.Statistic
          if(Categorical == FALSE) {
               out[i+1,4] <- Summ$L.criterion
               S.L <- c(S.L, Summ$S.L)}}
     if(Categorical == FALSE) cat("\n\nS.L:", S.L)
     colnames(out) <- c("BPIC","Concordance", "Discrep", "L-criterion")
     rownames(out) <- c("Full", paste("X[,-", 1:ncol(X.orig), "]", sep=""))
     attr(out, "S.L") <- S.L
     class(out) <- "importance"
     cat("\n\n")
     return(out)
     }

#End
