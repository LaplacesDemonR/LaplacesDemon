###########################################################################
# Validate                                                                #
#                                                                         #
# The purpose of the Validate function is to perform hold-out validation  #
# with BPIC.                                                              #
###########################################################################

Validate <- function(object, Model, Data, plot=FALSE, PDF=FALSE)
     {
     ### Initial Checks
     if(missing(object)) stop("The object argument is required.")
     if(missing(Model)) stop("The Model argument is required.")
     if(missing(Data)) stop("The Data argument is required.")
     if(!identical(class(object), "demonoid") &
        !identical(class(object), "pmc"))
          stop("object must be of class demonoid or pmc.")
     if(!is.list(Data)) stop("Data must be a list.")
     if(length(Data) != 2) stop("Data must have length 2.")
     if(is.null(Data[[1]][["y"]]) & is.null(Data[[1]][["Y"]]))
          stop("Data must have y or Y.")
     if(!is.null(Data[[1]][["y"]])) y1 <- as.vector(Data[[1]][["y"]])
     if(!is.null(Data[[1]][["Y"]])) y1 <- as.vector(Data[[1]][["Y"]])
     if(!is.null(Data[[2]][["y"]])) y2 <- as.vector(Data[[2]][["y"]])
     if(!is.null(Data[[2]][["Y"]])) y2 <- as.vector(Data[[2]][["Y"]])
     ### p(y[rep] | y)
     if(identical(class(object), "demonoid")) {
          post <- as.matrix(object$Posterior1)
          if(is.matrix(object$Posterior2) == TRUE) {
               post <- as.matrix(object$Posterior2)}}
     else {post <- as.matrix(object$Posterior2)}
     dev1 <- dev2 <- rep(NA, nrow(post))
     yhat1 <- matrix(NA, length(y1), nrow(post))
     yhat2 <- matrix(NA, length(y2), nrow(post))
     lengthcomp <- as.vector(Model(post[1,], Data[[1]])[["yhat"]])
     if(!identical(length(lengthcomp), length(y1)))
          stop("y and yhat differ in length.")
     for (i in 1:nrow(post)) {
          temp1 <- Model(post[i,], Data[[1]])
          temp2 <- Model(post[i,], Data[[2]])
          dev1[i] <- temp1[["Dev"]]
          dev2[i] <- temp2[["Dev"]]
          yhat1[,i] <- as.vector(temp1[["yhat"]])
          yhat2[,i] <- as.vector(temp2[["yhat"]])}
     ### BPIC
     Dbar.M <- round(mean(dev1, na.rm=TRUE),3)
     pD.M <- round(var(dev1, na.rm=TRUE)/2,3)
     BPIC.M <- Dbar.M + 2*pD.M
     Dbar.V <- round(mean(dev2, na.rm=TRUE),3)
     pD.V <- round(var(dev2, na.rm=TRUE)/2,3)
     BPIC.V <- Dbar.V + 2*pD.V
     bpic <- matrix(c(Dbar.M, pD.M, BPIC.M, Dbar.V, pD.V, BPIC.V), 3, 2)
     rownames(bpic) <- c("Dbar","pD","BPIC")
     colnames(bpic) <- c("Modeled","Validation")
     cat("\n")
     print(bpic)
     shorter <- min(length(dev1),length(dev2))
     cat("\np(Deviance.V > Deviance.M):",
          round(mean(dev2[1:shorter] > dev1[1:shorter], na.rm=TRUE),3))
     cat("\nE(Change in Deviance):",
          round(mean(dev2[1:shorter] - dev1[1:shorter], na.rm=TRUE),3),
          "\n\n")
     ### Warnings
     if(is.matrix(object$Posterior2) == FALSE) {
          warning("Non-stationary samples were used.")}
     if(any(is.na(yhat1))) cat("\nWARNING: Output matrix yhat.M has ",
          sum(is.na(yhat1)), " missing values.\n")
     if(any(is.na(yhat2))) cat("\nWARNING: Output matrix yhat.V has ",
          sum(is.na(yhat2)), " missing values.\n")
     if(any(is.nan(yhat1))) cat("\nWARNING: Output matrix yhat.M has ",
          sum(is.nan(yhat1)), " non-numeric (NaN) values.\n")
     if(any(is.nan(yhat2))) cat("\nWARNING: Output matrix yhat.V has ",
          sum(is.nan(yhat2)), " non-numeric (NaN) values.\n")
     if(any(is.infinite(yhat1))) cat("\nWARNING: Output matrix yhat.M has ",
          sum(is.infinite(yhat1)), " infinite values.\n")
     if(any(is.infinite(yhat2))) cat("\nWARNING: Output matrix yhat.V has ",
          sum(is.infinite(yhat2)), " infinite values.\n")
     ### Plot
     if(plot == TRUE) {
          if(PDF == TRUE) pdf("Validation.Plot.pdf")
          par(mfrow=c(2,1))
          kde1 <- density(dev1)
          kde2 <- density(dev2)
          plot(kde1, xlim=c(min(kde1$x,kde2$x),max(kde1$x,kde2$x)),
               ylim=c(min(kde1$y,kde2$y),max(kde1$y,kde2$y)), main="",
               xlab="Deviance: Modeled (Black), Validation (Red)",
               ylab="Density")
          polygon(kde1, col="black", border="black")
          lines(kde2, col="red")
          polygon(kde2, col="red", border="red")
          lines(kde1, col="black")
          dens.diff <- dev2[1:shorter] - dev1[1:shorter]
          kde <- density(dens.diff, na.rm=TRUE)
          plot(kde, col="gray", xlab="Change in Deviance", ylab="Density",
               main="", sub=paste("(D.V-D.M)=",
               round(mean(dens.diff, na.rm=TRUE),2), " (",
               round(as.vector(quantile(dens.diff, probs=0.025,
               na.rm=TRUE)),2), ",", round(as.vector(quantile(dens.diff,
               probs=0.975, na.rm=TRUE)),2), "), p(D.V > D.M) = ",
               round(mean(dev2[1:shorter] > dev1[1:shorter],
               na.rm=TRUE),2), sep=""))
          polygon(kde, col="gray", border="gray")
          abline(v=0, col="red", lty=2)
          if(PDF == TRUE) dev.off()}
     ### Create Output
     predicted <- list(list(y=y1, yhat=yhat1, Deviance=dev1),
          list(y=y2, yhat=yhat2, Deviance=dev2),
          BPIC=bpic)
     names(predicted) <- c("Modeled","Validation","BPIC")
     if(class(object) == "demonoid") class(predicted) <- "demonoid.val"
     else class(predicted) <- "pmc.val"
     return(predicted)
     }

#End
