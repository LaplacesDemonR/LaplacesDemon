###########################################################################
# Levene.Test                                                             #
#                                                                         #
# The purpose of the Levene.Test function is to apply Levene's test to    #
# residuals as a test of homoegeneity.                                    #
###########################################################################

Levene.Test <- function(x, Method="U", G=NULL, Data=NULL)
     {
     ### Initial Checks
     if(missing(x)) stop("The x argument is required.")
     if({class(x) != "demonoid.ppc"} & {class(x) != "iterquad.ppc"} &
        {class(x) != "laplace.ppc"} & {class(x) != "pmc.ppc"} &
        {class(x) != "vb.ppc"})
          stop("x is not of class demonoid.ppc, iterquad.ppc, laplace.ppc, pmc.ppc, or vb.ppc.")
     if({Method == "C"} & is.null(Data))
          stop("Data is required for Method C.")
     if({Method == "R"} & is.null(Data))
          stop("Data is required for Method R.")
     if({Method == "C" | Method == "R"} & is.null(Data[["Y"]]))
          stop("Y is required in Data.")
     y <- x[["y"]]
     yhat <- x[["yhat"]]
     if(is.null(Data) & is.null(G)) {
          G <- rep(1:4, each=round(nrow(yhat)/4), len=nrow(yhat))
          if(length(G) != length(y))
               stop("Lengths of G and y differ.")}
     if(!is.null(Data) & is.null(G) & {Method == "C" || Method == "R"}) {
          if(Method == "C") {
               G <- matrix(rep(1:4, each=round(nrow(Data[["Y"]])/4),
                    len=nrow(Data[["Y"]])), nrow(Data[["Y"]]),
                    ncol(Data[["Y"]]))}
          if(Method == "R") {
               G <- matrix(rep(1:4, each=round(ncol(Data[["Y"]])/4),
                    len=ncol(Data[["Y"]])), nrow(Data[["Y"]]),
                    ncol(Data[["Y"]]), byrow=TRUE)}}
     if(Method == "U") K <- length(unique(G))
     if(Method == "C") K <- length(unique(G[,1]))
     if(Method == "R") K <- length(unique(G[1,]))
     N <- nrow(yhat)
     S <- ncol(yhat)
     epsilon.obs <- y - yhat
     if(Method == "U") {
          epsilon.rep <- matrix(rnorm(N*S, mean(epsilon.obs, na.rm=TRUE),
               sd(as.vector(epsilon.obs), na.rm=TRUE)), N, S)}
     if(Method == "C") {
          epsilon.rep <- epsilon.obs
          for (j in 1:ncol(Data[["Y"]])) {
               point <- matrix(FALSE, nrow(Data[["Y"]]), ncol(Data[["Y"]]))
               point[,j] <- TRUE
               point <- as.vector(point)
               epsilon.rep[point,] <- rnorm(nrow(Data[["Y"]])*S,
                    mean(epsilon.obs[point,], na.rm=TRUE),
                    sd(as.vector(epsilon.obs[point,]), na.rm=TRUE))}}
     if(Method == "R") {
          epsilon.rep <- epsilon.obs
          for (i in 1:nrow(Data[["Y"]])) {
               point <- matrix(FALSE, nrow(Data[["Y"]]), ncol(Data[["Y"]]))
               point[i,] <- TRUE
               point <- as.vector(point)
               epsilon.rep[point,] <- rnorm(ncol(Data[["Y"]])*S,
                    mean(epsilon.obs[point,], na.rm=TRUE),
                    sd(as.vector(epsilon.obs[point,]), na.rm=TRUE))}}
     ### Levene Test
     if(Method == "U") {
          W.obs <- W.rep <- rep(0, S)
          for (s in 1:S) {
               epsilon.G.obs <- as.vector(by(epsilon.obs[,s], G, mean,
                    na.rm=TRUE))
               epsilon.G.rep <- as.vector(by(epsilon.rep[,s], G, mean,
                    na.rm=TRUE))
               Z.obs <- abs(epsilon.obs[,s] - epsilon.G.obs[G])
               Z.rep <- abs(epsilon.rep[,s] - epsilon.G.rep[G])
               Zbar.obs <- mean(Z.obs, na.rm=TRUE)
               Zbar.rep <- mean(Z.rep, na.rm=TRUE)
               Zbar.G.obs <- as.vector(by(Z.obs, G, mean, na.rm=TRUE))
               Zbar.G.rep <- as.vector(by(Z.rep, G, mean, na.rm=TRUE))
               W.obs[s] <- ((N-K)*sum((Zbar.G.obs[G] - Zbar.obs)^2,
                    na.rm=TRUE)) / ((K-1)*sum((Z.obs - Zbar.G.obs[G])^2,
                    na.rm=TRUE))
               W.rep[s] <- ((N-K)*sum((Zbar.G.rep[G] - Zbar.rep)^2,
                    na.rm=TRUE)) / ((K-1)*sum((Z.rep - Zbar.G.rep[G])^2,
                    na.rm=TRUE))}
          p <- round(mean(W.obs > W.rep, na.rm=TRUE),3)
          result <- "Homoskedastic"
          if((p < 0.025) | (p > 0.975)) result <- "Heteroskedastic"
          par(mfrow=c(1,1))
          d.W.obs <- density(W.obs)
          d.W.rep <- density(W.rep)
          plot(d.W.obs, xlim=c(min(d.W.obs$x, d.W.rep$x),
                    max(d.W.obs$x, d.W.rep$x)),
               ylim=c(0, max(d.W.obs$y, d.W.rep$y)),
               col=rgb(0,0,0,50,maxColorValue=255),
               main="Levene's Test",
               xlab="W",
               sub=paste("W.obs=", round(mean(W.obs, na.rm=TRUE),2),
                    " (", round(as.vector(quantile(W.obs, probs=0.025,
                    na.rm=TRUE)),2), ",", round(as.vector(quantile(W.obs,
                    probs=0.975, na.rm=TRUE)),2), "), p(W.obs > W.rep) = ",
                    p, " = ", result, sep=""),
               ylab="Density")
          polygon(d.W.obs, col=rgb(0,0,0,50,maxColorValue=255), border=NA)
          polygon(d.W.rep, col=rgb(255,0,0,50,maxColorValue=255),
               border=NA)}
     if(Method == "C") {
          par(mfrow=c(1,1), ask=TRUE)
          p <- rep(0, ncol(Data[["Y"]]))
          for (j in 1:ncol(Data[["Y"]])) {
               point <- matrix(FALSE, nrow(Data[["Y"]]), ncol(Data[["Y"]]))
               point[,j] <- TRUE
               point <- as.vector(point)
               W.obs <- W.rep <- rep(0, S)
               for (s in 1:S) {
                    epsilon.G.obs <- as.vector(by(epsilon.obs[point,s],
                         G[,j], mean, na.rm=TRUE))
                    epsilon.G.rep <- as.vector(by(epsilon.rep[point,s],
                         G[,j], mean, na.rm=TRUE))
                    Z.obs <- abs(epsilon.obs[point,s] -
                         epsilon.G.obs[G[,j]])
                    Z.rep <- abs(epsilon.rep[point,s] -
                         epsilon.G.rep[G[,j]])
                    Zbar.obs <- mean(Z.obs, na.rm=TRUE)
                    Zbar.rep <- mean(Z.rep, na.rm=TRUE)
                    Zbar.G.obs <- as.vector(by(Z.obs, G[,j], mean,
                         na.rm=TRUE))
                    Zbar.G.rep <- as.vector(by(Z.rep, G[,j], mean,
                         na.rm=TRUE))
                    W.obs[s] <- ((N-K)*sum((Zbar.G.obs[G[,j]] -
                         Zbar.obs)^2, na.rm=TRUE)) / ((K-1)*sum((Z.obs -
                         Zbar.G.obs[G[,j]])^2, na.rm=TRUE))
                    W.rep[s] <- ((N-K)*sum((Zbar.G.rep[G[,j]] -
                         Zbar.rep)^2, na.rm=TRUE)) / ((K-1)*sum((Z.rep -
                         Zbar.G.rep[G[,j]])^2, na.rm=TRUE))}
               p[j] <- round(mean(W.obs > W.rep, na.rm=TRUE),3)
               result <- "Homoskedastic"
               if((p[j] < 0.025) | (p[j] > 0.975))
                    result <- "Heteroskedastic"
               d.W.obs <- density(W.obs)
               d.W.rep <- density(W.rep)
               plot(d.W.obs, xlim=c(min(d.W.obs$x, d.W.rep$x),
                         max(d.W.obs$x, d.W.rep$x)),
                    ylim=c(0, max(d.W.obs$y, d.W.rep$y)),
                    col=rgb(col2rgb(j)[1], col2rgb(j)[2], col2rgb(j)[3],
                         50, maxColorValue=255),
                    main="Levene's Test",
                    xlab=paste("W.obs=", round(mean(W.obs, na.rm=TRUE),2),
                         " (", round(as.vector(quantile(W.obs, probs=0.025,
                         na.rm=TRUE)),2), ",",
                         round(as.vector(quantile(W.obs, probs=0.975,
                         na.rm=TRUE)),2), "), p(W.obs > W.rep) = ",
                         p[j], " = ", result, sep=""),
                    sub=paste("Y[,", j, "]", sep=""),
                    ylab="Density")
               polygon(d.W.obs, col=rgb(0,0,0,50,maxColorValue=255),
                    border=NA)
               polygon(d.W.rep, col=rgb(255,0,0,50,maxColorValue=255),
                    border=NA)}}
     if(Method == "R") {
          par(mfrow=c(1,1), ask=TRUE)
          p <- rep(0, nrow(Data[["Y"]]))
          for (i in 1:nrow(Data[["Y"]])) {
               point <- matrix(FALSE, nrow(Data[["Y"]]), ncol(Data[["Y"]]))
               point[i,] <- TRUE
               point <- as.vector(point)
               W.obs <- W.rep <- rep(0, S)
               for (s in 1:S) {
                    epsilon.G.obs <- as.vector(by(epsilon.obs[point,s],
                         G[i,], mean, na.rm=TRUE))
                    epsilon.G.rep <- as.vector(by(epsilon.rep[point,s],
                         G[i,], mean, na.rm=TRUE))
                    Z.obs <- abs(epsilon.obs[point,s] -
                         epsilon.G.obs[G[i,]])
                    Z.rep <- abs(epsilon.rep[point,s] -
                         epsilon.G.rep[G[i,]])
                    Zbar.obs <- mean(Z.obs, na.rm=TRUE)
                    Zbar.rep <- mean(Z.rep, na.rm=TRUE)
                    Zbar.G.obs <- as.vector(by(Z.obs, G[i,], mean,
                         na.rm=TRUE))
                    Zbar.G.rep <- as.vector(by(Z.rep, G[i,], mean,
                         na.rm=TRUE))
                    W.obs[s] <- ((N-K)*sum((Zbar.G.obs[G[i,]] -
                         Zbar.obs)^2, na.rm=TRUE)) /
                         ((K-1)*sum((Z.obs - Zbar.G.obs[G[i,]])^2,
                         na.rm=TRUE))
                    W.rep[s] <- ((N-K)*sum((Zbar.G.rep[G[i,]] -
                         Zbar.rep)^2, na.rm=TRUE)) /
                         ((K-1)*sum((Z.rep -
                         Zbar.G.rep[G[i,]])^2, na.rm=TRUE))}
               p[i] <- round(mean(W.obs > W.rep, na.rm=TRUE),3)
               result <- "Homoskedastic"
               if((p[i] < 0.025) | (p[i] > 0.975))
                    result <- "Heteroskedastic"
               d.W.obs <- density(W.obs)
               d.W.rep <- density(W.rep)
               plot(d.W.obs, xlim=c(min(d.W.obs$x, d.W.rep$x),
                         max(d.W.obs$x, d.W.rep$x)),
                    ylim=c(0, max(d.W.obs$y, d.W.rep$y)),
                    col=rgb(col2rgb(j)[1], col2rgb(j)[2], col2rgb(j)[3],
                         50, maxColorValue=255),
                    main="Levene's Test",
                    xlab=paste("W.obs=", round(mean(W.obs, na.rm=TRUE),2),
                         " (", round(as.vector(quantile(W.obs, probs=0.025,
                         na.rm=TRUE)),2), ",",
                         round(as.vector(quantile(W.obs, probs=0.975,
                         na.rm=TRUE)),2), "), p(W.obs > W.rep) = ",
                         p[i], " = ", result, sep=""),
                    sub=paste("Y[", i, ",]", sep=""),
                    ylab="Density")
               polygon(d.W.obs, col=rgb(0,0,0,50,maxColorValue=255),
                    border=NA)
               polygon(d.W.rep, col=rgb(255,0,0,50,maxColorValue=255),
                    border=NA)}}
     return(p)
     }

#End
