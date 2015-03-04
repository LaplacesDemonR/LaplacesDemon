###########################################################################
# Combine                                                                 #
#                                                                         #
# The purpose of the Combine function is to combine multiple objects of   #
# class demonoid.                                                         #
###########################################################################

Combine <- function(x, Data, Thinning=1)
     {
     ### Initial Checks
     if(missing(x)) stop("x is a required argument.")
     if(missing(Data)) stop("Data is a required argument.")
     Thinning <- abs(round(Thinning))
     len.x <- length(x)
     if(!all(sapply(x, class) == "demonoid")) {
          stop("At least one item in list x is not of class demonoid.")}
     Acceptance.Rate <- round(sum(sapply(x, with, Acceptance.Rate) *
               sapply(x, with, Iterations)) /
               sum(sapply(x, with, Iterations)),5)
     Algorithm <- x[[1]]$Algorithm
     Covar <- x[[len.x]]$Covar
     Iterations <- sum(sapply(x, with, Iterations))
     Model <- x[[1]]$Model
     Minutes <- sum(sapply(x, with, Minutes))
     Status <- round(mean(sapply(x, with, Status)))
     LIV <- x[[1]]$Parameters
     ### Combine
     thinned <- x[[1]]$Posterior1
     Dev <- matrix(x[[1]]$Deviance)
     Mon <- x[[1]]$Monitor
     for (i in 2:len.x) {
          thinned <- rbind(thinned, x[[i]]$Posterior1)
          Dev <- rbind(Dev, matrix(x[[i]]$Deviance))
          Mon <- rbind(Mon, x[[i]]$Monitor)
          }
     ### Thinning
     if(Thinning > 1) {
          thinned <- Thin(thinned, By=Thinning)
          Dev <- matrix(Thin(Dev, By=Thinning))
          Mon <- Thin(Mon, By=Thinning)}
     thinned.rows <- nrow(thinned)
     ### Assess Stationarity
     cat("\nAssessing Stationarity\n")
     if(thinned.rows %% 10 == 0) thinned2 <- thinned
     if(thinned.rows %% 10 != 0) thinned2 <- thinned[1:(10*trunc(thinned.rows/10)),]
     HD <- BMK.Diagnostic(thinned2, batches=10)
     Ind <- 1 * (HD > 0.5)
     BurnIn <- thinned.rows
     batch.list <- seq(from=1, to=nrow(thinned2), by=floor(nrow(thinned2)/10))
     for (i in 1:9) {
          if(sum(Ind[,i:9]) == 0) {
               BurnIn <- batch.list[i] - 1
               break
               }
          }
     Stat.at <- BurnIn + 1
     rm(batch.list, HD, Ind, thinned2)
     ### Assess Thinning and ESS Size for all parameter samples
     cat("Assessing Thinning and ESS\n")
     acf.rows <- trunc(10*log10(thinned.rows))
     acf.temp <- matrix(1, acf.rows, LIV)
     ESS1 <- Rec.Thin <- rep(1, LIV)
     for (j in 1:LIV) {
          temp0 <- acf(thinned[,j], lag.max=acf.rows, plot=FALSE)
          if(length(temp0$acf[-1,1,1]) == acf.rows)
               acf.temp[,j] <- abs(temp0$acf[-1,1,1])
          ESS1[j] <- ESS(thinned[,j])
          Rec.Thin[j] <- which(acf.temp[,j] <= 0.1)[1]*Thinning}
     Rec.Thin[which(is.na(Rec.Thin))] <- nrow(acf.temp)
     ### Assess ESS for all deviance and monitor samples
     ESS2 <- ESS(Dev)
     ESS3 <- ESS(Mon)
     ### Assess ESS for stationary samples
     if(Stat.at < thinned.rows) {
          ESS4 <- ESS(thinned[Stat.at:thinned.rows,])
          ESS5 <- ESS(Dev[Stat.at:thinned.rows,])
          ESS6 <- ESS(Mon[Stat.at:thinned.rows,]) }
     ### Posterior Summary Table 1: All Thinned Samples
     cat("Creating Summaries\n")
     Num.Mon <- ncol(Mon)
     Summ1 <- matrix(NA, LIV, 7, dimnames=list(Data[["parm.names"]],
          c("Mean","SD","MCSE","ESS","LB","Median","UB")))
     Summ1[,1] <- colMeans(thinned)
     Summ1[,2] <- sqrt(.colVars(thinned))
     Summ1[,3] <- 0
     Summ1[,4] <- ESS1
     Summ1[,5] <- apply(thinned, 2, quantile, c(0.025), na.rm=TRUE)
     Summ1[,6] <- apply(thinned, 2, quantile, c(0.500), na.rm=TRUE)
     Summ1[,7] <- apply(thinned, 2, quantile, c(0.975), na.rm=TRUE)
     for (i in 1:ncol(thinned)) {
          temp <- try(MCSE(thinned[,i]), silent=TRUE)
          if(!inherits(temp, "try-error")) Summ1[i,3] <- temp
          else Summ1[i,3] <- MCSE(thinned[,i], method="sample.variance")}
     Deviance <- rep(NA,7)
     Deviance[1] <- mean(Dev)
     Deviance[2] <- sd(as.vector(Dev))
     temp <- try(MCSE(as.vector(Dev)), silent=TRUE)
     if(inherits(temp, "try-error"))
          temp <- MCSE(as.vector(Dev), method="sample.variance")
     Deviance[3] <- temp
     Deviance[4] <- ESS2
     Deviance[5] <- as.numeric(quantile(Dev, probs=0.025, na.rm=TRUE))
     Deviance[6] <- as.numeric(quantile(Dev, probs=0.500, na.rm=TRUE))
     Deviance[7] <- as.numeric(quantile(Dev, probs=0.975, na.rm=TRUE))
     Summ1 <- rbind(Summ1, Deviance)
     for (j in 1:Num.Mon) {
          Monitor <- rep(NA,7)
          Monitor[1] <- mean(Mon[,j])
          Monitor[2] <- sd(as.vector(Mon[,j]))
          temp <- try(MCSE(as.vector(Mon[,j])), silent=TRUE)
          if(inherits(temp, "try-error"))
               temp <- MCSE(Mon[,j], method="sample.variance")
          Monitor[3] <- temp
          Monitor[4] <- ESS3[j]
          Monitor[5] <- as.numeric(quantile(Mon[,j], probs=0.025,
               na.rm=TRUE))
          Monitor[6] <- as.numeric(quantile(Mon[,j], probs=0.500,
               na.rm=TRUE))
          Monitor[7] <- as.numeric(quantile(Mon[,j], probs=0.975,
               na.rm=TRUE))
          Summ1 <- rbind(Summ1, Monitor)
          rownames(Summ1)[nrow(Summ1)] <- Data[["mon.names"]][j]
          }
     ### Posterior Summary Table 2: Stationary Samples
     Summ2 <- matrix(NA, LIV, 7, dimnames=list(Data[["parm.names"]],
          c("Mean","SD","MCSE","ESS","LB","Median","UB")))
     if(Stat.at < thinned.rows) {
          thinned2 <- matrix(thinned[Stat.at:thinned.rows,],
               thinned.rows-Stat.at+1, ncol(thinned))
          Dev2 <- matrix(Dev[Stat.at:thinned.rows,],
               thinned.rows-Stat.at+1, ncol(Dev))
          Mon2 <- matrix(Mon[Stat.at:thinned.rows,],
               thinned.rows-Stat.at+1, ncol(Mon))
          Summ2[,1] <- colMeans(thinned2)
          Summ2[,2] <- sqrt(.colVars(thinned2))
          Summ2[,3] <- 0
          Summ2[,4] <- ESS4
          Summ2[,5] <- apply(thinned2, 2, quantile, c(0.025), na.rm=TRUE)
          Summ2[,6] <- apply(thinned2, 2, quantile, c(0.500), na.rm=TRUE)
          Summ2[,7] <- apply(thinned2, 2, quantile, c(0.975), na.rm=TRUE)
          for (i in 1:ncol(thinned2)) {
               temp <- try(MCSE(thinned2[,i]), silent=TRUE)
               if(!inherits(temp, "try-error")) Summ2[i,3] <- temp
               else Summ2[i,3] <- MCSE(thinned2[,i],
                    method="sample.variance")}
          Deviance <- rep(NA,7)
          Deviance[1] <- mean(Dev2)
          Deviance[2] <- sd(as.vector(Dev2))
          temp <- try(MCSE(as.vector(Dev2)), silent=TRUE)
          if(inherits(temp, "try-error"))
               temp <- MCSE(as.vector(Dev2), method="sample.variance")
          Deviance[3] <- temp
          Deviance[4] <- ESS5
          Deviance[5] <- as.numeric(quantile(Dev2, probs=0.025,
               na.rm=TRUE))
          Deviance[6] <- as.numeric(quantile(Dev2, probs=0.500,
               na.rm=TRUE))
          Deviance[7] <- as.numeric(quantile(Dev2, probs=0.975,
               na.rm=TRUE))
          Summ2 <- rbind(Summ2, Deviance)
          for (j in 1:Num.Mon) {
               Monitor <- rep(NA,7)
               Monitor[1] <- mean(Mon2[,j])
               Monitor[2] <- sd(as.vector(Mon2[,j]))
               temp <- try(MCSE(as.vector(Mon[,j])), silent=TRUE)
               if(inherits(temp, "try-error"))
                    temp <- MCSE(as.vector(Mon[,j]),
                    method="sample.variance")
               Monitor[3] <- temp
               Monitor[4] <- ESS6[j]
               Monitor[5] <- as.numeric(quantile(Mon2[,j],
                    probs=0.025, na.rm=TRUE))
               Monitor[6] <- as.numeric(quantile(Mon2[,j],
                    probs=0.500, na.rm=TRUE))
               Monitor[7] <- as.numeric(quantile(Mon2[,j],
                    probs=0.975, na.rm=TRUE))
               Summ2 <- rbind(Summ2, Monitor)
               rownames(Summ2)[nrow(Summ2)] <- Data[["mon.names"]][j]}
          }
     ### Column names to samples
     if(identical(ncol(Mon), length(Data[["mon.names"]])))
          colnames(Mon) <- Data[["mon.names"]]
     if(identical(ncol(thinned), length(Data[["parm.names"]]))) {
          colnames(thinned) <- Data[["parm.names"]]}
     ### Logarithm of the Marginal Likelihood
     LML <- list(LML=NA, VarCov=NA)
     if(Algorithm %in% c("Adaptive Griddy-Gibbs",
          "Affine-Invariant Ensemble Sampler",
          "Automated Factor Slice Sampler",
          "Componentwise Hit-And-Run Metropolis",
          "Delayed Rejection Metropolis",
          "Elliptical Slice Sampler",
          "Gibbs Sampler",
          "Griddy-Gibbs",
          "Hamiltonian Monte Carlo",
          "Hit-And-Run Metropolis",
          "Independence Metropolis",
          "Metropolis-Adjusted Langevin Algorithm",
          "Metropolis-Coupled Markov Chain Monte Carlo",
          "Metropolis-within-Gibbs",
          "Multiple-Try Metropolis",
          "No-U-Turn Sampler",
          "Oblique Hyperrectangle Slice Sampler",
          "Preconditioned Crank-Nicolson",
          "Random Dive Metropolis-Hastings",
          "Random-Walk Metropolis",
          "Reflective Slice Sampler",
          "Refractive Sampler",
          "Reversible-Jump",
          "Sequential Metropolis-within-Gibbs",
          "Slice Sampler",
          "Stochastic Gradient Langevin Dynamics",
          "Tempered Hamiltonian Monte Carlo",
          "t-walk",
          "Univariate Eigenvector Slice Sampler") &
          {Stat.at < thinned.rows}) {
          cat("Estimating Log of the Marginal Likelihood\n")
          LML <- LML(theta=thinned2, LL=as.vector(Dev2)*(-1/2),
               method="NSIS")}
     ### Compile Output
     cat("Creating Output\n")
     LaplacesDemon.out <- list(Acceptance.Rate=Acceptance.Rate,
          Algorithm=Algorithm,
          Call=x[[1]]$Call,
          Covar=Covar,
          CovarDHis=x[[len.x]]$CovarDHis,
          Deviance=as.vector(Dev),
          DIC1=c(mean(as.vector(Dev)),
               var(as.vector(Dev))/2,
               mean(as.vector(Dev)) + var(as.vector(Dev))/2),
          DIC2=if(Stat.at < thinned.rows) {
               c(mean(as.vector(Dev2)),
               var(as.vector(Dev2))/2,
               mean(as.vector(Dev2)) + 
               var(as.vector(Dev2))/2)}
               else rep(NA,3),
          Initial.Values=x[[len.x]]$Initial.Values,
          Iterations=Iterations,
          LML=LML[[1]],
          Minutes=Minutes,
          Model=Model,
          Monitor=Mon,
          Parameters=LIV,
          Posterior1=thinned,
          Posterior2=if(Stat.at < thinned.rows) {
               thinned[Stat.at:thinned.rows,]}
               else thinned[thinned.rows,],
          Rec.BurnIn.Thinned=BurnIn,
          Rec.BurnIn.UnThinned=BurnIn*Thinning,
          Rec.Thinning=min(1000, max(Rec.Thin)),
          Specs=x[[1]]$Specs,
          Status=Status,
          Summary1=Summ1,
          Summary2=Summ2,
          Thinned.Samples=thinned.rows,
          Thinning=Thinning)
     class(LaplacesDemon.out) <- "demonoid"
     cat("\nLaplace's Demon has finished.\n")
     return(LaplacesDemon.out)
     }

#End
