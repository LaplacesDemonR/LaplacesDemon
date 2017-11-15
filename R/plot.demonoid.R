###########################################################################
# plot.demonoid                                                           #
#                                                                         #
# The purpose of the plot.demonoid function is to plot an object of class #
# demonoid.                                                               #
###########################################################################

plot.demonoid <- function(x, BurnIn=0, Data=NULL, PDF=FALSE, Parms=NULL,
                          FileName = paste0("laplacesDemon-plot_", format(Sys.time(), "%Y-%m-%d_%T"), ".pdf"), ...)
{
  ### Initial Checks
  if(missing(x)) stop("The x argument is required.")
  if(is.null(Data)) stop("The Data argument is NULL.")
  if(BurnIn >= nrow(x$Posterior1)) BurnIn <- 0
  Stat.at <- BurnIn + 1
  ### Selecting Parms
  if(is.null(Parms)) {Posterior <- x$Posterior1}
  else {
    Parms <- sub("\\[","\\\\[",Parms)
    Parms <- sub("\\]","\\\\]",Parms)
    Parms <- sub("\\.","\\\\.",Parms)
    if(length(grep(Parms[1], colnames(x$Posterior1))) == 0)
      stop("Parameter in Parms does not exist.")
    keepcols <- grep(Parms[1], colnames(x$Posterior1))
    if(length(Parms) > 1) {
      for (i in 2:length(Parms)) {
        if(length(grep(Parms[i],
                       colnames(x$Posterior1))) == 0)
          stop("Parameter in Parms does not exist.")
        keepcols <- c(keepcols,
                      grep(Parms[i], colnames(x$Posterior1)))}}
    Posterior <- as.matrix(x$Posterior1[,keepcols])
    colnames(Posterior) <- colnames(x$Posterior1)[keepcols]}
  if(PDF == TRUE) {
    pdf(FileName)
    par(mfrow=c(3,3))
  }
  else {par(mfrow=c(3,3), ask=TRUE)}
  ### Plot Parameters
  for (j in 1:ncol(Posterior))
  {
    plot(Stat.at:x$Thinned.Samples,
         Posterior[Stat.at:x$Thinned.Samples,j],
         type="l", xlab="Thinned Samples", ylab="Value",
         main=colnames(Posterior)[j])
    panel.smooth(Stat.at:x$Thinned.Samples,
                 Posterior[Stat.at:x$Thinned.Samples,j], pch="")
    plot(density(Posterior[Stat.at:x$Thinned.Samples,j]),
         xlab="Value", main=colnames(Posterior)[j])
    polygon(density(Posterior[Stat.at:x$Thinned.Samples,j]),
            col="black", border="black")
    abline(v=0, col="red", lty=2)
    ### Only plot an ACF if there's > 1 unique values
    if(!is.constant(Posterior[Stat.at:x$Thinned.Samples,j])) {
      z <- acf(Posterior[Stat.at:x$Thinned.Samples,j], plot=FALSE)
      se <- 1/sqrt(length(Posterior[Stat.at:x$Thinned.Samples,j]))
      plot(z$lag, z$acf, ylim=c(min(z$acf,-2*se),1), type="h",
           main=colnames(Posterior)[j], xlab="Lag",
           ylab="Correlation")
      abline(h=(2*se), col="red", lty=2)
      abline(h=(-2*se), col="red", lty=2)
    }
    else {plot(0,0, main=paste(colnames(Posterior)[j],
                               "is a constant."))}
  }
  ### Plot Deviance
  plot(Stat.at:length(x$Deviance),
       x$Deviance[Stat.at:length(x$Deviance)],
       type="l", xlab="Thinned Samples", ylab="Value", main="Deviance")
  panel.smooth(Stat.at:length(x$Deviance),
               x$Deviance[Stat.at:length(x$Deviance)], pch="")
  plot(density(x$Deviance[Stat.at:length(x$Deviance)]),
       xlab="Value", main="Deviance")
  polygon(density(x$Deviance[Stat.at:length(x$Deviance)]), col="black",
          border="black")
  abline(v=0, col="red", lty=2)
  ### Only plot an ACF if there's > 1 unique values
  if(!is.constant(x$Deviance[Stat.at:length(x$Deviance)])) {
    z <- acf(x$Deviance[Stat.at:length(x$Deviance)], plot=FALSE)
    se <- 1/sqrt(length(x$Deviance[Stat.at:length(x$Deviance)]))
    plot(z$lag, z$acf, ylim=c(min(z$acf,-2*se),1), type="h",
         main="Deviance", xlab="Lag", ylab="Correlation")
    abline(h=(2*se), col="red", lty=2)
    abline(h=(-2*se), col="red", lty=2)
  }
  else {plot(0,0, main="Deviance is a constant.")}
  ### Plot Monitored Variables
  if(is.vector(x$Monitor)) {J <- 1; nn <- length(x$Monitor)}
  else if(is.matrix(x$Monitor)) {
    J <- ncol(x$Monitor); nn <- nrow(x$Monitor)}
  for (j in 1:J)
  {
    plot(Stat.at:nn, x$Monitor[Stat.at:nn,j],
         type="l", xlab="Thinned Samples", ylab="Value",
         main=Data[["mon.names"]][j])
    panel.smooth(Stat.at:nn, x$Monitor[Stat.at:nn,j], pch="")
    plot(density(x$Monitor[Stat.at:nn,j]),
         xlab="Value", main=Data[["mon.names"]][j])
    polygon(density(x$Monitor[Stat.at:nn,j]), col="black",
            border="black")
    abline(v=0, col="red", lty=2)
    ### Only plot an ACF if there's > 1 unique values
    if(!is.constant(x$Monitor[Stat.at:nn,j])) {
      z <- acf(x$Monitor[Stat.at:nn,j], plot=FALSE)
      se <- 1/sqrt(length(x$Monitor[Stat.at:nn,j]))
      plot(z$lag, z$acf, ylim=c(min(z$acf,-2*se),1), type="h",
           main=Data[["mon.names"]][j], xlab="Lag",
           ylab="Correlation")
      abline(h=(2*se), col="red", lty=2)
      abline(h=(-2*se), col="red", lty=2)
    }
    else {plot(0,0, main=paste(Data[["mon.names"]][j],
                               "is a constant."))}
  }
  ### Diminishing Adaptation
  if(nrow(x$CovarDHis) > 1) {
    if(x$Algorithm %in% c("Adaptive Hamiltonian Monte Carlo",
                          "Hamiltonian Monte Carlo with Dual-Averaging",
                          "No-U-Turn Sampler")) {
      plot(x$CovarDHis[,1], type="l", xlab="Adaptations",
           main="Step-Size", ylab=expression(epsilon))}
    else {
      if(x$Algorithm %in% c("Oblique Hyperrectangle Slice Sampler",
                            "Univariate Eigenvector Slice Sampler"))
        title <- "Eigenvectors"
      else if(x$Algorithm %in% c("Metropolis-Adjusted Langevin Algorithm"))
        title <- "Lambda"
      else if(x$Algorithm %in% c("Componentwise Hit-And-Run Metropolis",
                                 "Hit-And-Run Metropolis"))
        title <- "Proposal Distance"
      else if(x$Algorithm %in% c("Adaptive Griddy-Gibbs",
                                 "Adaptive Metropolis-within-Gibbs",
                                 "Sequential Adaptive Metropolis-within-Gibbs",
                                 "Updating Sequential Adaptive Metropolis-within-Gibbs"))
        title <- "Proposal S.D."
      else if(x$Algorithm %in% c("Differential Evolution Markov Chain"))
        title <- "Z"
      else if(x$Algorithm %in% c("Adaptive Factor Slice Sampler",
                                 "Refractive Sampler"))
        title <- "Step-Size"
      else title <- "Proposal Variance"
      Diff <- abs(diff(x$CovarDHis))
      adaptchange <- matrix(NA, nrow(Diff), 3)
      for (i in 1:nrow(Diff)) {
        adaptchange[i,1:3] <- as.vector(quantile(Diff[i,],
                                                 probs=c(0.025, 0.500, 0.975)))}
      plot(1:nrow(Diff), adaptchange[,2],
           ylim=c(min(adaptchange), max(adaptchange)),
           type="l", col="red", xlab="Adaptations",
           ylab="Absolute Difference", main=title,
           sub="Median=Red, Interval=Transparent Red")
      polygon(c(1:nrow(Diff),rev(1:nrow(Diff))),
              c(adaptchange[,1], rev(adaptchange[,3])),
              col=rgb(255, 0, 0, 50, maxColorValue=255),
              border=FALSE)
      lines(adaptchange[,2], col="red")}
  }
  if(PDF == TRUE) dev.off()
}

# End
