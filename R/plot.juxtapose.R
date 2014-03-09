###########################################################################
# plot.juxtapose                                                          #
#                                                                         #
# The purpose of the plot.juxtapose function is to plot a comparison of   #
# MCMC algorithms according either to IAT or ISM.                         #
###########################################################################

plot.juxtapose <- function(x, Style="ISM", ...)
     {
     ### Initial Checks
     if(missing(x)) stop("The x argument is required.")
     if(!identical(class(x), "juxtapose"))
          stop("The x argument is not of class juxtapose.")
     if((Style != "IAT") & (Style != "ISM"))
          stop("Style must be IAT or ISM")
     Title <- "MCMC Juxtaposition"
     if(identical(Style, "IAT")) {
          ### Basic Plot
          plot(0, 0, ylim=c(0, ncol(x) + 1), xlim=c(1, max(x[6,])),
               main=Title, sub="", xlab="Integrated Autocorrelation Time",
               ylab="", type="n", ann=TRUE, yaxt="n")
          abline(v=1, col="gray")
          ### Add Medians
          points(x[5,], ncol(x):1, pch=20)
          ### Add Horizontal Lines for 2.5%-97.5% Quantiles
          for (i in 1:ncol(x)) {
               lines(x[c(4,6),i], c(ncol(x)-i+1, ncol(x)-i+1))}
          ### Add y-axis labels
          yy <- ncol(x):1
          cex.labels <- 1 / {log(ncol(x))/5 + 1}
          axis(2, labels=colnames(x), tick=FALSE, las=1, at=yy,
               cex.axis=cex.labels)}
     else {
          ### Basic Plot
          plot(0, 0, ylim=c(0, ncol(x) + 1), xlim=c(0, max(x[9,])),
               main=Title, sub="", xlab="Independent Samples per Minute",
               ylab="", type="n", ann=TRUE, yaxt="n")
          abline(v=0, col="gray")
          ### Add Medians
          points(x[8,], ncol(x):1, pch=20)
          ### Add Horizontal Lines for 2.5%-97.5% Quantiles
          for (i in 1:ncol(x)) {
               lines(x[c(7,9),i], c(ncol(x)-i+1, ncol(x)-i+1))}
          ### Add y-axis labels
          yy <- ncol(x):1
          cex.labels <- 1 / {log(ncol(x))/5 + 1}
          axis(2, labels=colnames(x), tick=FALSE, las=1, at=yy,
               cex.axis=cex.labels)}
     }

#End
