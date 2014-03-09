###########################################################################
# CenterScale                                                             #
#                                                                         #
# The purpose of the CenterScale function is to center and scale a        #
# continuous variable. Options are also provided for binary variables.    #
# This function is very similar to Gelman's rescale function in his arm   #
# package.                                                                #
###########################################################################

CenterScale <- function(x, Binary="none", Inverse=FALSE, mu, sigma, Range,
     Min)
     {
     if(identical(Inverse, FALSE)) {
          ### Initial Checks
          if(!is.numeric(x)){
               x <- as.numeric(factor(x))
               x.obs <- x[is.finite(x)]}
          x.obs <- x[is.finite(x)]
          ### Binary Variables
          if(identical(length(unique(x.obs)), 2)){
               if(identical(Binary, "none")){
                    return((x-min(x.obs)) / (max(x.obs)-min(x.obs)))}
               else if(identical(Binary, "center")) {
                    return(x-mean(x.obs))}
               else if(identical(Binary, "center0")) {
                    x <- (x-min(x.obs)) / (max(x.obs)-min(x.obs))
                    return(x-0.5)}
               else if(identical(Binary, "centerscale")) {
                    return({x-mean(x.obs)} / {2*sd(x.obs)})}
               }
          ### Continuous Variables
          else {return({x-mean(x.obs)} / {2*sd(x.obs)})}}
     else {
          ### Initial Checks
          if(!is.numeric(x)){
               x <- as.numeric(factor(x))
               x.obs <- x[is.finite(x)]}
          x.obs <- x[is.finite(x)]
          ### Binary Variables
          if(identical(length(unique(x.obs)), 2)){
               if(identical(Binary, "none")) {
                    return(x * Range + Min)}
               else if(identical(Binary, "center")) {
                    return(x + mu)}
               else if(identical(Binary, "center0")) {
                    return(x * Range + Min)}
               else if(identical(Binary, "centerscale")) {
                    return(x * (2*sigma) + mu)}
               }
          ### Continuous Variables
          else {return(x * (2*sigma) + mu)}
          }
     }

#End
