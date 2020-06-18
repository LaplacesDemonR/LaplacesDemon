###########################################################################
# Mode                                                                    #
#                                                                         #
# The purpose of these functions is to return the mode or modes of a      #
# vector, or test for the presence or number of modes.                    #
###########################################################################

is.amodal <- function(x, min.size=0.1)
     {
     if(any(is.na(Modes(x, min.size)[[1]]))) return(TRUE)
     else return(FALSE)
     }
is.bimodal <- function(x, min.size=0.1)
     {
     if(length(Modes(x, min.size)[[1]]) == 2) return(TRUE)
     else return(FALSE)
     }
is.multimodal <- function(x, min.size=0.1)
     {
     if(length(Modes(x, min.size)[[1]]) > 1) return(TRUE)
     else return(FALSE)
     }
is.trimodal <- function(x, min.size=0.1)
     {
     if(length(Modes(x, min.size)[[1]]) == 3) return(TRUE)
     else return(FALSE)
     }
is.unimodal <- function(x, min.size=0.1)
     {
     if((length(Modes(x, min.size)[[1]]) == 1) &
          (!any(is.na(Modes(x, min.size)[[1]])))) return(TRUE)
     else return(FALSE)
     }
Mode <- function(x)
     {
     ### Initial Checks
     if(missing(x)) stop("The x argument is required.")
     if(!is.vector(x)) x <- as.vector(x)
     x <- x[is.finite(x)]
     ### Amodal
     if(is.constant(x)) return(NA)
     ### Discrete
     if(all(x == round(x))) {
          Mode <- as.numeric(names(which.max(table(x))))}
     ### Continuous (using kernel density)
     else {
          x <- as.vector(as.numeric(as.character(x)))
          kde <- density(x)
          Mode <- kde$x[kde$y == max(kde$y)][1]
          }
     return(Mode)
     }
Modes <- function(x, min.size=0.1) {
     ### Initial Checks
     if(missing(x)) stop("The x argument is required.")
     x <- as.vector(as.numeric(as.character(x)))
     x <- x[is.finite(x)]
     ### Amodal
     if(is.constant(x))
          return(list(modes=NA, mode.dens=NA, size=1))
     ### Differentiate kernel density by x
     length(density(x)$y)
     dens.y.diff <- density(x)$y[-1] - density(x)$y[-length(density(x)$y)]
     incr <- dens.y.diff
     incr[which(dens.y.diff > 0)] <- 1
     incr[which(dens.y.diff <= 0)] <- 0
     ### Kernel density by increasing/decreasing density regions
     begin <- 1; count <- 1
     for (i in 2:length(incr)) {
          if(incr[i] != incr[i-1]) {
               count <- count + 1
               begin <- c(begin, i)}
          }
     begin <- c(begin, length(incr))
     size <- modes <- mode.dens <- rep(0, count/2)
     init <- 1
     dens <- density(x); sumdens <- sum(dens$y)
     if(incr[1] == 0) {
          size[1] <- sum(dens$y[1:begin[2]]) / sumdens
          init <- 2}
     j <- init
     for (i in init:length(size)) {
          size[i] <- sum(dens$y[begin[j]:begin[j+2]]) / sumdens
          kde <- dens
          kde$x <- kde$x[begin[j]:begin[j+2]]
          kde$y <- kde$y[begin[j]:begin[j+2]]
          modes[i] <- kde$x[kde$y == max(kde$y)][1]
          mode.dens[i] <- kde$y[kde$y == max(kde$y)][1]
          j <- j + 2
          }
     ### Order everything by density
     size <- size[order(mode.dens, decreasing=TRUE)]
     modes <- modes[order(mode.dens, decreasing=TRUE)]
     mode.dens <- mode.dens[order(mode.dens, decreasing=TRUE)]
     ### Remove modes with size < 10%
     if(any(size < 0.1)) {
         modes <- modes[-which(size < min.size)]
         mode.dens <- mode.dens[-which(size < min.size)]
         size <- size[-which(size < min.size)]
         }
     if(sum(size) > 1) size <- size / sum(size)
     #Output
     return(list(modes=modes, mode.dens=mode.dens, size=size))
     }


#End
