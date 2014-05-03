###########################################################################
# de.Finetti.Game                                                         #
#                                                                         #
# The purpose of the de.Finetti.Game function is to elicit the interval   #
# of a subjective probability about a possible event in the near future.  #
###########################################################################

de.Finetti.Game <- function(width)
     {
     if(missing(width)) stop("The width argument is required.")
     if((width <= 0) | (width > 1))
          stop("The width argument is incoherent.")
     ques <- paste("\nDescribe a possible event in the near",
          "future (such as ``rain tomorrow''): ")
     event <- readline(ques)
     region <- c(0,1)
     while((region[2] - region[1]) > width) {
          x <- round(mean(region) * 100)
          y <- 100 - x
          cat("\nYou have two options:")
          cat("\n\n1.Wait and receive $1 if the event happens.")
          cat("\n2.Draw a marble from an urn with", x,
               "black marbles and", y, "white marbles. Drawing a black",
               "marble results in receiving $1.")
          ans <- readline("\n\nChoose 1 or 2: ")
          if(ans == 1) region[1] <- x / 100
          else if(ans == 2) region[2] <- x / 100
          else region <- c(0,0)
          }
     if(sum(region) == 0) cat("\nTry again. Valid answers are 1 or 2.")
     else {cat("\n\nYour subjective probability is in the interval [",
          region[1], ",", region[2], "] regarding ", event, ".\n\n",
          sep="")}
     return(region)
     }

#End
