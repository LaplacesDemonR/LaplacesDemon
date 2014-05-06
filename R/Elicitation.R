###########################################################################
# Elicitation                                                             #
#                                                                         #
# The purpose of these functions is to facilitate prior elicitation and   #
# its use in model specification.                                         #
###########################################################################

delicit <- function(theta, x, a=-Inf, b=Inf, log=FALSE)
     {
     ### Initial Checks
     if(missing(theta)) stop("The theta argument is required.")
     if(missing(x)) stop("The x argument is required.")
     if(a >= b) stop("Lower bound a is not less than upper bound b.")
     if(is.finite(a) & is.infinite(b)) {
          if(a > 0) {theta <- log(theta); x <- log(x)}
          else if(a == 0) {
               theta <- log(theta + 1e-04); x <- log(x + 1e-04)
               }
          else { #(a < 0)
               theta <- log(theta - a + 1e-04)
               x <- log(x - a + 1e-04)
               }
          }
     if(is.infinite(a) & is.finite(b)) {
          a <- .Machine$double.xmin
          theta[which(theta <= a)] <- a + 1e-04
          theta[which(theta >= b)] <- b - 1e-04
          x[which(x <= a)] <- a + 1e-04
          x[which(x >= b)] <- b - 1e-04
          theta <- log((theta-a) / (b-theta))
          x <- log((x-a) / (b-x))
          }
     if(is.finite(a) & is.finite(b)) {
          theta[which(theta <= a)] <- a + 1e-04
          theta[which(theta >= b)] <- b - 1e-04
          x[which(x <= a)] <- a + 1e-04
          x[which(x >= b)] <- b - 1e-04
          theta <- log((theta-a) / (b-theta))
          x <- log((x-a) / (b-x))
          }
     ### Estimate Density
     kde <- density(x)
     dens <- approx(kde$x, kde$y, theta)$y
     if(log == TRUE) dens <- log(dens)
     return(dens)
     }
elicit <- function(n, cats, cat.names, show.plot=FALSE)
     {
     ### Initial Checks
     if(missing(n)) stop("The n argument is required.")
     if(missing(cats)) stop("The cats argument is required.")
     if(missing(cat.names)) stop("The cat.names argument is required.")
     if(!identical(length(cats),length(cat.names)))
          stop("Different lengths found for cats and cat.names.")
     cat.labels <- letters[1:length(cats)]
     ### Introduction
     cat("\nYou have", n, "chips.")
     cat("\nEach chip must be allocated to a category.")
     cat("\nThe categories are:")
     cat("\n\n", cat.names)
     cat("\n\nYou will be asked two questions until all chips are allocated:")
     cat("\n\n1. How many chips would you like to allocate now?")
     cat("\n2. To which category do you allocate these chips?\n")
     cat("\nCategories:", cat.names)
     cat("\nCategory Entry:", cat.labels, "\n\n")
     readline("Press Enter or Return when ready to begin: ")
     ### Elicitation
     while(n > 0) {
          cat("\n\nYou have", n ,"chips remaining.\n")
          N <- 0
          while ((N <= 0) | (N > n))
               N <- readline("How many chips would you like to allocate now? ")
          N <- as.numeric(N)
          cat("\nTo which category do you allocate these chips?\n")
          cat("\nCategories:", cat.names)
           cat("\nCategory Entry:", cat.labels, "\n\n")
          answer <- "LaplacesDemon"
          while (all(cat.labels != answer))
               answer <- readline("Category: ")
          pos <- which(cat.labels == answer)
          if(!exists("out")) out <- rep(cats[pos], N)
          else out <- c(out, rep(cats[pos], N))
          n <- n - N
          ### Barplot
          if(exists("out")) {if(show.plot == TRUE) {
               out.table <- table(out)
               count <- rep(0,length(cats))
               count[as.numeric(names(out.table))] <- as.vector(out.table)
               barplot(count, names.arg=cat.names, xlab="Category",
                    ylab="Chips", col="red")}}
          }
     cat("\n\nThank you for participating.\n")
     #Output
     return(out)
     }


#End

