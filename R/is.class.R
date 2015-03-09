###########################################################################
# is.class                                                                #
#                                                                         #
# The purpose of the is.class functions is to provide logical tests       #
# regarding the classes of objects.                                       #
###########################################################################

is.bayesfactor <- function(x)
     {
     bayesfactor <- FALSE
     if(identical(class(x), "bayesfactor")) bayesfactor <- TRUE
     return(bayesfactor)
     }
is.blocks <- function(x)
     {
     if(identical(class(x), "blocks")) blocks <- TRUE
     return(blocks)
     }
is.bmk <- function(x)
     {
     bmk <- FALSE
     if(identical(class(x), "bmk")) bmk <- TRUE
     return(bmk)
     }
is.demonoid <- function(x)
     {
     demonoid <- FALSE
     if(identical(class(x), "demonoid")) demonoid <- TRUE
     return(demonoid)
     }
is.demonoid.hpc <- function(x)
     {
     demonoid.hpc <- FALSE
     if(identical(class(x), "demonoid.hpc")) demonoid.hpc <- TRUE
     return(demonoid.hpc)
     }
is.demonoid.ppc <- function(x)
     {
     demonoid.ppc <- FALSE
     if(identical(class(x), "demonoid.ppc")) demonoid.ppc <- TRUE
     return(demonoid.ppc)
     }
is.demonoid.val <- function(x)
     {
     demonoid.val <- FALSE
     if(identical(class(x), "demonoid.val")) demonoid.val <- TRUE
     return(demonoid.val)
     }
is.hangartner <- function(x)
     {
     hangartner <- FALSE
     if(identical(class(x), "hangartner")) hangartner <- TRUE
     return(hangartner)
     }
is.heidelberger <- function(x)
     {
     heidelberger <- FALSE
     if(identical(class(x), "heidelberger")) heidelberger <- TRUE
     return(heidelberger)
     }
is.importance <- function(x)
     {
     importance <- FALSE
     if(identical(class(x), "importance")) importance <- TRUE
     return(importance)
     }
is.iterquad <- function(x)
     {
     iterquad <- FALSE
     if(identical(class(x), "iterquad")) iterquad <- TRUE
     return(iterquad)
     }
is.iterquad.ppc <- function(x)
     {
     iterquad.ppc <- FALSE
     if(identical(class(x), "iterquad.ppc")) iterquad.ppc <- TRUE
     return(iterquad.ppc)
     }
is.juxtapose <- function(x)
     {
     juxtapose <- FALSE
     if(identical(class(x), "juxtapose")) juxtapose <- TRUE
     return(juxtapose)
     }
is.laplace <- function(x)
     {
     laplace <- FALSE
     if(identical(class(x), "laplace")) laplace <- TRUE
     return(laplace)
     }
is.laplace.ppc <- function(x)
     {
     laplace.ppc <- FALSE
     if(identical(class(x), "laplace.ppc")) laplace.ppc <- TRUE
     return(laplace.ppc)
     }
is.miss <- function(x)
     {
     miss <- FALSE
     if(identical(class(x), "miss")) miss <- TRUE
     return(miss)
     }
is.pmc <- function(x)
     {
     pmc <- FALSE
     if(identical(class(x), "pmc")) pmc <- TRUE
     return(pmc)
     }
is.pmc.ppc <- function(x)
     {
     pmc.ppc <- FALSE
     if(identical(class(x), "pmc.ppc")) pmc.ppc <- TRUE
     return(pmc.ppc)
     }

is.pmc.val <- function(x)
     {
     pmc.val <- FALSE
     if(identical(class(x), "pmc.val")) pmc.val <- TRUE
     return(pmc.val)
     }
is.posteriorchecks <- function(x)
     {
     posteriorchecks <- FALSE
     if(identical(class(x), "posteriorchecks")) posteriorchecks <- TRUE
     return(posteriorchecks)
     }
is.raftery <- function(x)
     {
     raftery <- FALSE
     if(identical(class(x), "raftery")) raftery <- TRUE
     return(raftery)
     }
is.rejection <- function(x)
     {
     rejection <- FALSE
     if(identical(class(x), "rejection")) rejection <- TRUE
     return(rejection)
     }
is.sensitivity <- function(x)
     {
     sensitivity <- FALSE
     if(identical(class(x), "sensitivity")) sensitivity <- TRUE
     return(sensitivity)
     }
is.vb <- function(x)
     {
     vb <- FALSE
     if(identical(class(x), "vb")) vb <- TRUE
     return(vb)
     }
is.vb.ppc <- function(x)
     {
     vb.ppc <- FALSE
     if(identical(class(x), "vb.ppc")) vb.ppc <- TRUE
     return(vb.ppc)
     }

#End
