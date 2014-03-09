###########################################################################
# server_Listening                                                        #
#                                                                         #
# The server_Listening function is not intended to be called directly by  #
# the user. It is an internal-only function that is intended to prevent   #
# cluster problems while using the INCA algorithm through the             #
# LaplacesDemon.hpc function.                                             #
###########################################################################

server_Listening <- function(n=2, port=19009)
     {
     slist <- vector('list', n)
     for (i in 1:n) {
          slist[[i]] <- socketConnection("localhost", port, server=TRUE,
               open="r+")
          cat("\nClient", i, "Connected")}
     tmp <- NULL
     trow <- 0
     stop_server <- FALSE
     cat("\nStart listening...")
     repeat
          {
          ready <- which(socketSelect(slist, TRUE))
          for (i in ready) {
               #print(paste("Socket", i, "ready to write"))
               con <- slist[[i]]
               #print("Write message...")
               if(is.null(tmp)) serialize(tmp, con)
               else serialize(tmp[-(((i-1)*trow+1):(i*trow))], con)
               #print("Read message...")
               buf <- try(unserialize(con), silent=TRUE)
               if(is.matrix(buf)) {
                    if(is.null(tmp)) {
                         tmp <- matrix(0, nrow=n*nrow(buf), ncol=ncol(buf))
                         trow <- nrow(buf)
                         }
                    tmp[((i-1)*trow+1):(i*trow),] <- buf
                    }
               else {
                    stop_server <- TRUE
                    break
                    }
               }
          if(stop_server == TRUE) break
          }
     for (i in 1:n) {
          close(slist[[i]])
          cat("\nClose connection", i)
          }
     cat("\n")
     }

#End
