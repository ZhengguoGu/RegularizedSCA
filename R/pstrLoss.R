######################################################
# a function for calculating the objective function
# associated to the pstr.R file
#
# author: Katrijn Van Deun
#####################################################

pstrLoss <- function(B, Tmat, Target, W){

  DEV <- Tmat %*% B - Target
  wDEV <- W * DEV
  Loss <- sum(wDEV^2)

  return(Loss)
}
