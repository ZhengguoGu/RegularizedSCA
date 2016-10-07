# This is to set up a simulation to test whether the algorithms works
# The same simulation can be found in Gu & Van Deun 2016.

set.seed(112)
testalgorithm <- 2 # 1: for testing VarSelectComDistPre.R
                    # 2: for testing VarSelectFriedman.R
I <- 28
J1 <- 144
J2 <-44
Jk <- c(J1, J2)
sumJk <- sum(J1 + J2)
R <- 5
Comm <- 1

PropNoise <- 0.05
Perc0 <- 0.3     # this is for generating sparseness for testing VarSelectFriedman.R
Perc0Com <- 0.9  # this is for generating sparseness for testing VarSelectComDistPre.R

NRSTARTS <- 5
Ndataset <- 1
MAXITER <- 400

LASSO <- 0.005
GROUPLASSO <- .01

Tucker <- array()
ProportionComm <- array()
ProportionDist <- array()
Proportion <- array()

PoutBest <- list()
ToutBest <- list()
TuckerValues <- array()
PoutBestPermu <- list()

for (Nd in 1:Ndataset){

  if (testalgorithm == 1){
    ############## For testing VarSelectComDistPre##############
    DATA1 <- matrix(rnorm(I*J1, mean = 0, sd = 1), I, J1)
    DATA2 <- matrix(rnorm(I*J2, mean = 0, sd = 1), I, J2)
    DATA <- cbind(DATA1, DATA2)

    svddata <- svd(DATA, R, R)
    Ttrue <- svddata$u
    PTrueC <- as.matrix(svddata$v) %*% diag(svddata$d[1:R])   #note that only the first R eigen values are needed.

    PTrueCBlock1 <- PTrueC[1:J1,]
    PTrueCBlock2 <- PTrueC[(J1+1):(J1+J2),]

    v1 <- c(2, 3)
    PTrueCBlock1[, v1] <- 0
    v2 <- c(4, 5)
    PTrueCBlock2[, v2] <- 0
    dist <- c(v1, v2)

    PTrueCnew <- rbind(PTrueCBlock1, PTrueCBlock2)
    Pcommon <- PTrueCnew[, Comm]
    lengthPtrue <- length(Pcommon)
    v <- sample(1:lengthPtrue, size = round(Perc0Com*(J1+J2)), replace=F)
    Pcommon[v] <- 0
    PTrueCnew[ , Comm] <- Pcommon

    XTrue <- Ttrue %*% t(PTrueCnew)
    SSXtrue <- sum(XTrue ^ 2)

    Noise <- matrix(rnorm(I*(J1+J2), mean = 0, sd = 1), I, J1+J2)
    SSNoise <- sum(Noise ^ 2)
    g <- sqrt(PropNoise*SSXtrue/(SSNoise-PropNoise*SSNoise))
    NoiseNew <- g*Noise
    SSNoiseNew <- sum(NoiseNew ^ 2)
    Xgenerate <- XTrue + NoiseNew
    SSXgenerate <- sum(Xgenerate ^ 2)
    NoiseVSgenerate <- SSNoiseNew/SSXgenerate

    GroupStructure <- PTrueCnew
    GroupStructure[ , Comm] <- 1


    Pout3d <- list()
    Tout3d <- list()
    LOSS <- array()
    LOSSVEC <- list()
    IterVec <- list()

    for (i in 1:NRSTARTS){
      VarSelectResult <- VarSelectComDistPre(Xgenerate, Jk, R, Comm, GroupStructure, LASSO, MAXITER)
      Pout3d[[i]] <- VarSelectResult$Pmatrix
      Tout3d[[i]] <- VarSelectResult$Tmatrix
      LOSS[i] <- VarSelectResult$Loss
      LOSSVEC[[i]] <- VarSelectResult$Lossvec
      IterVec[[i]] <- VarSelectResult$iter
    }
    k <- which(LOSS == max(LOSS))
    if (length(k)>1){
      pos <- sample(1:length(k), 1)
      k <- k[pos]
    }
    PoutBest[[Nd]] <- Pout3d[[k]]
    ToutBest[[Nd]] <- Tout3d[[k]]

    TuckerResults <- TuckerCoef(Ttrue, Tout3d[[k]])
    TuckerValues[Nd] <- TuckerResults$tucker_value

    indSelectedC <- which(PoutBest[[Nd]][, Comm] != 0)
    indDropedC <- which(PoutBest[[Nd]][, Comm] == 0)
    PTrueCnewComm <- PTrueCnew[, Comm]
    ProportionComm[Nd] <- (sum(PTrueCnewComm[indSelectedC] != 0) + sum(PTrueCnewComm[indDropedC] == 0)) / sumJk

    ############## END : For testing VarSelectComDistPre##############
  }
else if (testalgorithm == 2){

  DATA1 <- matrix(rnorm(I*J1, mean = 0, sd = 1), I, J1)
  DATA2 <- matrix(rnorm(I*J2, mean = 0, sd = 1), I, J2)
  DATA <- cbind(DATA1, DATA2)

  svddata <- svd(DATA, R, R)
  Ttrue <- svddata$u
  PTrueC <- as.matrix(svddata$v) %*% diag(svddata$d[1:R])   #note that only the first R eigen values are needed.

  PTrueCBlock1 <- PTrueC[1:J1,]
  PTrueCBlock2 <- PTrueC[(J1+1):(J1+J2),]

  v1 <- c(1, 2, 3)
  PTrueCBlock1[, v1] <- 0
  v2 <- c(4, 5)
  PTrueCBlock2[, v2] <- 0



  PTrueCBlock1_vec <- as.vector(PTrueCBlock1[, v2])
  v <- sample(1:(J1*2), size = round(Perc0*(J1*2)), replace=F)
  PTrueCBlock1_vec[v] <- 0
  PTrueCBlock1[, v2] <- matrix(PTrueCBlock1_vec, nrow = J1, ncol = 2)

  PTrueCBlock2_vec <- as.vector(PTrueCBlock2[, v1])
  v <- sample(1:(J2*3), size = round(Perc0*(J2*3)), replace=F)
  PTrueCBlock2_vec[v] <- 0
  PTrueCBlock2[, v1] <- matrix(PTrueCBlock2_vec, nrow = J2, ncol = 3)

  PTrueCnew <- rbind(PTrueCBlock1, PTrueCBlock2)

  XTrue <- Ttrue %*% t(PTrueCnew)
  SSXtrue <- sum(XTrue ^ 2)

  Noise <- matrix(rnorm(I*(J1+J2), mean = 0, sd = 1), I, J1+J2)
  SSNoise <- sum(Noise ^ 2)
  g <- sqrt(PropNoise*SSXtrue/(SSNoise-PropNoise*SSNoise))
  NoiseNew <- g*Noise
  SSNoiseNew <- sum(NoiseNew ^ 2)
  Xgenerate <- XTrue + NoiseNew
  SSXgenerate <- sum(Xgenerate ^ 2)
  NoiseVSgenerate <- SSNoiseNew/SSXgenerate


    Pout3d <- list()
    Tout3d <- list()
    LOSS <- array()
    LOSSVEC <- list()
    IterVec <- list()

    for (i in 1:NRSTARTS){
      VarSelectResult <- VarSelectFriedman(Xgenerate, Jk, R, LASSO, GROUPLASSO, MAXITER)
      Pout3d[[i]] <- VarSelectResult$Pmatrix
      Tout3d[[i]] <- VarSelectResult$Tmatrix
      LOSS[i] <- VarSelectResult$Loss
      LOSSVEC[[i]] <- VarSelectResult$Lossvec
      IterVec[[i]] <- VarSelectResult$iter
    }
    k <- which(LOSS == max(LOSS))
    if (length(k)>1){
      pos <- sample(1:length(k), 1)
      k <- k[pos]
    }
    PoutBest[[Nd]] <- Pout3d[[k]]
    ToutBest[[Nd]] <- Tout3d[[k]]

    TuckerResults <- TuckerCoef(Ttrue, Tout3d[[k]])
    TuckerValues[Nd] <- TuckerResults$tucker_value

    indSelectedC <- which(PoutBest[[Nd]] != 0)
    indDropedC <- which(PoutBest[[Nd]] == 0)
    Proportion[Nd] <- (sum(PTrueCnew[indSelectedC] != 0) + sum(PTrueCnew[indDropedC] == 0)) / (sumJk*R)
  }


}
