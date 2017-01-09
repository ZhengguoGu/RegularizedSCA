## This is to test the algorithm CDfriedman.R

#----------------------- Situation CD.1 -----------------------------------
# 2 blocks, 5 components:
# 0 0 0 1 1
# 1 1 1 0 0
# Furthermore, the distinctive components are sparse.
# (note: sparse distinctive component here means some of the loadings in the distictive
# component are 0's. See the code below.)

set.seed(112)

I <- 28
J1 <- 144
J2 <-44
Jk <- c(J1, J2)
sumJk <- sum(J1 + J2)
R <- 5
Comm <- 1 #not used

PropNoise <- 0.05
Perc0 <- 0.3

NRSTARTS <- 20
Ndataset <- 50
MAXITER <- 400

LASSO <- 3.8  #.3
GROUPLASSO <- .6  #.1

Tucker <- array()
ProportionComm <- array()
ProportionDist <- array()
Proportion <- array()

PoutBest <- list()
ToutBest <- list()
TuckerValues <- array()
PoutBestPermu <- list()

for (Nd in 1:Ndataset){

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
    VarSelectResult <- CDfriedman(Xgenerate, Jk, R, LASSO, GROUPLASSO, MAXITER)
    Pout3d[[i]] <- VarSelectResult$Pmatrix
    Tout3d[[i]] <- VarSelectResult$Tmatrix
    LOSS[i] <- VarSelectResult$Loss
    LOSSVEC[[i]] <- VarSelectResult$Lossvec
    IterVec[[i]] <- VarSelectResult$iter
  }
  k <- which(LOSS == min(LOSS))
  if (length(k)>1){
    pos <- sample(1:length(k), 1)
    k <- k[pos]
  }
  PoutBest[[Nd]] <- Pout3d[[k]]
  ToutBest[[Nd]] <- Tout3d[[k]]

  TuckerResults <- TuckerCoef(Ttrue, Tout3d[[k]])
  TuckerValues[Nd] <- TuckerResults$tucker_value
  PoutBest[[Nd]] <- PoutBest[[Nd]][, TuckerResults$perm]

  indSelectedC <- which(PoutBest[[Nd]] != 0)
  indDropedC <- which(PoutBest[[Nd]] == 0)
  Proportion[Nd] <- (sum(PTrueCnew[indSelectedC] != 0) + sum(PTrueCnew[indDropedC] == 0)) / (sumJk*R)
}

#---------------------- Situation CD.2 ---------------------------------------
# 3 blocks, 3 components:
# 1 1 0
# 1 0 1
# 1 1 0
# Furthermore, the distinctive components are sparse.
# (note: sparse distinctive component here means some of the loadings in the distictive
# component are 0's. See the code below.)

set.seed(112)

I <- 28
J1 <- 44
J2 <-44
J3 <- 44
Jk <- c(J1, J2, J3)
sumJk <- sum(J1 + J2 + J3)
R <- 3
Comm <- 1

PropNoise <- 0.05
Perc0 <- 0.3

NRSTARTS <- 20
Ndataset <- 50
MAXITER <- 400

LASSO <- 0.2
GROUPLASSO <- 0.1

Tucker <- array()
ProportionComm <- array()
ProportionDist <- array()
Proportion <- array()

PoutBest <- list()
ToutBest <- list()
TuckerValues <- array()
PoutBestPermu <- list()

for (Nd in 1:Ndataset){

  DATA1 <- matrix(rnorm(I*J1, mean = 0, sd = 1), I, J1)
  DATA2 <- matrix(rnorm(I*J2, mean = 0, sd = 1), I, J2)
  DATA3 <- matrix(rnorm(I*J2, mean = 0, sd = 1), I, J3)
  DATA <- cbind(DATA1, DATA2, DATA3)

  svddata <- svd(DATA, R, R)
  Ttrue <- svddata$u
  PTrueC <- as.matrix(svddata$v) %*% diag(svddata$d[1:R])   #note that only the first R eigen values are needed.

  PTrueCBlock1 <- PTrueC[1:J1,]
  PTrueCBlock2 <- PTrueC[(J1+1):(J1+J2),]
  PTrueCBlock3 <- PTrueC[(J1+J2+1):(J1+J2+J3),]

  # 3 blocks, 3 components:
  # 1 1 0
  # 1 0 1
  # 1 1 0
  PTrueCBlock1[, 3] <- 0
  PTrueCBlock2[, 2] <- 0
  PTrueCBlock3[, 3] <- 0
  PTrueCnew <- rbind(PTrueCBlock1, PTrueCBlock2, PTrueCBlock3)

  v1 <- sample(1:sumJk, size = round(Perc0*sumJk), replace=F)
  PTrueCnew[v1, 1] <- 0
  v2 <- sample(1:J1, size = round(Perc0*J1), replace=F)
  PTrueCnew[v2, 2] <- 0
  v2 <- sample(89:132, size = round(Perc0*J3), replace=F)
  PTrueCnew[v2, 2] <- 0
  v3 <- sample(45:88, size = round(Perc0*44), replace=F)
  PTrueCnew[v3, 3] <- 0

  XTrue <- Ttrue %*% t(PTrueCnew)
  SSXtrue <- sum(XTrue ^ 2)

  Noise <- matrix(rnorm(I*(J1+J2+J3), mean = 0, sd = 1), I, J1+J2+J3)
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
    VarSelectResult <- CDfriedman(Xgenerate, Jk, R, LASSO, GROUPLASSO, MAXITER)
    Pout3d[[i]] <- VarSelectResult$Pmatrix
    Tout3d[[i]] <- VarSelectResult$Tmatrix
    LOSS[i] <- VarSelectResult$Loss
    LOSSVEC[[i]] <- VarSelectResult$Lossvec
    IterVec[[i]] <- VarSelectResult$iter
  }
  k <- which(LOSS == min(LOSS))
  if (length(k)>1){
    pos <- sample(1:length(k), 1)
    k <- k[pos]
  }
  PoutBest[[Nd]] <- Pout3d[[k]]
  ToutBest[[Nd]] <- Tout3d[[k]]

  TuckerResults <- TuckerCoef(Ttrue, Tout3d[[k]])
  TuckerValues[Nd] <- TuckerResults$tucker_value
  PoutBest[[Nd]] <- PoutBest[[Nd]][, TuckerResults$perm]

  indSelectedC <- which(PoutBest[[Nd]] != 0)
  indDropedC <- which(PoutBest[[Nd]] == 0)
  Proportion[Nd] <- (sum(PTrueCnew[indSelectedC] != 0) + sum(PTrueCnew[indDropedC] == 0)) / (sumJk*R)
}
