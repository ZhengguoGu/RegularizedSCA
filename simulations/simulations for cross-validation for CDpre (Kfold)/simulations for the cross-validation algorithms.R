#-------------- CV-Situation 1: two block 5 components:
# 1 0 0 1 1
# 1 1 1 0 0
# Note that the first component is common component but very sparse.

set.seed(112)

I <- 28
J1 <- 144
J2 <-44
Jk <- c(J1, J2)
sumJk <- sum(J1 + J2)
R <- 5
Comm <- 1

PropNoise <- 0.05
Perc0Com <- 0.9  # this is for generating sparseness for testing CDpre.R

NRSTARTS <- 20
Ndataset <- 50
MAXITER <- 400

LASSO <- 0.3

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

#----
target <- matrix(c(1,1,0,1,0,1,1,0,1,0), 2, 5)
#      [,1] [,2] [,3] [,4] [,5]
#[1,]    1    0    0    1    1
#[2,]    1    1    1    0    0

comStr <- component_structure(Jk, R, target)
#cros_results <- cv_CDpreKf(Xgenerate, Jk, R, CommPosition=Comm, component_structure=comStr, MaxIter=MAXITER, NRSTARTS=NRSTARTS)

cros_results <- cv_CDpreKf(Xgenerate, Jk, R, CommPosition=Comm, component_structure=comStr, MaxIter=MAXITER, NRSTARTS=NRSTARTS, LassoSequence = seq(from=0.0001, to=1.5, length.out = 10))


Pout3d <- list()
Tout3d <- list()
LOSS <- array()
LOSSVEC <- list()
IterVec <- list()
for (i in 1:NRSTARTS){
  VarSelectResult <- CDpre(Xgenerate, Jk, R, Comm, GroupStructure=comStr, LASSO=.5, MAXITER)
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
PoutBest <- Pout3d[[k]]
ToutBest <- Tout3d[[k]]


#------- cv-Situation 2 ---------------------------------------------------------------------

# 1 0 0 1 1
# 1 1 1 0 0
# 1 1 0 1 0
# Note that the first component is common component but very sparse.

set.seed(112)

I <- 28
J1 <- 44
J2 <- 44
J3 <- 44
Jk <- c(J1, J2, J3)
sumJk <- sum(J1 + J2 + J3)
R <- 5
Comm <- 1

PropNoise <- 0.05
Perc0Com <- 0.9  # this is for generating sparseness for testing CDpre.R

NRSTARTS <- 20
Ndataset <- 50
MAXITER <- 400

LASSO <- 0.3

DATA1 <- matrix(rnorm(I*J1, mean = 0, sd = 1), I, J1)
DATA2 <- matrix(rnorm(I*J2, mean = 0, sd = 1), I, J2)
DATA3 <- matrix(rnorm(I*J3, mean = 0, sd = 1), I, J3)
DATA <- cbind(DATA1, DATA2, DATA3)

svddata <- svd(DATA, R, R)
Ttrue <- svddata$u
PTrueC <- as.matrix(svddata$v) %*% diag(svddata$d[1:R])   #note that only the first R eigen values are needed.

PTrueCBlock1 <- PTrueC[1:J1,]
PTrueCBlock2 <- PTrueC[(J1+1):(J1+J2),]
PTrueCBlock3 <- PTrueC[(J1+J2+1):(J1+J2+J3),]

# 1 0 0 1 1
# 1 1 1 0 0
# 1 1 0 1 0
v1 <- c(2, 3)
PTrueCBlock1[, v1] <- 0
v2 <- c(4, 5)
PTrueCBlock2[, v2] <- 0
v3 <- c(3, 5)
PTrueCBlock3[, v3] <- 0

dist <- c(v1, v2) #dont need v3 here

PTrueCnew <- rbind(PTrueCBlock1, PTrueCBlock2, PTrueCBlock3)
Pcommon <- PTrueCnew[, Comm]
lengthPtrue <- length(Pcommon)
v <- sample(1:lengthPtrue, size = round(Perc0Com*(J1+J2+J3)), replace=F)
Pcommon[v] <- 0
PTrueCnew[ , Comm] <- Pcommon

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

#----
target <- matrix(c(1,1,1,0,1,1,0,1,0,1,0,1,1,0,0), 3, 5)
#      [,1] [,2] [,3] [,4] [,5]
#[1,]    1    0    0    1    1
#[2,]    1    1    1    0    0
#[3,]    1    1    0    1    0

comStr <- component_structure(Jk, R, target)
#cros_results <- cv_CDpreKf(Xgenerate, Jk, R, CommPosition=Comm, component_structure=comStr, MaxIter=MAXITER, NRSTARTS=NRSTARTS)

cros_results <- cv_CDpreKf(Xgenerate, Jk, R, CommPosition=Comm, component_structure=comStr, MaxIter=MAXITER, NRSTARTS=NRSTARTS, LassoSequence = seq(from=0.0001, to=1, length.out = 10))


Pout3d <- list()
Tout3d <- list()
LOSS <- array()
LOSSVEC <- list()
IterVec <- list()
for (i in 1:NRSTARTS){
  VarSelectResult <- CDpre(Xgenerate, Jk, R, Comm, GroupStructure=comStr, LASSO=.44, MAXITER)
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
PoutBest <- Pout3d[[k]]
ToutBest <- Tout3d[[k]]

TuckerResults <- TuckerCoef(Ttrue, ToutBest)
TuckerValues <- TuckerResults$tucker_value
PoutBest <- PoutBest[, TuckerResults$perm]

indSelectedC <- which(PoutBest != 0)
indDropedC <- which(PoutBest == 0)
Proportion <- (sum(PTrueCnew[indSelectedC] != 0) + sum(PTrueCnew[indDropedC] == 0)) / (sumJk*R)


