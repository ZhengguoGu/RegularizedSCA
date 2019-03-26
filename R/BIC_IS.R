#'Variable selection for regularized SCA using BIC and Index of Sparseness.
#'
#'\code{BIC_IS} is used to find the proper combination of Lasso and Group Lasso tuning parameters for regularized SCA based on BIC and Index of Sparseness.
#'
#'
#'@param DATA The concatenated data block, with rows representing subjects.
#'@param Jk A vector. Each element of this vector is the number of columns of a data block.
#'@param R The number of components (R>=2).
#'@param NRSTARTS The number of multistarts for this algorithm. The default value is 5.
#'@param MaxIter Maximum number of iterations for this algorithm. The default value is 400.
#'@param LassoSequence The range of Lasso tuning parameters. The default value is a sequence of 20 numbers from 0.00000001
#'to the smallest Lasso tuning parameter value that makes all the component loadings equal to zero. Note that by default the 20 numbers are equally spaced on the log scale. 
#'@param GLassoSequence The range of Group Lasso tuning parameters. The default value is a sequence of 20 numbers from 0.00000001
#'to the smallest Group Lasso tuning parameter value that makes all the component loadings equal to zero. Note that by default the 20 numbers are equally spaced (but not on the log scale). 
#'@param method "datablock" or "component". These are two options with respect to the grouping of the loadings as used in the Group Lasso penalty. 
#'If \code{method="component"}, the block-grouping of the coefficients is applied per component separately. If \code{method = "datablock"}, the grouping
#'is applied on the concatenated data block, with loadings of all components together. If \code{method} is missing, then the "component" method is used 
#'by default.  
#'
#'@return
#'\item{BIC}{A matrix of BIC values.}
#'\item{IS}{A matrix of IS values.}
#'\item{BIC_tuning}{Recommended tuning parameters for Lasso and Group Lasso based on BIC}
#'\item{IS_tuning}{Recommended tuning parameters for Lasso and Group Lasso based on IS}

#'@examples
#'\dontrun{
#'DATA1 <- matrix(rnorm(50), nrow=5)
#'DATA2 <- matrix(rnorm(100), nrow=5)  
#'DATA <- cbind(DATA1, DATA2)
#'Jk <- c(10, 20) 
#'BIC_IS(DATA, Jk, R=5, NRSTARTS = 40, MaxIter = 100)
#'}
#'
#'@references Gajjar, S., Kulahci, M., & Palazoglu, A. (2017). Selection of non-zero loadings in sparse principal component analysis. \emph{Chemometrics and Intelligent Laboratory Systems}, \emph{162}, 160-171.
#'@references Trendafilov, N. T. (2014). From simple structure to sparse components: a review. \emph{Computational Statistics}, \emph{29}(3-4), 431-454.
#'@references Zou, H., Hastie, T., & Tibshirani, R. (2006). Sparse principal component analysis. \emph{Journal of computational and graphical statistics}, \emph{15}(2), 265-286.
#'@references Croux, C., Filzmoser, P., & Fritz, H. (2013). Robust sparse principal component analysis. \emph{Technometrics}, \emph{55}(2), 202-214.
#'@references Guo, J., James, G., Levina, E., Michailidis, G., & Zhu, J. (2010). Principal component analysis with sparse fused loadings. \emph{Journal of Computational and Graphical Statistics}, \emph{19}(4), 930-946.
#'@export
BIC_IS <- function(DATA, Jk, R, NRSTARTS, MaxIter, LassoSequence, GLassoSequence, method){
  
  DATA <- data.matrix(DATA)
  n_sub <- dim(DATA)[1]
  
  if(missing(MaxIter)){
    MaxIter = 400
  }
  
  if(missing(NRSTARTS)){
    NRSTARTS = 5
  }
  
  if(missing(LassoSequence)){
    LassoSequence = exp(seq(from = log(0.00000001), to = log(maxLGlasso(DATA, Jk, R)$Lasso), length.out = 20))
  }
  
  if(missing(GLassoSequence)){
    GLassoSequence = seq(from = 0.00000001, to = maxLGlasso(DATA, Jk, R)$Glasso, length.out = 20)  #note that Glasso is not on the log scale, because it is not helpful.
  }
  
    if(missing(method)){
    method = "component"
  }
  VarSelect0 <- sparseSCA(DATA, Jk, R, LASSO = 0, GROUPLASSO = 0, MaxIter, NRSTARTS, method = "component")
  P_hat0 <- VarSelect0$Pmatrix
  T_hat0 <- VarSelect0$Tmatrix
  
  V_0 <- sum((DATA - T_hat0%*%t(P_hat0))^2)  # this is for BIC_Croux and BIC_GUO
  #error_var <- V_0 / n_sub  #this is for BIC_Guo
  
  V_oo <- sum(DATA^2)  # this is for Index of sparseness (IS)
  V_s <- sum((T_hat0%*%t(P_hat0))^2)  # this is for IS
  
  BIC_Croux <- matrix(NA, length(LassoSequence), length(GLassoSequence))
  #BIC_Guo <- matrix(NA, length(LassoSequence), length(GLassoSequence))
  IS <- matrix(NA, length(LassoSequence), length(GLassoSequence))
  for(i in 1:length(LassoSequence)){
    for(j in 1:length(GLassoSequence)){
      
      VarSelect <- sparseSCA(DATA, Jk, R, LASSO = LassoSequence[i], GROUPLASSO = GLassoSequence[j], MaxIter, NRSTARTS, method)
      P_hat <- VarSelect$Pmatrix
      T_hat <- VarSelect$Tmatrix
      V_tilde <- sum((DATA - T_hat %*% t(P_hat))^2)
      
      BIC_Croux[i, j] <- V_tilde / V_0 + sum(P_hat != 0) * log(n_sub) / n_sub  # this is the BIC index prop0sed by Croux et al.
      #BIC_Guo[i, j] <- V_tilde / error_var + sum(P_hat != 0) * log(n_sub)      # this is the BIC index proposed by Guo et al.
      
      V_a <- sum((T_hat %*% t(P_hat))^2)  # this is for IS
      IS[i, j] <- V_a * V_s / V_oo^2 * sum(P_hat == 0) /(sum(Jk) * R)          # this is index of sparseness
    }
  }
  
  Croux_index <- which(BIC_Croux == min(BIC_Croux), arr.ind = T)
  Lasso_croux <- max(LassoSequence[Croux_index[1]])  #max() is used in case multiple lasso values are chosen. 
  GLasso_croux <- max(GLassoSequence[Croux_index[2]])
  
  IS_index <- which(IS == max(IS), arr.ind = T)
  Lasso_IS <- max(LassoSequence[IS_index[1]])
  Glasso_IS <- max(GLassoSequence[IS_index[2]])
  
  BIC_IS_list <- list()
  BIC_IS_list$BIC <- BIC_Croux
  BIC_IS_list$IS <- IS
  BIC_IS_list$BIC_tuning <- c(Lasso_croux, GLasso_croux)
  names(BIC_IS_list$BIC_tuning) <- c("Lasso", "Group Lasso")
  BIC_IS_list$IS_tuning <- c(Lasso_IS, Glasso_IS)
  names(BIC_IS_list$IS_tuning) <- c("Lasso", "Group Lasso")
  attr(BIC_IS_list, "class") <- "BIC_IS"
  return(BIC_IS_list)
}