#' Standardize the given data matrix per column, over the rows, with multiple imputation for missing data. 
#'
#' @param DATA A data matrix
#' @return a standardized matrix
#' @examples
#' \dontrun{
#' mySTD(matrix(1:12, nrow = 3, ncol = 4))
#' }
#' @note
#' More details regarding data pre-processing, please see:
#'
#' Van Deun, K., Smilde, A.K., van der Werf, M.J., Kiers, H.A.L., & Mechelen, I.V. (2009). A structured overview of simultaneous component based data integration. \emph{BMC Bioinformatics}, 10:246.
#' 
#' The missing values are handled by means of Multivariate Imputation by Chained Equations (MICE). The number of multiple imputation is 5. More details see:
#' 
#' Buuren, S. V., & Groothuis-Oudshoorn, K. (2010). mice: Multivariate imputation by chained equations in R. \emph{Journal of statistical software}, 1-68.
#'@export
mySTD <- function(DATA) {
  
  DATA <- data.matrix(DATA)
  if(sum(is.na(DATA))>0){
    
    print("The data contain missing values, which will be imputed.")
    DATA <- mice::complete(mice::mice(DATA, seed=1, printFlag = F))
    
    DATA <- data.matrix(DATA)
    
  }
  nrow_data <- dim(DATA)[1]
  ncol_data <- dim(DATA)[2]
  v <- matrix(1, nrow = nrow_data, ncol = nrow_data)  # note here: it's a square matrix
  DATAc <- DATA - v %*% DATA/nrow_data
  CP <- v %*% (DATAc ^ 2)
  STDDATA <- sqrt(nrow_data-1) * DATAc / (CP ^ 0.5)

  return(STDDATA)
}
