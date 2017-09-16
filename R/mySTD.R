#' Standardize the given data matrix per column, over the rows. (This function will no longer be exported
#' and it will eventually removed from the package.)
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

mySTD <- function(DATA) {
  
  DATA <- data.matrix(DATA)
  nrow_data <- dim(DATA)[1]
  ncol_data <- dim(DATA)[2]
  v <- matrix(1, nrow = nrow_data, ncol = nrow_data)  # note here: it's a square matrix
  DATAc <- DATA - v %*% DATA/nrow_data
  CP <- v %*% (DATAc ^ 2)
  STDDATA <- sqrt(nrow_data-1) * DATAc / (CP ^ 0.5)

  return(STDDATA)
}
