# mySTD function standardizes the Data per column, over the rows

mySTD <- function(DATA) {
  nrow_data <- dim(DATA)[1]
  ncol_data <- dim(DATA)[2]
  v <- matrix(1, nrow = nrow_data, ncol = nrow_data)  # note here: it's a square matrix
  DATAc <- DATA - v %*% DATA/nrow_data
  CP <- v %*% (DATAc ^ 2)
  STDDATA <- sqrt(nrow_data) * DATAc / (CP ^ 0.5)

  return(STDDATA)
}
