#'Display a summary of the results of \code{sparseSCA()}.
#'
#'@param object Object of class inheriting from 'sparseSCA'.
#'@param ...  Argument to be passed to or from other methods. 
#'@examples
#'\dontrun{
#'## S3 method for class 'sparseSCA'
#'summary(object)
#'}
#'
#'@export
summary.sparseSCA <- function(object, ...){
  
  cat(sprintf("\nThe estimated P matrix is\n"))
  print(object$Pmatrix)
  
  cat(sprintf("\nThe estimated T matrix is\n"))
  print(object$Tmatrix)
}