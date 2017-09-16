#'Display a summary of the results of \code{structuredSCA()}.
#'
#'@param object Object of class inheriting from 'structuredSCA'.
#'@param ...  Argument to be passed to or from other methods. 
#'@examples
#'\dontrun{
#'## S3 method for class 'structuredSCA'
#'summary(object)
#'}
#'
#'@export
summary.structuredSCA <- function(object, ...){
  
  cat(sprintf("\nThe estimated P matrix is\n"))
  print(object$Pmatrix)
  
  cat(sprintf("\nThe estimated T matrix is\n"))
  print(object$Tmatrix)
}