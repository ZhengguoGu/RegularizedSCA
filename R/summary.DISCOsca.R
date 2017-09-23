#'Display a summary of the results of \code{DISCOsca()}.
#'
#'@param object Object of class inheriting from 'DISCOsca'.
#'@param disp The default is \code{simple}; in this case, the best-fitted common/distinctive 
#'            structure is displayed. 
#'            If \code{full}, then information is displayed regarding 1) the best-fitted 
#'            common/distinctive structure, 2) Estimated component score matrix (i.e., T),
#'            3) Estimated component loading matrix (i.e., P), and 4) Proportion of variance 
#'            per component.
#'@param ...  Argument to be passed to or from other methods. 
#'@examples
#'\dontrun{
#'## S3 method for class 'DISCOsca'
#'summary(object, disp="full")
#'}
#'
#'@export
summary.DISCOsca <- function(object, disp, ...){

if(missing(disp)){
  disp <- "simple"
}

if(disp == "simple"){
  
  cat(sprintf("\nThe best-fitted common/distinctive structure is\n"))
  print(object$comdist)
  cat(sprintf("\nNote: 0 indicates that the loadings of the entire column should be zero's.\n"))
}else if(disp == "full"){
  
  cat(sprintf("\nThe best-fitted common/distinctive structure is\n"))
  print(object$comdist)
  cat(sprintf("\nNote: 0 indicates that the loadings of the entire column should be zero's.\n"))
  
  cat(sprintf("\nThe estimated component score matrix (i.e., T) is\n"))
  print(object$Trot_best)
  
  cat(sprintf("\nThe estimated component loading matrix (i.e., P) is\n"))
  print(object$Prot_best)
  
  cat(sprintf("\nroportion of variance per component:\n"))
  print(object$propExp_component)
  
}else{
  stop("either simple or full")
}

}