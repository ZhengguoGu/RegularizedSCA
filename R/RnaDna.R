#' RNA and DNA data of 44 samples
#'
#' This dataset contains the gene expression and DNA copy number data
#' of 44 samples.
#'
#' @format The dataset contains the following list:
#' \describe{
#'  \item{dna}{a 44x269 matrix of CGH samples x spots.}
#'  \item{rna}{a 2459x44 matrix of samples x genes.}
#'  \item{chrom}{a 269-vector of chromosomal location of each CGH spot.}
#'  \item{nuc}{a 269-vector of nucleotide position for each CGH spot.}
#'  \item{gene}{a 2459-vector with an accession number for each gene.}
#'  \item{genenames}{a 2459-vector with a name for each gene.}
#'  \item{genechr}{a 2459-vector with a chromosomal location for each gene.}
#'  \item{genedesc}{a 2459-vector with a description for each gene.}
#'  \item{genepos}{a 2459-vector with a nucleotide position for each gene.}
#' }
#'
#' @source This dataset is part of a much larger dataset, namely \code{breastdata}, in the \code{R} package 'PMA' by Witten, D.,
#' Tibshirani, R., Gross, S., and Narasimhan, B. The \code{breastdata} dataset was published in the paper
#' Chin, K., DeVries, S., Fridlyand, J., Spellman, P., Roydasgupta, R., Kuo,W.-L., Lapuk, A., Neve,
#' R., Qian, Z., Ryder, T., Chen, F., Feiler, H., Tokuyasu, T., Kingsley, C., Dairkee, S., Meng, Z., Chew,
#' K., Pinkel, D., Jain, A., Ljung, B., Esserman, L., Albertson, D.,Waldman, F. & Gray, J. (2006),
#' 'Genomic and transcriptional aberrations linked to breast cancer pathophysiologies', Cancer Cell
#' 10, 529-541. The original raw data are publicly available at ArrayExpress \url{http://www.ebi.ac.uk/arrayexpress/},
#' with accession number E-TABM-158.
"RnaDna"




