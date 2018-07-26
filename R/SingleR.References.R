#' SingleR: A package for reference-based single-cell RNA-seq annotation
#' 
#' @docType package
#' @name SingleR
NULL

#' Immgen reference dataset for mouse
#'
#' @format A list with the following items:
#' \describe{
#'   \item{data}{Gene expression data matrix}
#'   \item{types}{vector of annotation per column in the matrix}
#'   \item{main_types}{vector of broad annotations per column in the matrix}
#'   \item{name}{name if the reference set}
#'   \item{sd.thres}{threshold of the standard deviation for genes to use in SingleR}
#'   \item{de.genes.main}{list of lists of differentially expressed genes between every two cell types in main_types}
#'   \item{de.genes}{list of lists of differentially expressed genes between every two cell types in types}
#' }
"immgen"

#' Blueprint+Encode reference dataset for human
#'
#' @format A list with the following items:
#' \describe{
#'   \item{data}{Gene expression data matrix}
#'   \item{types}{vector of annotation per column in the matrix}
#'   \item{main_types}{vector of broad annotations per column in the matrix}
#'   \item{name}{name if the reference set}
#'   \item{sd.thres}{threshold of the standard deviation for genes to use in SingleR}
#'   \item{de.genes.main}{list of lists of differentially expressed genes between every two cell types in main_types}
#'   \item{de.genes}{list of lists of differentially expressed genes between every two cell types in types}
#' }
"blueprint_encode"

#' Mouse-RNAseq reference dataset for mouse
#'
#' @format A list with the following items:
#' \describe{
#'   \item{data}{Gene expression data matrix}
#'   \item{types}{vector of annotation per column in the matrix}
#'   \item{main_types}{vector of broad annotations per column in the matrix}
#'   \item{name}{name if the reference set}
#'   \item{sd.thres}{threshold of the standard deviation for genes to use in SingleR}
#'   \item{de.genes.main}{list of lists of differentially expressed genes between every two cell types in main_types}
#'   \item{de.genes}{list of lists of differentially expressed genes between every two cell types in types}
#' }
"mouse.rnaseq"

#' Human Primary Cell Atlas (HPCA) reference dataset for human
#'
#' @format A list with the following items:
#' \describe{
#'   \item{data}{Gene expression data matrix}
#'   \item{types}{vector of annotation per column in the matrix}
#'   \item{main_types}{vector of broad annotations per column in the matrix}
#'   \item{name}{name if the reference set}
#'   \item{sd.thres}{threshold of the standard deviation for genes to use in SingleR}
#'   \item{de.genes.main}{list of lists of differentially expressed genes between every two cell types in main_types}
#'   \item{de.genes}{list of lists of differentially expressed genes between every two cell types in types}
#' }
"hpca"

#' Reference matrix for Kang et al. classification
#'
#' @format A list with the reference matrix
"cell.type.classification"

#' Human signatures
#'
#' @format A GeneSetCollection of 5 signatures
"human.egc"

#' Length of human genes for TPM calculation
#'
#' @format A named vector of gene lengths
"human_lengths"

#' Length of mouse genes for TPM calculation
#'
#' @format A named vector of gene lengths
"mouse_lengths"

#' Human signatures
#'
#' @format A GeneSetCollection of 5 signatures
"human.egc"

#' Mouse signatures
#'
#' @format A GeneSetCollection of 5 signatures
"mouse.egc"

SingleR.numCores = min(detectCores(all.tests = FALSE,logical = TRUE)-1,16)

# Colors
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
singler.colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, 
                               rownames(qual_col_pals)))
singler.colors = singler.colors[c(-4,-27)]
singler.colors = c(singler.colors,singler.colors,singler.colors)