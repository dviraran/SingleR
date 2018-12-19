
#' An S4 class to represent SingleR objects
#' This object can be upload to the SingleR browser (http)
#'
#' @slot project.name a string. Must be the same as the file name. 
#' @slot xy matrix of the xy coordinates. Row names must be single-cell names. 
#' @slot labels data.frame of the SingleR annotations. Rows are single-cell, a column for each 
#' annotation by reference. Column names are the reference name (e.g. Immgen, Immgen.main, Blueprint+Encode, etc.)
#' @slot labels.NFT Same as labels. Annotation before fine-tuning. 
#' @slot labels.clusters data.frame. Same as labels, but a row for each cluster, not single-cell.
#' @slot labels.clusters.NFT Same as labels.clusters. Annotation before fine-tuning.
#' @slot clusters data.frame. Rows are single-cell, columns are clusters. Multiple clsutering are possible. 
#' @slot ident data.frame. Rows are single-cell, columns are caterogical meta.data (e.g. orig.ident). 
#' @slot other data.frame. Rows are single-cell, columns are continous variables (e.g. cell cycle score).
#' @slot scores list. Each item in the list corresponds to the column labels in labels. Items are the SingleR score matrices.
#' @slot meta.data character array. Each item is a named charachter (e.g. Citation, Organism, Technology) 
#' @slot expr a matrix. The single-cell expression matrix.
#' 
setClass("SingleR", representation(project.name = "character", 
                                   xy = "matrix", 
                                   labels='data.frame',
                                   labels.NFT='data.frame',
                                   labels.clusters='data.frame',
                                   labels.clusters.NFT='data.frame',
                                   clusters='data.frame',
                                   ident='data.frame',
                                   other = 'data.frame',
                                   scores='list',
                                   meta.data='character',
                                   expr = 'ANY'))

#' Convert the SingleR output object to a SingleR S4 object that can be uploaded to the SingleR browser.
#'
#' @param singler the SingleR output object
#'
#' @return S4 SingleR object that can be uploaded to the SingleR browser
convertSingleR2Browser = function(singler) { 
  singler.small = new('SingleR',project.name=singler$meta.data$project.name,
                      xy=singler$meta.data$xy,
                      labels = data.frame('Immgen'=singler$singler[[1]]$SingleR.single$labels,
                                          'Immgen.main'=singler$singler[[1]]$SingleR.single.main$labels,
                                          'MouseRNAseq'=singler$singler[[2]]$SingleR.single$labels,
                                          'MouseRNAseq.main'=singler$singler[[2]]$SingleR.single.main$labels
                      ),        
                      labels.NFT = data.frame('Immgen'=singler$singler[[1]]$SingleR.single$labels1,
                                              'Immgen.main'=singler$singler[[1]]$SingleR.single.main$labels1,
                                              'MouseRNAseq'=singler$singler[[2]]$SingleR.single$labels1,
                                              'MouseRNAseq.main'=singler$singler[[2]]$SingleR.single.main$labels1
                      ),
                      labels.clusters = data.frame('Immgen'=singler$singler[[1]]$SingleR.clusters$labels,
                                                   'Immgen.main'=singler$singler[[1]]$SingleR.clusters.main$labels,
                                                   'MouseRNAseq'=singler$singler[[2]]$SingleR.clusters$labels,
                                                   'MouseRNAseq.main'=singler$singler[[2]]$SingleR.clusters.main$labels
                      ),        
                      labels.clusters.NFT = data.frame('Immgen'=singler$singler[[1]]$SingleR.clusters$labels1,
                                                       'Immgen.main'=singler$singler[[1]]$SingleR.clusters.main$labels1,
                                                       'MouseRNAseq'=singler$singler[[2]]$SingleR.clusters$labels1,
                                                       'MouseRNAseq.main'=singler$singler[[2]]$SingleR.clusters.main$labels1
                      ),
                      scores = list('Immgen' = singler$singler[[1]]$SingleR.single$scores,
                                    'Immgen.main' = singler$singler[[1]]$SingleR.single.main$scores,
                                    'MouseRNAseq' = singler$singler[[2]]$SingleR.single$scores,
                                    'MouseRNAseq.main' = singler$singler[[2]]$SingleR.single.main$scores
                      ),
                      clusters = data.frame(clusters=singler$meta.data$clusters),
                      ident = data.frame(orig.ident = singler$meta.data$orig.ident),
                      other = data.frame(singler$signatures),
                      expr = singler$seurat@data,
                      meta.data = c('Citation' = singler$singler[[2]]$about$Citation, 
                                    'Organism' = singler$singler[[2]]$about$Organism,
                                    'Technology' = singler$singler[[2]]$about$Technology)
  )
  singler.small
}

if (FALSE) {
for (i in 43:length(files)) { 
  load(paste0('~/Documents/SingleR/SingleRbrowseR/data/',files[i]))
  singler.small = convertSingleR2Browser(singler)
  saveRDS(singler.small,file=paste0('~/Documents/SingleR/SingleRbrowseR/data/',gsub('RData','rds',files[i])))
}
}
