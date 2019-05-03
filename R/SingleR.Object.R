
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
#' @param use.singler.cluster.annot if FALSE, do not add the SingleR clusters annotations. Default is TRUE.
#' This is used for objects that where subseted from the full object, but the SingleR clusters annotations
#' where not regenerated.
#' 
#' @return S4 SingleR object that can be uploaded to the SingleR browser
convertSingleR2Browser = function(singler,use.singler.cluster.annot=T) { 
  
  ref.names = unlist(lapply(singler$singler,FUN=function(x) x$about$RefData))
  
  cell.names = rownames(singler$singler[[1]]$SingleR.single$labels)
  
  labels = as.data.frame(sapply(singler$singler,FUN=function(x) x$SingleR.single$labels))
  if (!is.null(singler$singler[[1]]$SingleR.single.main)) {
    labels.main = as.data.frame(sapply(singler$singler,FUN=function(x) x$SingleR.single.main$labels))
    labels = cbind(labels,labels.main)
    colnames(labels) = c(ref.names,paste0(ref.names,'.main'))
  } else {
    colnames(labels) = c(ref.names)
  }
  rownames(labels) = cell.names
  
  labels1 = data.frame()
  if (!is.null(singler$singler[[1]]$SingleR.single$labels1)) {
    labels1 = as.data.frame(sapply(singler$singler,FUN=function(x) x$SingleR.single$labels1))
    if (!is.null(singler$singler[[1]]$SingleR.single.main)) {
      labels1.main = as.data.frame(sapply(singler$singler,FUN=function(x) x$SingleR.single.main$labels1))
      labels1 = cbind(labels1,labels1.main)
      colnames(labels1) = c(ref.names,paste0(ref.names,'.main'))
    } else {
      colnames(labels1) = c(ref.names)
    }
    rownames(labels1) = cell.names
  }
  
  labels.clusters = data.frame()
  labels.clusters1 = data.frame()
  clusters = data.frame()
  
  if (use.singler.cluster.annot==T) {
    if (length(levels(singler$meta.data$clusters))>1) {
      if (!is.null(singler$singler[[1]]$SingleR.clusters)) {
        labels.clusters = as.data.frame(sapply(singler$singler,FUN=function(x) x$SingleR.clusters$labels))
        if (!is.null(singler$singler[[1]]$SingleR.clusters.main)) {
          labels.clusters.main = as.data.frame(sapply(singler$singler,FUN=function(x) x$SingleR.clusters.main$labels))
          labels.clusters = cbind(labels.clusters,labels.clusters.main)
          colnames(labels.clusters) = c(ref.names,paste0(ref.names,'.main'))
        } else {
          colnames(labels.clusters) = c(ref.names)
        }
        print(labels.clusters)
        rownames(labels.clusters) = levels(singler$meta.data$clusters)
      }
      
      if (!is.null(singler$singler[[1]]$SingleR.cluster$labels1)) {
        
        if (!is.null(singler$singler[[1]]$SingleR.clusters)) {
          labels.clusters1 = as.data.frame(sapply(singler$singler,FUN=function(x) x$SingleR.clusters$labels1))
          if (!is.null(singler$singler[[1]]$SingleR.clusters.main)) {
            labels.clusters.main = as.data.frame(sapply(singler$singler,FUN=function(x) x$SingleR.clusters.main$labels1))
            labels.clusters1 = cbind(labels.clusters1,labels.clusters.main)
            colnames(labels.clusters1) = c(ref.names,paste0(ref.names,'.main'))
          } else {
            colnames(labels.clusters1) = c(ref.names)
          }
          rownames(labels.clusters1) = levels(singler$meta.data$clusters)
          
        }
      }
    }
    
    clusters = data.frame(clusters=singler$meta.data$clusters)
    rownames(clusters) = cell.names
  } 
  scores = lapply(singler$singler,FUN=function(x) x$SingleR.single$scores)
  if (!is.null(singler$singler[[1]]$SingleR.single.main)) {
    scores.main = lapply(singler$singler,FUN=function(x) x$SingleR.single.main$scores)
    scores = c(scores,scores.main)
    names(scores) = c(ref.names,paste0(ref.names,'.main'))
  } else {
    names(scores) = c(ref.names)
  }
  
  ident = data.frame(orig.ident = singler$meta.data$orig.ident)
  rownames(ident) = cell.names
  
  expr = NULL
  if (!is.null(singler$seurat)) {
    if (packageVersion('Seurat')>=3) {
      expr = singler$seurat@assays$RNA@data
    } else {
      expr = singler$seurat@data
    }
  } 
  singler.small = new('SingleR',project.name=singler$meta.data$project.name,
                      xy=singler$meta.data$xy,
                      labels = labels,
                      labels.NFT = labels1,
                      labels.clusters = labels.clusters,        
                      labels.clusters.NFT = labels.clusters1,
                      scores = scores,
                      clusters = clusters,
                      ident = ident,
                      other = data.frame(singler$signatures),
                      expr = expr,
                      meta.data = c('Citation' = singler$singler[[1]]$about$Citation, 
                                    'Organism' = singler$singler[[1]]$about$Organism,
                                    'Technology' = singler$singler[[1]]$about$Technology)
  )
  singler.small
}

if (FALSE) {
  files = dir('~/Documents/SingleR/SingleR/data/',pattern='RData')
  for (i in 1:length(files)) { 
    print(files[i])
    load(paste0('~/Documents/SingleR/SingleR/data/',files[i]))
    singler.small = convertSingleR2Browser(singler)
    saveRDS(singler.small,file=paste0('~/Documents/SingleR/SingleRbrowseR/data/',gsub('RData','rds',files[i])))
  }
}

SingleR.SubsetS4 = function(singler,subsetdata) {
  s = singler
  s@xy = s@xy[subsetdata,]
  s@labels = s@labels[subsetdata,]
  s@labels.NFT = s@labels.NFT[subsetdata,]
  s@clusters = as.data.frame(s@clusters[subsetdata,])
  s@ident = as.data.frame(s@ident[subsetdata,])
  s@other = as.data.frame(s@other[subsetdata,])
  for (k in names(s@scores)) {
    s@scores[[k]] = s@scores[[k]][subsetdata,]
  }
  s@expr = s@expr[,subsetdata]
  s
}
