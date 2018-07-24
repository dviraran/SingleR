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


# Colors
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
singler.colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, 
                               rownames(qual_col_pals)))
singler.colors = singler.colors[c(-4,-27)]
singler.colors = c(singler.colors,singler.colors,singler.colors)

#' Calculate median values by groups
#'
#' @param mat a matrix of numeric values NxM
#' @param groups a vector of size M with annotation of columns in mat
#'
#' @return a matrix with one column per group
medianMatrix = function(mat,groups) {
  fgroups = levels(factor(groups))
  mat.group <- do.call(cbind, lapply(fgroups, function(g) {
    A = groups==g
    if(sum(A)==1) {
      mat[,A]
    } else {
      rowMedians(mat[,A],na.rm=T)
    }
  }))
  colnames(mat.group) = fgroups
  rownames(mat.group) = rownames(mat)
  mat.group
}

#' Calculate a quantile value by groups
#'
#' @param mat a matrix of numeric values NxM
#' @param groups a vector of size M with annotation of columns in mat
#' @param q the quantile to calculate per group
#'
#' @return a matrix with one column per group
quantileMatrix = function(mat,groups,q) {
  fgroups = levels(factor(groups))
  mat.group <- do.call(cbind, lapply(fgroups, function(g) {
    A = groups==g
    if (nrow(mat)==1) {
      quantile(mat[A],na.rm=T,probs=q)
    } else {
      if(sum(A)==1) {
        mat[,A]
      } else {
        rowQuantiles(mat[,A],na.rm=T,probs=q)
      }
    }
  }))
  colnames(mat.group) = fgroups
  rownames(mat.group) = rownames(mat)
  mat.group
}

#' SingleR fine-tuning function. Takes top scores per cell type, finds variable genes in the
#' reference dataset for those cell types and recalculates the scores. This is performed until
#' only one cell type is left.
#'
#' @param sc_data  the single-cell RNA-seq data set as a matrix with genes as rownames. If the data if from a full-length platform, counts must be normalized to gene length (TPM, RPKM, FPKM, etc.).
#' @param ref_data the reference dataset with genes as rownames. Gene names must be in the same format as the single cell data (if sc_data uses genes symbols, ref_data must have the same)
#' @param types a list of cell type names corresponding to ref_data. Number of elements in types must be equal to number of columns in ref_data
#' @param scores the calculated scores per cell type produced by SingleR.ScoreData function
#' @param quantile.use correlation coefficients are aggregated for multiple cell types in the reference data set. This parameter allows to choose how to sort the cell types scores, by median (0.5) or any other number between 0 and 1. The default is 0.9.
#' @param fine.tune.thres the fine tuning step performs the scoring procedure for the top scoring labels using only genes that vary between those cell types in the reference data set. The top labels are those with a score lower than the top score by less than fine.tune.thres
#' @param genes list of genes to use for the annotations, or a method for extracting the genes from the data. Available methods are standard deviation and gene dispersion (\code{"sd"} or \code{"de"}). default is "de".
#' @param sd.thres if genes=='sd' then this is the threshold for defining a variable gene.
#' @param mean_mat if genes='de' then this is a matrix with one column per cell type, each column is the
#' average of all samples for the corresponding cell type.
#'
#' @return the top labels for each single cell
SingleR.FineTune <- function(sc_data,ref_data,types,scores,quantile.use,
                             fine.tune.thres,genes,sd.thres,mean_mat) {
  N = dim(sc_data)[2]
  numCores = min(detectCores(all.tests = FALSE, logical = TRUE)-1,16)
  print(paste("Fine-tuning round on top cell types (using", numCores, 
              "CPU cores):"))
  labels = pbmclapply(1:N,FUN=function(i){
    max_score = max(scores[i,])
    topLabels = names(scores[i,scores[i,]>=max_score-fine.tune.thres])
    if (length(topLabels)==0) {
      return (names(which.max(scores[i,])))
    } else {
      while(length(topLabels)>1) {
        topLabels = fineTuningRound(topLabels,types,ref_data,genes,
                                    mean_mat[,topLabels],sd.thres,
                                    sc_data[,i],quantile.use,fine.tune.thres)
      }
      return (topLabels)
    }
  },mc.cores=numCores)
  labels = as.matrix(unlist(labels))
  
  if (dim(sc_data)[2]>1) {
    rownames(labels)=t(colnames(sc_data))
  }
  return(labels)
  
}

#' Internal function for a step in the fine-tuning process
#'
#' @param topLabels list of the cell type to use in the fine-tuning round
#' @param types list of all cell types
#' @param ref_data the reference expression data
#' @param genes list of genes to use for the analysis, or 'de' or 'sd' 
#' @param mat the median expression matrix for the reference
#' @param sd.thres the standard deviation threshold to use for choosing genes. Only used if genes == 'sd'
#' @param sc_data the single cell data
#' @param quantile.use the quantile for the aggregating scores
#' @param fine.tune.thres a threshold of top cell types to use for the next round.
#'
#' @return the cell types to use for the next round or one cell type which is the final annotation.
fineTuningRound = function(topLabels,types,ref_data,genes,mat,sd.thres,
                           sc_data,quantile.use,fine.tune.thres) {
  labels.use = is.element(types,topLabels)
  ref_data.filtered = as.matrix(ref_data[,labels.use])
  types.filtered = types[labels.use]
  if (typeof(genes)=='list') {
    n = round(1000*(2/3)^(log2(c(ncol(mat)))))
    utypes = colnames(mat)
    genes.filtered = unique(unlist(unlist(lapply(utypes,function(j) 
      lapply(utypes, function(i) genes[[i]][[j]][1:n])))))
    genes.filtered = intersect(genes.filtered,rownames(mat))
  } else if (genes[1] == "de") {
    n = round(500*(2/3)^(log2(c(ncol(mat)))))
    genes.filtered = unique(unlist(unlist(
      lapply(1:ncol(mat), function(j) {
        lapply(1:ncol(mat), function(i) {
          s=sort(mat[,j]-mat[,i],decreasing=T);
          s=s[s>0];names(s)[1:min(n,length(s))]})
      }))))[-1]
  } else if (genes[1] == "sd") {
    sd =  rowSds(mat)
    thres = min(sort(sd,decreasing = TRUE)[500],sd.thres)
    genes.filtered = intersect(rownames(ref_data)[sd>=thres],names(sc_data))
  } else {
    genes.filtered=intersect(genes,intersect(rownames(sc_data),
                                             (rownames(ref_data))))
  }
  
  ref_data.filtered = ref_data.filtered[genes.filtered,]
  sc_data.filtered = as.matrix(sc_data[genes.filtered])
  #data = sc_data.filtered[,i]
  if (sd(sc_data.filtered)>0) {
    r=cor(sc_data.filtered,ref_data.filtered,method='spearman')
    agg_scores = quantileMatrix(r,types.filtered,quantile.use);
    max_score = max(agg_scores)
    agg_scores = agg_scores[,-which.min(agg_scores)]
    topLabels = names(agg_scores)[agg_scores>=max_score-fine.tune.thres]
  } else {
    topLabels = topLabels[1]
  }
  topLabels
}

#' Scoring single cells using reference data set
#'
#' @param sc_data  the single-cell RNA-seq data set as a matrix with genes as rownames. If the data if from a full-length platform, counts must be normalized to gene length (TPM, RPKM, FPKM, etc.).
#' @param ref_data the reference dataset with genes as rownames. Gene names must be in the same format as the single cell data (if sc_data uses genes symbols, ref_data must have the same)
#' @param genes the list of genes to use.
#' @param types a list of cell type names corresponding to ref_data. Number of elements in types must be equal to number of columns in ref_data
#' @param quantile.use correlation coefficients are aggregated for multiple cell types in the reference data set. This parameter allows to choose how to sort the cell types scores, by median (0.5) or any other number between 0 and 1. The default is 0.9.
#'
#' @return a list with the scores, the raw correlation coefficients and the top labels
SingleR.ScoreData <- function(sc_data,ref_data,genes,types,quantile.use) {
  sc_data = as.matrix(sc_data[genes,])
  ref_data = as.matrix(ref_data[genes,])
  r=cor(sc_data,ref_data,method='spearman')
  
  agg_scores = quantileMatrix(r,types,quantile.use);
  #agg_scores = aggregate(t(r)~types,FUN = quantile, probs  = quantile.use)
  labels = colnames(agg_scores)[max.col(agg_scores)]
  output = list()
  
  
  if (dim(sc_data)[2]>1) {
    names(labels)=t(colnames(sc_data))
  }
  
  output$scores = as.matrix(t(agg_scores))
  
  output$labels = as.matrix(labels)
  output$r = r
  output$scores = t(output$scores)
  
  return(output)
}

#' The main SingleR function
#'
#' Given single-cell RNAseq data and reference dataset the function returns the best annotation for each single-cell.
#'
#' @param method annotating each single-cell or as a group by cluster: (\code{"single"} or \code{"cluster"})
#' @param sc_data the single-cell RNA-seq data set as a matrix with genes as rownames. If the data if from a full-length platform, counts must be normalized to gene length (TPM, RPKM, FPKM, etc.).
#' @param ref_data the reference dataset with genes as rownames. Gene names must be in the same format as the single cell data (if sc_data uses genes symbols, ref_data must have the same)
#' @param types a list of cell type names corresponding to ref_data. Number of elements in types must be equal to number of columns in ref_data
#' @param clusters only if using the "cluster" method. Please provide grouping variables as a factor. The number of elements of clusters must be equal to the number of columns in sc_data
#' @param genes list of genes to use for the annotations, or a method for extracting the genes from the data. Available methods are standard deviation and gene dispersion (\code{"sd"} or \code{"de"}). default is "de".
#' @param quantile.use correlation coefficients are aggregated for multiple cell types in the reference data set. This parameter allows to choose how to sort the cell types scores, by median (0.5) or any other number between 0 and 1. The default is 0.9.
#' @param p.threshold Chi-square outlier detection is used to assess the significance power of the top correlation. Single-cell with an annotation of p-value > p.threshold are designated as "X". Only applies for non fine-tuned annotations.
#' @param fine.tune perform the fine tuning step? default is TRUE.
#' @param fine.tune.thres the fine tuning step performs the scoring procedure for the top scoring labels using only genes that vary between those cell types in the reference data set. The top labels are those with a score lower than the top score by less than fine.tune.thres
#' @param sd.thres if genes=='sd' then this is the threshold for defining a variable gene.
#' @param do.pvals compute chi-squared outlier test p-values
#'
#' @return a list with the labels and scores
SingleR <- function(method = "single", sc_data, ref_data, types, 
                    clusters = NULL, genes = "de", quantile.use = 0.8, 
                    p.threshold = 0.05, fine.tune = TRUE, 
                    fine.tune.thres = 0.05,sd.thres=1, do.pvals = T) {
  rownames(ref_data) = tolower(rownames(ref_data))
  rownames(sc_data) = tolower(rownames(sc_data))
  A = intersect(rownames(ref_data),rownames(sc_data))
  sc_data = as.matrix(sc_data[A,])
  ref_data = ref_data[A,]
  if (ncol(sc_data)>1) {
    not.use = rowSums(is.na(ref_data))>0 | rowSums(is.na(sc_data))>0 | 
      rowSums(ref_data)==0
    ref_data = ref_data[!not.use,]
    sc_data = sc_data[!not.use,]
  }
  
  mat = medianMatrix(ref_data,types)
  
  if (typeof(genes)=='list') {
    utypes = unique(types)
    n = round(1000*(2/3)^(log2(c(ncol(mat)))))
    genes.filtered = unique(unlist(unlist(lapply(utypes,function(j) 
      lapply(utypes, function(i) genes[[i]][[j]][1:n])))))
    genes.filtered = intersect(genes.filtered,rownames(mat))
    print(paste0("Number of DE genes:", length(genes.filtered)))
  } else if (genes[1] == "de") {
    n = round(500*(2/3)^(log2(c(ncol(mat)))))
    genes.filtered = unique(unlist(unlist(lapply(1:ncol(mat), function(j) {
      lapply(1:ncol(mat), function(i) {
        s=sort(mat[,j]-mat[,i],decreasing=T);
        s=s[s>0];
        names(s)[1:min(n,length(s))]
      })}))))[-1]
    print(paste0("Number of DE genes:", length(genes.filtered)))
  } else if (genes[1] == "sd") {
    sd =  rowSds(as.matrix(mat))
    genes.filtered=intersect(rownames(mat)[sd>sd.thres],rownames(sc_data))
    print(paste0("Number of genes with SD>",sd.thres,": ",
                 length(genes.filtered)))
  } else {
    print(paste("Number of genes using in analysis:",length(genes.filtered)))
    genes.filtered=intersect(genes,intersect(rownames(sc_data),
                                             (rownames(ref_data))))
  }
  
  cell.names = colnames(sc_data)
  
  if (method == "single") {
    print(paste("Number of cells:",dim(sc_data)[2]))
  } else if (method == "cluster") {
    n = length(levels(clusters))
    print(paste("Number of clusters:",n))
    data = matrix(nrow=dim(sc_data)[1],ncol=n)
    for (i in 1:n) {
      data[,i] = rowSums(as.matrix(sc_data[,is.element(clusters,
                                                       levels(clusters)[i])]))
    }
    colnames(data) = levels(clusters)
    rownames(data) = rownames(sc_data)
    sc_data = data
  } else {
    print ("Error: method must be 'single' or 'cluster'")
    return(0)
  }
  output = SingleR.ScoreData(sc_data,ref_data,genes.filtered,types,quantile.use)
  if (do.pvals == T)
    output$pval = apply(output$scores, 1,function(x) chisq.out.test(x)$p.value)
  
  # second round with top labels
  if (fine.tune==TRUE & length(unique(types)) > 2) {
    labels = SingleR.FineTune(sc_data,ref_data,types,output$scores,
                              quantile.use,fine.tune.thres,genes = genes,
                              sd.thres,mat)
    output$labels1 = as.matrix(output$labels)
    output$labels = as.matrix(labels)
    output$labels1.thres = c(output$labels)
    if (do.pvals == T)
      output$labels1.thres[output$pval>p.threshold] = "X"
  } else {
    labels = as.matrix(output$labels)
    output$labels.thres = c(output$labels)
    if (do.pvals == T)
      output$labels.thres[output$pval>p.threshold] = "X"
  }
  output$cell.names = cell.names
  output$quantile.use = quantile.use
  output$types = types
  output$method = method
  
  return (output)
}

#' Plot a scatter plot of a single cell vs. a reference sample
#'
#' @param sc_data the single-cell RNA-seq data set as a matrix with genes as rownames.
#' @param cell_id a number of the single cell to use
#' @param ref the reference dataset with genes as rownames. Gene names must be in the same format as the single cell data (if sc_data uses genes symbols, ref_data must have the same)
#' @param sample_id a number of the sample to use
#'
#' @return a ggplot
SingleR.DrawScatter = function(sc_data, cell_id, ref,sample_id) {
  rownames(sc_data) = tolower(rownames(sc_data))
  rownames(ref$data) = tolower(rownames(ref$data))
  A = intersect(rownames(sc_data),rownames(ref$data))
  df = data.frame(sc_data[A,cell_id],ref$data[A,sample_id])
  colnames(df) = c('x','y')
  ggplot(df,aes(x=x, y=y)) + geom_point(size=0.5,alpha=0.5,color='blue') +
    geom_smooth(method='lm',color='red')+
    theme(legend.position="none") + xlab('Single cell') + ylab('Reference sample') +
    ggtitle(paste('R =', round(1000*cor(df$x,df$y,method='spearman'))/1000)) + 
    theme_classic()
}

#' Plot boxplots for each label for a given single cell.
#'
#' @param sc_data the single-cell RNA-seq data set as a matrix with genes as rownames.
#' @param cell_id a number of the single cell to use
#' @param ref the reference dataset with genes as rownames. Gene names must be in the same format as the single cell data (if sc_data uses genes symbols, ref_data must have the same)
#' @param labels.use a list of labels to use. If NULL uses all labels.
#' @param quantile.order same a quantile.use - by which percentile to order the boxplots.
#' @param main_types aggregate labels by main types or by all types
#' @param top.n number of boxplots to present (starting from top)
#' @param tit title for the boxplot
#' @param colors colors to use. Default is singler.colors
#'
#' @return a list with a ggplot and the scores for the single cell
SingleR.DrawBoxPlot = function(sc_data, cell_id, ref, labels.use=NULL, 
                               quantile.order = 0.8, main_types=F, top.n=50, tit = NULL, colors=singler.colors) {
  names(singler.colors) = levels(factor(ref$main_types))
  main_colors = singler.colors[levels(factor(ref$main_types))]
  sub_colors = singler.colors[unique(cbind(ref$main_types,ref$types))[,1]]
  names(sub_colors) = unique(cbind(ref$main_types,ref$types))[,2]
  
  if (main_types==T) {
    types = ref$main_types
    sub_colors = main_colors
  } else {
    types = ref$types
  }
  
  if (!is.null(labels.use)) {
    types.use = types %in% labels.use
    types = types[types.use]
  } else {
    types.use = rep(TRUE,length(types))
  }
  
  res = SingleR(sc_data=as.matrix(sc_data[,c(cell_id,cell_id)]),
                ref_data = ref$data[,types.use],types=types,fine.tune=F,
                sd.thres=ref$sd.thres,do.pvals=F)
  if (is.null(top.n)) {
    top.n = length(unique(types))
  }
  A = types %in% names(sort(res$scores[1,],decreasing = T)[1:top.n])
  df = data.frame(Spearman=res$r[1,A],Types=types[A])
  fac <- with(df, reorder(Types, Spearman, 
                          function(x) quantile(x, probs  = quantile.order), 
                          order = TRUE))
  df$Types <- factor(df$Types, levels = levels(fac))
  
  p = ggplot(df, aes(x = Types, y = Spearman,  color = Types)) +
    geom_boxplot(alpha = 0.2)+theme_classic() +
    geom_point(alpha = 0.5, position = "jitter",shape=16) +
    xlab('') + ggtitle(paste(colnames(sc_data)[cell_id])) +
    
    theme(legend.position="none"
          , axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5)
          , panel.background = element_rect(fill = "transparent") # bg of the panel
          , plot.background = element_rect(fill = "transparent") # bg of the plot
          , legend.background = element_rect(fill = "transparent") # get rid of legend bg
          , legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
    )
  
  list(plot=p,res=res)
  
}

#' Plot a heatmap of the scores for all the single cells
#'
#' @param SingleR the output from the SingleR function
#' @param cells.use single cells to present, if NULL all single cells presented
#' @param types.use cell types to present, if NULL all cell types presented
#' @param clusters a clustering to present as annotation in the heatmap
#' @param top.n number of cell types to presents. Default is 40. This can have an effect on the clustering which is performed only on the cell types presented.
#' @param normalize if TRUE scores are normalized to a 0-1 scale.
#' @param order.by.clusters if TRUE columns are ordered by the input clusters, and are not clustered again
#' @param cells_order an input order for the column
#' @param silent if TRUE do not draw the plot  
SingleR.DrawHeatmap = function(SingleR,cells.use = NULL, types.use = NULL,
                               clusters=NULL,top.n=40,normalize=T,
                               order.by.clusters=F,cells_order=NULL,silent=F,
                               ...) {
  scores = SingleR$scores
  if (!is.null(cells.use)) {
    scores = scores[cells.use,]
  }
  if (!is.null(types.use)) {
    scores = scores[,types.use]
  }
  
  m = apply(t(scale(t(scores))),2,max)
  
  thres = sort(m,decreasing=TRUE)[min(top.n,length(m))]
  
  data = as.matrix(scores)
  
  if (normalize==T) {
    mmax = rowMaxs(data)
    mmin = rowMins(data)
    data = (data-mmin)/(mmax-mmin)
    data = data^3
  }
  data = data[,m>(thres-1e-6)]
  
  
  data = t(data)
  
  if (!is.null(clusters)) {
    clusters = as.data.frame(clusters)
    colnames(clusters) = 'Clusters'
    rownames(clusters) = colnames(data)
    
  }
  additional_params = list(...)
  if (is.null(additional_params$annotation_colors)) {
    annotation_colors = NA
  } else {
    annotation_colors = additional_params$annotation_colors
  }
  clustering_method = 'ward.D2'
  if (order.by.clusters==T) {
    data = data[,order(clusters$Clusters)]
    clusters = clusters[order(clusters$Clusters),,drop=F]
    pheatmap(data,border_color=NA,show_colnames=FALSE,
             clustering_method=clustering_method,fontsize_row=9,
             annotation_col = clusters,cluster_cols = F,silent=silent, 
             annotation_colors=annotation_colors)
  } else if (!is.null(cells_order)) {
    data = data[,cells_order]
    clusters = clusters[cells_order,,drop=F]
    pheatmap(data,border_color=NA,show_colnames=FALSE,
             clustering_method=clustering_method,fontsize_row=9,
             annotation_col = clusters,cluster_cols = F,silent=silent, 
             annotation_colors=annotation_colors)
  } else {
    if (!is.null(clusters)) {
      pheatmap(data,border_color=NA,show_colnames=FALSE,
               clustering_method=clustering_method,fontsize_row=9,
               annotation_col = clusters,silent=silent, 
               annotation_colors=annotation_colors)
    } else {
      pheatmap(data[,sample(ncol(data))],border_color=NA,show_colnames=FALSE,
               clustering_method=clustering_method,fontsize_row=9,
               silent=silent, annotation_colors=annotation_colors)
      
    }
  }
}

#' Plot a colored tSNE plot according to a set of labels
#'
#' @param SingleR the output from the SingleR function
#' @param xy a matrix with the coordinates of the single cells
#' @param labels labels for the single cells, Default is the SingleR labels
#' @param score.thres single-cells with a score lower than this number are marked as 'X'.
#' @param clusters if the SingleR method is 'clusters' then this vector is the cluster id for each single cell
#' @param do.letters if TRUE shows letters the first letter of the annotation
#' @param dot.size the size of the dot. Default is 1
#' @param do.labels only relevant in a 'cluster' method. If TRUE The cluster annotation is shown on top of the cluster.
#' @param do.legend if TRUE legend is presented.
#' @param label.size if do.labels is TRUE then this determines the font size of the labels
#' @param title title for the plot
#' @param colors colors to use
#' @param font.size size of fonts in the plot
#' @param alpha an alpha for the transparency of the dots
#'
#' @return a list with ggplot and a data frame with the coordinates and annotations
SingleR.PlotTsne = function(SingleR, xy, labels=SingleR$labels, score.thres=0, 
                            clusters = NULL, do.letters = TRUE, dot.size = 1, 
                            do.labels = FALSE, do.legend = TRUE, label.size=3, 
                            title = "",colors=singler.colors,font.size=NULL,
                            alpha=0.5) {
  if (do.labels == TRUE)
    do.letters = FALSE
  
  df = data.frame(row.names = SingleR$cell.names)
  df$x = xy[,1]
  df$y = xy[,2]
  
  if (SingleR$method == "cluster") {
    df$ident = clusters.map.values(clusters,labels)
  } else {
    df$ident = labels
  }
  
  if (score.thres>0) {
    max.score = apply(SingleR$scores,1,max)  
    df$ident[max.score<score.thres] = 'X'
  }
  df$ident = factor(df$ident)
  SYMBOLS = c(LETTERS,tolower(letters),c(0:9))
  df$initIdent = SYMBOLS[as.numeric(df$ident)]
  
  num.levels = length(levels(df$ident))
  
  p = ggplot(df, aes(x = x, y = y))
  
  p = p + geom_point(aes(color=ident), size=dot.size,alpha=alpha,stroke=0)
  
  if( do.letters == TRUE) {
    symbols = SYMBOLS[1:num.levels]
    names(symbols) =   lev = levels(df$ident)
    p = p + geom_point(aes(shape=ident), size=2*dot.size/5,color='black')
    p = p + scale_shape_manual(values=symbols)
    #p = p + geom_point(aes(shape=initIdent), size=dot.size/2)
    #p = p + scale_shape_identity()
  }
  
  if (do.labels == TRUE) {
    df %>% dplyr::group_by(ident) %>% summarize(x = median(x), 
                                                y = median(y)) -> centers
    p = p + geom_point(data = centers, aes(x=x, y=y), size=0, alpha=0) + 
      geom_text(data=centers, aes(label=ident), size = label.size,color='black')
    p = p + guides(colour=guide_legend(override.aes = list(size=3,alpha=1)))
    x.range = layer_scales(p)$x$range$range
    add_to_x = sum(abs(x.range))*0.03
    p = p + xlim(x.range[1]-add_to_x,x.range[2]+add_to_x)
  } else {
    if (is.null(font.size)) {
      font.size = 250*(1/num.levels)
      font.size = max(font.size,5)
      font.size = min(font.size,10)
    }
    if (num.levels>35 & num.levels<60) {
      p = p + theme(legend.position="bottom",legend.direction="vertical",
                    legend.text=element_text(size=6),legend.title = element_blank()) + guides(col=guide_legend(ncol=5,override.aes = list(size=2,alpha=1)))
    } else if (num.levels>60){
      p = p + theme(legend.position="bottom",legend.direction="vertical",
                    legend.text=element_text(size=6),legend.title = element_blank()) + guides(col=guide_legend(ncol=9,override.aes = list(size=2,alpha=1)))
    } else {
      p = p + theme(legend.text=element_text(size=font.size),legend.title = 
                      element_blank())+ 
        guides(color=guide_legend(ncol=1,override.aes = list(size=3,alpha=1)))
    }
  }
  
  lev = levels(df$ident)
  cols = colors[1:length(lev)]
  names(cols) = lev
  cols[names(cols)=='X']='black'
  
  p = p + scale_color_manual(values = cols)
  p = p + xlab("tSNE 1") + ylab("tSNE 2") + ggtitle(title)
  
  if (do.legend==FALSE) {
    p = p + theme(legend.position="none")
  }
  p = p + theme(panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                panel.background = element_blank(), 
                axis.line = element_line(colour = "black"))
  
  out = list(p=p,df=df,num.levels=num.levels)
  
}

#' Plot a feature on the t-SNE plot
#'
#' @param SingleR the output from the SingleR function
#' @param seurat the seurat object
#' @param plot.feature the feature to plot - 'MaxScore' for coloring according to the top score per cell, 'nGene' for number of non-zero genes, 'nUMI' for number of UMIs, or a gene name.
#' @param dot.size size of the dot in the plot.
#'
#' @return ggplot2 object
SingleR.PlotFeature = function(SingleR,seurat, plot.feature='MaxScore', 
                               dot.size=1,title=NULL) {
  df = data.frame(row.names = rownames(seurat@cell.names),
                  seurat@dr$tsne@cell.embeddings)
  if (length(plot.feature)==nrow(df)) {
    df$Feature=plot.feature
    tit = 'Feature'
  } else if (plot.feature=='MaxScore') {
    df$Feature = apply(SingleR$scores,1,max) 
    tit = 'Max Score'
  } else if (plot.feature=='nGene') {
    df$Feature = seurat@meta.data$nGene
    tit = 'nGene'
  } else if (plot.feature=='nUMI') {
    df$Feature = seurat@meta.data$nUMI
    tit = 'nUMI'
  } else {
    df$Feature = seurat@data[plot.feature,]
    tit = plot.feature
  }
  if (is.null(title)) {
    title = tit
  }
  ggplot(df,aes(x=tSNE_1,y=tSNE_2)) + 
    geom_point(aes(color=Feature), size=dot.size)+
    scale_colour_gradient(low='gray',high='blue')+
    ggtitle(title) + theme_classic()
}

#' Calculate single-sample gene set enrichment (ssGSEA) for each single cell
#'
#' @param sc_data  the single-cell RNA-seq data set as a matrix with genes as rownames. If the data if from a full-length platform, counts must be normalized to gene length (TPM, RPKM, FPKM, etc.).
#' @param species (\code{"Mouse"} or \code{"Human"})
#' @param signatures a GeneSetCollection object, or NULL to use default signatures
#' @param n.break run ssGSEA for n.break at a time.

#'
#' @return scores for each signature and single cell
calculateSignatures = function(sc_data,species='Human',signatures=NULL, 
                               n.break=1000) {
  #data('signatures')
  if(is.null(signatures)) {
    if (species=="Human") {
      egc = human.egc
    } else if (species=="Mouse") {
      egc = mouse.egc
    }
  } else {
    egc = signatures
  }
  
  sc_data = as.matrix(sc_data)
  rownames(sc_data) = tolower(rownames(sc_data))
  numClusters = min(detectCores(all.tests = FALSE, logical = TRUE)-1,4)
  # break to groups of n.break cells
  scores = matrix(NA,length(egc),ncol(sc_data))
  wind = seq(1,ncol(sc_data),by=n.break)
  print(paste('Using sets of',n.break, 'cells. Running',length(wind),'times.'))
  for (i in wind) {
    last = min(ncol(sc_data),i+n.break-1)
    a = gsva(sc_data[,i:last],egc,method='ssgsea',ssgsea.norm=F,
             parallel.sz=numClusters,parallel.type='FORK')
    scores[,i:last] = a
  }
  mmin = rowMins(scores)
  mmax = rowMaxs(scores)
  scores = scores/(mmax-mmin)
  rownames(scores) = rownames(a)
  #rownames(scores) = c('G1/S','G2/M')
  output = data.frame(t(scores))
  
  if (is.null(signatures) && species=="Human") {
    output$Cell_Cycle = rowMeans(t(scores[c('G1S','G2M'),]))
    output[,c('G1S','G2M')] = c()
  } else {
    colnames(output) = c('G1/S','G2/M','M','M/G1','S')
  }
  output
}

#' Map clusters to values
#'
#' @param cluster_ids the cluster ids of each single cell
#' @param singler_clusters the cluster annotation for each cluster
#'
#' @return annotation for each single cell
clusters.map.values = function(cluster_ids,singler_clusters) {
  singler_clusters = as.matrix(singler_clusters)
  cluster_ids = as.matrix(cluster_ids)
  singler_clusters = as.matrix(singler_clusters[sort.int(as.numeric(
    rownames(singler_clusters)),index.return=TRUE)$ix,])
  A = rownames(singler_clusters) %in% unique(cluster_ids)
  map = paste0(rownames(singler_clusters)[A],": ",singler_clusters[A])
  names(map) = sort(as.numeric(unique(cluster_ids)))
  clusters = recode(cluster_ids,!!!map)
}

#' SingleR.Cluster - using the scores data to cluster single cells.
#'
#' @param SingleR the output from the SingleR function
#' @param num.clusts number of clusters to output
#' @param normalize_rows if TRUE the rows of the SingleR score matrix (single-cells) are scaled. Default is TRUE. 
#' @param normalize_cols if TRUE the columns of the SingleR score matrix (cell types) are scaled. Default is FALSE. 
#'
#' @return cluster id for each single cell
SingleR.Cluster = function(SingleR,num.clusts=10,normalize_rows=F,
                           normalize_cols=T) {
  if (normalize_rows==T) {
    SingleR$scores = scale(SingleR$scores)
  }
  if (normalize_cols==T) {
    scores = t(scale(t(SingleR$scores^5)))
  } else {
    scores = SingleR$scores
  }
  r <- cor(t(scores), method="pearson")
  d <- as.dist(1-r)
  hc = hclust(d,method='ward.D2')
  #hc = hclust(dist(scores,method='euclidean'),method='ward.D2')
  cl = cutree(hc,k=num.clusts)
  list(hc=hc,cl=factor(cl))
}

#' Subseting a SingleR object. This function subsets all the SingleR data and the Seurat object if included.
#'
#' @param singler as SingleR object
#' @param subsetdata a logical vector of single-cells to include in the subset object
#'
#' @return a subset of the original SingleR vector
SingleR.Subset = function(singler,subsetdata) {
  s = singler
  
  if (!is.null(s$seurat)) {
    s$seurat = SubsetData(s$seurat,colnames(s$seurat@data)[subsetdata])
    subsetdata = unlist(lapply(s$seurat@cell.names,FUN=function(x) 
      which(singler$singler[[1]]$SingleR.single$cell.names==x)))
  }
  
  for (i in 1:length(s$singler)) {
    s$singler[[i]]$SingleR.single$cell.names = 
      s$singler[[i]]$SingleR.single$cell.names[subsetdata]
    s$singler[[i]]$SingleR.clusters$cell.names = 
      s$singler[[i]]$SingleR.clusters$cell.names[subsetdata]
    s$singler[[i]]$SingleR.single$scores = 
      s$singler[[i]]$SingleR.single$scores[subsetdata,]
    s$singler[[i]]$SingleR.single$labels = 
      as.matrix(s$singler[[i]]$SingleR.single$labels[subsetdata,])
    s$singler[[i]]$SingleR.single$labels1 = 
      as.matrix(s$singler[[i]]$SingleR.single$labels1[subsetdata,])
    s$singler[[i]]$SingleR.single$clusters$cl = 
      s$singler[[i]]$SingleR.single$clusters$cl[subsetdata]
    
    if(!is.null(s$singler[[i]]$SingleR.single.main)) {
      s$singler[[i]]$SingleR.single.main$cell.names = 
        s$singler[[i]]$SingleR.single.main$cell.names[subsetdata]
      s$singler[[i]]$SingleR.clusters.main$cell.names = 
        s$singler[[i]]$SingleR.clusters.main$cell.names[subsetdata]
      s$singler[[i]]$SingleR.single.main$scores = 
        s$singler[[i]]$SingleR.single.main$scores[subsetdata,]
      s$singler[[i]]$SingleR.single.main$labels = 
        as.matrix(s$singler[[i]]$SingleR.single.main$labels[subsetdata,])
      s$singler[[i]]$SingleR.single.main$labels1 = 
        as.matrix(s$singler[[i]]$SingleR.single.main$labels1[subsetdata,])
      s$singler[[i]]$SingleR.single.main$clusters$cl = 
        s$singler[[i]]$SingleR.single.main$clusters$cl[subsetdata]
      
    }
  }
  if (!is.null(s[["signatures"]])) {
    s$signatures = s$signatures[subsetdata,]
  }
  if(!is.null(s[['other']])) {
    s$other = s$other[subsetdata]
  }
  
  if (!is.null(s$meta.data)) {
    s$meta.data$orig.ident = factor(as.character(
      s$meta.data$orig.ident[subsetdata]))
    s$meta.data$xy = s$meta.data$xy[subsetdata,]
    s$meta.data$clusters = factor(as.character(
      s$meta.data$clusters[subsetdata]))
  }
  s
}

#' Remove data from a SingleR object to make it smaller
#'
#' @param singler.data a SingleR object
#'
#' @return a smaller SingleR object
remove.Unnecessary.Data.single = function(singler.data) {
  for (j in 1:length(singler.data$singler)) {
    singler.data$singler[[j]]$SingleR.single = 
      singler.data$singler[[j]]$SingleR.single[
        -c(which(names(singler.data$singler[[j]]$SingleR.single) 
                 %in% c('r','labels1.thres','types')))]
    singler.data$singler[[j]]$SingleR.clusters = 
      singler.data$singler[[j]]$SingleR.clusters[-c(
        which(names(singler.data$singler[[j]]$SingleR.clusters) 
              %in% c('r','labels1.thres','types')))]
    if (!is.null(singler.data$singler[[j]]$SingleR.single.main)) {
      singler.data$singler[[j]]$SingleR.single.main = 
        singler.data$singler[[j]]$SingleR.single.main[-c(
          which(names(singler.data$singler[[j]]$SingleR.single.main) 
                %in% c('r','labels1.thres','types')))]
      singler.data$singler[[j]]$SingleR.clusters.main = 
        singler.data$singler[[j]]$SingleR.clusters.main[-c(
          which(names(singler.data$singler[[j]]$SingleR.clusters.main) 
                %in% c('r','labels1.thres','types')))]
    }
  }
  singler.data
}



#' Normalize counts to gene length
#'
#' @param counts counts data
#' @param lengths genes' lengths. If null use a precalculated gene lengths
#'
#' @return TPM values
TPM = function(counts,lengths=NULL) {
  #if (is.null(lengths)) {
  #data('gene_lengths')
  #}
  rownames(counts) = tolower(rownames(counts))
  names(lengths) = tolower(names(lengths))
  
  A = intersect(rownames(counts),names(lengths))
  counts = counts[A,]
  lengths = lengths[A]
  rate = counts / lengths
  apply(rate,2,function(x) 1e6*x/sum(x))
}

capitalize <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

#' Create a list for differential genes for all pairwise cell types in the a reference data set.
#' Using this significantly reduces computation time for the 'de' variable gene method.
#'
#' @param ref_data the reference dataset gene expression matrix
#' @param types labels for each column in ref_data
#' @param n number of genes per pairwise comparison
#'
#' @return list of lists. Each list contains a list for each cell types with n differential genes.
CreateVariableGeneSet = function(ref_data,types,n) {
  mat = medianMatrix(ref_data,types)
  genes = lapply(1:ncol(mat), function(j) {
    lapply(1:ncol(mat), function(i) {
      s=sort(mat[,j]-mat[,i],decreasing=T)
      s=s[s>0]
      tolower(names(s)[1:min(n,length(s))])
    })
  } )
  names(genes) = colnames(mat)
  for (i in 1:length(genes)) {
    names(genes[[i]]) = colnames(mat)
  }
  genes
}

#' Wrapper function to create a SingleR object
#'
#' @param sc.data a matrix of single cell expression data
#' @param ref a reference set object. 
#' This object must be a list containing: data - log2 normalized expression data;
#' types - annotations for each of the samples; 
#' main_types - annotations for each of the samples, but less detailed; 
#' name - name for the reference set; 
#' sd.thres - a threshold for sd (used in 'sd' mode); 
#' de.genes - lists of lists of differentially expressed genes. Can be created using the CreateVariableGeneSet function.
#' de.genes.main - lists of lists of differentially expressed genes. Can be created using the CreateVariableGeneSet function.
#' @param clusters a numeric vector of cluster ids for each single cell. If NULL uses SingleR clustering.
#' @param species The species of the sample ('Human' or 'Mouse').
#' @param citation a citation for the project.
#' @param technology The technology used for creating the single-cell data.
#' @param variable.genes variable gene method to use - 'sd' or 'de'. Default is 'de'.
#' @param fine.tune perform fine tuning. Default is TRUE. Fine-tuning may take long to run.
#' @param do.main.types if TRUE runs a main cell type annotation using the main_types annotation.
#'
#' @return a SingleR object object
SingleR.CreateObject <- function(sc.data,ref,clusters=NULL,species='Human',
                                 citation='-',technology='-',variable.genes='sd',
                                 fine.tune=T,do.main.types=T) {
  types = ref$types
  
  print(paste0('Annotating data with ',ref$name,'...'))
  
  print(paste('Variable genes method:',variable.genes))
  
  if (variable.genes=='de') {
    if (!is.null(ref$de.genes)) {
      variable.genes = ref$de.genes
      variable.genes.main = ref$de.genes.main
    } else {
      variable.genes = CreateVariableGeneSet(ref$data,ref$types,200)
      variable.genes.main = CreateVariableGeneSet(ref$data,ref$main_types,300)
    }
  } else {
    variable.genes.main = variable.genes
  }
  
  SingleR.single = SingleR("single",sc.data,ref$data,types=types,
                           sd.thres = ref$sd.thres,genes = variable.genes,
                           fine.tune = fine.tune)
  
  SingleR.single$clusters = SingleR.Cluster(SingleR.single,10)
  
  if (is.null(clusters)) {
    clusters = SingleR.single$clusters$cl
  }
  
  SingleR.clusters = SingleR("cluster",sc.data,ref$data,types=types, 
                             clusters = factor(clusters),
                             sd.thres = ref$sd.thres,
                             genes = variable.genes,
                             fine.tune = fine.tune)
  
  about = list(Organism = capitalize(species),Citation=citation,
               Technology = technology,RefData=ref$name)
  
  
  singler = list(SingleR.single = SingleR.single, 
                 SingleR.clusters = SingleR.clusters,about=about)
  
  if (do.main.types==T) {
    print(paste0('Annotating data with ',ref$name,' (Main types)...'))
    types = ref$main_types
    singler$SingleR.single.main = SingleR("single",sc.data,ref$data,
                                          types=types,sd.thres = ref$sd.thres, 
                                          quantile.use = 0.8, 
                                          genes = variable.genes.main,
                                          fine.tune = fine.tune)
    singler$SingleR.single.main$clusters = 
      SingleR.Cluster(singler$SingleR.single.main,10)
    singler$SingleR.clusters.main = 
      SingleR("cluster",sc.data,ref$data,types=types, 
              clusters=factor(clusters),sd.thres = ref$sd.thres, 
              quantile.use = 0.8,genes = variable.genes.main,
              fine.tune = fine.tune)
  }
  
  if (!(ref$name %in% c('Immgen','RNAseq','HPCA','Blueprint_Encode',
                        'Fantom','GSE43005'))) {
    singler$about$refernce = ref
  }
  
  singler
}

#' Wrapper function to create a Seurat object
#'
#' @param project.name the project name
#' @param sc.data a matrix of single cell expression data
#' @param min.genes Include cells where at least this many genes are detected.
#' @param min.cells include genes with detected expression in at least this many cells. Will subset the raw.data matrix as well. To reintroduce excluded genes, create a new object with a lower cutoff.
#' @param regress.out variables to regress out (previously latent.vars in RegressOut). For example, nUMI, or percent.mito.
#' @param npca a vector of the dimensions to use in construction of the clustering and tSNE plot
#' @param resolution clustering resolution. See Seurat manual for more details.
#' @param temp.dir used by FindClusters function.
#'
#' @return a Seurat object
SingleR.CreateSeurat <- function(project.name,sc.data,min.genes = 500,
                                 min.cells = 2,regress.out = 'nUMI',
                                 npca = 10,resolution=0.8,temp.dir=NULL) {
  sc = CreateSeuratObject(raw.data = sc.data, min.cells = min.cells, 
                          min.genes = min.genes, project = project.name)
  mtgenes = '^mt-'
  mito.genes <- grep(pattern = mtgenes, x = rownames(x = sc@data), 
                     value = TRUE,ignore.case=TRUE)
  percent.mito <- colSums((sc.data[mito.genes, ]))/colSums(sc.data)
  sc <- AddMetaData(object = sc, metadata = percent.mito, 
                    col.name = "percent.mito")
  
  sc <- NormalizeData(object = sc, 
                      normalization.method = "LogNormalize", 
                      scale.factor = 10000)
  sc <- FindVariableGenes(object = sc, mean.function = ExpMean, 
                          dispersion.function = LogVMR, 
                          x.low.cutoff = 0.0125, x.high.cutoff = 3, 
                          y.cutoff = 0.5, do.contour = F, do.plot = F)
  
  sc <- ScaleData(object = sc, vars.to.regress = regress.out)
  
  sc <- RunPCA(object = sc, pc.genes = sc@var.genes, do.print = FALSE)
  
  PCElbowPlot(object = sc)
  
  sc <- FindClusters(object = sc, reduction.type = "pca", 
                     dims.use = 1:npca,resolution = resolution, 
                     print.output = 0, save.SNN = F, 
                     temp.file.location = temp.dir)
  
  if (ncol(sc@data)<100) {
    sc <- RunTSNE(sc, dims.use = 1:npca, do.fast = T,perplexity=10  )
  } else {
    sc <- RunTSNE(sc, dims.use = 1:npca, do.fast = T)
    
  }
  
  sc
}

#' Run annotation based on Kang et al. Nature Biotechnology 2017.
#'
#' @param sc.data a matrix of single cell expression data.
#'
#' @return list of the correlation matrix and the annotations.
SingleR.CreateKangAnnotations = function(sc.data) {
  #if (!exists('cell.type.classification'))
  #  data('cell.type.cor.classification')
  rownames(cell.type.classification$cell.types.avg) = 
    tolower(rownames(cell.type.classification$cell.types.avg))
  A = intersect(rownames(sc.data),
                rownames(cell.type.classification$cell.types.avg))
  r = cor(cell.type.classification$cell.types.avg[A,],
          as.matrix(sc.data[A,]),method='spearman')
  kang_annotation = max.col(t(r))
  
  map = colnames(cell.type.classification$cell.types.avg)
  names(map) = seq(1,ncol(cell.type.classification$cell.types.avg))
  kang_annotation = recode(kang_annotation,!!!map)
  names(kang_annotation) = colnames(r)
  
  list(r=r,kang_annotation=kang_annotation)
}


#' Internal function to read single-cell data
#'
#' @param counts a tab delimited text file containing the counts matrix, a 10X directory name or a matrix with the counts.
#' @param annot a tab delimited text file or a data.frame. Rownames correspond to column names in the counts data
#'
#' @return list with sc.data - the counts matrix, orig.ident - a named vector of the originial identities.
ReadSingleCellData = function(counts,annot) {
  if (typeof(counts) == 'character') {
    if (file.info(counts)$isdir==T) {
      counts = as.matrix(Read10X(counts))
    } else if (file.info(counts)$isdir==F) {
      counts <- as.matrix(read.table(counts, header=TRUE, sep="\t", 
                                     row.names=1, as.is=TRUE,comment.char='!'))
    } else {
      stop('Cannot find file or directory.')
    }
  } 
  
  A = tolower(rownames(counts))
  dupA = duplicated(A)
  if (sum(dupA)>0) {
    counts = counts[!dupA,]
  }
  
  #  if (species=='Mouse') {
  #    rownames(counts) = paste(toupper(substr(rownames(counts), 1, 1)), tolower(substr(rownames(counts), 2, nchar(rownames(counts)))), sep="")
  #  }
  
  if (!is.null(annot)) {
    if (length(annot) == 1) {
      types <- read.table(annot, header=TRUE, sep="\t", row.names=1, as.is=TRUE)
      orig.ident = types[,1]
      names(orig.ident) = rownames(types)
    } else {
      orig.ident = annot
      if (is.null(names(orig.ident))) {
        names(orig.ident) = colnames(counts)
      }
    }
  } else {
    orig.ident = rep('NA',ncol(counts))
    names(orig.ident)=colnames(counts)
  }
  
  list(counts=counts,orig.ident=orig.ident)
}

#' Wrapper function to create a SingleR object + Seurat object
#'
#' @param counts a tab delimited text file containing the counts matrix, a 10X directory name or a matrix with the counts.
#' @param annot a tab delimited text file or a data.frame. Rownames correspond to column names in the counts data
#' @param project.name the project name
#' @param min.genes Include cells where at least this many genes are detected.
#' @param technology The technology used for creating the single-cell data.
#' @param species The species of the sample ('Human' or 'Mouse')
#' @param citation a citation for the project
#' @param ref.list a list of reference objects. If NULL uses the predefined reference objects - Mouse: ImmGen and Mouse.RNAseq, Human: HPCA and Blueprint+Encode. 
#' @param normalize.gene.length if a full-length method set to TRUE, if a 3' method set to FALSE.
#' @param variable.genes variable gene method to use - 'sd' or 'de'. Default is 'de'.
#' @param fine.tune perform fine tuning. Default is TRUE. Fine-tuning may take long to run.
#' @param reduce.file.size remove less used SingleR fields that increase the object size.
#' @param do.signatures create signatures data
#' @param min.cells include genes with detected expression in at least this many cells. Will subset the raw.data matrix as well. To reintroduce excluded genes, create a new object with a lower cutoff.
#' @param regress.out variables to regress out (previously latent.vars in RegressOut). For example, nUMI, or percent.mito.
#' @param npca a vector of the dimensions to use in construction of the clustering and tSNE plot
#' @param do.main.types run the SingleR pipeline for main cell types (cell types grouped together) as well.
#' @param reduce.seurat.object if TRUE removes the raw data and other high memory objects from the Seurat object.
#' @param temp.dir used by the SingleR web app.
#'
#' @return a SingleR object containing a Seurat object
CreateSinglerSeuratObject = function(counts,annot=NULL,project.name,
                                     min.genes=500,technology='10X',
                                     species='Human',citation='',
                                     ref.list=list(),normalize.gene.length=F,
                                     variable.genes='de',fine.tune=T,
                                     reduce.file.size=T,do.signatures=T,
                                     min.cells=2,npca=10,regress.out='nUMI',
                                     do.main.types=T,reduce.seurat.object=T,
                                     temp.dir=NULL) {
  print(project.name)
  print('Reading single-cell data...')
  sc.data = ReadSingleCellData(counts,annot)
  print('Create Seurat object...')
  seurat = SingleR.CreateSeurat(project.name,sc.data$counts,min.genes=
                                  min.genes,min.cells=min.cells,
                                regress.out=regress.out,npca=npca,
                                temp.dir=temp.dir)
  
  orig.ident = sc.data$orig.ident[colnames(seurat@data)]
  counts = sc.data$counts[,colnames(seurat@data)]
  
  seurat@meta.data$orig.ident = factor(orig.ident)
  
  clusters = seurat@ident
  
  if(reduce.seurat.object==T) {
    seurat@raw.data = c()
    seurat@scale.data = c()
    seurat@calc.params = list()
  }
  
  print('Creat SingleR object...')
  
  singler = CreateSinglerObject(counts,orig.ident,project.name,
                                min.genes=0,technology,species,
                                citation,ref.list,
                                normalize.gene.length,variable.genes,
                                fine.tune,do.signatures,
                                seurat@ident,do.main.types,
                                reduce.file.size,temp.dir)
  
  singler$seurat = seurat 
  singler$meta.data$xy = seurat@dr$tsne@cell.embeddings
  singler$meta.data$clusters = seurat@ident
  
  singler
}

#' Wrapper function to create a SingleR object
#'
#' @param counts a tab delimited text file containing the counts matrix, a 10X directory name or a matrix with the counts.
#' @param annot a tab delimited text file or a data.frame. Rownames correspond to column names in the counts data
#' @param project.name the project name
#' @param min.genes Include cells where at least this many genes are detected (number non-zero genes).
#' @param technology The technology used for creating the single-cell data.
#' @param species The species of the sample ('Human' or 'Mouse').
#' @param citation a citation for the project.
#' @param ref.list a list of reference objects. If NULL uses the predefined reference objects - Mouse: ImmGen and Mouse.RNAseq, Human: HPCA and Blueprint+Encode. 
#' @param normalize.gene.length if a full-length method set to TRUE, if a 3' method set to FALSE.
#' @param variable.genes variable gene method to use - 'sd' or 'de'. Default is 'de'.
#' @param fine.tune perform fine tuning. Default is TRUE. Fine-tuning may take long to run.
#' @param do.signatures create signatures data
#' @param clusters input cluster id for each of the cells with at least min.genes, if NULL uses SingleR clusterings.
#' @param do.main.types run the SingleR pipeline for main cell types (cell types grouped together) as well.
#' @param reduce.file.size remove less used SingleR fields that increase the object size.
#' @param temp.dir used by the SingleR webtool.
#'
#' @return a SingleR object
CreateSinglerObject = function(counts,annot=NULL,project.name,
                               min.genes=500,technology='10X',
                               species='Human',citation='',
                               ref.list=list(),normalize.gene.length=F,
                               variable.genes='de',fine.tune=T,
                               do.signatures=T,clusters=NULL,
                               do.main.types=T,reduce.file.size=T,
                               temp.dir=NULL) {
  
  sc.data = ReadSingleCellData(counts,annot)
  
  print(paste0('Dimensions of counts data: ',
               nrow(sc.data$counts),'x',ncol(sc.data$counts)))
  
  singler = list()
  
  
  N = colSums(counts>0)
  orig.ident = sc.data$orig.ident[N>=min.genes]
  counts = sc.data$counts[,N>=min.genes]
  
  if (normalize.gene.length == F) {
    sc.data.gl = counts
    rownames(sc.data.gl) = tolower(rownames(sc.data.gl))
  } else {
    if (species == 'Human') {
      sc.data.gl = TPM(counts,human_lengths)
    } else if (species == 'Mouse') {
      sc.data.gl = TPM(counts,mouse_lengths)
    }
  }
  
  if (length(ref.list)==0) {
    if (species == 'Mouse') {
      #if (!exists('immgen'))
      #  data('Immgen')
      #if (!exists('mouse.rnaseq'))
      #  data('Mouse-RNAseq')
      res = list(SingleR.CreateObject(sc.data.gl,immgen,clusters,species,
                                      citation,technology,
                                      do.main.types=do.main.types,
                                      variable.genes=variable.genes,
                                      fine.tune=fine.tune),
                 SingleR.CreateObject(sc.data.gl,mouse.rnaseq,clusters,
                                      species,citation,technology,
                                      do.main.types=do.main.types,
                                      variable.genes=variable.genes,
                                      fine.tune=fine.tune)
      )
    } else if (species == 'Human') {
      #if(!exists('hpca'))
      #  data ('HPCA')
      #if (!exists('blueprint_encode'))
      #  data('Blueprint_Encode')
      res = list(SingleR.CreateObject(sc.data.gl,hpca,clusters,species,
                                      citation,technology,
                                      do.main.types = do.main.types,
                                      variable.genes=variable.genes,
                                      fine.tune=fine.tune),
                 SingleR.CreateObject(sc.data.gl,blueprint_encode,
                                      clusters,species,citation,technology,
                                      do.main.types = do.main.types,
                                      variable.genes=variable.genes,
                                      fine.tune=fine.tune))
    }
  } else {
    res = lapply(ref.list, FUN=function(x) {
      SingleR.CreateObject(sc.data.gl,x,clusters,species,citation,technology,
                           do.main.types=do.main.types,
                           variable.genes=variable.genes,fine.tune=fine.tune)
    })
  }
  
  singler$singler = res
  
  if (do.signatures==TRUE) {
    signatures = calculateSignatures(sc.data.gl,species=species)
    singler$signatures = signatures
    
  }
  
  if (species == 'Human') {
    kang = SingleR.CreateKangAnnotations(sc.data.gl)
    singler$other = kang$kang_annotation
  }
  
  singler$meta.data = list(project.name=project.name,orig.ident=orig.ident)
  
  if (reduce.file.size==T) {
    singler = remove.Unnecessary.Data.single(singler)
  }
  
  singler
  
}

#' Helper function to combine multiple 10X datasets
#'
#' @param dirs a list of the directories of the 10X samples. dirs should be the absolute path. dirs = list.dirs(dir.path,full.names = T)
#' @param random.sample number of cells to choose randomly from each sample
#' @param min.genes a threshold on the number of non-zero genes.
#'
#' @return a list with sc.data - the count matrix, and orig.ident - the directory name of the sample
Combine.Multiple.10X.Datasets = function(dirs,random.sample=0,min.genes=500) {
  print(paste(1,basename(dirs[1])))
  sc.data = as.matrix(Read10X(dirs[1]))
  orig.ident = rep(basename(dirs[1]),ncol(sc.data))
  
  if (random.sample>0) {
    A = colSums(as.matrix(sc.data)>0)
    use = sample(which(A>min.genes),random.sample)
    sc.data = sc.data[,use]
    orig.ident = orig.ident[use]
  }
  
  for (j in 2:length(dirs)) {
    print(paste(j,basename(dirs[j])))
    data <- Read10X(dirs[j])
    ident = rep(basename(dirs[j]),ncol(data))
    if (random.sample>0) {
      A = colSums(as.matrix(data)>0)
      use = sample(which(A>min.genes),random.sample)
      data = data[,use]
      ident = ident[use]
    }
    
    orig.ident = c(orig.ident,ident)
    sc.data = cbind(sc.data,data)
  }
  
  colnames(sc.data) = make.unique(colnames(sc.data))
  names(orig.ident) = colnames(sc.data)
  
  list(sc.data=sc.data,orig.ident=orig.ident)
}

#' Combining SingleR objects together to one object. 
#' 
#' This function can be useful when the memory is insufficient to analyze the whole counts matrix. We can analyze subsets of the full data and combine them together afterwards.
#'
#' @param singler.list a list of SingleR objects
#' @param order a vector of the cells name to order the SingleR object by. Important if there is a pre-built single-cell object, and the order of cells in SingleR needs to be in the same order.
#' @param clusters a vector of clusters  per cell. Optional.
#' @param expr the counts matrix. Optional.
#' @param xy the t-SNE coordinates of all cells.
#'
#' @return a combined SingleR object
SingleR.Combine = function(singler.list,order=NULL,clusters=NULL,expr=NULL,
                           xy=NULL) {
  
  singler=c()
  singler$singler = singler.list[[1]]$singler
  for (j in 1:length(singler.list[[1]]$singler)) {
    singler$singler[[j]]$SingleR.cluster = c()
    singler$singler[[j]]$SingleR.cluster.main = c()
    singler$singler[[j]]$SingleR.single$clusters=c()
  }
  singler$meta.data = singler.list[[1]]$meta.data
  singler$meta.data$clusters=c()
  singler$meta.data$xy=c()
  singler$meta.data$data.sets = rep(singler$meta.data$project.name,
                                    length(singler$meta.data$orig.ident))
  
  for (i in 2:length(singler.list)) {
    for (j in 1:length(singler$singler)) {
      if (singler.list[[i]]$singler[[j]]$about$RefData!=
          singler.list[[1]]$singler[[j]]$about$RefData) {
        stop('The objects are not ordered by the same reference data.')
      }
      
      singler$singler[[j]]$about$Organism = 
        c(singler$singler[[j]]$about$Organism,singler.list[[i]]$singler[[j]]$about$Organism)
      singler$singler[[j]]$about$Citation = 
        c(singler$singler[[j]]$about$Citation,singler.list[[i]]$singler[[j]]$about$Citation)
      singler$singler[[j]]$about$Technology = 
        c(singler$singler[[j]]$about$Technology,singler.list[[i]]$singler[[j]]$about$Technology)
      #singler$singler[[j]]$about$RefData = c(singler$singler[[j]]$about$RefData,singler.list[[i]]$singler[[j]]$about$RefData)
      
      singler$singler[[j]]$SingleR.single$labels =
        rbind(singler$singler[[j]]$SingleR.single$labels,singler.list[[i]]$singler[[j]]$SingleR.single$labels)
      singler$singler[[j]]$SingleR.single$labels1 =  
        rbind(singler$singler[[j]]$SingleR.single$labels1,singler.list[[i]]$singler[[j]]$SingleR.single$labels1)
      singler$singler[[j]]$SingleR.single$scores =  
        rbind(singler$singler[[j]]$SingleR.single$scores,singler.list[[i]]$singler[[j]]$SingleR.single$scores)
      
      singler$singler[[j]]$SingleR.single.main$labels =  
        rbind(singler$singler[[j]]$SingleR.single.main$labels,singler.list[[i]]$singler[[j]]$SingleR.single.main$labels)
      singler$singler[[j]]$SingleR.single.main$labels1 =  
        c(singler$singler[[j]]$SingleR.single.main$labels1,singler.list[[i]]$singler[[j]]$SingleR.single.main$labels1)
      singler$singler[[j]]$SingleR.single.main$scores =  
        rbind(singler$singler[[j]]$SingleR.single.main$scores,singler.list[[i]]$singler[[j]]$SingleR.single.main$scores)
      
      singler$singler[[j]]$SingleR.single$cell.names =  
        c(singler$singler[[j]]$SingleR.single$cell.names,singler.list[[i]]$singler[[j]]$SingleR.single$cell.names)
      singler$singler[[j]]$SingleR.single.main$cell.names =  
        c(singler$singler[[j]]$SingleR.single$cell.names,singler.list[[i]]$singler[[j]]$SingleR.single.main$cell.names)
      
    }
    singler$meta.data$project.name = paste(singler$meta.data$project.name,singler.list[[i]]$meta.data$project.name,sep='+')
    singler$meta.data$orig.ident = c(singler$meta.data$orig.ident,singler.list[[i]]$meta.data$orig.ident)
    singler$meta.data$data.sets = c(singler$meta.data$data.sets,rep(singler.list[[i]]$meta.data$project.name,length(singler.list[[i]]$meta.data$orig.ident)))
    
  }
  
  for (j in 1:length(singler$singler)) {
    if (!is.null(order)) {
      singler$singler[[j]]$SingleR.single$labels = 
        singler$singler[[j]]$SingleR.single$labels[order,]
      singler$singler[[j]]$SingleR.single$labels1 = 
        singler$singler[[j]]$SingleR.single$labels1[order,]
      singler$singler[[j]]$SingleR.single$scores = 
        singler$singler[[j]]$SingleR.single$scores[order,]
      
      singler$singler[[j]]$SingleR.single.main$labels = 
        singler$singler[[j]]$SingleR.single.main$labels[order,]
      singler$singler[[j]]$SingleR.single.main$labels1 = 
        singler$singler[[j]]$SingleR.single.main$labels1[order,]
      singler$singler[[j]]$SingleR.single.main$scores = 
        singler$singler[[j]]$SingleR.single.main$scores[order,]
      
    }
    
    #  singler$singler[[j]]$SingleR.single$clusters = SingleR.Cluster(singler$singler[[j]]$SingleR.single,10)
    #  singler$singler[[j]]$SingleR.single.main$clusters = SingleR.Cluster(singler$singler[[j]]$SingleR.single.main,10)
  }
  
  if (!is.null(clusters) && !is.null(expr)) {
    #data('Immgen')
    #data('Mouse-RNAseq')
    #data ('HPCA')
    #data('Blueprint_Encode')
    for (j in 1:length(singler$singler)) {
      if (is.character(singler$singler[[j]]$about$RefData)) {
        ref = get(singler$singler[[j]]$about$RefData )
      } else {
        ref = singler$singler[[j]]$about$RefData
      }
      singler$singler[[j]]$SingleR.clusters = 
        SingleR("cluster",expr,ref$data,types=ref$types, 
                clusters = factor(clusters),sd.thres = ref$sd.thres,
                genes = 'de',fine.tune = T)
      singler$singler[[j]]$SingleR.clusters.main = 
        SingleR("cluster",expr,ref$data,types=ref$main_types, 
                clusters = factor(clusters),sd.thres = ref$sd.thres,
                genes = 'de',fine.tune = T)
    }
    singler$meta.data$clusters = clusters
    if (!is.null(xy)) {
      singler$meta.data$xy = xy
    }
  }
  
  singler
}

#' Creating a subset of the super reference dataset
#'
#' @param ref.name name of the reference dataset
#' @param ref the super reference dataset
#' @param ids numeric vector of indexes to use
#'
#' @return a specfic super reference dataset
superRef = function(ref.name,ref,ids) {
  ids = which(ids)
  dup = duplicated(names(ref$titles)[ids])
  ids = ids[!dup]
  de.genes = CreateVariableGeneSet(ref$data[,ids],ref$titles[ids],200)
  #de.genes.main = CreateVariableGeneSet(ref$data,main_types,300)
  de.genes.main = de.genes 
  list(name=ref.name,data=ref$data[,ids],types=ref$titles[ids],
       main_types=ref$titles[ids],de.genes=de.genes,
       de.genes.main=de.genes.main,sd.thres=1)
}
