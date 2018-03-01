# install.packages(c('Seurat',devtools','dplyr','Matrix','rowr','RColorBrewer','outliers','matrixStats','plyr','made4','ggplot2','grid','gridExtra','statmod','foreach','doSNOW','parallel'))
#source("https://bioconductor.org/biocLite.R")
#biocLite(c('made4','GSEABase','GSVA','GEOquery'))

library(RColorBrewer)
library(Seurat)
library(outliers)
library(matrixStats)
library(plyr)
library(made4)
library(ggplot2)
library(grid)
library(gridExtra)
require(statmod)
library(foreach)
library(doSNOW)
library(parallel)
library(Matrix)
library(GSEABase)
library(GSVA)
library(pheatmap)
library(limma)

# Colors
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
singler.colors = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
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

#' SingleR fine-tuning function. Takes top scores per cell type, finds variable genes in the
#' reference dataset for those cell types and recalculates the scores. This is performed until
#' only one cell type is left.
#'
#' @param sc_data  the single-cell RNA-seq data set as a matrix with genes as rownames. If the data if from a full-length platform, counts must be normalized to gene length (TPM, RPKM, FPKM, etc.).
#' @param ref_data the reference dataset with genes as rownames. Gene names must be in the same format as the single cell data (if sc_data uses genes symbols, ref_data must have the same)
#' @param types a list of cell type names corresponding to ref_data. Number of elements in types must be equal to number of columns in ref_data
#' @param scores the calculated scores per cell type produced by SingleR.ScoreData function
#' @param quantile.use correlation coefficients are aggregated for multiple cell types in the reference data set. This parameter allows to choose how to sort the cell types scores, by median (0.5) or any other number between 0 and 1. The deafult is 0.9.
#' @param fine.tune.thres the fine tunning step performs the scoring procedure for the top scoring labels using only genes that vary between those cell types in the reference data set. The top labels are those with a score lower than the top score by less than fine.tune.thres
#' @param genes list of genes to use for the annotations, or a method for extracting the genes from the data. Available methods are standard deviation and gene dispersion (\code{"sd"} or \code{"de"}). Deafult is "de".
#' @param sd.thres if genes=='sd' then this is the threshold for defining a variable gene.
#' @param mean_mat if genes='de' then this is a matrix with one column per cell type, each column is the
#' average of all samples for the corresponding cell type.
#'
#' @return the top labels for each single cell
SingleR.FineTune <- function(sc_data,ref_data,types,scores,quantile.use,fine.tune.thres,genes,sd.thres,mean_mat) {
  n = dim(sc_data)[2]
  numClusters = min(detectCores(all.tests = FALSE, logical = TRUE),16)
  print(paste("Fine-tunning round on top cell types (using", numClusters, "CPU cores):"))
  cl = makeCluster(numClusters)
  registerDoSNOW(cl)
  pb = txtProgressBar(max = n, style = 3)

  progress = function(n) setTxtProgressBar(pb, n)
  opts = list(progress = progress)

  labels = foreach(i = 1:n, .combine = rbind, .options.snow = opts, .export=c("Matrix","rowSds")) %dopar% {
    max_score = max(scores[i,])
    topLabels = names(scores[i,scores[i,]>=max_score-fine.tune.thres])
    if (length(topLabels)==0) {
      return (names(which.max(scores[i,])))
    } else {
      while(length(topLabels)>1) {
        labels.use = is.element(types,topLabels)
        ref_data.filtered = as.matrix(ref_data[,labels.use])
        types.filtered = types[labels.use]
        mat = mean_mat[,topLabels]

        if (genes[1] == "de") {
          n = round(500*(2/3)^(log2(c(ncol(mat)))))
          genes.filtered = unique(unlist(unlist(lapply(1:ncol(mat), function(j) {lapply(1:ncol(mat), function(i) {s=sort(mat[,j]-mat[,i],decreasing=T);s=s[s>0];names(s)[1:min(n,length(s))]})}))))[-1]
        } else if (genes[1] == "sd") {
          sd =  rowSds(mat)
          thres = min(sort(sd,decreasing = TRUE)[500],sd.thres)
          genes.filtered=intersect(rownames(ref_data)[sd>=thres],rownames(sc_data))
        } else {
          genes.filtered=intersect(genes,intersect(rownames(sc_data),(rownames(ref_data))))
        }

        ref_data.filtered = ref_data.filtered[genes.filtered,]
        sc_data.filtered = as.matrix(sc_data[genes.filtered,])
        data = sc_data.filtered[,i]
        if (sd(data)>0) {
          r=cor(data,ref_data.filtered,method='spearman')
          agg_scores = aggregate(t(r)~types.filtered,FUN = quantile, probs  = quantile.use)
          max_score = max(agg_scores[,-1])
          agg_scores = agg_scores[-which.min(agg_scores[,-1]),]
          topLabels = agg_scores[agg_scores[,-1]>=max_score-fine.tune.thres,1]
        } else {
          topLabels = agg_scores[which.max(agg_scores[,-1]),1]
        }
      }
      return (topLabels)
    }
  }
  close(pb)
  stopCluster(cl)

  if (dim(sc_data)[2]>1) {
    rownames(labels)=t(colnames(sc_data))
  }
  return(labels)

}

#' Scoring single cells using reference data set
#'
#' @param sc_data  the single-cell RNA-seq data set as a matrix with genes as rownames. If the data if from a full-length platform, counts must be normalized to gene length (TPM, RPKM, FPKM, etc.).
#' @param ref_data the reference dataset with genes as rownames. Gene names must be in the same format as the single cell data (if sc_data uses genes symbols, ref_data must have the same)
#' @param genes the list of genes to use.
#' @param types a list of cell type names corresponding to ref_data. Number of elements in types must be equal to number of columns in ref_data
#' @param quantile.use correlation coefficients are aggregated for multiple cell types in the reference data set. This parameter allows to choose how to sort the cell types scores, by median (0.5) or any other number between 0 and 1. The deafult is 0.9.
#'
#' @return a list with the scores, the raw correlation coefficients and the top labels
SingleR.ScoreData <- function(sc_data,ref_data,genes,types,quantile.use) {
  sc_data = as.matrix(sc_data[genes,])
  ref_data = as.matrix(ref_data[genes,])
  r=cor(sc_data,ref_data,method='spearman')

  agg_scores = aggregate(t(r)~types,FUN = quantile, probs  = quantile.use)
  labels = agg_scores[max.col(t(agg_scores[,-1])),1]
  output = list()


  if (dim(sc_data)[2]>1) {
    names(labels)=t(colnames(sc_data))
  }

  output$scores = as.matrix(agg_scores[,-1])
  rownames(output$scores) = agg_scores[,1]

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
#' @param genes list of genes to use for the annotations, or a method for extracting the genes from the data. Available methods are standard deviation and gene dispersion (\code{"sd"} or \code{"de"}). Deafult is "de".
#' @param quantile.use correlation coefficients are aggregated for multiple cell types in the reference data set. This parameter allows to choose how to sort the cell types scores, by median (0.5) or any other number between 0 and 1. The deafult is 0.9.
#' @param p.threshold Chi-square outlier detection is used to assess the signifcance power of the top correlation. Single-cell with an annotation of p-value > p.threshold are designated as "X". Only applies for non fine-tuned annotations.
#' @param fine.tune perform the fine tunning step? Deafult is TRUE.
#' @param fine.tune.thres the fine tunning step performs the scoring procedure for the top scoring labels using only genes that vary between those cell types in the reference data set. The top labels are those with a score lower than the top score by less than fine.tune.thres
#' @param sd.thres if genes=='sd' then this is the threshold for defining a variable gene.
#'
#' @return a list with the labels and scores
SingleR <- function(method = "single", sc_data, ref_data, types, clusters = NULL, genes = "sd", quantile.use = 0.9, p.threshold = 0.05, fine.tune = TRUE, fine.tune.thres = 0.05,sd.thres=1) {
  rownames(ref_data) = tolower(rownames(ref_data))
  rownames(sc_data) = tolower(rownames(sc_data))
  A = intersect(rownames(ref_data),rownames(sc_data))
  sc_data = as.matrix(sc_data[A,])
  ref_data = ref_data[A,]
  if (ncol(sc_data)>1) {
    not.use = rowSums(is.na(ref_data))>0 | rowSums(is.na(sc_data))>0 | rowSums(ref_data)==0
    ref_data = ref_data[!not.use,]
    sc_data = sc_data[!not.use,]
  }
  mat = medianMatrix(ref_data,types)
  if (genes[1] == "de") {
    n = round(500*(2/3)^(log2(c(ncol(mat)))))
    genes.filtered = unique(unlist(unlist(lapply(1:ncol(mat), function(j) {lapply(1:ncol(mat), function(i) {s=sort(mat[,j]-mat[,i],decreasing=T);s=s[s>0];names(s)[1:min(n,length(s))]})}))))[-1]
    print(paste0("Number of DE genes:", length(genes.filtered)))
  } else if (genes[1] == "sd") {
    sd =  rowSds(as.matrix(mat))
    genes.filtered=intersect(rownames(mat)[sd>sd.thres],rownames(sc_data))
    print(paste0("Number of genes with SD>",sd.thres,": ",length(genes.filtered)))
  } else {
    print(paste("Number of genes using in analysis:",length(genes.filtered)))
    genes.filtered=intersect(genes,intersect(rownames(sc_data),(rownames(ref_data))))
  }

  cell.names = colnames(sc_data)

  if (method == "single") {
    print(paste("Number of cells:",dim(sc_data)[2]))
  } else if (method == "cluster") {
    n = length(levels(clusters))
    print(paste("Number of clusters:",n))
    data = matrix(nrow=dim(sc_data)[1],ncol=n)
    for (i in 1:n) {
      data[,i] = rowSums(as.matrix(sc_data[,is.element(clusters,levels(clusters)[i])]))
    }
    colnames(data) = levels(clusters)
    rownames(data) = rownames(sc_data)
    sc_data = data
  } else {
    print ("Error: method must be 'single' or 'cluster'")
    return(0)
  }
  output = SingleR.ScoreData(sc_data,ref_data,genes.filtered,types,quantile.use)
  output$pval = apply(output$scores, 1,function(x) chisq.out.test(x)$p.value)

  # second round with top labels
  if (fine.tune==TRUE & length(unique(types)) > 2) {
    labels = SingleR.FineTune(sc_data,ref_data,types,output$scores,quantile.use,fine.tune.thres,genes = genes[1],sd.thres,mat)
    output$labels1 = as.matrix(output$labels)
    output$labels = as.matrix(labels)
    output$labels1.thres = c(output$labels)
    output$labels1.thres[output$pval>p.threshold] = "X"
  } else {
    labels = as.matrix(output$labels)
    output$labels.thres = c(output$labels)
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
    ggtitle(paste('R =', round(1000*cor(df$x,df$y,method='spearman'))/1000))
}

#' Plot boxplots for each lablel for a given single cell.
#'
#' @param sc_data the single-cell RNA-seq data set as a matrix with genes as rownames.
#' @param cell_id a number of the single cell to use
#' @param ref the reference dataset with genes as rownames. Gene names must be in the same format as the single cell data (if sc_data uses genes symbols, ref_data must have the same)
#' @param labels.use a list of labels to use. If NULL uses all labels.
#' @param quantile.order same a quantile.use - by which percentile to order the boxplots.
#' @param main_types aggregate labels by main types or by all types
#' @param top.n number of boxplots to present (starting from top)
#' @param tit title for the boxplot
#' @param colors colors to use. Defualt is singler.colors
#'
#' @return a list with a ggplot and the scores for the single cell
SingleR.DrawBoxPlot = function(sc_data, cell_id, ref, labels.use=NULL, quantile.order = 0.9, main_types=F, top.n=50, tit = NULL, colors=singler.colors) {
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

  res = SingleR(sc_data=as.matrix(sc_data[,c(cell_id,cell_id)]),ref_data=ref$data[,types.use],types=types,fine.tune=F,sd.thre=ref$sd.thres)
  if (is.null(top.n)) {
    top.n = length(unique(types))
  }
  A = types %in% names(sort(res$scores[1,],decreasing = T)[1:top.n])
  df = data.frame(Spearman=res$r[1,A],Types=types[A])
  fac <- with(df, reorder(Types, Spearman, function(x) quantile(x, probs  = quantile.order), order = TRUE))
  df$Types <- factor(df$Types, levels = levels(fac))

  p = ggplot(df, aes(x = Types, y = Spearman,  color = Types)) +
    geom_point(alpha = 0.5, position = "jitter",shape=16) +
    geom_boxplot(alpha = 0, colour = "black")+theme_classic() +
    xlab('') + ggtitle(paste(colnames(sc_data)[cell_id])) +

    theme(legend.position="none"
          , axis.text.x = element_text(angle = 45, hjust = 1),plot.title = element_text(hjust = 0.5)
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
#' @param labels labels to present, if NULL all labels presented
#' @param clusters a clustering to present as annotation in the heatmap
#' @param top.n number of cell types to presents. Defualt is 40. This can have an effect on the clustering which is performed only on the cell types presented.
#' @param normalize if TRUE scores are normalized to a 0-1 scale.
#' @param order.by.clusters if TRUE columns are ordered by the input clusters, and are not clustered again
#' @param cells_order an input order for the column
SingleR.DrawHeatmap = function(SingleR,labels = NULL,clusters=NULL,top.n=40,normalize=T,order.by.clusters=F,cells_order=NULL) {
  if (is.null(labels)) {
    labels = rownames(SingleR$scores)
  }
  m = apply(SingleR$scores[labels,],2,max)

  thres = sort(m,decreasing=TRUE)[min(top.n,length(m))]

  if (normalize==T) {
    mmax = rowMaxs(data)
    mmin = rowMins(data)
    data = (data-mmin)/(mmax-mmin)
  }

  if (SingleR$method == "cluster") {
    data = as.matrix(SingleR$scores^3)
  } else if (SingleR$method == "single") {
    data = as.matrix(SingleR$scores[labels,m>(thres-1e-6)]^3)
  }
  data = t(data)
  if (!is.null(clusters)) {
    clusters = as.data.frame(clusters)
    colnames(clusters) = 'Clusters'
    rownames(clusters) = colnames(data)

  }
  if (order.by.clusters==T) {
    data = data[,order(clusters$Clusters)]
    clusters = clusters[order(clusters$Clusters),,drop=F]
    pheatmap(data,border_color=NA,show_colnames=FALSE,clustering_method='ward.D',fontsize_row=6,annotation_col = clusters,cluster_cols = F)
  } else if (!is.null(cells_order)) {
    data = data[,order(cells_order)]
    clusters = clusters[order(cells_order),,drop=F]
    pheatmap(data,border_color=NA,show_colnames=FALSE,clustering_method='ward.D',fontsize_row=6,annotation_col = clusters,cluster_cols = F)
  } else {
    if (!is.null(clusters)) {
      pheatmap(data,border_color=NA,show_colnames=FALSE,clustering_method='ward.D',fontsize_row=6,annotation_col = clusters)
    } else {
      pheatmap(data,border_color=NA,show_colnames=FALSE,clustering_method='ward.D',fontsize_row=6)

    }
  }
}

#' Plot a colored tSNE plot according to a set of labels
#'
#' @param SingleR the output from the SingleR function
#' @param xy a matrix with the coordinates of the single cells
#' @param labels labels for the single cells, defualt is the SingleR labels
#' @param clusters if the SingleR method is 'clusters' then this vector is the cluster id for each single cell
#' @param do.letters if TRUE shows letters the first letter of the annotation
#' @param dot.size the size of the dot. Default is 1
#' @param do.labels only relevant in a 'cluster' method. If TRUE The cluster annotation is shown on top of the cluster.
#' @param do.legend if TRUE legend is presened.
#' @param label.size if do.labels is TRUE then this determines the font size of the labels
#' @param title title for the plot
#' @param colors colors to use
#' @param font.size size of fonts in the plot
#' @param alpha an alpha for the transperancy of the dots
#'
#' @return a list with ggplot and a data frame with the coordinates and annotations
SingleR.PlotTsne = function(SingleR, xy, labels=SingleR$labels, clusters = NULL, do.letters = TRUE, dot.size = 1, do.labels = FALSE, do.legend = TRUE, label.size=3, title = "",colors=singler.colors,font.size=NULL,alpha=0.5) {
  if (do.labels == TRUE)
    do.letters = FALSE

  df = data.frame(row.names = SingleR$cell.names)
  df$x = xy[,1]
  df$y = xy[,2]
  if (SingleR$method == "cluster") {
    df$ident = clusters.map.values(clusters,labels)
    #df$ident <- getClusterLabels(clusters,labels)
    #as.factor(plyr::mapvalues(clusters, from = current.cluster.ids, to = new.cluster.ids))
  } else {
    df$ident = as.factor(labels)
  }
  if (sum(levels(df$ident) %in% 'Other')) {
    lev = levels(df$ident)
    df$ident = factor(df$ident,levels = c(lev[-which(lev %in% 'Other')],'Other'))
  }

  num.levels = length(levels(df$ident))
  if (num.levels<0)
    colors = getcol(c(1:length(unique(labels))))

  p = ggplot(df, aes(x = x, y = y))

  p = p + geom_point(aes(color=ident), size=dot.size,alpha=alpha)

  if( do.letters == TRUE) {
    p = p + geom_point(aes(shape=as.character(ident)), size=dot.size/2)
    p = p + scale_shape_identity()
  } else {
  }

  if (do.labels == TRUE) {
    ctx <- ddply(df, .(ident), summarise, x = median(x))
    cty <- ddply(df, .(ident), summarise, y = median(y))
    centers = data.frame(ident = ctx$ident,x = ctx$x,y=cty$y)
    p = p + geom_point(data = centers, aes(x=x, y=y), size=0, alpha=0) + geom_text(data=centers, aes(label=ident), size = label.size, fontface="bold")
    p = p + guides(colour=FALSE)
    #x.range = ggplot_build(p)$panel$ranges[[1]]$x.range
    x.range = layer_scales(p)$x$range$range
    add_to_x = sum(abs(x.range))*0.03
    p = p + xlim(x.range[1]-add_to_x,x.range[2]+add_to_x)
    #p = p + + xlim(-5000, 5000)
  } else {
    #  n = 1
    # if (num.levels>35) {
    #   n = round(0.5+num.levels/35)
    #   p = p + theme(legend.position="bottom",legend.direction="vertical",legend.text=element_text(size=8),legend.title = element_blank()) + guides(col=guide_legend(ncol=5))
    # } else {
    if (is.null(font.size)) {
      font.size = 250*(1/num.levels)
      font.size = max(font.size,5)
      font.size = min(font.size,10)
    }
    p = p + theme(legend.text=element_text(size=font.size),legend.title = element_blank())+ guides(col=guide_legend(ncol=1))
    #  }
  }
  p = p + scale_color_manual(name="Type",values =colors)
  p = p + xlab("tSNE 1") + ylab("tSNE 2") + ggtitle(title)

  if (do.legend==FALSE) {
    p = p + theme(legend.position="none")
  }
  p + theme_bw()
  out = list(p=p,df=df,num.levels=num.levels)

}

#' Calculate single-sample gene set enrichment (ssGSEA) for each single cell
#'
#' @param sc_data  the single-cell RNA-seq data set as a matrix with genes as rownames. If the data if from a full-length platform, counts must be normalized to gene length (TPM, RPKM, FPKM, etc.).
#' @param species (\code{"Mouse"} or \code{"Human"})
#' @param signatures a GeneSetCollection oject, or NULL to use default signatures
#'
#' @return scores for each signature and single cell
calculateSignatures = function(sc_data,species='Human',signatures=NULL) {
  data('signatures')
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
  numClusters = min(detectCores(all.tests = FALSE, logical = TRUE),16)
  scores = gsva(sc_data,egc,method='ssgsea',parallel.sz=numClusters)
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
#' @param seurat_clusters the clusters for each single cell
#' @param singler_clusters the cluster annotation for each cluster
#'
#' @return annotation for each single cell
clusters.map.values = function(seurat_clusters,singler_clusters) {
  singler_clusters = as.matrix(singler_clusters)
  seurat_clusters = as.matrix(seurat_clusters)
  singler_clusters = as.matrix(singler_clusters[sort.int(as.numeric(rownames(singler_clusters)),index.return=TRUE)$ix,])
  A = rownames(singler_clusters) %in% unique(seurat_clusters)
  clusters = mapvalues(seurat_clusters, from = sort(as.numeric(unique(seurat_clusters))), to = paste0(rownames(singler_clusters)[A],": ",singler_clusters[A]))
}

#' SingleR.Cluster - using the scores data to cluster single cells.
#'
#' @param SingleR the output from the SingleR function
#' @param num.clusts number of clusters to output
#'
#' @return cluster id for each single cell
SingleR.Cluster = function(SingleR,num.clusts=10) {
  hc = hclust(dist(SingleR$scores),method='ward.D')
  cl = cutree(hc,k=num.clusts)
  cl
}

#' Subseting a SingleR object
#'
#' @param singler as SingleR object which includes a Seurat object.
#' @param subsetdata
#'
#' @return
#' @export
#'
#' @examples
SingleR.Subset = function(singler,subsetdata) {
  singler$seurat@data = singler$seurat@data[,subsetdata]
  singler$seurat@dr$tsne@cell.embeddings = singler$seurat@dr$tsne@cell.embeddings[subsetdata,]
  singler$seurat@meta.data = subset(singler$seurat@meta.data,subsetdata)
  singler$seurat@ident = singler$seurat@ident[subsetdata]
  singler$seurat@cell.names = singler$seurat@cell.names[subsetdata]

  for (i in 1:length(singler$singler)) {
    singler$singler[[i]]$SingleR.single$cell.names = singler$singler[[i]]$SingleR.single$cell.names[subsetdata]
    singler$singler[[i]]$SingleR.clusters$cell.names = singler$singler[[i]]$SingleR.clusters$cell.names[subsetdata]
    singler$singler[[i]]$SingleR.single$scores = singler$singler[[i]]$SingleR.single$scores[subsetdata,]
    singler$singler[[i]]$SingleR.single$labels = as.matrix(singler$singler[[i]]$SingleR.single$labels[subsetdata,])
  }
  if (!is.null(singler[["signatures"]])) {
    singler$signatures = singler$signatures[subsetdata,]
  }
  if(!is.null(singler[['other']])) {
    singler$other = singler$other[subsetdata]
  }
  singler
}

#' Remove data from a SingleR object to make it smaller
#'
#' @param singler.data a SingleR object
#'
#' @return a smaller SingleR object
remove.Unnecessary.Data.single = function(singler.data) {
  for (j in 1:length(singler.data$singler)) {
    singler.data$singler[[j]]$SingleR.single = singler.data$singler[[j]]$SingleR.single[-c(which(names(singler.data$singler[[j]]$SingleR.single) %in% c('r','pval','labels1.thres','types')))]
    singler.data$singler[[j]]$SingleR.clusters = singler.data$singler[[j]]$SingleR.clusters[-c(which(names(singler.data$singler[[j]]$SingleR.clusters) %in% c('r','pval','labels1.thres','types')))]
    if (!is.null(singler.data$singler[[j]]$SingleR.single.main)) {
      singler.data$singler[[j]]$SingleR.single.main = singler.data$singler[[j]]$SingleR.single.main[-c(which(names(singler.data$singler[[j]]$SingleR.single.main) %in% c('r','pval','labels1.thres','types')))]
      singler.data$singler[[j]]$SingleR.clusters.main = singler.data$singler[[j]]$SingleR.clusters.main[-c(which(names(singler.data$singler[[j]]$SingleR.clusters.main) %in% c('r','pval','labels1.thres','types')))]
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
  if (is.null(lengths)) {
    data('gene_lengths')
  }
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

#' Wrapper function to create a SingleR object
#'
#' @param sc.data
#' @param ref
#' @param clusters
#' @param do.main.types
#' @param species
#' @param citation
#' @param technology
#'
#' @return
#' @export
#'
#' @examples
SingleR.CreateObject <- function(sc.data,ref,clusters,do.main.types=T,species='Human',citation='-',technology='-') {

  types = ref$types

  SingleR.single = SingleR("single",sc.data,ref$data,types=types,sd.thres = ref$sd.thres)
  SingleR.clusters = SingleR("cluster",sc.data,ref$data,types=types, clusters=factor(clusters),sd.thres = ref$sd.thres)

  about = list(Organism = capitalize(species),Citation=citation,Technology=technology,RefData=ref$name)
  singler = list(SingleR.single = SingleR.single, SingleR.clusters = SingleR.clusters,about=about)

  if (do.main.types==T) {
    types = ref$main_types
    singler$SingleR.single.main = SingleR("single",sc.data,ref$data,types=types,sd.thres = ref$sd.thres, quantile.use = 0.9)
    singler$SingleR.clusters.main = SingleR("cluster",sc.data,ref$data,types=types, clusters=factor(clusters),sd.thres = ref$sd.thres, quantile.use = 0.9)
  }

  if (!(ref$name %in% c('Immgen','RNAseq','HPCA','Blueprint_Encode','Fantom','GSE43005'))) {
    singler$about$refernce = ref
  }
  singler
}

#' Wrapper function to create a Seurat object
#'
#' @param project
#' @param sc.data
#' @param min.genes
#' @param min.cells
#' @param regress.out
#' @param npca
#' @param resolution
#' @param species
#' @param temp.dir
#'
#' @return
#' @export
#'
#' @examples
SingleR.CreateSeurat <- function(project,sc.data,min.genes = 500,min.cells = 2,regress.out = 'nUMI',npca = 10,resolution=0.8,species='Human',temp.dir=NULL) {
  sc = CreateSeuratObject(raw.data = sc.data, min.cells = min.cells, min.genes = min.genes, project = project)
  if (species == 'Human') {
    mtgenes = '^MT-'

  } else {
    mtgenes = '^Mt-'
  }
  mito.genes <- grep(pattern = mtgenes, x = rownames(x = sc@data), value = TRUE)
  percent.mito <- colSums((sc.data[mito.genes, ]))/colSums(sc.data)
  sc <- AddMetaData(object = sc, metadata = percent.mito, col.name = "percent.mito")
  #sc <- FilterCells(object = sc, subset.names = c("nGene", "percent.mito"), low.thresholds = c(200, -Inf), high.thresholds = c(2500, 0.05))

  sc <- NormalizeData(object = sc, normalization.method = "LogNormalize", scale.factor = 10000)
  sc <- FindVariableGenes(object = sc, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.contour = F, do.plot = F)

  sc <- ScaleData(object = sc, vars.to.regress = regress.out)

  sc <- RunPCA(object = sc, pc.genes = sc@var.genes, do.print = FALSE)

  sc <- FindClusters(object = sc, reduction.type = "pca", dims.use = 1:npca,resolution = resolution, print.output = 0, save.SNN = F, temp.file.location = temp.dir)

  if (ncol(sc@data)<100) {
    sc <- RunTSNE(sc, dims.use = 1:npca, do.fast = T,perplexity=10  )
  } else {
    sc <- RunTSNE(sc, dims.use = 1:npca, do.fast = T)

  }

  sc
}

#' Wrapper function to create a SingleR+Seurat object
#'
#' @param counts
#' @param annot
#' @param project.name
#' @param min.genes
#' @param min.cells
#' @param regress.out
#' @param npca
#' @param technology
#' @param species
#' @param citation
#' @param ref.list
#' @param reduce.file.size
#' @param normalize.gene.length
#' @param do.signatures
#' @param temp.dir
#'
#' @return
#' @export
#'
#' @examples
SingleR.CreateFromCounts = function(counts,annot=NULL,project.name,min.genes=500,min.cells=2,regress.out='nUMI',npca=10,technology='10X',species='Human',citation='',ref.list=list(),reduce.file.size=T,normalize.gene.length=F,do.signatures=T,temp.dir=NULL) {
  if (typeof(counts) == 'character') {
    counts <- as.matrix(read.table(counts, header=TRUE, sep="\t", row.names=1, as.is=TRUE,comment.char='!'))
  }
  print(dim(counts))

  A = tolower(rownames(counts))
  dupA = duplicated(A)
  if (sum(dupA)>0) {
    counts = counts[!dupA,]
  }

  if (species=='Mouse') {
    rownames(counts) = paste(toupper(substr(rownames(counts), 1, 1)), tolower(substr(rownames(counts), 2, nchar(rownames(counts)))), sep="")
  }

  if (!is.null(annot)) {
    if (length(annot) == 1) {
      types <- read.table(annot, header=TRUE, sep="\t", row.names=1, as.is=TRUE)
      orig.ident = types[,1]
      names(orig.ident) = rownames(types)
    } else {
      orig.ident = annot
    }
  } else {
    orig.ident = rep('NA',ncol(counts))
    names(orig.ident)=colnames(counts)
  }

  print(project.name)
  sc = SingleR.CreateSeurat(project.name,counts,min.genes=min.genes,min.cells=min.cells,regress.out=regress.out,npca=npca,species=species,temp.dir=temp.dir)

  orig.ident = orig.ident[colnames(sc@data)]
  sc.data = counts[,colnames(sc@data)]

  if (normalize.gene.length == F) {
    sc.data.gl = sc.data
    rownames(sc.data.gl) = tolower(rownames(sc.data.gl))
  } else {
    if (species == 'Human') {
      sc.data.gl = TPM(sc.data,human_lengths)
    } else if (species == 'Mouse') {
      sc.data.gl = TPM(sc.data,mouse_lengths)
    }
  }

  clusters = sc@ident

  if (length(ref.list)==0) {
    if (species == 'Mouse') {
      if (!exists('immgen'))
        data('Immgen')
      if (!exists('mouse.rnaseq'))
        data('Mouse-RNAseq')
      res = list(SingleR.CreateObject(sc.data.gl,immgen,clusters,species,citation,technology,do.main.types=T),
                 SingleR.CreateObject(sc.data.gl,mouse.rnaseq,clusters,species,citation,technology,do.main.types=T)
      )
    } else if (species == 'Human') {
      if(!exists('hpca'))
        data ('HPCA')
      if (!exists('blueprint_encode'))
        data('Blueprint_Encode')
      res = list(SingleR.CreateObject(sc.data.gl,hpca,clusters,species,citation,technology,do.main.types = T),
                 SingleR.CreateObject(sc.data.gl,blueprint_encode,clusters,species,citation,technology,do.main.types = T))
    }
  } else {
    res = lapply(ref.list, FUN=function(x) {
      SingleR.CreateObject(sc.data.gl,x,clusters,species,citation,technology,do.main.types=T)
    })
  }
  sc@meta.data$orig.ident = factor(orig.ident)

  if (do.signatures==TRUE) {
    signatures = calculateSignatures(sc.data.gl,species=species)
    singler = list(seurat = sc,singler = res)

  } else {
    singler = list(seurat = sc,signatures = signatures,singler = res)
  }


  if (species == 'Human') {
    if (!exists('cell.type.classification'))
      data('cell.type.cor.classification')
    rownames(cell.type.classification$cell.types.avg) = tolower(rownames(cell.type.classification$cell.types.avg))
    A = intersect(rownames(sc.data.gl),rownames(cell.type.classification$cell.types.avg))
    r = cor(cell.type.classification$cell.types.avg[A,],as.matrix(sc.data.gl[A,]),method='spearman')
    other_annotation = max.col(t(r))
    other_annotation=mapvalues(other_annotation,from=seq(1,ncol(cell.type.classification$cell.types.avg)),to=colnames(cell.type.classification$cell.types.avg))
    names(other_annotation) = colnames(r)
    singler$other = other_annotation
  }

  if(reduce.file.size==T) {
    singler$seurat@raw.data = c()
    singler$seurat@scale.data = c()
    singler$seurat@calc.params = list()
  }

  singler = remove.Unnecessary.Data.single(singler)

  singler

}


