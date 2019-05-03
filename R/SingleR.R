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
#' @param numCores Number of cores to use.
#'
#' @return the top labels for each single cell
SingleR.FineTune <- function(sc_data,ref_data,types,scores,quantile.use,
                             fine.tune.thres,genes,sd.thres,mean_mat,
                             numCores = SingleR.numCores) {
  N = dim(sc_data)[2]
  print(paste("Fine-tuning round on top cell types (using", numCores, 
              "CPU cores):"))
  labels = pbmclapply(1:N,FUN=function(i){
    max_score = max(scores[i,])
    topLabels = names(scores[i,scores[i,]>=max_score-fine.tune.thres])
    if (length(topLabels)==0) {
      return (names(which.max(scores[i,])))
    } else {
      k=1
      while(length(topLabels)>1) {
        topLabels = fineTuningRound(topLabels,types,ref_data,genes,
                                    mean_mat[,topLabels],sd.thres,
                                    sc_data[,i],quantile.use,fine.tune.thres)
        k=k+1
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
  if (length(genes.filtered)<20) {
    return (topLabels[1])
  }
  if (sd(sc_data.filtered)>0) {
    r=cor(sc_data.filtered,ref_data.filtered,method='spearman')
    agg_scores = quantileMatrix(r,types.filtered,quantile.use)[1,];
    agg_scores = sort(agg_scores,decreasing = T)
    agg_scores = agg_scores[-length(agg_scores)]
    topLabels = names(agg_scores)[agg_scores>=agg_scores[1]-fine.tune.thres]
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
#' @param numCores Number of cores to use.
#' @param step number of cells in each correlation analysis. The correlation analysis memory requirements may be too high, thus it can be split to smaller sets.
#'
#' @return a list with the scores, the raw correlation coefficients and the top labels
SingleR.ScoreData <- function(sc_data,ref_data,genes,types,quantile.use,numCores=1,step=10000) {
  sc_data = as.matrix(sc_data[genes,])
  ref_data = as.matrix(ref_data[genes,])
  
  if (ncol(sc_data)>step) {
    n = ncol(sc_data)
    s = seq(step+1,n,by=step)
    if (FALSE) {
      cl <- makeCluster(numCores)
      doParallel::registerDoParallel(cl)
      tmpr = foreach (i = 0:length(s)) %dopar% {
        if(i == 0){
          res = data.table::data.table(cor(sc_data[,1:step],ref_data,method='spearman'))
        } else {
          A = seq(s[i],min(s[i]+step-1,n))
          # r=rbind(r,cor(sc_data[,A],ref_data,method='spearman'))
          res = data.table::data.table(cor(sc_data[,A],ref_data,method='spearman'))
        }
        res
      }
      r = data.table::rbindlist(tmpr, use.names = F)
      r = as.matrix(r)
      rownames(r) = colnames(sc_data)
      on.exit(stopCluster(cl))
    } else {
      s = seq(step+1,n,by=step)
      r=cor(sc_data[,1:step],ref_data,method='spearman')
      for (i in 1:length(s)) {
        A = seq(s[i],min(s[i]+step-1,n))
        r = rbind(r,cor(sc_data[,A],ref_data,method='spearman'))
      }
    }
    
  } else {
    r=cor(sc_data,ref_data,method='spearman')
  }
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
#' @param numCores Number of cores to use.
#'
#' @return a list with the labels and scores
SingleR <- function(method = "single", sc_data, ref_data, types, 
                    clusters = NULL, genes = "de", quantile.use = 0.8, 
                    p.threshold = 0.05, fine.tune = TRUE, 
                    fine.tune.thres = 0.05,sd.thres=1, do.pvals = T, 
                    numCores = SingleR.numCores, ...) {
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
    n = round(500*(2/3)^(log2(c(ncol(mat)))))
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
    genes.filtered=intersect(tolower(genes),intersect(tolower(rownames(sc_data)),
                                                      tolower(rownames(ref_data))))
    print(paste("Number of genes using in analysis:",length(genes.filtered)))
    
  }
  
  cell.names = colnames(sc_data)
  
  if (method == "single") {
    if (dim(sc_data)[2]>2) {
      print(paste("Number of cells:",dim(sc_data)[2]))
    }
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
  output = SingleR.ScoreData(sc_data,ref_data,genes.filtered,types,quantile.use,numCores=numCores,...)
  if (do.pvals == T) {
    output$pval = SingleR.ConfidenceTest(output$scores)
  }
  
  # second round with top labels
  if (fine.tune==TRUE & length(unique(types)) > 2) {
    labels = SingleR.FineTune(sc_data,ref_data,types,output$scores,
                              quantile.use,fine.tune.thres,genes = genes,
                              sd.thres,mat,numCores = numCores)
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

#' Test confidence of SingleR annotation.
#'
#' Calculates outliers chisquare-test for the cell typres scores for each cell.
#'    
#' @param scores a SingleR scores matrix
#' 
#' @return vector of p-values
SingleR.ConfidenceTest = function(scores) {
  apply(scores, 1,function(x) chisq.out.test(x)$p.value)
}

#' Calculate single-sample gene set enrichment (ssGSEA) for each single cell
#'
#' @param sc_data  the single-cell RNA-seq data set as a matrix with genes as rownames. If the data if from a full-length platform, counts must be normalized to gene length (TPM, RPKM, FPKM, etc.).
#' @param species (\code{"Mouse"} or \code{"Human"})
#' @param signatures a GeneSetCollection object, or NULL to use default signatures
#' @param n.break run ssGSEA for n.break at a time.
#' @param numCores Number of cores to use.
#'
#' @return scores for each signature and single cell
calculateSignatures = function(sc_data,species='Human',signatures=NULL, 
                               n.break=1000, numCores = SingleR.numCores) {
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
  # break to groups of n.break cells
  scores = matrix(NA,length(egc),ncol(sc_data))
  wind = seq(1,ncol(sc_data),by=n.break)
  print(paste('Using sets of',n.break, 'cells. Running',length(wind),'times.'))
  for (i in wind) {
    last = min(ncol(sc_data),i+n.break-1)
    a = gsva(sc_data[,i:last],egc,method='ssgsea',ssgsea.norm=F,
             parallel.sz=numCores,parallel.type='FORK')
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

#' Calculate single-sample gene set enrichment (using singscore) for each single cell
#'
#' @param sc_data  the single-cell RNA-seq data set as a matrix with genes as rownames. If the data if from a full-length platform, counts must be normalized to gene length (TPM, RPKM, FPKM, etc.).
#' @param species (\code{"Mouse"} or \code{"Human"})
#' @param signatures a GeneSetCollection object, or NULL to use default signatures
#'
#' @return scores for each signature and single cell
calculateSingScores = function(sc_data,species='Human',signatures=NULL) {
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
  
  rankedData <- rankGenes(sc_data)
  scores = matrix(NA,ncol(rankedData),length(egc))
  options(warn=-1)
  for (i in 1:length(egc))
    scores[,i] <- simpleScore(rankedData, geneIds(egc[[i]]), centerScore = TRUE)$TotalScore
  options(warn=0)
  scores = t(scores)
  rownames(scores) = names(egc)
  colnames(scores) = colnames(rankedData)
  scores[is.na(scores)] = 0
  
  mmin = rowMins(scores)
  mmax = rowMaxs(scores)
  scores = scores/(mmax-mmin)
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
    scores = t(scale(t(SingleR$scores^3)))
  } else {
    scores = SingleR$scores
  }
  #r <- cor(t(scores), method="pearson")
  #d <- as.dist(1-r)
  #hc = hclust(d,method='ward.D2')
  hc = hclust(dist(scores,method='euclidean'),method='ward.D2')
  
  cl = cutree(hc,k=num.clusts)
  list(hc=hc,cl=factor(cl))
}

#' Subseting a SingleR object. This function subsets all the SingleR data and the Seurat object if included.
#'
#' @param singler as SingleR object
#' @param subsetdata a logical vector of single-cells to include in the subset object
#' @param rerun.seurat if TRUE reruns the Seurat analyses for the subset data.
#'
#' @return a subset of the original SingleR vector
SingleR.Subset = function(singler,subsetdata,rerun.seurat=F,rerun.singler.clustering=F) {
  s = singler
  
  if (!is.null(s$seurat)) {
    if (packageVersion('Seurat')>=3) {
      s$seurat = SubsetData(s$seurat,'RNA',colnames(s$seurat@assays$RNA)[subsetdata])
      subsetdata = unlist(lapply(colnames(s$seurat@assays$RNA),FUN=function(x) 
        which(singler$singler[[1]]$SingleR.single$cell.names==x)))
    } else {
      s$seurat = SubsetData(s$seurat,colnames(s$seurat@data)[subsetdata])
      subsetdata = unlist(lapply(s$seurat@cell.names,FUN=function(x) 
        which(singler$singler[[1]]$SingleR.single$cell.names==x)))
    }
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
    if (!is.null(s$singler[[i]]$SingleR.single$labels1)) {
      s$singler[[i]]$SingleR.single$labels1 = 
        as.matrix(s$singler[[i]]$SingleR.single$labels1[subsetdata,])
    }
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
      if (!is.null(s$singler[[i]]$SingleR.single.main$labels1)) {
        s$singler[[i]]$SingleR.single.main$labels1 = 
          as.matrix(s$singler[[i]]$SingleR.single.main$labels1[subsetdata,])
      }
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
  
  if (rerun.seurat==T) {
    if (packageVersion('Seurat')>=3) {
      s$seurat <- SCTransform(object = s$seurat, vars.to.regress = "percent.mito", verbose = FALSE,
                              do.correct.umi=T)
      s$seurat <- RunPCA(object = s$seurat,verbose = FALSE)
      s$seurat <- FindNeighbors(object = s$seurat, dims = 1:30)
      s$seurat <- FindClusters(object = s$seurat)
      if (ncol(s$seurat@assays$RNA@data)<100) {
        s$seurat <- RunTSNE(s$seurat,perplexity=10,dims = 1:30)
      } else {
        s$seurat <- RunTSNE(s$seurat,dims = 1:30)
      }
      s$seurat <- RunUMAP(s$seurat,dims = 1:30, verbose = FALSE)
      
      s$meta.data$xy=s$seurat@reductions$umap@cell.embeddings
      s$meta.data$clusters=s$seurat@active.ident
    } else {
      s$seurat <- FindVariableGenes(object = s$seurat, mean.function = ExpMean,
                                    dispersion.function = LogVMR,
                                    x.low.cutoff = 0.0125, x.high.cutoff = 3,
                                    y.cutoff = 0.5, do.contour = F, do.plot = F)
      regress.out='nUMI'
      s$seurat <- ScaleData(object = s$seurat, vars.to.regress = regress.out)
      s$seurat <- RunPCA(object = s$seurat, pc.genes = s$seurat@var.genes, do.print = FALSE)
      PCElbowPlot(object = s$seurat)
      npca=10
      resolution=0.8
      s$seurat <- FindClusters(object = s$seurat, reduction.type = "pca",
                               dims.use = 1:npca,resolution = resolution,
                               print.output = 0, save.SNN = F)
      s$seurat <- RunTSNE(s$seurat, dims.use = 1:npca, do.fast = T)
      s$meta.data$xy=s$seurat@dr$tsne@cell.embeddings
      s$meta.data$clusters=s$seurat@ident
    }
    if (rerun.singler.clustering == T) {
      for (i in 1:length(s$singler)) {
        ref = get(gsub('-','\\.',tolower(s$singler[[i]]$about$RefData)))
        s$singler[[i]]$SingleR.clusters = 
          SingleR(method = "cluster", s$seurat@data, ref$data,ref$types,
                  clusters = s$meta.data$clusters, genes=ref$de.genes)
        s$singler[[i]]$SingleR.clusters.main = 
          SingleR(method = "cluster", s$seurat@data, ref$data,ref$main_types,
                  clusters = s$meta.data$clusters, genes=ref$de.genes.main)
      }
    }
  }
  s
}
