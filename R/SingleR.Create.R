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
#' @param numCores Number of cores to use.
#'
#' @return a SingleR objects
SingleR.CreateObject <- function(sc.data,ref,clusters=NULL,species='Human',
                                 citation='-',technology='-',variable.genes='de',
                                 fine.tune=T,do.main.types=T,
                                 numCores = SingleR.numCores) {
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
                           fine.tune = fine.tune,numCores = numCores)
  
  if (is.null(clusters)) {
    SingleR.single$clusters = SingleR.Cluster(SingleR.single,10)
    clusters = SingleR.single$clusters$cl
  }
  
  SingleR.clusters = SingleR("cluster",sc.data,ref$data,types=types, 
                             clusters = factor(clusters),
                             sd.thres = ref$sd.thres,
                             genes = variable.genes,
                             fine.tune = fine.tune,numCores = numCores)
  
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
                                          fine.tune = fine.tune,
                                          numCores = numCores)
    if (is.null(clusters)) {
      singler$SingleR.single.main$clusters = 
        SingleR.Cluster(singler$SingleR.single.main,10)
    }
    singler$SingleR.clusters.main = 
      SingleR("cluster",sc.data,ref$data,types=types, 
              clusters=factor(clusters),sd.thres = ref$sd.thres, 
              quantile.use = 0.8,genes = variable.genes.main,
              fine.tune = fine.tune,numCores = numCores)
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
#' @param npca a vector of the dimensions to use in construction of the clustering and tSNE plot (not used in Seurat v3)
#' @param resolution clustering resolution. See Seurat manual for more details.
#' @param temp.dir used by FindClusters function.
#'
#' @return a Seurat object
SingleR.CreateSeurat <- function(project.name,sc.data,min.genes = 200,
                                 min.cells = 2,regress.out = 'nUMI',
                                 npca = 10,resolution=0.8,temp.dir=NULL) {
  mtgenes = '^mt-'
  
  if (packageVersion('Seurat')>=3) {
    sc = CreateSeuratObject(sc.data, min.cells = min.cells, 
                            min.features = min.genes, project = project.name)
    percent.mito <- PercentageFeatureSet(object = sc, pattern = "^(?i)mt-")
    # mito.features <- grep(pattern = mtgenes, x = rownames(x = sc), value = TRUE,ignore.case=TRUE)
    #  percent.mito <- Matrix::colSums(x = GetAssayData(object = sc, slot = 'counts')[mito.features, ]) / Matrix::colSums(x = GetAssayData(object = sc, slot = 'counts'))
    sc <- AddMetaData(object = sc, metadata = percent.mito, 
                      col.name = "percent.mito")
  } else {
    sc = CreateSeuratObject(sc.data, min.cells = min.cells, 
                            min.genes = min.genes, project = project.name)
    mito.genes <- grep(pattern = mtgenes, x = rownames(x = sc@data), 
                       value = TRUE,ignore.case=TRUE)
    percent.mito <- colSums((sc.data[mito.genes, ]))/colSums(sc.data)
    sc <- AddMetaData(object = sc, metadata = percent.mito, 
                      col.name = "percent.mito")
    
    sc <- NormalizeData(object = sc, 
                        normalization.method = "LogNormalize", 
                        scale.factor = 10000)
  }
  
  if (packageVersion('Seurat')>=3) {
    sc <- SCTransform(object = sc, vars.to.regress = "percent.mito", verbose = FALSE,
                      do.correct.umi=T)
    
    # sc <- FindVariableFeatures(object = sc, selection.method = 'mean.var.plot', 
    #                             mean.cutoff = c(0.0125, 3), 
    #                            dispersion.cutoff = c(0.5, Inf) ,
    #                           do.contour = F, do.plot = F)
    
    #sc <- ScaleData(object = sc,use.umi=T)
    sc <- RunPCA(object = sc,verbose = FALSE)
    sc <- FindNeighbors(object = sc, dims = 1:30)
    sc <- FindClusters(object = sc)
    if (ncol(sc@assays$RNA@data)<100) {
      sc <- RunTSNE(sc,perplexity=10,dims = 1:npca)
    } else {
      sc <- RunTSNE(sc,dims = 1:30)
    }
    sc <- RunUMAP(sc,dims = 1:30, verbose = FALSE)
  } else {
    sc <- FindVariableGenes(object = sc, mean.function = ExpMean, 
                            dispersion.function = LogVMR, 
                            x.low.cutoff = 0.0125, x.high.cutoff = 3, 
                            y.cutoff = 0.5, do.contour = F, do.plot = F)
    if (!is.null(regress.out)) {
      sc <- ScaleData(object = sc, vars.to.regress = regress.out)
    } else {
      sc <- ScaleData(object = sc)
    }
    sc <- RunPCA(object = sc, pc.genes = sc@var.genes, do.print = FALSE)
    #PCElbowPlot(object = sc)
    sc <- FindClusters(object = sc, reduction.type = "pca", 
                       dims.use = 1:npca,resolution = resolution, 
                       print.output = 0, save.SNN = F, 
                       temp.file.location = temp.dir)
    if (ncol(sc@data)<100) {
      sc <- RunTSNE(sc, dims.use = 1:npca, do.fast = T,perplexity=10  )
    } else {
      sc <- RunTSNE(sc, dims.use = 1:npca, do.fast = T,check_duplicates = FALSE)
      
    }
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
  
  colnames(counts) = make.unique(colnames(counts))
  names(orig.ident) = colnames(counts)
  
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
#' @param numCores Number of cores to use.
#'
#' @return a SingleR object containing a Seurat object
CreateSinglerSeuratObject = function(counts,annot=NULL,project.name,
                                     min.genes=200,technology='10X',
                                     species='Human',citation='',
                                     ref.list=list(),normalize.gene.length=F,
                                     variable.genes='de',fine.tune=T,
                                     reduce.file.size=T,do.signatures=F,
                                     min.cells=2,npca=10,regress.out='nUMI',
                                     do.main.types=T,reduce.seurat.object=T,
                                     temp.dir=NULL, numCores = SingleR.numCores) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat needed for this function to work. Please install it.",
         call. = FALSE)
  }
  
  print(project.name)
  print('Reading single-cell data...')
  sc.data = ReadSingleCellData(counts,annot)
  print('Create Seurat object...')
  seurat = SingleR.CreateSeurat(project.name,sc.data$counts,min.genes=
                                  min.genes,min.cells=min.cells,
                                regress.out=regress.out,npca=npca,
                                temp.dir=temp.dir)
  if (packageVersion('Seurat')>=3) {
    data = seurat@assays$RNA@data
    clusters = seurat@active.ident
  } else {
    data = seurat@data
    clusters = seurat@ident
    
  }
  
  orig.ident = sc.data$orig.ident[colnames(data)]
  counts = as.matrix(sc.data$counts[,colnames(data)])
  
  seurat@meta.data$orig.ident = factor(orig.ident)
  
  if(reduce.seurat.object==T) {
    if (packageVersion('Seurat')>=3) {
      seurat@assays$RNA@counts = matrix()
      seurat@assays$RNA@scale.data = matrix()
    } else {
      seurat@raw.data = c()
      seurat@scale.data = c()
      seurat@calc.params = list()
    }
  }
  
  print('Creat SingleR object...')
  
  singler = CreateSinglerObject(counts,orig.ident,project.name,
                                min.genes=min.genes,technology,species,
                                citation,ref.list,
                                normalize.gene.length,variable.genes,
                                fine.tune,do.signatures,
                                clusters,do.main.types,
                                reduce.file.size,temp.dir,numCores = numCores)
  
  singler$seurat = seurat 
  
  if (packageVersion('Seurat')>=3) {
    singler$meta.data$xy = seurat@reductions$tsne@cell.embeddings
    singler$meta.data$clusters = seurat@active.ident
    
  } else {
    singler$meta.data$xy = seurat@dr$tsne@cell.embeddings
    singler$meta.data$clusters = seurat@ident
  }
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
#' @param numCores Number of cores to use.
#'
#' @return a SingleR object
CreateSinglerObject = function(counts,annot=NULL,project.name,
                               min.genes=0,technology='10X',
                               species='Human',citation='',
                               ref.list=list(),normalize.gene.length=F,
                               variable.genes='de',fine.tune=T,
                               do.signatures=F,clusters=NULL,
                               do.main.types=T,reduce.file.size=T,
                               temp.dir=NULL,numCores = SingleR.numCores) {
  
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
                                      fine.tune=fine.tune,numCores = numCores),
                 SingleR.CreateObject(sc.data.gl,mouse.rnaseq,clusters,
                                      species,citation,technology,
                                      do.main.types=do.main.types,
                                      variable.genes=variable.genes,
                                      fine.tune=fine.tune,numCores = numCores)
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
                                      fine.tune=fine.tune,numCores = numCores),
                 SingleR.CreateObject(sc.data.gl,blueprint_encode,
                                      clusters,species,citation,technology,
                                      do.main.types = do.main.types,
                                      variable.genes=variable.genes,
                                      fine.tune=fine.tune,numCores = numCores))
    }
  } else {
    res = lapply(ref.list, FUN=function(x) {
      SingleR.CreateObject(sc.data.gl,x,clusters,species,citation,technology,
                           do.main.types=do.main.types,
                           variable.genes=variable.genes,fine.tune=fine.tune,
                           numCores = numCores)
    })
  }
  
  singler$singler = res
  
  if (do.signatures==TRUE) {
    signatures = calculateSingScores(sc.data.gl,species=species)
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
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Seurat needed for this function to work. Please install it.",
         call. = FALSE)
  }
  print(paste(1,basename(dirs[1])))
  sc.data = as.matrix(Read10X(dirs[1]))
  
  A = colSums(as.matrix(sc.data)>0)>min.genes
  if(sum(A)>0) {
    sc.data = sc.data[,A]
    if (random.sample>0) {
      use = sample(ncol(sc.data),min(random.sample,ncol(sc.data)))
      sc.data = sc.data[,use]
    } 
    orig.ident = rep(basename(dirs[1]),ncol(sc.data))
  } else {
    print('No cells with nGene > min.genes')
    orig.ident = c()
    sc.data = c()
  }
  for (j in 2:length(dirs)) {
    print(paste(j,basename(dirs[j])))
    data <- Read10X(dirs[j])
    
    A = colSums(as.matrix(data)>0)>min.genes
    if(sum(A)>0) {
      
      data = data[,A]
      
      if (random.sample>0) {
        use = sample(ncol(data),min(random.sample,ncol(data)))
        data = data[,use]
      }
      orig.ident = c(orig.ident,rep(basename(dirs[j]),ncol(data)))
      g = intersect(rownames(sc.data),rownames(data))
      sc.data = cbind(sc.data[g,],data[g,])
    } else {
      print('No cells with nGene > min.genes')
    }
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
      if (!is.null(singler$singler[[j]]$SingleR.single$labels1)) {
        singler$singler[[j]]$SingleR.single$labels1 =  
          rbind(singler$singler[[j]]$SingleR.single$labels1,singler.list[[i]]$singler[[j]]$SingleR.single$labels1)
      }
      singler$singler[[j]]$SingleR.single$scores =  
        rbind(singler$singler[[j]]$SingleR.single$scores,singler.list[[i]]$singler[[j]]$SingleR.single$scores)
      
      singler$singler[[j]]$SingleR.single.main$labels =  
        rbind(singler$singler[[j]]$SingleR.single.main$labels,singler.list[[i]]$singler[[j]]$SingleR.single.main$labels)
      if (!is.null(singler$singler[[j]]$SingleR.single.main$labels1)) {
        
        singler$singler[[j]]$SingleR.single.main$labels1 =  
          rbind(singler$singler[[j]]$SingleR.single.main$labels1,singler.list[[i]]$singler[[j]]$SingleR.single.main$labels1)
      }
      singler$singler[[j]]$SingleR.single.main$scores =  
        rbind(singler$singler[[j]]$SingleR.single.main$scores,singler.list[[i]]$singler[[j]]$SingleR.single.main$scores)
      
      singler$singler[[j]]$SingleR.single$cell.names =  
        c(singler$singler[[j]]$SingleR.single$cell.names,singler.list[[i]]$singler[[j]]$SingleR.single$cell.names)
      singler$singler[[j]]$SingleR.single.main$cell.names =  
        c(singler$singler[[j]]$SingleR.single.main$cell.names,singler.list[[i]]$singler[[j]]$SingleR.single.main$cell.names)
      
      if (!is.null(singler$singler[[j]]$SingleR.single.main$pval)) {
        singler$singler[[j]]$SingleR.single.main$pval =  
          c(singler$singler[[j]]$SingleR.single.main$pval,singler.list[[i]]$singler[[j]]$SingleR.single.main$pval)
      }
      
      if (!is.null(singler$singler[[j]]$SingleR.single$pval)) {
        singler$singler[[j]]$SingleR.single$pval =  
          c(singler$singler[[j]]$SingleR.single$pval,singler.list[[i]]$singler[[j]]$SingleR.single$pval)
      }
    }
    singler$meta.data$project.name = paste(singler$meta.data$project.name,singler.list[[i]]$meta.data$project.name,sep='+')
    singler$meta.data$orig.ident = c(singler$meta.data$orig.ident,singler.list[[i]]$meta.data$orig.ident)
    singler$meta.data$data.sets = c(singler$meta.data$data.sets,rep(singler.list[[i]]$meta.data$project.name,length(singler.list[[i]]$meta.data$orig.ident)))
    
  }
  
  for (j in 1:length(singler$singler)) {
    if (!is.null(order)) {
      singler$singler[[j]]$SingleR.single$labels = 
        singler$singler[[j]]$SingleR.single$labels[order,]
      if (!is.null(singler$singler[[j]]$SingleR.single$labels1)) {
        singler$singler[[j]]$SingleR.single$labels1 = 
          singler$singler[[j]]$SingleR.single$labels1[order,]
      }
      singler$singler[[j]]$SingleR.single$scores = 
        singler$singler[[j]]$SingleR.single$scores[order,]
      singler$singler[[j]]$SingleR.single$cell.names = 
        singler$singler[[j]]$SingleR.single$cell.names[order]
      
      singler$singler[[j]]$SingleR.single.main$labels = 
        singler$singler[[j]]$SingleR.single.main$labels[order,]
      if (!is.null(singler$singler[[j]]$SingleR.single.main$labels1)) {
        singler$singler[[j]]$SingleR.single.main$labels1 = 
          singler$singler[[j]]$SingleR.single.main$labels1[order,]
      }
      singler$singler[[j]]$SingleR.single.main$scores = 
        singler$singler[[j]]$SingleR.single.main$scores[order,]
      singler$singler[[j]]$SingleR.single.main$cell.names = 
        singler$singler[[j]]$SingleR.single.main$cell.names[order]
      
      if (!is.null(singler$singler[[j]]$SingleR.single$pval)) {
        singler$singler[[j]]$SingleR.single$pval = 
          singler$singler[[j]]$SingleR.single$pval[order]
        singler$singler[[j]]$SingleR.single.main$pval = 
          singler$singler[[j]]$SingleR.single.main$pval[order]
      }
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

#' Analyze very big data sets. Runs SingleR on small chunks (10,000 cells per run) and then combines them together.
#' 
#' @param counts a tab delimited text file containing the counts matrix, a 10X directory name or a matrix with the counts.
#' @param annot a tab delimited text file or a data.frame. Rownames correspond to column names in the counts data
#' @param project.name the project name
#' @param xy a matrix with the xy coordinates. From the original single-cell object.
#' @param clusters the clusters identities.From the original single-cell object.
#' @param N number of cells in each iteration. Default is 10000.
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
#' @param temp.dir used by the SingleR webtool.
#' @param numCores Number of cores to use.
CreateBigSingleRObject = function(counts,annot=NULL,project.name,xy,clusters,N=10000,
                                  min.genes=200,technology='10X',
                                  species='Human',citation='',
                                  ref.list=list(),normalize.gene.length=F,
                                  variable.genes='de',fine.tune=T,
                                  reduce.file.size=T,do.signatures=F,
                                  do.main.types=T,
                                  temp.dir=getwd(), numCores = SingleR.numCores) {
  
  n = ncol(counts)
  s = seq(1,n,by=N)
  dir.create(paste0(temp.dir,'/singler.temp/'), showWarnings = FALSE)
  for (i in s) {
    print(i)
    A = seq(i,min(i+N-1,n))
    singler = CreateSinglerObject(counts[,A], annot = annot[A], project.name=project.name, 
                                  min.genes = min.genes,  technology = technology, 
                                  species = species, citation = citation,
                                  do.signatures = do.signatures, clusters = NULL,
                                  numCores = numCores)
    
    save(singler,file=paste0(temp.dir,'/singler.temp/',project.name,'.',i,'.RData'))
  }
  
  singler.objects.file <- list.files(paste0(temp.dir,'/singler.temp/'), 
                                     pattern='RData',full.names=T)
  
  singler.objects = list()
  for (i in 1:length(singler.objects.file)) {
    load(singler.objects.file[[i]])
    singler.objects[[i]] = singler
  }
  
  singler = SingleR.Combine(singler.objects,order = colnames(counts), 
                            clusters=clusters,xy=xy)
  
  singler
}
