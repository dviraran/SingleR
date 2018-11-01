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
    ggtitle(paste('R =', round(1000*cor(df$x,df$y,method='spearman',use='pairwise'))/1000)) + 
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
                               fontsize_row=9,...) {
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
             clustering_method=clustering_method,fontsize_row=fontsize_row,
             annotation_col = clusters,cluster_cols = F,silent=silent, 
             annotation_colors=annotation_colors)
  } else if (!is.null(cells_order)) {
    data = data[,cells_order]
    clusters = clusters[cells_order,,drop=F]
    pheatmap(data,border_color=NA,show_colnames=FALSE,
             clustering_method=clustering_method,fontsize_row=fontsize_row,
             annotation_col = clusters,cluster_cols = F,silent=silent, 
             annotation_colors=annotation_colors)
  } else {
    if (!is.null(clusters)) {
      pheatmap(data,border_color=NA,show_colnames=FALSE,
               clustering_method=clustering_method,fontsize_row=fontsize_row,
               annotation_col = clusters,silent=silent, 
               annotation_colors=annotation_colors)
    } else {
      pheatmap(data[,sample(ncol(data))],border_color=NA,show_colnames=FALSE,
               clustering_method=clustering_method,fontsize_row=fontsize_row,
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
  
  SYMBOLS = c(LETTERS,letters,c(0:9),
              c("!","@","#","$","%","^","&","*","(",")",")","-",
                "+","_","=",";","/","|","{","}","~"))

  if (score.thres>0) {
    max.score = apply(SingleR$scores,1,max)  
    df$ident[max.score<score.thres] = 'X'
  }
  
  if (do.letters==T && length(unique(df$ident))>length(SYMBOLS)) {
    n = table(df$ident)
    thres = sort(n,decreasing = T)[length(SYMBOLS)]
    df$ident[df$ident %in% names(n)[n <= thres]] = 'X'
  }
  
  df$ident = factor(df$ident)

  num.levels = length(levels(df$ident))
  
  df$initIdent = SYMBOLS[as.numeric(df$ident)]
  
  
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
    p = p + guides(shape=guide_legend(override.aes = list(size=3)))
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
                    legend.text=element_text(size=6),legend.title = element_blank()) + guides(shape=guide_legend(ncol=5,override.aes = list(size=2.5,alpha=1)))
    } else if (num.levels>60){
      p = p + theme(legend.position="bottom",legend.direction="vertical",
                    legend.text=element_text(size=6),legend.title = element_blank()) + guides(shape=guide_legend(ncol=9,override.aes = list(size=2,alpha=1)))
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