library(shiny)
library(SingleR)
library(dendextend)
library(Matrix)
library(canvasXpress)
library(corrplot)
library(ggplot2)
library(matrixStats)
library(plyr)
library(reshape)
library(dplyr)
library(htmlwidgets)
library(pheatmap)
library(Seurat)


#library(psych)
#library(gplots)
#library(RColorBrewer)
#library(shinysky)
#library(rCharts)
#library(gridExtra)
#library(heatmap3)
#library(scatterD3)
#library(magrittr)

loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

analysis_sets = c('SingleR annotations','SingleR (no fine-tuning)','Meta.data','Other','Gene')

getReferenceList = function(singler) {
  unlist(lapply(singler,function(x) return(x$about$RefData)))
}

singlerDrawBoxPlot = function(sc_data,ref,cell_id,n.show = 50) {
  SingleR.DrawBoxPlot(sc_data,cell_id,ref,top.n = n.show)
  
}


singlerDrawBoxPlot2 = function(scores,ref,cell_id,n.show = 50,is.main.types=F) {
  Scores=scores[cell_id,]
  Scores = sort(Scores,decreasing = T)[1:min(n.show,length(Scores))]
  s = data.frame(Cell.Types=names(Scores),Scores)
  if (is.main.types==T) {
    s$Main_types = s$Cell.Types
  } else {
    s$Main_types = unlist(lapply(s$Cell.Types,FUN=function(x) ref$main_types[ref$types==x][[1]]))
  }
  ggplot(s, aes(x = reorder(Cell.Types, Scores), y = Scores, fill = Main_types)) + 
    geom_bar(stat = "identity")+theme(axis.text.x = element_text(angle = 45, hjust = 1))+
    coord_cartesian(ylim = c(min(Scores)*0.9, max(Scores)*1.1))+
    xlab('Note: this plot shows only the non fine-tuned scores.')+ggtitle(rownames(scores)[cell_id])
}
#singlerDrawBoxPlot = function(singler,ref,types,r,cell_id,labels.use = NULL,n.show = 50) {
# par(mar=c(10,6,4,2))
# if (!is.null(labels.use)) {
#   types.use = types %in% labels.use
# } else {
#   types.use = rep(TRUE,length(types))
# }
# 
# bymedian <- reorder(factor(types[types.use]), r, function(x) quantile(x, probs  = singler$quantile.use))
# tit = sprintf('%s: %s',names(singler$labels[cell_id,]),singler$labels[cell_id,])
# if (!is.null(labels.use)) {
#   n = length(labels.use)
# } else {
#   n = length(levels(bymedian))
# }
#  if (n>n.show) {
#     A = bymedian %in% levels(bymedian)[(n-(n.show-1)):n]
#   } else {
# A = rep(TRUE,length(r))
#   }
# names(colors) = levels(factor(ref$main_types[A]))
# main_colors = colors[levels(factor(ref$main_types[A]))]
# sub_colors = colors[unique(cbind(ref$main_types[A],types[A]))[,1]]
# names(sub_colors) = unique(cbind(ref$main_types[A],types[A]))[,2]
# bymedian = factor(bymedian[A])
# r = r[A]
# boxplot(r ~ bymedian,las=2,cex.axis=0.75,ylab="Spearman R",outer=TRUE,col=sub_colors[levels(bymedian)],main = tit,drop=TRUE)
#}


singlerDrawHeatmap = function(scores,labels = NULL,annot = NULL, cluster.by = NULL, n.show = 50) {
  if (!is.null(labels)) {
    scores = scores[labels,]
  }
  m = apply(t(scale(t(scores))),2,max)
  
  thres = sort(m,decreasing=TRUE)[min(n.show,length(m))]
  
  data = as.matrix(scores)
  mmax = rowMaxs(data)
  mmin = rowMins(data)
  data = (data-mmin)/(mmax-mmin)
  data = data^3
  data = data[,m>(thres-1e-6)]
  data = t(data)
  
  clustering_method = 'ward.D2'
  
  fontsize_row = 9
  cluster_cols = T
  if (!is.null(cluster.by)) {
    data = data[,order(cluster.by)]
    cluster_cols=F
  } 
  pheatmap(data,border_color=NA,show_colnames=FALSE,
           clustering_method=clustering_method,fontsize_row=fontsize_row,
           annotation_col = annot,cluster_cols=cluster_cols)
  
}

Identify.Diff.Markers = function(expr,cells_id1,cells_id2,test.use,min.pct,min.diff.pct) {
  min.pct = min.pct/100
  min.diff.pct = min.diff.pct/100
  thresh.min=1
  if (test.use=="Wilcoxon test") test.use='wilcox'
  else if (test.use=="Bimod test") test.use='bimod'
  else if (test.use=="ROC test") test.use='roc'
  else if (test.use=="t-Test") test.use='t'
  else if (test.use=="Tobit-censoring") test.use='tobit'
  else if (test.use=="Poisson test") test.use='poisson'
  else if (test.use=="Negative-binomial test") test.use='negbinom'
  else if (test.use=="MAST") test.use='MAST'
  
  
  cells.1 = colnames(expr)[cells_id1]
  if (sum(cells_id2)==0) {
    cells.2 = colnames(expr)[!cells_id1]
  } else {
    cells.2 = colnames(expr)[cells_id2 & !cells_id1]
  }
  
  seurat = CreateSeuratObject(expr[,c(cells.1,cells.2)])
  
  markers = FindMarkers(seurat, ident.1 = cells.1,ident.2 = cells.2,test.use=test.use,min.pct=min.pct,min.diff.pct=min.diff.pct)
  
  dm = list()
  dm$data = markers
  dm$cells1 = cells.1
  dm$cells2 = cells.2
  #dm$cells_id1 = cells_id1
  #dm$cells_id2 = cells_id2
  
  return(dm)
}

Make.HeatMap =function(expr,dm,n_genes,genes.use=NULL) {
  print('make.heatmap')
  disp.min=-2.5
  disp.max=2.5
  if (is.null(genes.use)) {
    dm$data = dm$data[!is.na(dm$data[,1]),]
    pos = dm$data[dm$data[,2]>0,]
    neg = dm$data[dm$data[,2]<0,]
    thres = sort(pos[,1])[min(dim(pos)[1],n_genes+1)]
    genes.use.pos=rownames(pos)[pos[,1]<thres]
    thres = sort(neg[,1])[min(dim(pos)[1],n_genes+1)]
    
    genes.use.neg=rownames(neg)[neg[,1]<thres]
    genes.use = c(genes.use.pos,genes.use.neg)
    genes.use = genes.use[!is.na(genes.use)]
    
    cells.use=c(dm$cells1,dm$cells2)
  } else {
    cells.use = 1:length(seurat@cell.names)
    cells.use = cells.use[order(seurat@ident)]
    cells.ident = seurat@ident[cells.use]
    names(cells.ident) = seurat@cell.names[cells.use]
  }
  
  
  data.use=t(scale(t(expr[genes.use,cells.use])))
  data.use=MinMax(data.use,min=disp.min,max=disp.max)
  print(length(cells.use))
  print(length(c(rep('Set1',length(dm$cells1)),rep('Set2',length(dm$cells2)))))
  cells.ident = data.frame(row.names=cells.use,Group=c(rep('Set1',length(dm$cells1)),rep('Set2',length(dm$cells2))))
  pheatmap(as.matrix(data.use),cluster_rows=FALSE,cluster_cols = FALSE,show_colnames = FALSE,
           annotation_col = cells.ident,fontsize = 12,fontsize_row = (90-n_genes)/5)
  
  #return(data.use)
  
}

FineTuningStep <- function(singler,expr,ref,types,topLabels,cell_id) {
  fine.tune.thres = 0.05
  if(length(topLabels)>1) {
    labels.use = types %in% topLabels
    types.filtered = types[labels.use]
    sd =  rowSds(as.matrix(ref$data[,labels.use]))
    thres = min(sort(sd,decreasing = TRUE)[250],ref$sd.thres)
    genes.filtered=intersect(rownames(expr)[sd>=thres],rownames(ref$data))
    genes.filtered = genes.filtered[!is.na(genes.filtered)]
    
    data = as.matrix(expr[genes.filtered,cell_id])
    out = list()
    if (sd(data)>0) {
      out$r=cor(data,ref$data[genes.filtered,],method="spearman")
      out$r = out$r[labels.use]
      agg_scores = aggregate(out$r~types.filtered,FUN = quantile, probs  = singler$quantile.use)
      max_score = max(agg_scores[,-1])
      agg_scores = agg_scores[-which.min(agg_scores[,-1]),]
      topLabels = agg_scores[agg_scores[,-1]>=max_score-fine.tune.thres,1]
    } else {
      topLabels = agg_scores[which.max(agg_scores[,-1]),1]
    }
  }
  out$labels=topLabels
  
  return(out)
}

createScatterPlot = function(obj=obj,input=input,gene.to.use) {
  xy = data.frame(x=obj$singler@xy[,1],y=obj$singler@xy[,2])
  col_continuous = FALSE
  colorSpectrum = c('#4575b4', '#91bfdb', '#e0f3f8', '#ffffbf', '#fee090', '#fc8d59', '#d73027')
  if (input$data2 == 'SingleR annotations') {
    if (input$by_clusters==T && !empty(obj$singler@labels.clusters)) {
      labels = clusters.map.values(obj$singler@clusters[,input$data_clusters],obj$singler@labels.clusters[,input$ref_data,drop=F])
    } else {
      labels = obj$singler@labels[,input$ref_data]
    }
    tit = paste(obj$singler@project.name,"SingleR fine-tuned annotations")
  } else if (input$data2 == 'SingleR (no fine-tuning)') {
    if (input$by_clusters==T && !empty(obj$singler@labels.clusters)) {
      labels = clusters.map.values(obj$singler@clusters[,input$data_clusters],obj$singler@labels.clusters.NFT[,input$ref_data,drop=F])
    } else {
      labels = obj$singler@labels.NFT[,input$ref_data]
    }
    tit = paste(obj$singler@project.name,"SingleR non fine-tuned annotations")
  } else if (input$data2 == 'Meta.data') {
    if (input$data7 %in% colnames(obj$singler@clusters)) {
      labels=obj$singler@clusters[,input$data7]
    } else {
      labels=obj$singler@ident[,input$data7]
    }
    tit = paste(paste(obj$singler@project.name),input$data7)
  } else if (input$data2 == 'Other') {
    if (input$by_clusters==T) {
      agg = aggregate(obj$singler@other[,input$data8],by=obj$singler@clusters,mean)
      labels = mapvalues(obj$singler@clusters[,input$data_clusters],from=agg[,1],to=agg[,2])
    } else {
      labels=obj$singler@other[,input$data8]
    }
    col_continuous = TRUE
    tit = paste(paste(obj$singler@project.name),input$data8)
  } else if (input$data2 == 'Gene') {
    if (is.null(input$data6)) {
      g = as.matrix(obj$singler@expr[obj$gene,])
      tit = paste(paste(obj$singler@project.name),'-',as.character(input$data6),'expression')
      
    } else if (length(input$data6)==1) {
      g = as.matrix(obj$singler@expr[as.character(input$data6[1]),])
      obj$gene = as.character(input$data6[1])
      tit = paste(paste(obj$singler@project.name),'-',as.character(input$data6),'expression')
      
    } else {
      g = as.matrix(obj$singler@expr[as.character(input$data6),])
      g = colSums(g>0)
      tit = paste0(input$data6,collapse='+')
      tit = paste(obj$singler@project.name,'-',tit,'expression')
      
      obj$gene = as.character(input$data6)
    }
    if (input$by_clusters==T) {
      agg = aggregate(g,by=list(obj$singler@clusters[,input$data_clusters]),mean)
      labels = mapvalues(obj$singler@clusters[,input$data_clusters],from=agg[,1],to=agg[,2])
    } else {
      labels = g
    }
    col_continuous = TRUE
  }
  obj$labels = labels
  if (col_continuous==F & length(unique(labels))>39) {
    tbl = table(labels)
    thres = sort(tbl,decreasing = T)[40]
    labels = as.character(labels)
    labels[labels %in% names(tbl)[tbl<=thres]] = paste0('zOther (N<=',thres,')')
  }
  annot = data.frame(labels,obj$singler@clusters[,1],obj$singler@ident[,1])
  colnames(annot) = c('labels','cluster','ident')
  rownames(annot) = rownames(xy)
  dot.size = input$dot.size
  xy$labels=as.character(labels)
  events = htmlwidgets::JS("{ 'mousemove' : function(o, e, t) {
                                                if (o.objectType == null) {
                                                  t.showInfoSpan(e, '<b>' + o.y.vars[0] + '</b><br/>' +
                                                  '<b>Value:</b> ' + o.z.labels + '<br/>' +
                                                  '<b>Cluster:</b> ' + o.z.cluster + '<br/>' +
                                                  '<b>Identity:</b> ' + o.z.ident);
                                                }
                                                else {
                                                  t.showInfoSpan(e, o.display);
                                                };
                                   },",
              "'select': function(o, e, t){var g = $('#selected');g.val(o.y.vars);g.trigger('change')}",
              "}")
  
  # "'click': function(o, e, t){var g = $('#selected');g.val(o.y.vars[0]);g.trigger('change')},",
  
  stringVariableFactors = list()
  if (col_continuous == F) {
    stringVariableFactors = list('labels')
  }
  if (col_continuous==T) {
    colorSpectrum=c('#DCDCDC','#4575b4')
  }
  setMaxX = max(xy[,1])
  setMinX = min(xy[,1])
  setMaxY = max(xy[,2])
  setMinY = min(xy[,2])
  
  canvasXpress(xy,varAnnot=annot,graphType='Scatter2D',colorSpectrum=colorSpectrum,
               colorBy='labels',dataPointSize=dot.size,outlineWidth=dot.size/40,
               sizeByShowLegend=F,xAxisMinorTicks=F,yAxisMinorTicks=F,
               xAxisMajorTicks=F,yAxisMajorTicks=F,stringVariableFactors=stringVariableFactors,
               title=tit,titleScaleFontFactor=0.5,xAxisTitle='tSNE1',yAxisTitle='tSNE2',
               events=events,showLegendTitle=F,axisAlgorithm="wilkinsonExtended",plotBox=FALSE,
               showDecorations=TRUE,
               setMaxX=setMaxX,setMinX=setMinX,setMaxY=setMaxY,setMinY=setMinY)
}

save_pheatmap_pdf <- function(x, filename, width=7, height=7) {
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
