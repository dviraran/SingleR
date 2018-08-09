library(RColorBrewer)
library(reshape2)
library(SingleR)
library(pheatmap)
library(ggplot2)
library(ggpubr)
library(grid)
library(Seurat)
library(corrplot)
library(DeconRNASeq)
library(RColorBrewer)


path.data = '~/Documents/SingleR/package/SingleR/manuscript_figures/FiguresData/'
path.out = '~/Documents/SingleR/package/SingleR/manuscript_figures/FiguresOut'

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

### Figure 1b - Fibroblasts/BMDC example

load (file.path(path.data,'GSE78779.RData'))
xy = singler$meta.data$xy
a = xy[,1]
xy[,1]=-xy[,2]
xy[,2]= a

# Original identities
pdf(file.path(path.out,'Fig1b_orig.ident.pdf'),width=5,height=3.75)
out = SingleR.PlotTsne(singler$singler[[1]]$SingleR.single,xy,do.label=FALSE,
                       do.letters =F,labels=singler$meta.data$orig.ident, 
                       dot.size = 3)
out$p
dev.off()

# SingleR annotations
cl.set2 =brewer.pal(8,'Set2')
cl.set2[3] = cl.set2[4];cl.set2[4] = cl.set2[1];cl.set2[1] = cl.set2[2];cl.set2[2] = cl.set2[6]

pdf(file.path(path.out,'Fig1b_singler.annot.pdf'),width=5.33,height=3.75)
out = SingleR.PlotTsne(singler$singler[[1]]$SingleR.single.main,do.letters = F,colors=cl.set2,
                       xy,do.label=FALSE, dot.size = 3)
out$p
dev.off()

# SingleR heatmap
pdf(file.path(path.out,'Fig1b_singler.heatmap.pdf'),width=4,height=3)
SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single.main,top.n=8)
dev.off()

# Validation of GM-CSF macrophages
BMDC = singler$meta.data$orig.ident=='BMDCs'

counts.file = file.path(path.data,'GSE78779_series_matrix.txt')
counts <- as.matrix(read.table(counts.file, header=TRUE, sep="\t", row.names=1, as.is=TRUE,comment.char='!'))
M = apply(counts,2,function(x) 1e6*x/sum(x)) #normalize
rownames(M) = unlist(lapply(tolower(rownames(M)),capitalize))
gse62631 <- read.table(file.path(path.data,'GSE62361.txt'), header=TRUE, sep="\t", row.names=1, as.is=TRUE,comment.char='!')

A = intersect(gse62631$Gene.symbol,rownames(M))
gse62631 = gse62631[gse62631$Gene.symbol %in% A,]
thres = sort(gse62631$logFC)[c(51,length(gse62631$logFC)-50)]
A = gse62631$logFC < thres[1] | gse62631$logFC > thres[2]
genes.use = data.frame(row.names=gse62631$Gene.symbol[A],logFC=gse62631$logFC[A])
genes.use$Group = 'GM-Macs'
genes.use$Group[genes.use$logFC<0] = 'GM-DCs'
genes.use$logFC=c()
genes.use$Group= factor(genes.use$Group)
a = t(scale(t(M[rownames(genes.use),BMDC])))
a[a>2]=2;a[a< -2]=-2

annotation_col = data.frame(Annotation=singler$singler[[1]]$SingleR.single.main$labels[BMDC,])

pdf(file.path(path.out,'Fig1b_helft.validation.pdf'),width=7,height=4.5)
pheatmap(a[order(genes.use$Group),order(singler$singler[[1]]$SingleR.single.main$labels[BMDC])],
         cluster_cols = F,cluster_rows = F,clustering_method='ward.D2',border_color = NA,
         annotation_col=annotation_col,
         annotation_row = genes.use,show_colnames=F,show_rownames = F,
         annotation_colors=list(Annotation=c('DC'=cl.set2[1],'Macrophages'=cl.set2[3])),
         fontsize = 14)
dev.off()

### Figure 1c
load(file.path(path.data,'SingleR.PBMC.3K.4K.RData'))
singler = singler.4k

pdf(file.path(path.out,'Fig1c-seurat-clusters.pdf'),width=4.63,height=3.75)
out = SingleR.PlotTsne(singler$singler[[2]]$SingleR.single.main, singler$meta.data$xy, 
                       do.letters = F,labels = singler$seurat@ident,
                       dot.size = 0.7,colors = singler.colors)
out$p
dev.off()

labels = singler$singler[[2]]$SingleR.single$labels
A = labels=='Monocytes' & grep('Monocyte',singler$singler[[1]]$SingleR.single$labels)
labels[A] = singler$singler[[1]]$SingleR.single$labels[A]

pdf(file.path(path.out,'Fig1c-singler.pdf'),width=5.13,height=3.75)
out = SingleR.PlotTsne(singler$singler[[2]]$SingleR.single.main, singler$meta.data$xy, 
                       do.letters = F,labels = labels,
                       dot.size = 0.7,colors = singler.colors)
out$p
dev.off()
pdf(file.path(path.out,'Fig1c-singler-legend.pdf'),width=5.13,height=3.75)
out = SingleR.PlotTsne(singler$singler[[2]]$SingleR.single.main, singler$meta.data$xy, 
                       do.letters = F,
                       dot.size = 1.5,colors = singler.colors)
out$p
dev.off()

s = SingleR.Subset(singler,grepl('(CD4|CD8|Treg)',singler$singler[[2]]$SingleR.single$labels))

pdf(file.path(path.out,'Fig1c-singler-tcells.pdf'),width=5.13,height=3.75)
out = SingleR.PlotTsne(s$singler[[2]]$SingleR.single, s$meta.data$xy, 
                       do.letters = F,
                       dot.size = 1.5,colors = singler.colors)
out$p
dev.off()
genes.use = c('CD4','CD8A','CCR7','SELL','CXCR5','IL7R','GNLY', 'NKG7','FOXP3','CTLA4')
agg = aggregate(t(as.matrix(s$seurat@data[genes.use,])),by=list(s$singler[[2]]$SingleR.single$labels),mean)
rownames(agg) = agg[,1]; agg = agg[,-1]
pdf(file.path(path.out,'Fig1c-tcells-markers.pdf'),width=2.5,height=5)
corrplot(t(scale(as.matrix(agg))),is.corr=F,addgrid.col = NA,tl.cex=0.7,tl.col='black',cl.cex=0.00001)
dev.off()

s = SingleR.Subset(singler,grepl('B-cells',singler$singler[[2]]$SingleR.single$labels))
pdf(file.path(path.out,'Fig1c-singler-bcells.pdf'),width=5.13,height=3.75)
out = SingleR.PlotTsne(s$singler[[2]]$SingleR.single, s$meta.data$xy, 
                       do.letters = F,
                       dot.size = 1.5,colors = singler.colors)
out$p
dev.off()
K = as.numeric(factor(grepl('naive',s$singler[[2]]$SingleR.single$labels)))
K[K==1] = 'Memory'
K[K==2] = 'Naive'

df = data.frame(CD27 = s$seurat@data['CD27',],Type=K)
pdf(file.path(path.out,'Fig1c-bcells-cd27.pdf'),width=2.5,2)
ggplot(df,aes(x=Type,y=CD27,color=Type,fill=Type))+geom_violin(alpha=0.2)+
  geom_point(size=0.5,alpha=0.5,position=position_jitterdodge())+
  ylab('CD27 expression')+xlab("")+scale_color_manual(values=singler.colors[2:3])+scale_fill_manual(values=singler.colors[2:3])+
  theme(axis.title.x=element_blank())+ theme(legend.position="none")
dev.off()

### Figure 2a - Case/Control of lung sc data
load(file.path(path.data,'GSE111664.RData'))

singler.macs = SingleR.Subset(singler,singler$seurat@ident %in% c(0,1,2,7))

full.cond = as.character(singler$meta.data$orig.ident)
full.cond[grepl('ontrol',full.cond)]='Control'
full.cond[grepl('leo',full.cond)]='Bleomycin'

p = SingleR.PlotTsne(singler$singler[[1]]$SingleR.single,singler$meta.data$xy,labels=full.cond,
                     do.letters=F,colors=brewer.pal(3,'Set2'),do.labels = F,dot.size=1,alpha=0.35)
pdf(file.path(path.out,'Fig2a.pdf'),width=3.2,height=2.4)
p$p+theme_void()+theme(legend.position="none")
dev.off()

full.cond2 = as.character(singler$meta.data$orig.ident)
full.cond2[grepl('ontrol',full.cond)]='Healthy lungs'
full.cond2[grepl('leo',full.cond)]='Bleomycin-injured lungs'

p = SingleR.PlotTsne(singler$singler[[1]]$SingleR.single,singler$meta.data$xy,labels=full.cond2,do.letters=F,
                     colors=brewer.pal(3,'Set2'),do.labels = F,dot.size=3,alpha=0.35)
pdf(file.path(path.out,'Fig2a_legend.pdf'),width=3.5,height=3.5)
p$p+theme_void()
dev.off()

### Figure 2b - SingleR annotations

p = SingleR.PlotTsne(singler$singler[[1]]$SingleR.single.main,singler$meta.data$xy,do.letters=F,
                     colors=sample(c(brewer.pal(12,'Set3'),brewer.pal(8,'Set2'))),do.labels = F,
                     dot.size=1,alpha=0.35)
pdf(file.path(path.out,'Fig2b-la.pdf'),width=3.2,height=2.4)
p$p+theme_void()+theme(legend.position="none")
dev.off()

### Figure 2c - SingleR annotations with lung myeloid data

labels = singler.macs$singler[[2]]$SingleR.single$labels
labels[grepl('(CD8-)',labels)] = 'pDC (CD8-)'
labels[labels=='pDC (CD8+ Flu d3)'] = 'pDC (CD8+)'
labels[grepl('CD103-',labels)] = 'cDC (CD11b+)'
labels[grepl('(CD103+)',labels)] = 'cDC (CD11b-)'
labels[labels=='AM'] = 'Alveolar Macs'
labels[labels=='IM1'] = 'Interstitial Macs (CD11c-)'
labels[labels=='IM2'] = 'Interstitial Macs (CD11c-)'
labels[labels=='IM3'] = 'Interstitial Macs (CD11c+)'

p = SingleR.PlotTsne(singler.macs$singler[[2]]$SingleR.single,singler.macs$meta.data$xy,labels=labels,do.letters=F,
                     colors=c(brewer.pal(8,'Accent'),brewer.pal(3,'Set1')),do.labels = F,dot.size=0.9,alpha=0.5)
pdf(file.path(path.out,'Fig2c.pdf'),width=2.5,height=3)
p$p+theme_void()+theme(legend.position="none")
dev.off()

p = SingleR.PlotTsne(singler.macs$singler[[3]]$SingleR.single,singler.macs$meta.data$xy,labels=labels,do.letters=F,
                     colors=c(brewer.pal(8,'Accent'),brewer.pal(3,'Set1')),do.labels = F,dot.size=3,alpha=0.5)
pdf(file.path(path.out,'Fig2c_legend.pdf'),width=3.5,height=4.3)
p$p+theme_void()
dev.off()

### Figure 2d - deconvolution

# create decon data
load(file.path(path.data,'GSE94135_GSE49932.RData'))
amim = cbind(rowMeans(gse94135_gse49932$data[,1:3]),rowMeans(gse94135_gse49932$data[,10:12]))
expr = singler.macs$seurat@data
decon = DeconRNASeq(as.data.frame(as.matrix(expr)),as.data.frame(amim),checksig=FALSE,known.prop = FALSE)

#load('~/Documents/singler/decon.rds')

df = data.frame(x=singler.macs$meta.data$xy[,1],y=singler.macs$meta.data$xy[,2])

df$Scale = decon$out.all[,1]
cl = brewer.pal(3,'Set1')
pdf(file.path(path.out,'Fig2d.pdf'),width=2.5,height=3)
ggplot(df) + geom_point(aes(x=x, y=y,color=Scale),size=1.2)+
  scale_color_gradient(low=cl[1],high=cl[2])+theme_void()+theme(legend.position="none")
dev.off()

### Figure 3a - clustering

col = brewer.pal(8,'Set2')
p = SingleR.Cluster(singler.macs$singler[[1]]$SingleR.single,3)

K = p$cl
K = plyr::mapvalues(K,from=1:3,to=c('C1','C2','C3'))

singler.macs$seurat@ident = K

pdf(file.path(path.out,'Fig3a.pdf'),width=2.5,height=3)
ggplot(df) + geom_point(aes(x=x, y=y,color=K),size=0.7)+theme_void()+theme(legend.position="none")
dev.off()

library(dendextend)
dend = as.dendrogram(p$hc)
d = cut(dend,h=1)$upper
n = length(labels(d))
labels(d) = rep('',n)
pdf(file.path(path.out,'Fig3a-dendrogram.pdf'),width=3,height=1)
par(mar=c(0,0,0,0))
plot(d,axes=F)
dev.off()

### Figure 2d - deconvolution +

df = data.frame(x=singler.macs$seurat@dr$tsne@cell.embeddings[,1],
                y=singler.macs$seurat@dr$tsne@cell.embeddings[,2],Cluster=K)
df$Scale = decon$out.all[,1]
df$PCA = prcomp(cbind(df$x,df$y))$x[,1]
cl = brewer.pal(3,'Set1')
pdf(file.path(path.out,'Fig2d-order.pdf'),width=3.5,height=2.8)
ggplot(df) + geom_point(aes(y=Scale,x=PCA,color=Scale),size=0.5,alpha=0.5)+
  scale_color_gradient(low=cl[1],high=cl[2])+xlab('t-SNE ordering')+ylab('IM-AM similarity')+
  geom_vline(xintercept = 0,linetype = "dashed",color='grey')+geom_vline(xintercept = 25,linetype = "dashed",color='grey')
dev.off()

### Figure 3a, revisited

cond = as.character(singler.macs$meta.data$orig.ident)
cond[grepl('ontrol',cond)]='Control'
cond[grepl('leo',cond)]='Bleomycin'
singler.macs$meta.data$orig.ident=cond

annotation_col = data.frame(Cluster=K,Condition=factor(cond))
rownames(annotation_col) = singler.macs$seurat@cell.names
cluster_colors = gg_color_hue(3)
cond_color = brewer.pal(3, 'Set2')
annotation_colors = list(Condition=c(Bleomycin=cond_color[1],Control=cond_color[2]),
                         Clusters=c('C1'=cluster_colors[1],'C2'=cluster_colors[2],'C3'=cluster_colors[3]))

cond.cluster = table(annotation_col$Condition,K)/cbind(table(full.cond),table(full.cond),table(full.cond))
cond.cluster = melt(cond.cluster)

pdf(file.path(path.out,'Fig3a-heatmap.pdf'),width=6,height=3)
SingleR.DrawHeatmap(singler.macs$singler[[1]]$SingleR.single,top.n=20,
                    clusters=K,order.by.clusters = T,annotation_colors=annotation_colors)
dev.off()

pdf(file.path(path.out,'Fig3a_insert.pdf'),width=3,height=2.5)
ggplot(cond.cluster,aes(K,value*100)) + geom_bar(stat = "identity", aes(fill =Var1), position = "dodge")+
  scale_fill_manual(values=cond_color[1:2])+xlab("")+ylab('% of cells')
dev.off()

# Figure 3b - de

markers <- FindMarkers(object = singler.macs$seurat, ident.1 = 'C1',ident.2 = 'C3', min.pct = 0.33,logfc.threshold = 1.5)

am = as.matrix(singler.macs$seurat@data[rownames(markers[markers$avg_logFC>0,]),])
im = as.matrix(singler.macs$seurat@data[rownames(markers[markers$avg_logFC<0,]),])
c2.ct.am = rowSums(am[,K=='C2']>0)
c2.ct.im = rowSums(im[,K=='C2']>0)

am = am[order(c2.ct.am),]
im = im[order(-c2.ct.im),]
d = t(scale(t(as.matrix(rbind(am,im)))))
d[d>2] = 2;d[d< -2]=-2

num.genes = rbind(colSums(am>0)/nrow(am),colSums(im>0)/nrow(im))

pdf(file.path(path.out,'Fig3b.pdf'),width=6.5,height=6)
pheatmap(d[,order(annotation_col$Cluster,num.genes[2,]-num.genes[1,])],
         annotation_col=annotation_col,annotation_colors=annotation_colors,cluster_cols=F,
         cluster_rows=F,show_colnames=F,border_color=NA,fontsize=9)
dev.off()

### Figure 3c - expressed genes

rownames(num.genes) = c('C1 genes','C3 genes')
res = c()
for (i in seq(0,0.85,by=0.01)) {
  res = rbind(res,(table(num.genes[1,]>i & num.genes[2,]>i,annotation_col$Cluster)/rbind(table(annotation_col$Cluster),table(annotation_col$Cluster)))[2,])
}
colnames(res) = c('C1','C2','C3')
f = melt(res)
f$value = f$value*100
pdf(file.path(path.out,'Fig3c.pdf'),width=4,height=4)
ggplot(f)+geom_smooth(aes(x=Var1,y=value,color=Var2),size=2)+ylab('% of cells')+
  xlab('% of expressed genes from C1 & C3')+
  theme(text = element_text(size=12),axis.text.x = element_text(size=12),axis.text.y = element_text(size=12))
dev.off()

df = data.frame(x=singler.macs$meta.data$xy[,1],y=singler.macs$meta.data$xy[,2])
for (i in 1:nrow(markers)) {
  df$Scale = as.matrix(singler.macs$seurat@data[rownames(markers)[i],])
  ggplot(df) + geom_point(aes(x=x, y=y,color=Scale),size=0.3,alpha=0.5)+
    scale_color_gradient(low='gray',high='darkblue')+ggtitle(rownames(markers)[i])+theme_void()+
    theme(title =element_text(size=10, face='bold'),plot.title = element_text(hjust = 0.5,vjust=-0.5))
  ggsave(file.path(path.out,paste0('Fig3b_',rownames(markers)[i],'.pdf')),width=2,height=2)
}

### Figure 3d - bulk validation

de=read.table(file.path(path.data,'BulkRNAseq_SiglecF+Cd11c+.txt'),sep="\t",header=TRUE,row.names=1, as.is=TRUE)

M = de[,13:26]
de$Gene[is.na(de$Gene)] = 'NA'
rownames(M) = make.unique(de$Gene)

rownames(M) = unlist(lapply(tolower(rownames(M)),capitalize))
genes = intersect(rownames(M),rownames(markers)[markers$avg_logFC<0])
m = M[genes,c(13,14,2,4,6,1,3,5,8,10,12,7,9,11)]
#c(14) = CTRL_H_3.norm
annotation_col = data.frame(c('Baseline','Baseline','2 wk MHC lo','2 wk MHC lo','2 wk MHC lo',
                              '2 wk MHC high','2 wk MHC high','2 wk MHC high',
                              '4 wk MHC lo','4 wk MHC lo','4 wk MHC lo',
                              '4 wk MHC high','4 wk MHC high','4 wk MHC high'))
rownames(annotation_col) = colnames(m)
colnames(annotation_col) = 'Group'

m = t(scale(t(m)))
m[m>2] = 2; m[m< -2] = -2 
pdf(file.path(path.out,'Fig3d-bulk.pdf'),width=6,height=3.5)
pheatmap(m,scale='row',cluster_cols = F,annotation_col=annotation_col,show_colnames=F,cluster_rows = T,
         border_color=NA,clustering_method='ward.D2')
dev.off()

### Figure 6a - pdgfa in bulk

df = data.frame(Pdgfa = t(M['Pdgfa',c(13,14,2,4,6,1,3,5,8,10,12,7,9,11)]), Group = annotation_col$Group)
df$Group = factor(df$Group,levels=c('Baseline','2 wk MHC lo','2 wk MHC high','4 wk MHC lo','4 wk MHC high'))
pdf(file.path(path.out,'Fig6a-bulk-pdgfa.pdf'),width=3.5,height=2.7)
ggplot(df,aes(x=Group,y=Pdgfa,color=Group,fill=Group))+
  #geom_boxplot(alpha=0.2)+
  geom_point(size=2,alpha=0.5,position=position_jitterdodge(jitter.width=1.5))+
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
               geom = "crossbar", width = 0.5)+
  stat_compare_means(label.x.npc=0.4,label='p.signif',comparisons=list(c('2 wk MHC lo','2 wk MHC high')))+
  theme(axis.title.x=element_blank())+theme(legend.position="none")+
  ylab('Pdgfa normalized expression')
dev.off()


## Fig 4a - telome dysfunction
spc = read.table('~/Documents/SingleR/spc-trf1.txt',header=T,row.names = 1,as.is=T,sep='\t')
spc.nodup = createMatrixNoDups(spc)
types= c('CTRL','CTRL','CTRL','Disease-3M','Disease-3M','Disease-3M','Disease-9M','Disease-9M','Disease-9M')

sig = list()
sig[[1]] = GeneSet(rownames(am), setName='C1')
sig[[2]] = GeneSet(rownames(im), setName='C3')
egc <- GeneSetCollection(sig)
scores = gsva(as.matrix(spc.nodup),egc,method='ssgsea',kcdf='Poisson')

df = data.frame(C3=scale(scores[2,]),Mouse=types)
pdf(file.path(path.out,'Fig4a.other-mouse-model-C3.pdf'),width=2.6,height=4)
ggplot(df,aes(x=Mouse,y=C3,color=Mouse,fill=Mouse))+
  #geom_boxplot(alpha=0.2)+
  geom_point(size=5,alpha=0.5,position=position_jitterdodge(jitter.width=1.5))+
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
               geom = "crossbar", width = 0.5)+
  theme(axis.title.x=element_blank())+theme(legend.position="none")+
  ylab('C3 score')
dev.off()

annotation_col = data.frame(row.names = colnames(spc.nodup),Categories=types)
d = t(scale(t(spc.nodup[intersect(rownames(spc.nodup),rownames(im)),])))
d[d>2] = 2; d[d< - 2] = -2

pdf(file.path(path.out,'Fig4a.other-mouse-model-C3genes.pdf'),width=6,height=3)
pheatmap(log2(spc.nodup[intersect(rownames(spc.nodup),rownames(im)),]),scale = 'row',cluster_rows = T,cluster_cols = F,annotation_col=annotation_col,show_colnames = F,border_color = NA)
dev.off()

### Figure 4b - human gsea ipf

library(GSEABase)
library(GSVA)
library(biomaRt)

load(file.path(path.data,'GSE32537_expr.RData'))
clinical = as.matrix(read.table(file.path(path.data,'GSE32537_clinical.txt'), header=TRUE, sep="\t", row.names=1, as.is=TRUE,comment.char='!'))
IPF = clinical[,'final.diagnosis']=='control'

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                 values = rownames(am) , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
human.am <- unique(genesV2[, 2])
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                 values = rownames(im) , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
human.im <- intersect(rownames(expr),c(unique(genesV2[, 2]),'CLEC10A'))
human.mhcii = rownames(expr)[which(grepl('HLA-D',rownames(expr)))]

human.am = c('CTSK', 'CHIA', 'CDC42EP3', 'ATP6V0D2', 'FABP4', 'RNASE3', 'RNASE2', 'MARCO', 'CA4', 'PLET1', 'IL18', 'FABP1', 'KRT19')
human.im = c('SRGN','S100A6','S100A4','IL1B','C1QA','C1QB','C1QC','CCR2','CXCL16','SOCS3','MAFB','MMP12','CX3CR1','AIF1','PLBD1','CLEC10A')

sig = list()
sig[[1]] = GeneSet(human.am, setName='C1')
sig[[2]] = GeneSet(human.im, setName='C3')
sig[[3]] = GeneSet(human.mhcii, setName='MHCII')

egc <- GeneSetCollection(sig)

scores = gsva(expr,egc,method='ssgsea',kcdf='Gaussian')

Group = clinical[,'final.diagnosis']!='control'
Group[Group==F] = 'Control'
Group[Group==T] = 'Fibrosis'


df = data.frame('C1 genes'=scale(scores[1,]),'MHCII genes'=scale(scores[3,]),
                Group=Group)  
colnames(df) <- gsub("\\.", " ", colnames(df))
colnames(df)
d = melt(df,variable.name='Group')
colnames(d) = c('Group','Set','Score')
pdf(file.path(path.out,'Fig4b-human.ipf.pdf'),width=4.2,height=3)
ggplot(d,aes(x=Group,y=Score,color=Group,fill=Group))+geom_boxplot(alpha=0.2,outlier.shape = NA)+
  geom_point(size=0.3,alpha=0.5,position=position_jitterdodge())+
  facet_wrap(~Set,nrow=1)+ 
  stat_compare_means(label.x.npc=0.4,size=8,label='p.signif')+
  ylab('Scaled scores')+xlab("")+
  theme(axis.title.x=element_blank(),axis.ticks.x = element_blank())+ theme(legend.position="none")
dev.off()

### Figure 5c - human ipf cxcr31

df = data.frame(CX3CR1=(expr['CX3CR1',]),
                Group=Group)  
colnames(df) <- gsub("\\.", " ", colnames(df))
colnames(df)
d = melt(df,variable.name='Group')
colnames(d) = c('Group','Set','Score')
pdf(file.path(path.out,'Fig5c-cx3cr1.pdf'),width=2.9,height=3)
ggplot(d,aes(x=Group,y=Score,color=Group,fill=Group))+geom_boxplot(alpha=0.2,outlier.shape = NA)+
  geom_point(size=0.5,alpha=0.5,position=position_jitterdodge())+
  stat_compare_means(label.x.npc=0.4,size=8,label='p.signif')+
  ylab('CX3CR1 expression')+xlab("")+
  theme(axis.title.x=element_blank())+ theme(legend.position="none")
dev.off()


### Figure 5b - qunatitation of tdtomato

df = read.table(file.path(path.data,'fig5b-quantitation.txt'), header=TRUE, sep="\t", row.names=NULL, as.is=TRUE)
colnames(df) = c('SiglecF+Tdtomato-','SiglecF+Tdotomato+')
df = melt(df)
pdf(file.path(path.out,'Fig5b-tdtomata.pdf'),width=1.5,height=2.8)
ggplot(df,aes(x=variable,y=value,fill=variable,color=variable))+
  geom_boxplot(alpha=0.2,outlier.shape = NA)+
  geom_point(size=1,alpha=0.5,position=position_jitter(width = 0.2))+
  scale_color_manual(values=c("#00BA38","#E58700"))+
  scale_fill_manual(values=c("#00BA38","#E58700"))+
  #stat_compare_means(label.x.npc=0.4,label='p.signif')+
  ylab('% of cells within fibroblast clusters')+xlab('')+
  theme(legend.position="none",text = element_text(size=10),axis.text.x = element_text(size=8))
dev.off()

### Figure 4c - human ihc ipf
df = read.table(file.path(path.data,'quant_ihc_ipf.txt'), header=TRUE, sep="\t", row.names=NULL, as.is=TRUE)
df = melt(df)
df$Group= rep('Control',nrow(df))
df$Group[grepl('Fibrosis',df$variable)]= 'Fibrosis'

pdf(file.path(path.out,'Fig3d.pdf'),width=3.5,height=2.8)
ggplot(df,aes(x=variable,y=value,color=variable))+
  geom_boxplot(alpha=0.2)+
  scale_color_manual(values=c("#E69F00","#E69F00","#E69F00","#56B4E9","#56B4E9","#56B4E9"))+
  geom_point(aes(color=variable,fill=variable),size=2,alpha=0.5,position=position_jitterdodge())+
  scale_color_manual(values=c("#E69F00","#E69F00","#E69F00","#56B4E9","#56B4E9","#56B4E9"))+
  theme(axis.title.x=element_blank())+theme(legend.position="none")+theme_classic() + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())+
  ylab('% of Mafb+ cells in a field of view\n(among cd68+ cells)')+
  theme(legend.position="none")
dev.off()

### Figure 5a - time course 
CreateContour = function(file) { 
  a <- as.matrix(read.csv(file, header=TRUE, row.names=NULL, as.is=TRUE))
  df = data.frame(x=(a[,2]),y=(a[,1]))
  
  x <- log10(df[,1])
  y <- log10(df[,2])
  print(file)
  print(length(x))
  print(mean(x<log10(5200) & y>log10(310)))
  print(mean(x<log10(5200) & y<log10(310)))
  print(mean(x>log10(5200) & y>log10(310)))
  print(mean(x>log10(5200) & y<log10(310)))
  
  xax <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    #showticklabels = FALSE,
    showgrid = FALSE,
    range = c(2,5.3)
  )
  yax <- list(
    title = "",
    zeroline = FALSE,
    showline = FALSE,
    #showticklabels = FALSE,
    showgrid = FALSE,
    range = c(1,3.7)
  )
  
  outliers = smoothScatter(log10(df),ret.selection=T,nrpoints=round(nrow(df)/10))
  plot_ly(x = x, y = y) %>% 
    add_histogram2dcontour(autocolorscale=T,showlegend=F,bgcolor=rgb(0,0,0),
                           contours = list(coloring='heatmap',end = 85, size = 10, start = 10)) %>%
    add_markers(x = x[outliers], y = y[outliers], marker=list(size=2), color=I("black"), 
                opacity=.5) %>%
    layout(xaxis = xax, yaxis = yax,
           showlegend=FALSE,
           shapes=list(list(type='line', x0 = log10(5200), x1=log10(5200), y0=0.1, y1=4, line=list(dash='dot', width=1)),
                       list(type='line', x0 = 0.7, x1=5.6, y0=log10(310), y1=log10(310), line=list(dash='dot', width=1))))
  
}

p1 = CreateContour(file.path(path.data,'export_0-1_Data Source - 1_tdTomato-A subset_APC-Cy7-A, APC-A subset.csv'))
export(p1, file = file.path(path.out,'Fig5a_0days.pdf'))
p2 = CreateContour(file.path(path.data,'export_7-1_Data Source - 1_tdTomato-A subset_APC-Cy7-A, APC-A subset.csv'))
export(p2, file = file.path(path.out,'Fig5a_7days.pdf'))
p3 = CreateContour(file.path(path.data,'export_14-1_Data Source - 1_tdTomato-A subset_APC-Cy7-A, APC-A subset.csv'))
export(p3, file = file.path(path.out,'Fig5a_14days.pdf'))

df = read.table(file.path(path.data,'cx3 timecourse flow data.txt'), header=TRUE, sep="\t", row.names=NULL, as.is=TRUE)
df$Days = factor(df$Days,levels=c('0 days','7 days','14 days'))
pdf(file.path(path.out,'Fig5a.pdf'),width=6,height=1.2)
ggplot(df,aes(x=Days,y=Value,color=Days,fill=Days))+
  #geom_boxplot(alpha=0.2,stat='median')+
  geom_point(size=2.5,alpha=0.5,position=position_jitterdodge())+
  stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean,
               geom = "crossbar", width = 0.5)+
  facet_wrap(~Measurement,nrow=1,scales = "free", strip.position = "left")  +
  ylab(NULL)+
  expand_limits(y = 0)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  theme(strip.background = element_blank(),
        strip.placement = "outside")
dev.off()


### Figure 6a - pdgfa in single-cell

agg=aggregate(singler.macs$seurat@data['Pdgfa',]>0,by=list(K),mean);
colnames(agg) = c('Clusters','value')
pdf(file.path(path.out,'Fig6a-sc-pdgfa.pdf'),width=3.8,2.3)
ggplot(agg,aes(Clusters,value*100)) + geom_bar(stat = "identity", aes(fill =Clusters), position = "dodge")+
  scale_fill_manual(values=cond_color[1:3])+xlab("")+ylab('% of expressing cells')
dev.off()

### Figure 6b - qunatitation of pdgfa

df = read.table(file.path(path.data,'fig6b-quantitation.txt'), header=TRUE, sep="\t", row.names=NULL, as.is=TRUE)
colnames(df) = c('Tdtomato-','Tdotomato+')
df = melt(df)
pdf(file.path(path.out,'Fig6b-tdtomata_pdgfa.pdf'),width=1.5,height=3.8)
ggplot(df,aes(x=variable,y=value,fill=variable,color=variable))+
  geom_boxplot(alpha=0.2,outlier.shape = NA)+
  geom_point(size=1,alpha=0.5,position=position_jitter(width = 0.2))+
  scale_color_manual(values=c("#00BA38","yellow3"))+
  scale_fill_manual(values=c("#00BA38","yellow3"))+
  #stat_compare_means(label.x.npc=0.4,label='p.signif')+
  ylab('% of cells in a field of view (among PDGF-AA+ cells)')+xlab('')+
  theme(legend.position="none",text = element_text(size=10),axis.text.x = element_text(size=8))
dev.off()

### Figure 6c - migration area

pdf(file.path(path.out,'Fig6c-migration.pdf'),width=4,height=4)
df = read.table('~/Documents/SingleR/migration_area.txt', header=TRUE, sep="\t", row.names=NULL, as.is=TRUE)
df$Group = factor(df$Group,levels=c('MHCIIlow','MHCIIlow+blockingAB','MHCIIhigh','MHCIIhigh+blockingAB'))
df$MigrationArea = (df$MigrationArea/mean(df$MigrationArea[df$Group=='MHCIIlow']))
ggplot(df,aes(x=Group,y=MigrationArea,color=Group,fill=Group))+
  geom_boxplot(alpha=0.2,outlier.shape = NA)+
  geom_point(size=2,alpha=0.5,position=position_jitterdodge())+
  scale_y_continuous(limits = c(0.6,1.5)) +
  theme(axis.title.x=element_blank())+theme(legend.position="none")+
  stat_compare_means(label.x.npc=0.4,label='p.signif',comparisons=list(c('MHCIIlow','MHCIIhigh'),c('MHCIIhigh','MHCIIhigh+blockingAB')))+
  ylab('Migration Area (A.U.)')
dev.off()

### Figure 6d - edu

edu = c(11.2, 10.4,12,7.37,7.02,7.62)
types = c('IgG1','IgG2','IgG3','PDGF-Ab1','PDGF-Ab2','PDGF-Ab3')
df = data.frame(edu,types)
ggplot(df,aes(types,edu,fill=types))+geom_col()+
  scale_fill_manual(values=c(singler.colors[10],singler.colors[10],singler.colors[10],singler.colors[12],singler.colors[12],singler.colors[12]))+
  ylab('% EDU-positive fibroblasts')
ggsave(file.path(path.out,'Fig6d-edu.pdf'),width=4.5,height=4)

edu = c(11.2, 10.4,12,7.37,7.02,7.62)
types = c('IgG','IgG','IgG','PDGF-AA Ab','PDGF-AA Ab','PDGF-AA Ab')
df = data.frame(edu,types)
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

df = data_summary(df,'edu','types')

ggplot(df, aes(x=types, y=edu, fill=types)) + 
  geom_bar(position=position_dodge(), stat="identity") +
  geom_errorbar(aes(ymin=edu-sd, ymax=edu+sd),
                width=.2, 
                position=position_dodge(.9))+
  ylab('% EDU-positive fibroblasts')+xlab('')+ 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(path.out,'Fig6d-edu.pdf'),width=3.5,height=4)


### Figure 6e - sc fibroblasts
singler$meta.data$orig.ident=full.cond
singler.fibs = SingleR.Subset(singler,singler$seurat@ident %in% c(9,11) & singler$meta.data$xy[,1]<10 & singler$meta.data$xy[,1]> -20 & singler$meta.data$xy[,2]< -38)

g2m.genes = geneIds(mouse.egc[[2]])
s.genes=geneIds(mouse.egc[[5]])
s = singler.fibs$seurat
rownames(s@data) = tolower(rownames(s@data))
s = CellCycleScoring(s,s.genes = s.genes, g2m.genes = g2m.genes)
SingleR.PlotFeature(singler.fibs$singler[[1]]$SingleR.single,singler.fibs$seurat,plot.feature = s@meta.data$S.Score+s@meta.data$G2M.Score,dot.size = 3)

SingleR.PlotFeature(singler.fibs$singler[[1]]$SingleR.single,singler.fibs$seurat,plot.feature = 'Mki67',dot.size = 3)

A = singler.fibs$meta.data$xy[,1] > -5 & singler.fibs$meta.data$xy[,2] > -42
B = singler.fibs$meta.data$xy[,1] < -3 & singler.fibs$meta.data$xy[,2] > -47
K = as.numeric(B)
K[A] = 2
K = factor(K)
df = data.frame(tSNE1=singler.fibs$meta.data$xy[,1],tSNE2=singler.fibs$meta.data$xy[,2],Condition=singler.fibs$meta.data$orig.ident,Group=K,Cell.Cycle = s@meta.data$S.Score+s@meta.data$G2M.Score,Mki67=singler.fibs$seurat@data['Mki67',])
pdf(file.path(path.out,'Fig6e.sc.fibs.pdf'),width=5.2,height=4.2)
ggplot(df,aes(x=tSNE1,y=tSNE2))+geom_point(aes(color=Cell.Cycle,shape=Condition),size=2,alpha=0.5)+
  scale_colour_gradient(low='gray',high='blue')
dev.off()

pdf(file.path(path.out,'Fig6e.sc.fibs-boxplot.pdf'),width=3.8,height=4)
df = data.frame(Group=A,Condition=singler.fibs$meta.data$orig.ident,Cell.Cycle = s@meta.data$S.Score+s@meta.data$G2M.Score)
df$Group[A==TRUE] = 'A'
df$Group[A==FALSE] = 'B & C'
df=melt(df,id.vars = c('Group','Condition'))
ggplot(df,aes(y=value,x=Group))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(aes(color=Condition),size=1.2,alpha=0.6,position=position_jitter())+ylab('Cell-cycle score')+xlab('Fibroblasts cluster')
dev.off()

### Figure 7a  - human PDGFA

df = read.table(file.path(path.data,'Quantitation_PDGFAA _in_CD68+_and_SH.txt'), header=TRUE, sep="\t", row.names=NULL, as.is=TRUE)
df = df[df$Patients %in% c('18542','18545','18546'),]
df$Group[df$Group=='High']='Hi'
pdf(file.path(path.out,'Fig7c.Quantitation_PDGFAA _in_CD68+_and_SH.pdf'),width=3,height=2)
ggplot(df,aes(x=Group,y=PDGF.AA.CD68,fill=Group,color=Group))+
  geom_boxplot(alpha=0.2,outlier.shape = NA)+facet_grid(~Patients)+
  geom_point(size=1,alpha=0.5,position=position_jitter(width = 0.2))+
  stat_compare_means(label.x.npc=0.4,label='p.signif')+
  ylab('PDGF-AA in CD68+ (MFI)')+xlab('SH signal for collagen')+
  theme(legend.position="none",text = element_text(size=10),axis.text.x = element_text(size=8))
dev.off()

### Figure 7b - hydroproxyline

df = read.table(file.path(path.data,'hydroxoproline.txt'), header=TRUE, sep="\t", row.names=NULL, as.is=TRUE)
df = df[1:24,]
df$Group = factor(df$Group)
df$Type = as.factor(df$Type)
pdf(file.path(path.out,'Fig7b.hydroxyproline.pdf'),width=4,height=3.4)
ggplot(df,aes(x=Group,y=Value,color=Group,fill=Group))+
  geom_boxplot(alpha=0.2,outlier.shape = NA)+
  geom_point(aes(shape=Type),size=2,alpha=0.5,position=position_jitterdodge())+
  #geom_point(aes(shape=Type),size=2,alpha=0.5,position=position_jitter(width = 0.2))+
  scale_shape_manual(values=c(17,16)) + 
  theme(axis.title.x=element_blank())+theme(legend.position="none")+theme_classic2()+
  ylim(0,2.3)+
  stat_compare_means(label.x.npc=0.4,label='p.signif',comparisons=list(c('1','2'),c('2','3')))
dev.off()


## Supplementary Figure 7

quants = read.table('~/Documents/SingleR/manuscript/quantification_pdgfra_pdgfrb.txt', header=TRUE, sep="\t", row.names=NULL, as.is=TRUE)
df= quants
df$Set = factor(df$Set,levels=c('WT','Cre+DTA+'))
pdf(file.path(path.out,'Fig7a-Quantitation-pdgfra-pdgfrb.pdf'),width=5,height=2.1)
ggplot(df,aes(x=Set,y=Quant,color=Set,fill=Set))+
  geom_boxplot(alpha=0.2,outlier.shape = NA)+
  geom_point(aes(shape=Mice),size=2,alpha=0.5,position=position_jitterdodge(dodge.width = 0.05))+
  facet_wrap(~Marker,scales = "free") +
  stat_compare_means(label.x.npc=0.4,label='p.signif',comparisons=list(c('WT','Cre+DTA+')))+ ylab('MFI')+
  theme(axis.title.x=element_blank(),legend.position="none",text = element_text(size=10),axis.text.x = element_text(size=8))
dev.off()

