library(RColorBrewer)
library(reshape2)
library(SingleR)
library(pheatmap)
library(ggplot2)
library(ggpubr)
require("biomaRt")

path.data = '~/Documents/SingleR/FiguresData'
path.out = '~/Documents/SingleR/FiguresOut'

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

# Top-left
pdf(file.path(path.out,'Fig1b_topleft.pdf'),width=5,height=3.75)
out = SingleR.PlotTsne(singler$singler[[1]]$SingleR.single,xy,do.label=FALSE,
                       do.letters =F,labels=singler$meta.data$orig.ident, 
                       dot.size = 3)
out$p
dev.off()

# Bottom-left
cl.set2 =brewer.pal(8,'Set2')
cl.set2[3] = cl.set2[4];cl.set2[4] = cl.set2[1];cl.set2[1] = cl.set2[2];cl.set2[2] = cl.set2[6]

pdf(file.path(path.out,'Fig1b_bottomleft.pdf'),width=5.33,height=3.75)
out = SingleR.PlotTsne(singler$singler[[1]]$SingleR.single.main,do.letters = F,colors=cl.set2,
                       xy,do.label=FALSE, dot.size = 3)
out$p
dev.off()

# Top-right
pdf(file.path(path.out,'Fig1b_topright.pdf'),width=4,height=3)
SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single.main,top.n=8)
dev.off()

# Bottom-right
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

pdf(file.path(path.out,'Fig1b_bottompright.pdf'),width=7,height=4.5)
pheatmap(a[order(genes.use$Group),order(singler$singler[[1]]$SingleR.single.main$labels[BMDC])],
         cluster_cols = F,cluster_rows = F,clustering_method='ward.D2',border_color = NA,
         annotation_col=annotation_col,
         annotation_row = genes.use,show_colnames=F,show_rownames = F,
         annotation_colors=list(Annotation=c('DC'=cl.set2[1],'Macrophages'=cl.set2[3])),
         fontsize = 14)
dev.off()

### Figure 1c - Case/Control of lung sc data

load(file.path(path,'GSE111664.RData'))
singler.macs = SingleR.Subset(singler,singler$seurat@ident %in% c(0,1,2,7))

full.cond = as.character(singler$meta.data$orig.ident)
full.cond[grepl('ontrol',full.cond)]='Control'
full.cond[grepl('leo',full.cond)]='Bleomycin'

p = SingleR.PlotTsne(singler$singler[[1]]$SingleR.single,singler$meta.data$xy,labels=full.cond,
                     do.letters=F,colors=brewer.pal(3,'Set2'),do.labels = F,dot.size=0.5,alpha=0.35)
pdf(file.path(path.out,'Fig1c.pdf'),width=4,height=5)
p$p+theme_void()+theme(legend.position="none")
dev.off()

full.cond2 = as.character(singler$meta.data$orig.ident)
full.cond2[grepl('ontrol',full.cond)]='Healthy lungs'
full.cond2[grepl('leo',full.cond)]='Bleomycin-injured lungs'

p = SingleR.PlotTsne(singler$singler[[1]]$SingleR.single,singler$meta.data$xy,labels=full.cond2,do.letters=F,
                     colors=brewer.pal(3,'Set2'),do.labels = F,dot.size=3,alpha=0.35)
pdf(file.path(path.out,'Fig1c_legend.pdf'),width=3.5,height=3.5)
p$p+theme_void()
dev.off()

### Figure 1d - SingleR annotations

p = SingleR.PlotTsne(singler$singler[[1]]$SingleR.single.main,singler$meta.data$xy,do.letters=F,
                     colors=sample(c(brewer.pal(12,'Set3'),brewer.pal(8,'Set2'))),do.labels = F,
                     dot.size=0.5,alpha=0.35)
pdf(file.path(path.out,'Fig1d.pdf'),width=5,height=5)
p$p+theme_void()+theme(legend.position="none")
dev.off()

### Figure 1e - SingleR annotations with lung myeloid data

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
pdf(file.path(path.out,'Fig1d.pdf'),width=3.5,height=4.3)
p$p+theme_void()+theme(legend.position="none")
dev.off()

p = SingleR.PlotTsne(singler.macs$singler[[3]]$SingleR.single,singler.macs$meta.data$xy,labels=labels,do.letters=F,
                     colors=c(brewer.pal(8,'Accent'),brewer.pal(3,'Set1')),do.labels = F,dot.size=3,alpha=0.5)
pdf(file.path(path.out,'Fig1d_legend.pdf'),width=3.5,height=4.3)
p$p+theme_void()
dev.off()

### Figure 2a - deconvolution

library(DeconRNASeq)
library(RColorBrewer)

# create decon data
load(file.path(path.data,'GSE94135_GSE49932.RData'))
amim = cbind(rowMeans(gse94135_gse49932$data[,1:3]),rowMeans(gse94135_gse49932$data[,10:12]))
#load('~/Documents/singler/macs_expr.rds')
decon = DeconRNASeq(as.data.frame(as.matrix(expr)),as.data.frame(amim),checksig=FALSE,known.prop = FALSE)

#load('~/Documents/singler/decon.rds')

df = data.frame(x=singler.macs$meta.data$xy[,1],y=singler.macs$meta.data$xy[,2])

df$Scale = decon$out.all[,1]
cl = brewer.pal(3,'Set1')
pdf(file.path(path.out,'Fig2a.pdf'),width=3.5,height=3)
ggplot(df) + geom_point(aes(x=x, y=y,color=Scale),size=1.2)+
  scale_color_gradient(low=cl[1],high=cl[2])+theme_void()
dev.off()


### Figure 2b - clustering

col = brewer.pal(8,'Set2')
p = SingleR.Cluster(singler.macs$singler[[1]]$SingleR.single,3)

K = p$cl
K = plyr::mapvalues(K,from=1:3,to=c('C1','C2','C3'))

singler.macs$seurat@ident = K
ggplot(df) + geom_point(aes(x=x, y=y,color=K),size=0.7)+theme_void()+theme(legend.position="none")

pdf(file.path(path.out,'Fig2b.pdf'),width=3,height=3)
ggplot(df) + geom_point(aes(x=x, y=y,color=K),size=0.7)+theme_void()+theme(legend.position="none")
dev.off()

### Figure 2b insert

cond = as.character(singler.macs$meta.data$orig.ident)
cond[grepl('ontrol',cond)]='Control'
cond[grepl('leo',cond)]='Bleomycin'
singler.macs$meta.data$orig.ident=cond

annotation_col = data.frame(Cluster=K,Condition=factor(cond))
rownames(annotation_col) = singler.macs$seurat@cell.names
cluster_colors = gg_color_hue(3)
cond_color = brewer.pal(3, 'Set2')
annotation_colors = list(Condition=c(Bleomycin=cond_color[2],Control=cond_color[1]),
                         Cluster=c('C1'=cluster_colors[1],'C2'=cluster_colors[2],'C3'=cluster_colors[3]))

cond.cluster = table(annotation_col$Condition,K)/cbind(table(full.cond),table(full.cond),table(full.cond))
cond.cluster = melt(cond.cluster)

pdf(file.path(path.out,'Fig2b_insert.pdf'),width=3,height=2.5)
ggplot(cond.cluster,aes(K,value*100)) + geom_bar(stat = "identity", aes(fill =Var.1), position = "dodge")+
  scale_fill_manual(values=cond_color[1:2])+xlab("")+ylab('% of cells')
dev.off()

# Figure 2c - de

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

pdf(file.path(path.out,'Fig2c.pdf'),width=6,height=4)
pheatmap(d[,order(annotation_col$Cluster,num.genes[2,]-num.genes[1,])],
         annotation_col=annotation_col,annotation_colors=annotation_colors,cluster_cols=F,
         cluster_rows=F,show_colnames=F,border_color=NA,fontsize=7)
dev.off()

### Figure 2d - expressed genes

rownames(num.genes) = c('C1 genes','C3 genes')
res = c()
for (i in seq(0,0.85,by=0.01)) {
  res = rbind(res,(table(num.genes[1,]>i & num.genes[2,]>i,annotation_col$Cluster)/rbind(table(annotation_col$Cluster),table(annotation_col$Cluster)))[2,])
}
colnames(res) = c('C1','C2','C3')
f = melt(res)

pdf(file.path(path.out,'Fig2d.pdf'),width=2.5,height=2)
ggplot(f)+geom_smooth(aes(x=X1,y=f$value*100,color=X2))+ylab('% of cells')+
  xlab('% of expressed genes from C1 & C3')+
  theme(text = element_text(size=7),axis.text.x = element_text(size=7),axis.text.y = element_text(size=7))
dev.off()

df = data.frame(x=singler.macs$meta.data$xy[,1],y=singler.macs$meta.data$xy[,2])
for (i in 1:nrow(markers)) {
  df$Scale = as.matrix(singler.macs$seurat@data[rownames(markers)[i],])
  ggplot(df) + geom_point(aes(x=x, y=y,color=Scale),size=0.3,alpha=0.5)+
    scale_color_gradient(low='gray',high='darkblue')+ggtitle(rownames(markers)[i])+theme_void()+
    theme(title =element_text(size=1, face='bold'))
  ggsave(file.path(path.out,paste0('Fig2d_',rownames(markers)[i],'.pdf')),width=2,height=2)
}

### Figure 2e - bulk validation

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
pdf(file.path(path.out,'Fig2e.pdf'),width=7,height=3.5)
pheatmap(m,scale='row',cluster_cols = F,annotation_col=annotation_col,show_colnames=F,cluster_rows = T,
         border_color=NA,clustering_method='ward.D2')
dev.off()

### Figure 2f - time course 
pdf(file.path(path.out,'Fig2f.pdf'),width=7.5,height=2.5)
df = read.table(file.path(path.data,'cx3 timecourse flow data.txt'), header=TRUE, sep="\t", row.names=NULL, as.is=TRUE)
df$Days = factor(df$Days,levels=c('0 days','7 days','14 days'))
ggplot(df,aes(x=Days,y=Value,color=Days,fill=Days))+
  geom_boxplot(alpha=0.2)+
  geom_point(size=4,alpha=0.5,position=position_jitterdodge())+
  facet_wrap(~Measurement,nrow=1,scales = "free")+
  expand_limits(y = 0)+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
dev.off()

### Figure 3a = pdgfa
pdf(file.path(path.out,'Fig3a.pdf'),width=4,height=3)
df = data.frame(Pdgfa = t(M['Pdgfa',c(13,14,2,4,6,1,3,5,8,10,12,7,9,11)]), Group = annotation_col$Group)
df$Group = factor(df$Group,levels=c('Baseline','2 wk MHC lo','2 wk MHC high','4 wk MHC lo','4 wk MHC high'))
ggplot(df,aes(x=Group,y=Pdgfa,color=Group,fill=Group))+
  geom_boxplot(alpha=0.2)+
  geom_point(size=3,alpha=0.5,position=position_jitterdodge())+
  theme(axis.title.x=element_blank())+theme(legend.position="none")+
  ylab('Pdgfa normalized expression')
dev.off()

### Figure 3b - migration area

pdf(file.path(path.out,'Fig3b.pdf'),width=4,height=3)
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

### Figure 3c - hydroproxyline

pdf(file.path(path.out,'Fig3c.pdf'),width=3.5,height=3)
df = read.table(file.path(path.data,'hydroxoproline.txt'), header=TRUE, sep="\t", row.names=NULL, as.is=TRUE)
df$Group = factor(df$Group,levels=c('Control','Bleomycin','CX3CR1 Bleomycin'))
ggplot(df,aes(x=Group,y=Value,color=Group,fill=Group))+
  geom_boxplot(alpha=0.2)+
  geom_point(size=2,alpha=0.5,position=position_jitterdodge())+
  theme(axis.title.x=element_blank())+theme(legend.position="none")+
  ylim(0,2.3)+
  stat_compare_means(label.x.npc=0.4,label='p.signif',comparisons=list(c('Control','Bleomycin'),c('Bleomycin','CX3CR1 Bleomycin')))
dev.off()

### Figure 3e - human samples
expr <- as.matrix(read.table(file.path(path.data,'GSE32537_expr.txt'), header=TRUE, sep="\t", row.names=1, as.is=TRUE,comment.char='!'))
colnames(expr) = unlist(lapply(colnames(expr),FUN=function(x) strsplit(x,'_')[[1]][1]))
clinical = as.matrix(read.table(file.path(path.data,'GSE32537_clinical.txt'), header=TRUE, sep="\t", row.names=1, as.is=TRUE,comment.char='!'))
IPF = clinical[,'final.diagnosis']=='control'
expr = expr[,rownames(clinical)]

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                 values = rownames(am) , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
human.am <- unique(genesV2[, 2])
genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", 
                 values = rownames(im) , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
human.im <- unique(genesV2[, 2])
human.mhcii = rownames(expr)[which(grepl('HLA-D',rownames(expr)))]

sig = list()
sig[[1]] = GeneSet(human.am, setName='C1')
sig[[2]] = GeneSet(human.im, setName='C3')
sig[[3]] = GeneSet(human.mhcii, setName='MHCII')

egc <- GeneSetCollection(sig)

scores = gsva(expr,egc,method='ssgsea',rnaseq=F)

Group = clinical[,'final.diagnosis']!='control'
Group[Group==F] = 'Control'
Group[Group==T] = 'Lung Fibrosis'

df = data.frame('C1 genes'=scale(scores[1,]),'MHCII genes'=scale(scores[3,]),
                CX3CR1=scale(expr['CX3CR1',]),Group=Group)  
colnames(df) <- gsub("\\.", " ", colnames(df))
colnames(df)
d = melt(df,variable.name='Group')
colnames(d) = c('Group','Set','Score')
pdf(file.path(path.data,'Fig3e.pdf'),width=8,height=4)
ggplot(d,aes(x=Group,y=Score,color=Group,fill=Group))+geom_boxplot(alpha=0.2)+
  geom_point(size=0.5,alpha=0.5,position=position_jitterdodge())+
  facet_wrap(~Set,nrow=1)+ 
  stat_compare_means(label.x.npc=0.4,size=8,label='p.signif')+
  ylab('Scaled values')+xlab("")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),axis.ticks.x=element_blank())
dev.off()

