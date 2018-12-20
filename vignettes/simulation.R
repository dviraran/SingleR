library(SingleR)
library(Seurat)
library(reshape2)
library(ggplot2)
library(directlabels)


dirs = dir(paste0('~/Documents/SingleR/10X_sorted'),full.names=T)
dirs = dirs[-2]
tenx1 = Combine.Multiple.10X.Datasets(dirs[-c(3,5)],random.sample=50,min.genes=1000)
tenx2 = Combine.Multiple.10X.Datasets(dirs[c(3,5)],random.sample=100,min.genes=900)

tenx = list()
tenx$sc.data = cbind(tenx1$sc.data,tenx2$sc.data)
tenx$orig.ident = c(tenx1$orig.ident,tenx2$orig.ident)

data = trimExpression(tenx$sc.data,1000)

#s = SingleR(method = "single", tenx$sc.data, blueprint_encode$data, blueprint_encode$types)

s = SingleR(method = "single", data, blueprint_encode$data, blueprint_encode$types)

n = 10
ngene = colSums(tenx$sc.data>0)
cells.use = c(rownames(s$labels)[ngene>=1000 & grepl('B-cell',s$labels) & tenx$orig.ident=='B-cells'][1:n],
              rownames(s$labels)[ngene>=1000 & (s$labels=='CD4+ Tcm' | s$labels=='CD4+ Tem') & tenx$orig.ident=='CD4+ Tm'][1:n],
              rownames(s$labels)[ngene>=1000 & (s$labels=='CD8+ Tcm' | s$labels=='CD8+ Tem') & tenx$orig.ident=='CD8+ Tcyto'][1:n],
              rownames(s$labels)[ngene>=1000 & s$labels %in% c('CLP','MEP','GMP') & tenx$orig.ident=='HSC'][1:n],
              rownames(s$labels)[ngene>=1000 & s$labels=='Monocytes' & tenx$orig.ident=='Monocytes'][1:n],
              rownames(s$labels)[ngene>=1000 & s$labels=='NK cells' & tenx$orig.ident=='NK cells'][1:n],
              rownames(s$labels)[ngene>=1000 & s$labels=='Tregs' & tenx$orig.ident=='Tregs'][1:n],
              rownames(s$labels)[s$labels=='CD8+ T-cells' & tenx$orig.ident=='CD8+ Tn'][1:n],              
              rownames(s$labels)[s$labels=='CD4+ T-cells' & tenx$orig.ident=='CD4+ Tn'][1:n])

#m = trimExpression(m[,cells.use],1000)

orig.ident = tenx$orig.ident[cells.use]

seq.genes = seq(1000,100,by=-50)

tbl= list()
tbl.nft = list()
pvals = list()
correct = matrix(NA,length(seq.genes))
rownames(correct) = seq.genes
correct.nft = correct


for (k in 1:10) {
  m = data[,cells.use]
  j=1
  for (i in seq.genes) {
    print(i)
    m = trimExpression(m,i)
    s = SingleR(method = "single", m, blueprint_encode$data, blueprint_encode$types)
    tbl[[j]] = table(s$labels,orig.ident)
    tbl.nft[[j]] = table(s$labels1,orig.ident)
    pvals[[j]] = s$pval
    print(tbl[[j]])
    correct[j] = countCorrect(s$labels,orig.ident)
    correct.nft[j] = countCorrect(s$labels1,orig.ident)
    j=j+1
  }
  save(tbl,tbl.nft,correct,correct.nft,pvals,file=paste0('~/Documents/SingleR/simulations/simulation.',k,'.RData'))
}

df = data.frame(nGenes = seq(1000,100,by=-50), Correct = correct, Correct.NFT = correct.nft)
df = melt(df,id.vars = 'nGenes')
df$nGenes = as.numeric(df$nGenes)
df$value = df$value/90
ggplot(df,aes(x=nGenes,y=value,color=variable)) + geom_line()+theme_classic()



files = dir('~/Documents/SingleR/simulations/',full.names = T,pattern = 'RData')

res = list()
res$correct = matrix(NA,19,length(files))
res$correct.nft = matrix(NA,19,length(files))
for (i in 1:length(files)) {
  load(files[i])
  res$correct[,i] = correct
  res$correct.nft[,i] = correct.nft 
}

res$correct = 100*res$correct/90
res$correct.nft = 100*res$correct.nft/90
rownames(res$correct) = seq(1000,100,by=-50)
rownames(res$correct.nft) = seq(1000,100,by=-50)
correct = melt(rbind(res$correct,res$correct.nft))
correct$FineTune=rep(c(rep('TRUE',nrow(res$correct)),rep('FALSE',nrow(res$correct))),ncol(res$correct))
colnames(correct)[1]='nGenes'
correct <- summarySE(correct, measurevar="value", groupvars=c("nGenes",'FineTune'))


ggplot(correct, aes(x=nGenes, y=value,color=FineTune)) + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1) +
  geom_line() +
  geom_point()+theme_classic()+ylab('% correct')

ggsave(file='~/Documents/SingleR/simulations/simulations.all.cells.pdf',width=6,height=2.7)


B = c()
for (k in 1:length(files)) {
  
  load(files[k])
  ct=unique(names(unlist(lapply(tbl,FUN=function(x) x[,'B-cells']))))
  ident = colnames(tbl[[1]])
  b = matrix(0,length(ident),length(tbl))
  for (j in 1:length(ident)) {
    a=matrix(0,length(ct),length(tbl))
    rownames(a) = ct
    
    for (i in 1:length(tbl)) {
      x = tbl[[i]][,ident[j]]
      a[names(x),i] = x 
      
    }
    b[j,] = switch(ident[j], 
                   'B-cells' = colSums(a[grepl('B-cell',rownames(a)),]),
                   'CD4+ Tm' = colSums(a[rownames(a) %in% c('CD4+ Tcm','CD4+ Tem'),]),
                   'CD4+ Tn' = a[rownames(a)=='CD4+ T-cells',],
                   'CD8+ Tcyto' = colSums(a[rownames(a) %in% c('CD8+ Tcm','CD8+ Tem'),]),
                   'CD8+ Tn' = a[rownames(a)=='CD8+ T-cells',],
                   'HSC' = colSums(a[rownames(a) %in% c('CLP','MEP','GMP','MPP'),]),
                   'Monocytes' = a[rownames(a)=='Monocytes',],
                   'NK cells' = a[rownames(a)=='NK cells',],
                   'Tregs' = a[rownames(a)=='Tregs',]
    )
    
    
  }  
  rownames(b) = ident
  colnames(b) = seq(1000,100,by=-50) 
  
  if (length(B)==0) {
    B = b
  } else {
    B = rbind(B,b)
  }
}



df = melt(B)
colnames(df) = c('CellType','nGenes','value')
df <- summarySE(df, measurevar="value", groupvars=c("nGenes",'CellType'))

ggplot(df, aes(x=nGenes, y=value,color=CellType)) + 
  geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.1) +
  geom_line() +
  geom_point()+theme_classic()+ylab('# correct')+scale_color_manual(values=singler.colors)

ggsave(file='~/Documents/SingleR/simulations/simulations.per.cells.pdf',width=6,height=3)
