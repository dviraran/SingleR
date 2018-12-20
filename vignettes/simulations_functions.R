
trimExpression = function(expr,n.genes) {
  m = apply(expr,2,function(x) {
    n = min(n.genes,sum(x>0))
    A = sample(which(x>0),n)
    y = matrix(0,length(x))
    y[A] = x[A]
    y
  })
  rownames(m) = rownames(expr)
  colnames(m) = colnames(expr)
  m
}


summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  # Rename the "mean" column
  datac <- rename(datac, c("mean" = measurevar))
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  return(datac)
}


countCorrect = function(labels,orig.ident) {
  sum(grepl('B-cell',labels) & orig.ident=='B-cells')+
    sum((labels=='CD4+ Tcm' | labels=='CD4+ Tem') & orig.ident=='CD4+ Tm')+
    sum(labels=='CD4+ T-cells' & orig.ident=='CD4+ Tn')+
    sum((labels=='CD8+ Tcm' | labels=='CD8+ Tem') & orig.ident=='CD8+ Tcyto')+
    sum(labels=='CD8+ T-cells' & orig.ident=='CD8+ Tn')+
    sum(labels %in% c('CLP','MEP','GMP','MPP') & orig.ident=='HSC')+
    sum(labels=='Monocytes' & orig.ident=='Monocytes')+
    sum(labels=='NK cells' & orig.ident=='NK cells')+
    sum(labels=='Tregs' & orig.ident=='Tregs')
}
