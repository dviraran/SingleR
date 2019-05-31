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


#' Stable correlation analysis. Courtesy of Thomas Wu.
#'
#' @param x x
#' @param y y
#' @param method method
#' #'
#' @return r
cor.stable <- function (x, y, method="pearson", ...) {
  omit1 <- which(apply(x, 2, sd) == 0)
  omit2 <- which(apply(y, 2, sd) == 0)
  if (length(omit1) > 0 && length(omit2) > 0) {
    r <- matrix(0, ncol(x), ncol(y))
    r[-omit1,-omit2] = cor(x[,-omit1], y[,-omit2], method=method, ...)
  } else if (length(omit1) > 0) {
    r <- matrix(0, ncol(x), ncol(y))
    r[-omit1,] = cor(x[,-omit1], y, method=method, ...)
  } else if (length(omit2) > 0) {
    r <- matrix(0, ncol(x), ncol(y))
    r[,-omit2] = cor(x, y[,-omit2], method=method, ...)
  } else {
    r = cor(x, y, method=method, ...)
  }
}