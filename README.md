# SingleR - Single-cell Recognition

Recent advances in single cell RNA-seq (scRNA-seq) have enabled an unprecedented level of granularity in characterizing gene expression changes in disease models. Multiple single cell analysis methodologies have been developed to detect gene expression changes and to cluster cells by similarity of gene expression. However, the classification of clusters by cell type relies heavily on known marker genes, and the annotation of clusters is performed manually. This strategy suffers from subjectivity and limits adequate differentiation of closely related cell subsets. Here, we present SingleR, a novel computational method for unbiased cell type recognition of scRNA-seq. SingleR leverages reference transcriptomic datasets of pure cell types to infer the cell of origin of each of the single cells independently. SingleR’s annotations combined with Seurat, a processing and analysis package designed for scRNA-seq, provide a powerful tool for the investigation of scRNA-seq data. We developed an R package to generate annotated scRNA-seq objects, that can then use the SingleR web tool for visualization and further analysis of the data – <http://comphealth.ucsf.edu/SingleR>.

For more informations please refer to the manuscript: Aran, Looney, Liu et al. Single-cell RNA-seq reveals profibrotic macrophages in lung fibrosis. bioRxiv, 2018

# Install

```R
devtools::install_github('dviraran/SingleR')
# this might take long, though mostly because of the installation of Seurat.
```

# Usage

```R
library(SingleR)

# Simplest use is running the wrapper function that creates both a SingleR and Seurat object:

# counts.file maybe a tab delimited text file, 10X directory or a matrix. annot is a tab delimited text file or
# a data.frame with the original identities. normalize.gene.length should be true if the data comes from a
# full-length platform. min.genes, npca and regress.out are passed to Seurat to create a Seurat object object.
singler = CreateSinglerSeuratObject(counts.file,annot,project.name,species,citation,normalize.gene.length,min.genes,regress.out,npca,technology)

# The object can then be saved and uploaded to the SingleR web-app for further analysis and visualization
# or using functions available in the SingleR package (see vignette).
save(singler,file=paste0(project.name,'.RData')
```

For more details and examples see [SingleR specifications](http://comphealth.ucsf.edu/sample-apps/SingleR/SingleR_specifications.html).

# Contributors

SingleR was developed by Dvir Aran. Please contact Dvir Aran: dvir.aran at ucsf edu for any questions or suggestions.

