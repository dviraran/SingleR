# SingleR - Single-cell Recognition

Recent advances in single cell RNA-seq (scRNA-seq) have enabled an unprecedented level of granularity in characterizing gene expression changes in disease models. Multiple single cell analysis methodologies have been developed to detect gene expression changes and to cluster cells by similarity of gene expression. However, the classification of clusters by cell type relies heavily on known marker genes, and the annotation of clusters is performed manually. This strategy suffers from subjectivity and limits adequate differentiation of closely related cell subsets. Here, we present SingleR, a novel computational method for unbiased cell type recognition of scRNA-seq. SingleR leverages reference transcriptomic datasets of pure cell types to infer the cell of origin of each of the single cells independently. SingleR’s annotations combined with Seurat, a processing and analysis package designed for scRNA-seq, provide a powerful tool for the investigation of scRNA-seq data. We developed an R package to generate annotated scRNA-seq objects, that can then use the SingleR web tool for visualization and further analysis of the data – <http://comphealth.ucsf.edu/SingleR>.

For more informations please refer to the manuscript: Aran et al. Single-cell RNA-seq reveals profibrotic Cx3cr1+ macrophages in lung fibrosis. XXX, 2018

# Install

devtools::install_github('dviraran/SingleR')

# Usage

library(SingleR)

See [https://github.com/dviraran/SingleR/vignettes/SingleR_specifications.Rmd](SingleR specifications) for examples of usage.

# Contributors

SingleR was developed by the Dvir Aran. Please contact Dvir Aran: dvir.aran at ucsf edu for any questions or suggestions.

