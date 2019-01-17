# SingleR - Single-cell Recognition

Recent advances in single cell RNA-seq (scRNA-seq) have enabled an unprecedented level of granularity in characterizing gene expression changes in disease models. Multiple single cell analysis methodologies have been developed to detect gene expression changes and to cluster cells by similarity of gene expression. However, the classification of clusters by cell type relies heavily on known marker genes, and the annotation of clusters is performed manually. This strategy suffers from subjectivity and limits adequate differentiation of closely related cell subsets. Here, we present SingleR, a novel computational method for unbiased cell type recognition of scRNA-seq. SingleR leverages reference transcriptomic datasets of pure cell types to infer the cell of origin of each of the single cells independently. SingleR’s annotations combined with Seurat, a processing and analysis package designed for scRNA-seq, provide a powerful tool for the investigation of scRNA-seq data. We developed an R package to generate annotated scRNA-seq objects, that can then use the SingleR web tool for visualization and further analysis of the data – <http://comphealth.ucsf.edu/SingleR>.

For more informations please refer to the manuscript: Aran, Looney, Liu et al. Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage. Nature Immunology (2019)

# Install

```R
devtools::install_github('dviraran/SingleR')
# this might take long, though mostly because of the installation of Seurat.
```

# Updates

**12.18.2018**: The SingleR browser has been upgraded allowing to use datasets with >100K single-cells. In addition, to make it more user friendly, the object used for SingleR has been upgraded as well to an S4 object. The web browser only accepts the new S4 objects. 

However, the functions in the package (outside of the browser) have not been yet been modified to work with new object. Thus, to use the browser, please convert the SingleR object (the object created by the CreateSingleR functions) to the S4 object using the following code:

```R
singler.new = convertSingleR2Browser(singler)
saveRDS(singler.new,file=paste0(singler.new@project.name,'.rds')
```

The new SingleR S4 object simlifies the access SingleR annotations and allows multipe identity columns (orig.ident) and clustering columns. See ?'SingleR-class' for more details.

**11.6.2018**: The Seurat team has roled out Seurat version 3.0 with many changes to the function names and variables. It is important to note that SingleR is independent of Seurat, and only requires xy representation of the single-cells, which can be obtained by any function. However, SingleR provides wrapper functions to streamline the analysis with Seuart which are affected by these changes. SingleR is now updated to run with Seurat 3.0, but it is still in beta mode. 

Best practice is to create the single-cell object first. Then, create the SingleR object with the same cells used in the single-cell object (after filtering low quality cells). Finally, copy to relevant fields to the SingleR object. The new t-SNE representations are in 
`sc@reductions$tsne@cell.embeddings`.

# Usage

```R
library(SingleR)

# Simplest use is running the wrapper function that creates both a SingleR and Seurat object:

# counts.file maybe a tab delimited text file, 10X directory or a matrix. annot is a tab delimited 
# text file or a data.frame with the original identities. normalize.gene.length should be true if 
# the data comes from a full-length platform. min.genes, min.cells, npca and regress.out are passed 
# to Seurat to create a Seurat object object:
singler = CreateSinglerSeuratObject(counts.file, annot, project.name,
  min.genes = 500, technology, species = "Human" (or "Mouse"), citation,
  normalize.gene.length = F, min.cells = 2, npca = 10
  regress.out = "nUMI", reduce.seurat.object = T)

# The object can then be saved and uploaded to the SingleR web-app for further analysis and visualization or using functions available in the SingleR package (see vignette).
save(singler,file=paste0(project.name,'.RData')
```

For more details on creating a SingleR object see [SingleR - create object](http://comphealth.ucsf.edu/sample-apps/SingleR/SingleR_create.html).

For more details about the SingleR method see [SingleR Supplementary Information 1](http://comphealth.ucsf.edu/sample-apps/SingleR/SupplementaryInformation1.html).

For examples of SingleR usage see [SingleR Supplementary Information 2](http://comphealth.ucsf.edu/sample-apps/SingleR/SupplementaryInformation2.html).

The fine-tuning process of SingleR may take very long, and in the current setting is not feasible for large datasets. Thus, for datasets with tens of thousands of cells it is recommended to run SingleR on subsets of the full data and then combine them together. See an example of analyzing 242,533 from the [Microwell-Seq Mouse Cell Atlas](http://comphealth.ucsf.edu/sample-apps/SingleR/SingleR.MCA.html).

# Contributors

SingleR was developed by Dvir Aran. Please contact Dvir Aran: dvir.aran at ucsf edu for any questions or suggestions.

