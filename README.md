# SoupLadle
SoupLadle offers a robust method to assess patient pooling for single-nuclei RNA-seq and assign patients efficiently to demultiplexed single-nuclei clusters. Hence, it can be used to support the pooling and demultiplexing of large single-cell and single-nuclei RNA-seq experiments. SoupLadle aims to achieve optimal SNP diversity within large patient cohorts, ensuring enhanced precision in deconvoluting single nuclei during snRNA-seq experiments.

SoupLadle key features:
 - Enhanced SNP Diversity Detection: Strategically pooling patients with high SNP diversity ensures improved performance in subsequent deconvolution algorithms.
 - Optimal Patient Assignment: Utilizing SoupLadle, we match the best SNPs from bulk datasets to efficiently assign deconvoluted nuclei back to their respective patients.

SoupLadle enables more accurate patient attribution, streamlining the analysis pipeline for snRNA-seq experiments and providing valuable insights into complex cohort datasets.

<img src="https://github.com/ToreBle/SoupLadle/blob/main/SoupLadle_Cover.png" width="500">


## Workflow
Processing for pooling patients of a large cohorts:

1. bulkRNA sequencing for each individual
2. SNP calling and genotyping bulkRNA-seq
3. Evaluation and SNP diversity with SoupLadle
4. snRNA-seq pooling and sequencing
5. snRNA-seq deconvolution with existing methods (e.g. [SouporCell](https://github.com/wheaton5/souporcell), [vireo](https://github.com/single-cell-genetics/vireo))
6. Patient assignment with based matching SNPs
7. QC and downstream anaylsis with Seurat or Scanpy

## Installation
SoupLadle can be directly installed from Github using the following command:

```r
if (!require("devtools", quietly = TRUE)) {
  install.packages("devtools")
}
devtools::install_github("ToreBle/SoupLadle“, build_vignettes = T)
```

The installation takes only take a minute.
The dependencies of the package are listed in the DESCRIPTION file of the package.
Afterwards the vignette, demonstrating the usage of SoupLadle for our PBMC data can be followed by the vignette:

```r
vignette("SoupLadle")
```

## Data
- Raw sequencing data: [GSE247708](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247708)

## Vignette
The vignette for the PBMC data of the manuscript can be found [here](https://github.com/ToreBle/SoupLadle/blob/main/vignettes/SoupLadle.Rmd).

## Reference
The manuscript has been submitted for peer-review
