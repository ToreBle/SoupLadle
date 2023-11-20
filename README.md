# SoupLadle
Here we present an approach to evaluate the pooling of patients with high SNP diversity with SoupLadle. The approach ensures that for the pooling of large cohorts for snRNA-seq experiemtns a high SNP diversity is archieved that subsequently will help to optimally demultiplex the single nuclei with exisiting tools.
Next the demultiplexed nuclei can be assigned to the patients with SoupLadel based on the best matching SNPs from bulk.

<img src="https://github.com/ToreBle/SoupLadle/SoupLadle_Cover.png" width="100">

Processing for pooling patients of a large cohorts:

1. bulkRNA sequencing for each individual
2. SNP calling and genotyping bulkRNA-seq
3. Evaluation and SNP diversity with SoupLadle
4. snRNA-seq pooling and sequencing
5. snRNA-seq demultiplexing with existing methods (e.g. SouporCell, Vireo)
6. Patient assignement with based matching SNPs
7. QC and downstream anaylsis with Seurat or Scanpy

##Installation of SoupLadle

