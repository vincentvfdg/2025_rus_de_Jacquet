# 2025_rus_de_Jacquet

Open source code for the article : Parkinson's disease risk factors are expressed at brain barriers.


## Code architecture

All the repositories, including the code, represent an example of architecture. You can change the architecture with the one you want. 

In the R code, variable can be set for repositories. In this case, `plots/`, `results`, and `data` are empty, but will be filled when running the chunks of the code. 

## Steps 

Raw datasets from Kamath, Feleke, and Smajic were all processed using the in-house package SingleCell. 

Yang was then processed, but only with plexus choroid control samples. 

Quality control were performed on all datasets based on the distribution of features and counts before and after processing. 

A base assignation was made. 

Some assignation were not precised enough on endothelial cells. All but Smajic datasets were subsampled and endothelial cells were reclustered and reassigned. 

Truer endothelial cells were clustered and cells were assign to three groups: Veinous, capillary cells, and arterial cells. 

Graphics were made using multiple homemade functions (see below). 

Fold change were calculated between conditions for every cell types in all datasets. 

## Raw data

In the `raw_data` repository, most raw and process files are available. Provenance of these files are noted in the article. 

In this repository, you will find `gwas.tsv` wich is the raw gwas file from `https://www.ebi.ac.uk/gwas/docs/file-downloads` catalog from GWAS (v1.0.2.1).
### The matrices you need to fetch 

#### Feleke

Every matrices were to heavy to push. You can find them here: 

* C36: Here `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5380921` under `GSM5380921_C36_raw_feature_bc_matrix.tar.gz`. 
* C48: Here `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5380922` under `GSM5380922_C48_raw_feature_bc_matrix.tar.gz`.
* PDC05: Here `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5380944` under `GSM5380944_PDC05_raw_feature_bc_matrix.tar.gz`
* PDC34: Here `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5380946` under `GSM5380946_PDC34_raw_feature_bc_matrix.tar.gz`
* PDC87: Here `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5380947` under `GSM5380947_PDC87_raw_feature_bc_matrix.tar.gz`.
* PDC91: Here `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5380948` under `GSM5380948_PDC91_raw_feature_bc_matrix.tar.gz`.
* PD416: Here `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5380932` under `GSM5380932_PD416_raw_feature_bc_matrix.tar.gz`.
* PD523: Here `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5380934` under `GSM5380934_PD523_raw_feature_bc_matrix.tar.gz`.
* PD683: Here `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5380940` under `GSM5380940_PD683_raw_feature_bc_matrix.tar.gz`
* PDC22: Here `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5380945` under `GSM5380945_PDC22_raw_feature_bc_matrix.tar.gz`
* PD666: Here `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5380938` under `GSM5380938_PD666_raw_feature_bc_matrix.tar.gz`.
* PD747: Here `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5380943` under `GSM5380943_PD747_raw_feature_bc_matrix.tar.gz`
* PD732: Here `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM5380942` under `GSM5380942_PD732_raw_feature_bc_matrix.tar.gz`


#### Kamath

* **matrix.mtx.gz**: This was too heavy to be dropped here. You can find the matrix here: `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE178265` under `GSE178265_Homo_matrix.mtx.gz`.

#### Smajic

The matrix was too heavy as well. It is located in `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157783` under `GSE157783_IPDCO_hg_midbrain_UMI.tar.gz`.

### RDS that would be generated

#### Feleke
* **feleke_cluster_19.rds**: Seurat object containing the cluster 19 (base resolution of 1) of Feleke. This had both endothelial and pericyte cells. 
* **feleke_endothelial.rds**: Seurat object containing the clusterisation of endothelial cells for Feleke. 
* **feleke_preprocessed.rds**: Seurat object of Feleke after our homemade package SingleCell. 
* **feleke_raw.rds**: Raw Seurat object of Feleke. 
* **feleke_reassembled.rds**: Seurat object after the reassignation of barcodes to their proper cell types (pericytes or endothelial). Contains every cell from Feleke. 

#### Kamath
* **kamath_cluster_8.rds**: Seurat object containing the cluster 8 (base resolution of 1) of Kamath. This object contains both endothelial and pericyte cells. 
* **kamath_endothelial.rds**: Seurat object containing the clusterisation of endothelial cells for Kamath.
* **kamath_preprocessed.rds**: Seurat object of Kamath after our homemade package SingleCell. 
* **kamath_raw.rds**: Raw Seurat object of Kamath. 
* **kamath_reassembled.rds**: Seurat object after reassigning barcodes to their proper cell types (pericytes or endothelial). Contains every cell from Kamath. 

#### Smajic
* **smajic_endothelial.rds**: Seurat object with the clusterisation of endothelial cells from Smajic. Reassignation was not performed on Smajic for endothelial cells. 
* **smajic_preprocessed.rds**: Seurat object of Smajic after our homemade package SingleCell. 
* **smajic_raw.rds**: Raw Seurat object for Smajic. 

#### Yang
* **Yang/**: Subrepository containing all matrices, barcodes, and features from the plexus choroid control samples from Yang. 
* **yang_plexus_choroid_control_samples_endothelial.rds**: Seurat object containing the clusterisation of endothelial cells from Yang. 
* **yang_plexus_choroid_control_samples_ependymal.rds**: Seurat object containing the clusterisation of ependymal cells from Yang. 
* **yang_plexus_choroid_control_samples_mesenchymal.rds**: Seurat object containing cluster 15 (base resolution of 1) of Yang. This had both endothelial and mesenchymal cells. 
* **yang_plexus_choroid_control_samples_preprocessed.rds**: Seurat object of Yang after the clusterisation with our homemade package SingleCell. 
* **yang_plexus_choroid_control_samples_reassembled.rds**: Seurat object after reassigning barcodes to their proper cell type (mesenchymal or endothelial). Contains very cell from Yang. 
* **yang_plexus_choroid_raw.rds**: Seurat object containing every cell after CreateSeuratObject with the matrices in inputs (Yang/).

## Functions inside the code

Multiple in-house functions are written in the code. We did not provide the functions, as they are simply decalque of existing functions. 

* **SingleCell::analyze_integrated**: Homemade function from the in-house package. Details are provided in the Material and Method, but here are resumed steps (filtration, normalization (if TRUE), PCA, clusterisation, UMAP). 
* **cell_assignation_automatisation** Homemade function that wraps multiple graphical and table codes. With a list of marker in input and a seurat object, it can perform series_featureplot, a dot plot, a violin plot, and a barplot with the proportion of each cell type in input. 

#### Assignation
* **assignation_cluster**: Assigning clusters from a meta data column into cell assignation based on a list. Will form a new meta data column with the informations in input. 
* **assignation_names**: Takes a subsampled Seurat object and combined it with its original Seurat object. Used after reassignation endothelial cells properly and combining the info into the preprocessed files. 
* **barcodes_reassignation**: Same function as `assignation_names`, but uses barcodes of each cells instead of predetermined clusters. 

#### Graphical and tables output
* **series_featureplot**: An equivalent to Seurat::FeaturePlot() with numerical values. This function provides us more control on aesthetics. 
* **label_clusters_umap**: An equivalent to Seurat::FeaturePlot() with categorical variable. This function provides us more control on aesthetics. The outputs are UMAPs with and without labels. 
* **table_dotplot**: Finds the logarithm expression of each gene in input based on each cell type, with the relative expression. 
* **homemade_dotplot**: Takes the table from `table_dotplot` and draw a dotplot, either with relative expression or log1p(expression). Equivalent to Seurat::DotPlot(). 
* **subsetted_umaps**: Homemade function. This as no equivalent. This subsample a Seurat object on two defined groups in a meta data column and draw a UMAP for each condition for each gene in input. 
* **homemade_dotplot_2FC**: Homemade function. This as no equivalent. Takes the table from `findmarkers_subsetted` and draw a dotplot with the color gradient based on the fold change in log2.  

#### Statistical analyses
* **findmarkers_subsetted**: Calculates the fold change of each cell type in input from a Seurat object based on two treatments. A table is created with each cell type, and a combined table is created with all the concatenated results. 
* **average_expression_annotation**: An equivalent to Seurat::AverageExpression(), which calculates the mean of expression of each cell type in input for all genes in the Seurat object. This function can perform by subsetted Seurat object based on a meta data column and a character string. 


