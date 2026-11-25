# Single-cell RNA-seq analysis: organized and anonymized research draft
#
# All machine-specific absolute paths have been replaced with named configuration
# entries. Comments have been translated to English. Analysis order, object names,
# statistical settings, and exploratory alternatives are retained from the draft.
# This is a collection of analysis blocks, not a validated end-to-end pipeline.
# Run relevant blocks interactively after checking their required input objects.
# Package installation, removal, and R updates remain in the setup block; review
# that block before running it. No packages or analyses were executed during editing.
# One incomplete assignment is commented out and marked TODO.
#
# Path configuration: replace every placeholder below with a local path.
# Repeated references to the same original path share the same configuration entry.
paths <- list(
  # Directory: filter
  path_01 = "<PATH_01>",
  # File: CytoTRACE_0.3.3.tar.gz
  path_02 = "<PATH_02>",
  # Directory: OA1
  path_03 = "<PATH_03>",
  # Directory: OA2
  path_04 = "<PATH_04>",
  # Directory: OA3
  path_05 = "<PATH_05>",
  # Directory: OA4
  path_06 = "<PATH_06>",
  # File: normalizednormal.rds
  path_07 = "<PATH_07>",
  # Directory: normal1
  path_08 = "<PATH_08>",
  # Directory: normal2
  path_09 = "<PATH_09>",
  # Directory: normal3
  path_10 = "<PATH_10>",
  # File: all_normalizednormal.rds
  path_11 = "<PATH_11>",
  # File: allcombined.rds
  path_12 = "<PATH_12>",
  # File: Scissor-master.zip
  path_13 = "<PATH_13>",
  # File: GSE114007_normal_normalized.counts.txt
  path_14 = "<PATH_14>",
  # File: GSE114007_oa_normalized.counts.txt
  path_15 = "<PATH_15>",
  # File: GSEA_allbulk.txt
  path_16 = "<PATH_16>"
)

# =============================================================================
# 01. Environment setup and package dependencies
# =============================================================================

getwd()
setwd(paths$path_01)
samples=list.files()  # List samples in the working directory.
dir <- file.path(samples)  # Construct sample directory paths.
names(dir) <- samples
samples
install.packages("installr")
install.packages("viridis")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.14")
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'Matrix.utils',
                       'HDF5Array', 'terra', 'ggrastr'))
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3')
library(installr)
updateR()
install.packages("ggsignif")
install.packages("Seurat")
install.packages("SeuratData")
library(ggsignif)
library(Seurat)
library(SeuratData)
install.packages('dpylr')
library(patchwork)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c('dpylr'))
BiocManager::install(c('SeuratData'))
library(dplyr)
library(magrittr)
library(limma)
library(monocle3)
library("ggsci")
install.packages("psych")
library(psych)
library(pheatmap)
library(corrplot)
library(ggplot2)
library(cowplot)
devtools::install_github("Coolgenome/iTALK", build_vignettes = TRUE)
library(iTALK)
library(dplyr)
library(Matrix)
devtools::install_github("sqjin/CellChat")
library(CellChat)
library(ggpubr)

##
# Install required R packages.
install.packages("ggplot2")
install.packages("ggridges")
install.packages("reshape")
install.packages("ggprism")
# Load required R packages.
library(ggplot2)
library(ggridges)
library(reshape)
library(ggprism)
install.packages("IDPmisc")
library(IDPmisc)
# install packages from CRAN
cran.packages <- c("msigdbr", "dplyr", "purrr", "stringr","magrittr",
                   "RobustRankAggreg", "tibble", "reshape2", "ggsci",
                   "tidyr", "aplot", "ggfun", "ggplotify", "ggridges",
                   "gghalves", "Seurat", "SeuratObject", "methods",
                   "devtools", "BiocManager","data.table","doParallel",
                   "doRNG")
if (!requireNamespace(cran.packages, quietly = TRUE)) {
  install.packages(cran.packages, ask = F, update = F)
}

# install packages from Bioconductor
bioconductor.packages <- c("GSEABase", "AUCell", "SummarizedExperiment",
                           "singscore", "GSVA", "ComplexHeatmap", "ggtree",
                           "Nebulosa")
if (!requireNamespace(bioconductor.packages, quietly = TRUE)) {
  BiocManager::install(bioconductor.packages, ask = F, update = F)
}
# install packages from Github
if (!requireNamespace("UCell", quietly = TRUE)) {
  devtools::install_github("carmonalab/UCell")
}
if (!requireNamespace("irGSEA", quietly = TRUE)) {
  devtools::install_github("chuiqin/irGSEA")
}
library(UCell)
library(irGSEA)
library(ggcorrplot)  # Install this package before loading it if necessary.
install.packages('ggcorrplot')
devtools::install_local(paths$path_02)


# =============================================================================
# 02. OA samples: data import, integration, and clustering
# =============================================================================


counts <- Read10X(data.dir = dir)
scRNA1 = CreateSeuratObject(counts,min.cells = 3, min.features = 200)
unique(scRNA1$orig.ident)  # Inspect sample identities stored in the metadata.
## [1] Control1 Control2 Model1   Model2  
## Levels: Control1 Control2 Model1 Model2
scRNA1$sample <- scRNA1$orig.ident  # Copy the original sample identities into a dedicated metadata column.
scRNA1 <- FindVariableFeatures(scRNA1, selection.method = "vst", nfeatures = 2000)
scRNA1 <- ScaleData(scRNA1, verbose = FALSE)
scRNA1 <- RunPCA(scRNA1, npcs = 30, verbose = FALSE)
scRNA1 <- RunUMAP(scRNA1, reduction = "pca", dims = 1:30)
scRNA1 <- FindNeighbors(scRNA1, reduction = "pca", dims = 1:30)
scRNA1 <- FindClusters(scRNA1, resolution = 0.4)
p0 <- DimPlot(scRNA1, reduction = "umap", split.by = "sample")
p0
p1 <- DimPlot(scRNA1, reduction = "umap", group.by = "sample")
p1
p2 <- DimPlot(scRNA1, reduction = "umap", label = TRUE, repel = TRUE)
p2
saveRDS(scRNA1, file = "unnormalization.rds")
OA1 <-Read10X(data.dir = paths$path_03)
OA1 <- CreateSeuratObject(OA1,min.cells = 3, min.features = 200)
OA2 <-Read10X(data.dir = paths$path_04)
OA2 <- CreateSeuratObject(OA2,min.cells = 3, min.features = 200)
OA3 <-Read10X(data.dir = paths$path_05)
OA3 <- CreateSeuratObject(OA3,min.cells = 3, min.features = 200)
OA4 <-Read10X(data.dir = paths$path_06)
OA4 <- CreateSeuratObject(OA4,min.cells = 3, min.features = 200)
OA1$replicate <- sample(c("OA1"), size = ncol(OA1), replace = TRUE)
OA2$replicate <- sample(c("OA2"), size = ncol(OA2), replace = TRUE)
OA3$replicate <- sample(c("OA3"), size = ncol(OA3), replace = TRUE)
OA4$replicate <- sample(c("OA4"), size = ncol(OA4), replace = TRUE)
class(OA1)

list <-c(OA1,OA2,OA3,OA4)
# normalize and identify variable features for each dataset independently
list <- lapply(X = list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list)
immune.anchors <- FindIntegrationAnchors(object.list = list, anchor.features = features)
# this command creates an 'integrated' data assay
immune.combined <- IntegrateData(anchorset = immune.anchors)
saveRDS(immune.combined,file ="normalizednormal.rds")
immune.combined <-readRDS(file =paths$path_07)
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
ElbowPlot(immune.combined)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:15)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:15)
immune.combined <- FindClusters(immune.combined, resolution = 0.40)
p3 <- DimPlot(immune.combined, reduction = "umap", group.by = "replicate",pt.size= 1)
p4 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE,pt.size= 1)
p4
p4<-p4+scale_fill_aaas()
p3|p4
immune.combined <- RunTSNE(immune.combined,dims = 1:15)
DimPlot(immune.combined, reduction = "tsne")
DimPlot(immune.combined,reduction = "umap")
table(immune.combined$seurat_clusters1)  # Count cells in each cluster.
#BM1  BM2  BM3  GM1  GM2  GM3 
#2754  747 2158 1754 1528 1983
prop.table(table(Idents(immune.combined)))
table(Idents(immune.combined), immune.combined$replicate)  # Tabulate cell counts by cluster and sample.
Cellratio <- prop.table(table(Idents(immune.combined), immune.combined$group), margin = 2)  # Calculate cluster proportions within each group or sample.
Cellratio <- prop.table(table(Idents(immune.combined), immune.combined$replicate), margin = 2)  # Calculate cluster proportions within each group or sample.
Cellratio
Cellratio <- as.data.frame(Cellratio)
colourCount = length(unique(Cellratio$Var1))
library(ggplot2)
p5<-ggplot(Cellratio) +
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = '#222222')+
  theme_classic() +
  labs(x='Sample',y = 'Ratio')+
  coord_flip()+
  theme(panel.border = element_rect(fill=NA,color="black", size=0.5, linetype="solid"))
p5
oa.markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
oa.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.table(oa.markers,file = "oamalmarkersres0.4_n.txt",sep ='\t',quote = F)
oa.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(immune.combined, features = top10$gene)
DefaultAssay(immune.combined) <- "integrated"
##cell cycle
cc.genes.updated.2019
immune.combined <- CellCycleScoring(
  object = immune.combined,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes)
VlnPlot(immune.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB","G2M.Score","S.Score"), ncol = 6)  #+scale_color_npg()
####cell cycle
umapem<-immune.combined@reductions$umap@cell.embeddings
umapem<-as.data.frame(umapem)
metada= immune.combined@meta.data
dim(umapem);dim(metada)

metada$bar<-rownames(metada)
umapem$bar<-rownames(umapem)
ccdata<-merge(umapem,metada,by="bar")
head(ccdata)
library(ggplot2)
plot<-ggplot(ccdata, aes(UMAP_1, UMAP_2,label=Phase))+geom_point(aes(colour = factor(Phase)))+
#plot<-plot+scale_colour_manual(values=c("#CC33FF","Peru","#660000","#660099","#990033","black","red", "#666600", "green","#6699CC","#339900","#0000FF","#FFFF00","#808080"))+
  labs("@yunlai",x = "", y="")
plot=plot+scale_color_aaas()  +
  theme_bw()+theme(panel.grid=element_blank(),legend.title=element_blank(),legend.text = element_text(color="black", size = 10, face = "bold"))
plot<-plot+guides(colour = guide_legend(override.aes = list(size=5))) +theme(plot.title = element_text(hjust = 0.5))
plot
write.table(ccdata,file ="ccdata.txt",sep = "\t", quote = F)
result <- table(ccdata$integrated_snn_res.0.4,ccdata$Phase)
result

# Explore correlations in average gene expression.
AverageExp<-AverageExpression(immune.combined,features=unique(top10$gene))
typeof(AverageExp)
head(AverageExp$RNA)
coorda<-corr.test(AverageExp$RNA,AverageExp$RNA,method="spearman")
pheatmap(coorda$r)
corrplot(cor(AverageExp$RNA),method = "color", addCoef.col = "grey")
# Alternative correlation analysis.
av <- AverageExpression(immune.combined,features = unique(oa.markers$gene),assays = "RNA")
typeof(av)
av= av[[1]]
head(av)
View(av)  # Review row alignment before continuing.

###

# Explore correlations in average gene expression.
AverageExp<-AverageExpression(immune.combined,features=unique(top10$gene))
typeof(AverageExp)
View(AverageExp$RNA)
cor(AverageExp$RNA,method = "spearman")
pheatmap::pheatmap(cor(AverageExp$RNA,method = "spearman"))

####
DefaultAssay(immune.combined) <- "RNA"
FeaturePlot(immune.combined,features = c("ESR1","SPP1","STAT3","MMP13","SQSTM1"),pt.size= 1)
vlnplot


# =============================================================================
# 03. OA samples: Monocle 3 trajectories and gene modules
# =============================================================================

# Monocle 3 trajectory analysis using the RNA count matrix.
immune.combined1 <- GetAssayData(immune.combined, assay = 'RNA',slot ='counts')
cell_metadata <- immune.combined@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(immune.combined1))
rownames(gene_annotation) <- rownames(immune.combined1)

# Create a Monocle cell_data_set object.
cds <- new_cell_data_set(immune.combined1,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds@colData$replicate<- immune.combined$replicate
# Run the preprocessing steps required by Monocle.

# Perform dimensionality reduction and clustering.
#NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)  # Preprocess the count matrix and compute principal components.
plot_pc_variance_explained(cds)  # Inspect the variance explained by principal components.
cds <- align_cds(cds, alignment_group = "replicate")

# Compute a UMAP embedding.
cds<- reduce_dimension(cds,preprocess_method = "PCA")  # Use PCA as the preprocessing method.


# Compute a t-SNE embedding.
cds <- reduce_dimension(cds, reduction_method="tSNE")
plot_cells(cds, reduction_method="tSNE", color_cells_by="seurat_clusters")

# Cluster cells.


cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(immune.combined, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed  # Reuse the Seurat UMAP coordinates after aligning cell names.
cds@colData$integrated_snn_res.0.4 <- immune.combined$integrated_snn_res.0.4

# Visualize the embedding.
plot_cells(cds, reduction_method="UMAP", color_cells_by="integrated_snn_res.0.4")
cds <- cluster_cells(cds)
cds <- learn_graph(cds)  # Learn the principal graph for trajectory inference.

# Color cells using the selected metadata field.
head(colData(cds))
p_plot_cells<-plot_cells(cds,
                         color_cells_by = "integrated_snn_res.0.4",
                         label_groups_by_cluster=FALSE,
                         label_leaves=FALSE,
                         label_branch_points=TRUE,
                         group_label_size=4,
                         cell_size=0.8)
p_plot_cells <- p_plot_cells+scale_color_npg()
p_plot_cells
# Inspect the inferred cell-state trajectory.
# Black lines represent the principal graph.
# White circles label graph leaves or terminal states.
# Black circles label branch points.
# Control graph labels with label_leaves and label_branch_points.
# Order cells along the trajectory.
cds = order_cells(cds)  # Select a biologically justified root when ordering cells.
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1.5)

# Test for genes whose expression varies along the principal graph.
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=6)
write.csv(Track_genes, "OA_Trajectory_genes_6clusters.csv", row.names = F)
Track_genes <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)

# Select the top ten genes for visualization.
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()

# Plot gene expression trends along pseudotime.
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="integrated_snn_res.0.4",
                         min_expr=0.5, ncol = 2)
p3<-plot_genes_in_pseudotime(cds[c("SQSTM1","COL2A1","NQO1","ACAN","MMP3","ADAMTS5","MMP13","COL10A1"),], color_cells_by="integrated_snn_res.0.4",
                             min_expr=0.5, ncol = 2)
p3 <-p3+scale_color_npg()
p3
# Visualize selected genes on the cell embedding.
plot_cells(cds, genes=c("SQSTM1","COL2A1","NQO1","ACAN","MMP3","ADAMTS5","MMP13","COL10A1"), show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)
# Interpret trajectory-associated genes alongside genes of biological interest.
# Identify modules of coexpressed genes.
genelist <- pull(Track_genes, gene_short_name) %>% as.character()
gene_module <- find_gene_modules(cds[genelist,], resolution=1e-2, cores = 10)
cell_group <- tibble::tibble(cell=row.names(colData(cds)),
                             cell_group=colData(cds)$integrated_snn_res.0.4)
agg_mat <- aggregate_gene_expression(cds, gene_module, cell_group)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2",color = colorRampPalette(colors = c("blue","white","red"))(100))
write.csv(gene_module,"oa_pseuodolgenemodule.csv",row.names = F )
###show module genes
plot_cells(cds_subset2,
           genes=gene_module %>% filter(module %in% c(102, 46)),
           label_cell_groups=FALSE,

# =============================================================================
# 04. GSVA pathway scoring and cell-cycle subsets
# =============================================================================

           show_trajectory_graph=FALSE)


##GSVA

library(GSVA)
install.packages("msigdbr")
library(msigdbr)

subset_G2M <- immune.combined[,immune.combined@meta.data$Phase %in% c("G2M")]  # Subset cells in the selected cell-cycle phase.

Expr_G2M <- as.matrix(subset_G2M@assays$RNA@counts)  # Extract the raw count matrix.
write.csv(Expr_G2M,"Expr_G2M.csv")
library(msigdbr)
library(GSVA)
# Prepare the expression matrix.
expr=as.matrix(object@assays$RNA@data)
# Prepare pathway gene sets.
msgdC2 = msigdbr(species = "Homo sapiens", category = "C2",subcategory = "KEGG")
keggSet = msgdC2 %>% split(x = .$gene_symbol, f = .$gs_description)

kegg <- gsva(Expr_G2M, gset.idx.list = keggSet, kcdf="Poisson",
             parallel.sz=1)  # Prepare the expression matrix.
###diff gene
DefaultAssay(immune.combined) <- "RNA"
immune.combined <- RenameIdents(immune.combined,immune.combined$Phase)
SetIdent(immune.combined,immune.combined$Phase)


# =============================================================================
# 05. Combined normal and OA samples: integration and exploration
# =============================================================================

# Integrate the normal and OA samples.
normal1 <-Read10X(data.dir = paths$path_08)
normal1 <- CreateSeuratObject(normal1,min.cells = 3, min.features = 200)
normal2 <-Read10X(data.dir = paths$path_09)
normal2 <- CreateSeuratObject(normal2,min.cells = 3, min.features = 200)
normal3 <-Read10X(data.dir = paths$path_10)
normal3 <- CreateSeuratObject(normal3,min.cells = 3, min.features = 200)
normal1$replicate <- sample(c("normal1"), size = ncol(normal1), replace = TRUE)
normal2$replicate <- sample(c("normal2"), size = ncol(normal2), replace = TRUE)
normal3$replicate <- sample(c("normal3"), size = ncol(normal3), replace = TRUE)
OA1 <-Read10X(data.dir = paths$path_03)
OA1 <- CreateSeuratObject(OA1,min.cells = 3, min.features = 200)
OA2 <-Read10X(data.dir = paths$path_04)
OA2 <- CreateSeuratObject(OA2,min.cells = 3, min.features = 200)
OA3 <-Read10X(data.dir = paths$path_05)
OA3 <- CreateSeuratObject(OA3,min.cells = 3, min.features = 200)
OA4 <-Read10X(data.dir = paths$path_06)
OA4 <- CreateSeuratObject(OA4,min.cells = 3, min.features = 200)
OA1$replicate <- sample(c("OA1"), size = ncol(OA1), replace = TRUE)
OA2$replicate <- sample(c("OA2"), size = ncol(OA2), replace = TRUE)
OA3$replicate <- sample(c("OA3"), size = ncol(OA3), replace = TRUE)
OA4$replicate <- sample(c("OA4"), size = ncol(OA4), replace = TRUE)
class(OA1)

list <-c(normal1,normal2,normal3,OA1,OA2,OA3,OA4)
# normalize and identify variable features for each dataset independently
list <- lapply(X = list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list)
immune.anchors <- FindIntegrationAnchors(object.list = list, anchor.features = features)
# this command creates an 'integrated' data assay
DefaultAssay(immune.combined) <- "RNA"
immune.combined <- IntegrateData(anchorset = immune.anchors)
saveRDS(immune.combined,file ="all_normalizednormal.rds")
immune.combined <-readRDS(file =paths$path_11)
immune.combined[["percent.mt"]] <- PercentageFeatureSet(immune.combined, pattern ="^MT-")
immune.combined <- subset(immune.combined, subset = nCount_RNA<150000&percent.mt <5)
VlnPlot(immune.combined, features = c("nFeatures_RNA","nCount_RNA","percent.mt"),ncol = 3)
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
ElbowPlot(immune.combined)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:15)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:15)
immune.combined <- FindClusters(immune.combined, resolution = 0.40)
table(immune.combined$orig.ident,immune.combined$replicate)
sample_normal <- sample(c("normal"),size = 9042,replace = TRUE)
sample_oa <- sample(c("oa"),size = 57198,replace = TRUE)
immune.combined$group <- c(sample_normal,sample_oa)

##cell cycle
cc.genes.updated.2019
immune.combined <- CellCycleScoring(
  object = immune.combined,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes)
VlnPlot(immune.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB","G2M.Score","S.Score"), ncol = 6)  #+scale_color_npg()
####cell cycle
umapem<-immune.combined@reductions$umap@cell.embeddings
umapem<-as.data.frame(umapem)
metada= immune.combined@meta.data
dim(umapem);dim(metada)

metada$bar<-rownames(metada)
umapem$bar<-rownames(umapem)
ccdata<-merge(umapem,metada,by="bar")
head(ccdata)
library(ggplot2)
plot<-ggplot(ccdata, aes(UMAP_1, UMAP_2,label=Phase))+geom_point(aes(colour = factor(Phase)))+
#plot<-plot+scale_colour_manual(values=c("#CC33FF","Peru","#660000","#660099","#990033","black","red", "#666600", "green","#6699CC","#339900","#0000FF","#FFFF00","#808080"))+
  labs("@yunlai",x = "", y="")
plot=plot+scale_color_aaas()  +
  theme_bw()+theme(panel.grid=element_blank(),legend.title=element_blank(),legend.text = element_text(color="black", size = 10, face = "bold"))
plot<-plot+guides(colour = guide_legend(override.aes = list(size=5))) +theme(plot.title = element_text(hjust = 0.5))
plot
write.table(ccdata,file ="ccdata.txt",sep = "\t", quote = F)
table(immune.combined$group,immune.combined$Phase)
table(immune.combined$integrated_snn_res.0.4,immune.combined$Phase)
table(immune.combined$group,immune.combined$integrated_snn_res.0.4)
###
as.data.frame(prop.table(table(Idents(immune.combined), immune.combined@meta.data[,"group"]), margin = 2))-> pdf -> td
library(tidyverse)
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
allcolour -> colour1
plt<- ggplot(td,aes(x=td[,2],y=td[,3],fill=td[,1]))+
  geom_bar(position = 'stack',stat="identity")+
  labs(x="replicate",y="Cells Ratio")+
  theme(panel.background=element_rect(fill='transparent', color='black'),
        legend.key=element_rect(fill='transparent', color='transparent'),axis.text = element_text(color="black"))+
  scale_y_continuous(expand=c(0.001,0.001))+
  scale_fill_manual(values=colour1)+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1,ncol=1,title = 'Cluster'))

plt
#############
DefaultAssay(immune.combined) <- "integrated"
p3 <- DimPlot(immune.combined, reduction = "umap", group.by = "replicate",pt.size= 1)
p4 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE,pt.size= 1,cols = c("0" ="#DC143C","1"= "#0000FF","2"= "#20B2AA","3"= "#FFA500","4"= "#9370DB","5"= "#98FB98","6"= "#F08080","7"= "#1E90FF","8"= "#7CFC00","9"= "#FFFF00"))
p4 <- DimPlot(immune.combined, reduction = "umap", split.by = "replicate",pt.size= 1,cols = c("0" ="#DC143C","1"= "#0000FF","2"= "#20B2AA","3"= "#FFA500","4"= "#9370DB","5"= "#98FB98","6"= "#F08080","7"= "#1E90FF","8"= "#7CFC00","9"= "#FFFF00"))
p3|p4
p3
p4
table(immune.combined$integrated_snn_res.0.4,immune.combined$replicate)
####
markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.table(markers,file = "allmarkersres0.4.txt",sep ='\t',quote = F)
markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

DoHeatmap(immune.combined, features = top10$gene,group.colors = c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00"))

##conserved cell type
DefaultAssay(immune.combined) <- "RNA"
conserved_markers <- FindConservedMarkers(immune.combined,ident.1= 0,grouping.var = "group")
conserved_markers0 <-conserved_markers

###

# =============================================================================
# 06. Combined samples: trajectory analysis
# =============================================================================

# Monocle 3 trajectory analysis using the RNA count matrix.
immune.combined1 <- GetAssayData(immune.combined, assay = 'RNA',slot ='counts')
cell_metadata <- immune.combined@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(immune.combined1))
rownames(gene_annotation) <- rownames(immune.combined1)

# Create a Monocle cell_data_set object.
# Create a Monocle cell_data_set object.
cds <- new_cell_data_set(immune.combined1,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds@colData$replicate<- immune.combined$replicate
cds@colData$group<- immune.combined$group
# Run the preprocessing steps required by Monocle.

# Perform dimensionality reduction and clustering.
#NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)  # Preprocess the count matrix and compute principal components.
plot_pc_variance_explained(cds)  # Inspect the variance explained by principal components.
cds <- align_cds(cds, alignment_group = "replicate")

# Compute a UMAP embedding.
cds<- reduce_dimension(cds,preprocess_method = "PCA")  # Use PCA as the preprocessing method.


# Compute a t-SNE embedding.
#cds <- reduce_dimension(cds, reduction_method="tSNE")
#plot_cells(cds, reduction_method="tSNE", color_cells_by="seurat_clusters")

# Cluster cells.


cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(immune.combined, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed  # Reuse the Seurat UMAP coordinates after aligning cell names.
cds@colData$integrated_snn_res.0.4 <- immune.combined$integrated_snn_res.0.4

# Visualize the embedding.
plot_cells(cds, reduction_method="UMAP", color_cells_by="integrated_snn_res.0.4")
cds <- cluster_cells(cds)
cds <- learn_graph(cds)  # Learn the principal graph for trajectory inference.

# Color cells using the selected metadata field.
head(colData(cds))
p_plot_cells<-plot_cells(cds,
                         color_cells_by = "integrated_snn_res.0.4",
                         label_groups_by_cluster=FALSE,
                         label_leaves=FALSE,
                         label_branch_points=TRUE,
                         group_label_size=4,
                         cell_size=0.8)
#p_plot_cells <- p_plot_cells+scale_color_npg()
p_plot_cells
# Inspect the inferred cell-state trajectory.
# Black lines represent the principal graph.
# White circles label graph leaves or terminal states.
# Black circles label branch points.
# Control graph labels with label_leaves and label_branch_points.
# Order cells along the trajectory.
cds = order_cells(cds)  # Select a biologically justified root when ordering cells.
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1.5)

# Test for genes whose expression varies along the principal graph.
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=6)
write.csv(Track_genes, "Trajectory_genes_all.csv", row.names = F)
Track_genes <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)

# Select the top ten genes for visualization.
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()

# Plot gene expression trends along pseudotime.
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="integrated_snn_res.0.4",
                         min_expr=0.5, ncol = 2)
p3<-plot_genes_in_pseudotime(cds[c("SQSTM1","COL2A1","NQO1","ACAN","MMP3","ADAMTS5","MMP13","COL10A1"),], color_cells_by="integrated_snn_res.0.4",
                             min_expr=0.5, ncol = 2)
# Interactive help: ?plot_genes_in_pseudotime()
p3 <-p3+scale_color_npg()
p3
p3 <- p3+scale_color_manual(values = c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00"))
p3
# Visualize selected genes on the cell embedding.
plot_cells(cds, genes=c("SQSTM1","COL2A1","NQO1","ACAN","MMP3","ADAMTS5","MMP13","COL10A1"), show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)
# Interpret trajectory-associated genes alongside genes of biological interest.
# Identify modules of coexpressed genes.
genelist <- pull(Track_genes, gene_short_name) %>% as.character()
gene_module <- find_gene_modules(cds[genelist,], resolution=1e-2, cores = 10)
cell_group <- tibble::tibble(cell=row.names(colData(cds)),
                             cell_group=colData(cds)$integrated_snn_res.0.4)
agg_mat <- aggregate_gene_expression(cds, gene_module, cell_group)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2",color = colorRampPalette(colors = c("blue","white","red"))(100))
write.csv(gene_module,"all_pseuodolgenemodule.csv",row.names = F )
###show module genes
plot_cells(cds_subset2,
           genes=gene_module %>% filter(module %in% c(102, 46)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)


###
conserved_markers1 <- FindConservedMarkers(immune.combined,ident.1= 1,grouping.var = "group")
conserved_markers2 <- FindConservedMarkers(immune.combined,ident.1= 2,grouping.var = "group")
conserved_markers3 <- FindConservedMarkers(immune.combined,ident.1= 3,grouping.var = "group")
conserved_markers4 <- FindConservedMarkers(immune.combined,ident.1= 4,grouping.var = "group")
conserved_markers5<- FindConservedMarkers(immune.combined,ident.1= 5,grouping.var = "group")
conserved_markers6 <- FindConservedMarkers(immune.combined,ident.1= 6,grouping.var = "group")
conserved_markers7 <- FindConservedMarkers(immune.combined,ident.1= 7,grouping.var = "group")
conserved_markers8 <- FindConservedMarkers(immune.combined,ident.1= 8,grouping.var = "group")
conserved_markers9<- FindConservedMarkers(immune.combined,ident.1= 9,grouping.var = "group")
conserved_markersEC <- FindConservedMarkers(immune.combined,ident.1= 'EC',grouping.var = "group")
conserved_markersSPP1 <- FindConservedMarkers(immune.combined,ident.1= 'SPP1+C',grouping.var = "group")
conserved_markerspreHTC <- FindConservedMarkers(immune.combined,ident.1= 'PreHTC',grouping.var = "group")
write.csv(conserved_markersEC,"conserved_markersEC.csv",row.names = T )
write.csv(conserved_markerspreHTC,"conserved_markersPreHTC.csv",row.names = T )
write.csv(conserved_markers0,"conserved_markers0.csv",row.names = T )
write.csv(conserved_markers1,"conserved_markers1.csv",row.names = T )
write.csv(conserved_markers2,"conserved_markers2.csv",row.names = T )
write.csv(conserved_markers3,"conserved_markers3.csv",row.names = T )
write.csv(conserved_markers4,"conserved_markers4.csv",row.names = T )
write.csv(conserved_markers5,"conserved_markers5.csv",row.names = T )
write.csv(conserved_markers6,"conserved_markers6.csv",row.names = T )
write.csv(conserved_markers7,"conserved_markers7.csv",row.names = T )
write.csv(conserved_markers8,"conserved_markers8.csv",row.names = T )
write.csv(conserved_markers9,"conserved_markers9.csv",row.names = T )


####
immune.combined@active.assay
DefaultAssay(immune.combined) <- "integrated"
cluster1.markers <- FindMarkers(immune.combined, ident.1 = 'SPP1+ C', logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster2.markers <- FindMarkers(immune.combined, ident.1 = 'SPP1+ C_2', logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
clusterspp1.markers <- FindMarkers(immune.combined, ident.1 = 2, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
# Interactive help: ?FindMarkers()

####
immune.combined@active.assay
DefaultAssay(immune.combined) <- "RNA"
FeaturePlot(immune.combined,features = c("RHOC","RAN","FHL2","SERPINE2","P3H2","PGAM1","KCNMA1"),pt.size= 1)  ##conserved 0 marker
FeaturePlot(immune.combined,features = c("FST","KCNMA1","SNORC","S100A16","S100A2","S100A10","TNFRSF12A","FRZB"),pt.size= 1)  ##find markers 0 marker
FeaturePlot(immune.combined,features = c("FN1"),pt.size= 1)
VlnPlot(immune.combined,features = c("FN1"),pt.size= 0)
####


# =============================================================================
# 07. Cell communication, correlations, and marker exploration
# =============================================================================

Idents(immune.combined) <- factor(Idents(immune.combined), levels = c("HSPC", "Mono/Mk Doublets",
                                                                      "pDC", "Eryth", "Mk", "DC", "CD14 Mono", "CD16 Mono", "B Activated", "B", "CD8 T", "NK", "T activated",
                                                                      "CD4 Naive T", "CD4 Memory T"))
markers.to.plot <- c("RHOC", "KCNMA1", "CAV1", "CCL20", "MMP3", "INHBA", "SOD2", "STEAP1", "SQSTM1",
                     "FTH1", "RPS13", "FTL", "COL10A1", "IL17B", "SNHG12", "SH3BGRL3", "HBB", "PCOLCE2", "DIO3",
                     "CHI3L2")
DotPlot(immune.combined, features = markers.to.plot, cols = c("blue", "red"), dot.scale = 6, split.by = "group") +
  RotatedAxis()
DotPlot(immune.combined, features =c("S100A4","MMP13","MMP3","SQSTM1","ESR1"), cols = c("blue", "red"), dot.scale = 6, split.by = "group") +
  RotatedAxis()
VlnPlot(immune.combined,features = c("S100A4"),split.by ="group" ,pt.size =0, combine =FALSE)
##iTALK
# Transpose the expression matrix so rows are cells and columns are genes.
immune.combined$seurat_clusters <- immune.combine@active.ident
iTalk_data <- as.data.frame(t(immune.combined@assays$RNA@counts))
# Add the cell_type column required by iTALK.
iTalk_data$cell_type <- immune.combined@meta.data$seurat_clusters
# Add the compare_group column identifying the biological group of each cell.
iTalk_data$compare_group <- immune.combined@meta.data$replicate

unique(iTalk_data$cell_type)
# "cd56_nk" "cd14_monocytes" "b_cells" "cytotoxic_t" "regulatory_t" "memory_t" "naive_t"
unique(iTalk_data$compare_group)
# "group1" "group2" "group3"

my10colors <- my36colors <-c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282')

highly_exprs_genes <- rawParse(iTalk_data, top_genes=50, stats="mean")
# Analyze communication categories with iTALK.
comm_list<-c('growth factor','other','cytokine','checkpoint')
cell_types <- unique(iTalk_data$cell_type)
cell_col <- structure(my10colors[1:length(cell_types)], names=cell_types)

iTalk_res <- NULL
for(comm_type in comm_list){
  res_cat <- FindLR(highly_exprs_genes, datatype='mean count', comm_type=comm_type)
  iTalk_res <- rbind(iTalk_res, res_cat)
}
write.csv(iTalk_res,"iTalk.csv",row.names = F )

iTalk_res <- iTalk_res[order(iTalk_res$cell_from_mean_exprs*iTalk_res$cell_to_mean_exprs,decreasing=T),][1:100,]
NetView(iTalk_res,col=cell_col,vertex.label.cex=1,arrow.width=1,edge.max.width=5)
LRPlot(iTalk_res[1:100,],datatype='mean count',cell_col=cell_col,link.arr.lwd=iTalk_res$cell_from_mean_exprs[1:100],link.arr.width=iTalk_res$cell_to_mean_exprs[1:100])
# Analyze expression correlations.
immune.combined@active.assay
DefaultAssay(immune.combined) <- "RNA"
immune.combined$seurat_clusters1 <- immune.combined@active.ident
av<- AverageExpression(immune.combined,group.by = "seurat_clusters1",assays = "RNA")
av = av[[1]]
head(av)
# Select the 1,000 genes with the largest standard deviations.
cg = names(tail(sort(apply(av,1,sd)),1000))
# Analyze expression correlations.
View(cor(av[cg,],method ='spearman'))
pheatmap::pheatmap(cor(av[cg,],method ='spearman'))


###cell marker

##EC
VlnPlot(immune.combined, features = c("TF","CHRDL2","TPM4","FGF1","RARRES2","DSC2","CYSTM1","KIF9"),pt.size =0)  ##0
FeaturePlot(immune.combined,features = c("TF","CHRDL2","FST","CRYAB","KIF9"),pt.size = 1)
FeaturePlot(immune.combined,features = c("TF","CHRDL2","FST","CRYAB","KIF9"),pt.size = 1)
FeaturePlot(immune.combined,features = c("JUN"),pt.size = 1,min.cutoff = 0)
VlnPlot(immune.combined, features = c("TF","CHRDL2","FST","CRYAB"),pt.size =0,group.by = "group")
VlnPlot(immune.combined, features = c("KIF9"),pt.size =0)
###RegC
VlnPlot(immune.combined, features = c("CKS2","TUBB4B","HMOX1","TMEM176A"),pt.size =0)  # 7
FeaturePlot(immune.combined,features = c("TUBB4B","HMOX1"),pt.size = 1,min.cutoff =2)
FeaturePlot(immune.combined,features = c("HSPA1A","HSPA1B","CKS2","TUBB4B"),pt.size = 1,min.cutoff =0)
FeaturePlot(immune.combined,features = c("TUBB4B"),pt.size = 1,min.cutoff =1.5)
FeaturePlot(immune.combined,features = c("CKS2"),pt.size = 1)
###ProC
VlnPlot(immune.combined, features = c("P3H2","DDX21","GLIS3","ATP13A3","NOP16","MTHFD2","TXNRD1","ACTN1","ASNS","ERCC1","FAM101B","BMP6"),pt.size =0)  #5
FeaturePlot(immune.combined,features = c("P3H2"),pt.size = 1,min.cutoff =2)
FeaturePlot(immune.combined,features = c("UPP1"),pt.size = 1,min.cutoff =0)
FeaturePlot(immune.combined,features = c("TOP2A"),pt.size = 1,min.cutoff =0)
FeaturePlot(immune.combined,features = c("MKI67"),pt.size = 1,min.cutoff =0)
##PreHTC
VlnPlot(immune.combined, features = c("TGFBI","S100A4","TNC","ABI3BP","CRTAC1","TNFAIP6","CRLF1"),pt.size =0)  #4,3
FeaturePlot(immune.combined,features = c("TGFBI","S100A4","TNC","ABI3BP"),pt.size = 1)
FeaturePlot(immune.combined,features = c("CHI3L2","CHI3L1"),pt.size = 1)
FeaturePlot(immune.combined,features = c("MMP1"),pt.size = 1,min.cutoff =1)
###FC
###FC
VlnPlot(immune.combined, features = c('COL1A1',	'COL1A2',	'TMSB4X',	'MTRNR2L2',	'S100A4',	'MGST3',	'MTRNR2L8',	'SERPINF1',	'ID3',	'IFITM2',	'MMP2',	'THY1',	'LGALS3BP',	'IFI16',	'PSME1',	'IGFBP7',	'IFI27',	'COL5A1',	'MTRNR2L1'),pt.size =0)  #4
VlnPlot(immune.combined, features = c('OGN'),pt.size =0)
FeaturePlot(immune.combined,features = c('OGN',	'COL1A2',	'ZFP36L1'),pt.size = 4)
FeaturePlot(immune.combined,features =c('FN1'),pt.size = 1,min.cutoff = 4)
FeaturePlot(immune.combined,features =c('S100A4'),pt.size = 1,min.cutoff =1)
FeaturePlot(immune.combined,features =c('OGN'),pt.size = 1,min.cutoff =1)
FeaturePlot(immune.combined,features =c('RUNX2'),pt.size = 1,min.cutoff =0)
FeaturePlot(immune.combined,features =c('ASPN'),pt.size = 1,min.cutoff =0)
FeaturePlot(immune.combined,features =c('TMSB4X'),pt.size = 1,min.cutoff =0)
FeaturePlot(immune.combined,features =c('FN1'),pt.size = 1,min.cutoff =4.5)
markers <- FindAllMarkers(object = immune.combined,
                          only.pos = TRUE,
                          logfc.threshold = 0.25)

FeaturePlot(immune.combined,features =c('DKK3','LOXL4'),pt.size = 1,min.cutoff = 0)
##HTC
VlnPlot(immune.combined, features = c('RNU6-2',	'WWP2',	'COL10A1',	'BHLHE41',	'IBSP',	'RMRP',	'PLCG2',	'SCIN',	'S100A1',	'MEF2C',	'CRISPLD1',	'VTRNA1-1'),pt.size =0)  #6
FeaturePlot(immune.combined,features =c('WWP2',	'COL10A1',	'BHLHE41',	'IBSP',	'RMRP',	'PLCG2',	'SCIN',	'S100A1'),pt.size = 1)
FeaturePlot(immune.combined,features =c('COL10A1',	'S100A1'),pt.size = 1)
FeaturePlot(immune.combined,features =c('COL11A1'),pt.size = 1,min.cutoff = 1)
FeaturePlot(immune.combined,features =c('IBSP'),pt.size = 1,min.cutoff = 0)
##HomC
VlnPlot(immune.combined, features = c('JUN',	'PPP1R10',	'BRD2',	'ID3',	'ID1',	'RGS16',	'HES1',	'NR4A1',	'RSRP1',	'TRA2B',	'RSRC2',	'NR4A2',	'CLK1',	'RASD1',	'CCNL1',	'IER2',	'NFKBIA',	'SNAI1',	'HEXIM1',	'KLF2'),pt.size =0)  #7
FeaturePlot(immune.combined,features = c('JUN','FOS'),pt.size = 1,min.cutoff = 1)
FeaturePlot(immune.combined,features = c('JUN'),pt.size = 1,min.cutoff = 1)
FeaturePlot(immune.combined,features = c('RSRC2'),pt.size = 1,min.cutoff = 1)
##HBB HBA1+ C
VlnPlot(immune.combined, features = c('HBB','HBA1'),pt.size =0,y.max =10,sort = TRUE)  #Red Blood Cell 8
FeaturePlot(immune.combined,features = c("HBB","HBA1","HBA2"),min.cutoff =1,pt.size = 1)

FeaturePlot(immune.combined,features = c("MMP13"),min.cutoff =1,pt.size = 1)
FeaturePlot(immune.combined,features = c("ADAMTS5"),min.cutoff =1,pt.size = 1)
FeaturePlot(immune.combined,features = c("PTGS2"),min.cutoff =1,pt.size = 1)
# Interactive help: ?FeaturePlot()
##spp1+ cartilage
VlnPlot(immune.combined, features = c('SPP1','OGN'),pt.size =0)
FeaturePlot(immune.combined,features = c('SPP1','OGN'),pt.size = 1,min.cutoff =1)  #1,2   
#9 undecided
#rename
new.cluster.ids <- c("EC", "PreHTC", "SPP1+C", "PreHTC", "FC", "ProC",
                     "HTC", "RegC/HomoC", "RBC","undecidedC")
names(new.cluster.ids) <- levels(immune.combined)
immune.combined <- RenameIdents(immune.combined, new.cluster.ids)
immune.combined <- RenameIdents(immune.combined, `0` = "EC", `1` = "SPP1+ C_1", `2` = "SPP1+ C_2",
                                `3` = "PreHTC_1", `4` = "PreHTC_2/FC", `5` = "ProC", `6` = "HTC", `7` = "RegC/HomoC", `8` = "RBC", `9` = "Undecided")
immune.combined <- RenameIdents(immune.combined, "EC"="0" , "SPP1+ C_1"="1"  ,"SPP1+ C_2"= "2" , "PreHTC_1"="3", "PreHTC_2/FC"="4", "ProC"="5",  "HTC"="6", "RegC/HomoC"="7", "RBC"="8", "Undecided"="9")
immune.combined <- RenameIdents(immune.combined, "EC"="0", "SPP1+ C_1"=`1`  ,"SPP1+ C_2"= `2` , "PreHTC_1"=`3`, "PreHTC_2/FC"=`4`, "ProC"=`5`,  "HTC"=`6`, "RegC/HomoC"=`7`, "RBC"=`8`, "Undecided"=`9`)
immune.combined@active.ident <- immune.combined$seurat_clusters
DimPlot(immune.combined, label = TRUE,pt.size = 1)
P2<-DimPlot(immune.combined, reduction = "umap", pt.size= 1,cols = c("EC" ="#DC143C","PreHTC"= "#0000FF","SPP1+C"= "#20B2AA","FC"= "#FFA500","ProC"= "#9370DB","HTC"= "#98FB98","RegC/HomoC"= "#F08080","RBC"= "#1E90FF","undecidedC"= "#7CFC00"))
##COL5A2
VlnPlot(immune.combined, features = c('C1orf56','COL5A2'),pt.size =0,sort = TRUE)  #?
p1<-VlnPlot(immune.combined, features = 'senmayo1',pt.size =0,group.by = "group")
p1<-VlnPlot(immune.combined, features = 'Angio1',pt.size =0,sort = TRUE)  #?
p1<-p1+geom_boxplot()
p1
p2<-VlnPlot(immune.combined, features = 'senmayo1',pt.size =0,sort = TRUE)  #?
p2<-p2+geom_boxplot()
p2

VlnPlot(immune.combined, features = 'gene_upautp1',pt.size =0,sort = TRUE)  #
VlnPlot(immune.combined, features = 'senmayo1',pt.size =0,sort = TRUE)  #
# Interactive help: ?VlnPlot()
VlnPlot(immune.combined, features = c('SPP1','OGN'),pt.size =0)
VlnPlot(immune.combined, features = c('IBSP'),pt.size =0,sort = T)
VlnPlot(immune.combined, features = c('ACAN','COL2A1'),pt.size =0,sort = T)

VlnPlot(immune.combined, features = c('SPP1','OGN'),pt.size =0)
FeaturePlot(immune.combined,features = c('SPP1','CXCL12'),pt.size = 1)
FeaturePlot(immune.combined,features = c('SPP1','HBB','HBA1','HBA2'),pt.size = 1)
FeaturePlot(immune.combined,features = c('SPP1','OGN','COL10A1','ALPL'))
FeaturePlot(immune.combined,features = c('PECAM1', 'VWF'),pt.size = 1)
FeaturePlot(immune.combined,features = c('HBB','HBA1','SLC4A1'),pt.size = 1)
FeaturePlot(immune.combined,features = c('ESR1'),pt.size = 1)
VlnPlot(immune.combined,features = c('ESR1'),pt.size = 1,sort = T)
###CellChat
data.input <- GetAssayData(immune.combined, assay = "RNA", slot = "data")  # normalized data matrix
labels <- Idents(immune.combined)
identity <- data.frame(group = labels, row.names = names(labels))  # create a dataframe of the cell labels
cellchat <- createCellChat(object = data.input)

cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels")  # set "labels" as default cell identity
levels(cellchat@idents)  # show factor levels of the cell labels

##save
saveRDS(immune.combined,file ="allcombined.rds")
immune.combined <-readRDS(file ="allcombined.rds")

VlnPlot(immune.combined,features = c('COL10A1'),pt.size = 1,sort = T)
FeaturePlot(immune.combined,features = c('CXCL8'),pt.size = 1,min.cutoff = 3)
FeaturePlot(immune.combined,features = c('ESR1'),pt.size = 1,order = T)

##MARKERS
cluster8.markers <- FindMarkers(immune.combined,ident.1=8, logfc.threshold =  0.25,test.use = "roc",only.pos = TRUE)
cluster0.markers <- FindMarkers(immune.combined,ident.1=0, logfc.threshold =  0.25,test.use = "roc",only.pos = TRUE)
cluster7.markers <- FindMarkers(immune.combined,ident.1=7, logfc.threshold =  0.25,test.use = "roc",only.pos = TRUE)
write.csv(cluster8.markers,"cluster8.csv",row.names = T )
write.csv(markers,"markersn.csv",row.names = T )
###DEG
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
EC <- subset(immune.combined, idents = "EC")
Idents(EC) <- "group"
avg.t.EC <- as.data.frame(log1p(AverageExpression(EC, verbose = FALSE)$RNA))
avg.t.EC$gene <- rownames(EC)
write.csv(avg.t.EC,"EC_diff.csv",row.names = F )

SPP1_C_1 <- subset(immune.combined, idents = "SPP1+ C_1")
Idents(SPP1_C_1) <- "group"
avg.SPP1_C_1 <- as.data.frame(log1p(AverageExpression(SPP1_C_1, verbose = FALSE)$RNA))
avg.SPP1_C_1$gene <- rownames(avg.SPP1_C_1)
write.csv(avg.SPP1_C_1,"SPP1_C_1.csv",row.names = F )
FeaturePlot(SPP1_C_1,features = c('SPP1'),pt.size = 1)

SPP1_C_2 <- subset(immune.combined, idents = "SPP1+ C_2")
Idents(SPP1_C_2) <- "group"
avg.SPP1_C_2 <- as.data.frame(log1p(AverageExpression(SPP1_C_2, verbose = FALSE)$RNA))
avg.SPP1_C_2$gene <- rownames(avg.SPP1_C_2)
write.csv(avg.SPP1_C_2,"SPP1_C_2.csv",row.names = F )
FeaturePlot(SPP1_C_2,features=c('SPP1'),split.by = 'group',pt.size = 1)
VlnPlot(SPP1_C_2,features=c('SPP1'),split.by = 'group',pt.size = 0)

PreHTC <- subset(immune.combined, idents = "PreHTC")
Idents(PreHTC) <- "group"
avg.PreHTC <- as.data.frame(log1p(AverageExpression(PreHTC, verbose = FALSE)$RNA))
avg.PreHTC$gene <- rownames(avg.PreHTC)
write.csv(avg.PreHTC,"PreHTC.csv",row.names = F )

PreHTC_2_FC <- subset(immune.combined, idents = "PreHTC_2/FC")
Idents(PreHTC_2_FC) <- "group"
avg.PreHTC_2_FC <- as.data.frame(log1p(AverageExpression(PreHTC_2_FC, verbose = FALSE)$RNA))
avg.PreHTC_2_FC$gene <- rownames(avg.PreHTC_2_FC)
write.csv(avg.PreHTC_2_FC,"PreHTC_2_FC.csv",row.names = F )

ProC <- subset(immune.combined, idents = "ProC")
Idents(ProC) <- "group"
avg.ProC <- as.data.frame(log1p(AverageExpression(ProC, verbose = FALSE)$RNA))
avg.ProC$gene <- rownames(avg.ProC)
write.csv(avg.ProC,"ProC.csv",row.names = F )

HTC <- subset(immune.combined, idents = "HTC")
Idents(HTC) <- "group"
avg.HTC <- as.data.frame(log1p(AverageExpression(HTC, verbose = FALSE)$RNA))
avg.HTC$gene <- rownames(avg.HTC)
write.csv(avg.HTC,"HTC.csv",row.names = F )

RegC_HomoC <- subset(immune.combined, idents = "RegC/HomoC")
Idents(RegC_HomoC) <- "group"
avg.RegC_HomoC <- as.data.frame(log1p(AverageExpression(RegC_HomoC, verbose = FALSE)$RNA))
avg.RegC_HomoC$gene <- rownames(avg.RegC_HomoC)
write.csv(avg.RegC_HomoC,"RegC_HomoC.csv",row.names = F )
##
FeaturePlot(immune.combined,)
VlnPlot(immune.combined,features = c('HBB','HBA1'),slot = "counts",split.by = "group")
VlnPlot(immune.combined,features = c('VEGFA'),split.by = "group",pt.size = 0)
VlnPlot(immune.combined,features = c('FN1'),split.by = "group",pt.size = 0)
VlnPlot(immune.combined,features = c('THY1'),split.by = "group",pt.size = 0)
VlnPlot(immune.combined,features = c ('BIRC5', 'CENPU','UNE2C', 'DHFR','STMN1'),pt.size = 0)
FeaturePlot(immune.combined,features =c ("nFeature_RNA"))
FeaturePlot(immune.combined,features =c ("THY1"),split.by = "group",pt.size = 1)
FeaturePlot(immune.combined,features =c ("S100A4"),split.by = "group",pt.size = 1)
FeaturePlot(immune.combined,features =c ('BIRC5', 'CENPU','UNE2C', 'DHFR','STMN1'),pt.size = 1)
FeaturePlot(immune.combined,features =c ('CD44'),split.by = "group",pt.size = 1)

FeaturePlot(immune.combined,features =c ('PECAM1'),split.by = "group",pt.size = 1)

genes.to.label = c("ISG15", "LY6E", "IFI6", "ISG20", "MX1", "IFIT2", "IFIT1", "CXCL10", "CCL8")
p1 <- ggplot(avg.t.cells, aes(CTRL, STIM)) + geom_point() + ggtitle("CD4 Naive T Cells")
p1 <- LabelPoints(plot = p1, points = genes.to.label, repel = TRUE)
p2 <- ggplot(avg.cd14.mono, aes(CTRL, STIM)) + geom_point() + ggtitle("CD14 Monocytes")
p2 <- LabelPoints(plot = p2, points = genes.to.label, repel = TRUE)
p1 + p2

####subset new normal chondro

subset_cartilage <- subset(immune.combined,idents = c("ProC", "RegC","PreHTC",
                                                      "HomoC/RegC","EC"))
DefaultAssay(subset_cartilage) <- "integrated"

subset_cartilage <- ScaleData(subset_cartilage, verbose = FALSE)
subset_cartilage <- RunPCA(subset_cartilage, npcs = 30, verbose = FALSE)
ElbowPlot(subset_cartilage)
subset_cartilage <- RunUMAP(subset_cartilage, reduction = "pca", dims = 1:15)
subset_cartilage <- FindNeighbors(subset_cartilage, reduction = "pca", dims = 1:15)
subset_cartilage <- FindClusters(subset_cartilage, resolution = 0.4)
DimPlot(subset_cartilage,reduction = "umap")

DefaultAssay(subset_cartilage) <- "RNA"
FeaturePlot(subset_cartilage,features = c("TF","CHRDL2","FST","CRYAB","KIF9"),pt.size = 1)
FeaturePlot(subset_cartilage,features = c("FOS","JUN","CKS2","TUBB4B"),pt.size = 1)
FeaturePlot(subset_cartilage,features = c("TGFBI","S100A4"),pt.size = 1)
DefaultAssay(immune.combined2) <- "RNA"
FeaturePlot(immune.combined2,features =c ('COL1A2','OGN'),pt.size = 1,order = T)
FeaturePlot(immune.combined2,features =c ('HMOX1','KEAP1','NFE2L2','NQO1','SQSTM1'),pt.size = 1,)
new.cluster.ids <- c("ProC", "ProC", "EC", "PreHTC", "EC","HomoC/RegC")
names(new.cluster.ids) <- levels(subset_cartilage)
subset_cartilage <- RenameIdents(subset_cartilage, new.cluster.ids)
FeaturePlot(subset_cartilage,features =c ('HMOX1','KEAP1','NFE2L2','NQO1','SQSTM1'),pt.size = 1)
# Interactive help: ?FeaturePlot()

# =============================================================================
# 08. Additional trajectory exploration
# =============================================================================

# Monocle 3 trajectory analysis using the RNA count matrix.
immune.combined2 <- immune.combined1
immune.combined <- subset_cartilage
immune.combined1 <- GetAssayData(immune.combined, assay = 'RNA',slot ='counts')
cell_metadata <- immune.combined@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(immune.combined1))
rownames(gene_annotation) <- rownames(immune.combined1)

# Create a Monocle cell_data_set object.
# Create a Monocle cell_data_set object.
cds <- new_cell_data_set(immune.combined1,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds@colData$replicate<- immune.combined$replicate
cds@colData$group<- immune.combined$group
# Run the preprocessing steps required by Monocle.

# Perform dimensionality reduction and clustering.
#NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)  # Preprocess the count matrix and compute principal components.
plot_pc_variance_explained(cds)  # Inspect the variance explained by principal components.
cds <- align_cds(cds, alignment_group = "replicate")

# Compute a UMAP embedding.
cds<- reduce_dimension(cds,preprocess_method = "PCA")  # Use PCA as the preprocessing method.


# Compute a t-SNE embedding.
#cds <- reduce_dimension(cds, reduction_method="tSNE")
#plot_cells(cds, reduction_method="tSNE", color_cells_by="seurat_clusters")

# Cluster cells.


cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(immune.combined, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed  # Reuse the Seurat UMAP coordinates after aligning cell names.
cds@colData$seurat_clusters <- immune.combined@active.ident

# Visualize the embedding.
plot_cells(cds, reduction_method="UMAP", color_cells_by="seurat_clusters")
cds <- cluster_cells(cds)
cds <- learn_graph(cds)  # Learn the principal graph for trajectory inference.

# Color cells using the selected metadata field.
head(colData(cds))
p_plot_cells<-plot_cells(cds,
                         color_cells_by = "seurat_clusters",
                         label_groups_by_cluster=FALSE,
                         label_leaves=FALSE,
                         label_branch_points=TRUE,
                         group_label_size=4,
                         cell_size=0.8)
#p_plot_cells <- p_plot_cells+scale_color_npg()
p_plot_cells
# Inspect the inferred cell-state trajectory.
# Black lines represent the principal graph.
# White circles label graph leaves or terminal states.
# Black circles label branch points.
# Control graph labels with label_leaves and label_branch_points.
# Order cells along the trajectory.
cds = order_cells(cds)  # Select a biologically justified root when ordering cells.
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1)
###add module score
gene_angio<-read.table("genelist_angio.txt",header = T)
gene<- as.list(gene_angio)
immune.combined <- AddModuleScore(object = immune.combined,features = gene,name = "Angio")
p1<-VlnPlot(immune.combined,features = 'Angio1',split.by = "group",pt.size = 0)

VlnPlot(immune.combined,features = 'Angio1',group.by = "group",pt.size = 0)
my_comparisons <- list(c('normal','oa'))
p.AddModuleScore <- ggviolin(immune.combined@meta.data, x = "group", y = "Angio1",
                             color = "group",add = 'mean_sd',fill = 'group',
                             add.params = list(color = "black")) +
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") +
#scale_color_manual(values = my9color) + 
#scale_fill_manual(values = my9color) +
  theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)) +
  NoLegend() + labs(x = '')

p.AddModuleScore
p1<-VlnPlot(immune.combined,features = 'Angio1',pt.size = 0)
p1<-p1+

p8 <- ggviolin(immune.combined@meta.data, x = "group", y = "Angio1",fill = "group",
               palette = c("#3B52FF", "#FF4319"),
               add = "boxplot", add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") +
#scale_color_manual(values = my9color) + 
#scale_fill_manual(values = my9color) +
  theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)) +
  NoLegend() + labs(x = '')


p8
immune.combined$seurat_clusters1 <- immune.combined@active.ident

p9 <- ggviolin(immune.combined@meta.data, x = "seurat_clusters1", y = "senmayo1",fill = "group",
              palette = c("#3B52FF", "#FF4319"),
              add = "boxplot") + stat_compare_means(aes(group = group),method = "wilcox.test",label = "p.signif")
#stat_compare_means(comparisons = my_comparisons,label = "p.signif") + 
#scale_color_manual(values = my9color) + 
#scale_fill_manual(values = my9color) +
  theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)) +
  NoLegend() + labs(x = '')
p9

p9 <- ggviolin(immune.combined@meta.data, x = "seurat_clusters1", y = "Angio1",fill = "group",
               palette = c("#3B52FF", "#FF4319"),
               add = "boxplot") + stat_compare_means(aes(group = group),method = "wilcox.test",label = "p.signif")
#stat_compare_means(comparisons = my_comparisons,label = "p.signif") + 
#scale_color_manual(values = my9color) + 
#scale_fill_manual(values = my9color) +
theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)) +
  NoLegend() + labs(x = '')
p9

p10 <- ggviolin(immune.combined@meta.data, x = "seurat_clusters1", y = "gene_dnautp1",fill = "group",
               palette = c("#3B52FF", "#FF4319"),
               add = "boxplot") + stat_compare_means(aes(group = group),method = "wilcox.test",label = "p.signif")
p10


gene_dnangio<-read.table("genelist_downangio.txt",header = T)
gene<- as.list(gene_dnangio)
immune.combined <- AddModuleScore(object = immune.combined,features = gene,name = "dnangio")
VlnPlot(immune.combined,features = 'dnangio1',split.by = "group",pt.size = 0)
VlnPlot(immune.combined,features = 'dnangio1',pt.size = 0)
VlnPlot(immune.combined,features = 'dnangio1',group.by = "group",pt.size = 0)
p.AddModuleScore_d <- ggviolin(immune.combined@meta.data, x = "group", y = "dnangio1",
                             color = "group",add = 'mean_sd',fill = 'group',
                             add.params = list(color = "black")) +
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") +
#scale_color_manual(values = my9color) + 
#scale_fill_manual(values = my9color) +
  theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)) +
  NoLegend() + labs(x = '')

p.AddModuleScore_d


gene_angiore<-read.table("genelist_angiore.txt",header = T)
gene<- as.list(gene_angiore)
immune.combined <- AddModuleScore(object = immune.combined,features = gene,name = "Angiore")
VlnPlot(immune.combined,features = 'Angiore1',split.by = "group",pt.size = 0)
VlnPlot(immune.combined,features = 'Angiore1',group.by = "group",pt.size = 0)
my_comparisons <- list(c('normal','oa'))
p.AddModuleScore_re <- ggviolin(immune.combined@meta.data, x = "group", y = "Angiore1",
                             color = "group",add = 'mean_sd',fill = 'group',
                             add.params = list(color = "black")) +
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") +
#scale_color_manual(values = my9color) + 
#scale_fill_manual(values = my9color) +
  theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)) +
  NoLegend() + labs(x = '')
p.AddModuleScore_re


gene_angiopro<-read.table("genelist_angiopro.txt",header = T)
gene<- as.list(gene_angiopro)
immune.combined <- AddModuleScore(object = immune.combined,features = gene,name = "Angiopro")
VlnPlot(immune.combined,features = 'Angiopro1',split.by = "group",pt.size = 0)
P2<-VlnPlot(immune.combined,features = 'Angiopro1',group.by = "group",pt.size = 0)
P2
p.AddModuleScore_pro <- ggviolin(immune.combined@meta.data, x = "group", y = "Angiopro1",
                                color = "group",add = 'mean_sd',fill = 'group',
                                add.params = list(color = "black")) +
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") +
#scale_color_manual(values = my9color) + 
#scale_fill_manual(values = my9color) +
  theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)) +
  NoLegend() + labs(x = '')
p.AddModuleScore_pro


gene_angiode<-read.table("genelist_angiode.txt",header = T)
gene<- as.list(gene_angiode)
immune.combined <- AddModuleScore(object = immune.combined,features = gene,name = "Angiode")
VlnPlot(immune.combined,features = 'Angiode1',split.by = "group",pt.size = 0)
P1<-VlnPlot(immune.combined,features = 'Angiode1',group.by = "group",pt.size = 0)
p.AddModuleScore_de <- ggviolin(immune.combined@meta.data, x = "group", y = "Angiode1",
                                 color = "group",add = 'mean_sd',fill = 'group',
                                 add.params = list(color = "black")) +
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") +
#scale_color_manual(values = my9color) + 
#scale_fill_manual(values = my9color) +
  theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)) +
  NoLegend() + labs(x = '')

p.AddModuleScore_de

gene_senmayo<-read.table("genelist_senmayo.txt",header = T)
gene<- as.list(gene_senmayo)
immune.combined <- AddModuleScore(object = immune.combined,features = gene,name = "senmayo")
VlnPlot(immune.combined,features = 'senmayo1',split.by = "group",pt.size = 0)
VlnPlot(immune.combined,features = 'senmayo1',group.by = "group",pt.size = 0)
p.AddModuleScore_sen <- ggviolin(immune.combined@meta.data, x = "group", y = "senmayo1",
                                color = "group",add = 'mean_sd',fill = 'group',
                                add.params = list(color = "black")) +
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") +
#scale_color_manual(values = my9color) + 
#scale_fill_manual(values = my9color) +
  theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)) +
  NoLegend() + labs(x = '')

p.AddModuleScore_sen

gene_senmayo_n<-read.table("genelist_senmayo_n.txt",header = T)
gene<- as.list(gene_senmayo_n)
immune.combined <- AddModuleScore(object = immune.combined,features = gene,name = "senmayo_n")
VlnPlot(immune.combined,features = 'senmayo_n1',split.by = "group",pt.size = 0)
VlnPlot(immune.combined,features = 'senmayo_n1',pt.size = 0)
VlnPlot(immune.combined,features = 'senmayo_n1',group.by = "group",pt.size = 0)
p.AddModuleScore_sen <- ggviolin(immune.combined@meta.data, x = "group", y = "senmayo_n1",
                                 color = "group",add = 'mean_sd',fill = 'group',
                                 add.params = list(color = "black")) +
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") +
#scale_color_manual(values = my9color) + 
#scale_fill_manual(values = my9color) +
  theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)) +
  NoLegend() + labs(x = '')

p.AddModuleScore_sen

gene_upautp<-read.table("genelist_upautophagy.txt",header = T)
gene<- as.list(gene_upautp)
immune.combined <- AddModuleScore(object = immune.combined,features = gene,name = "gene_upautp")
VlnPlot(immune.combined,features = 'gene_upautp1',split.by = "group",pt.size = 0)
p1<-VlnPlot( immune.combined,features = 'gene_upautp1',group.by = "group",pt.size = 0)
p1 <- p1+geom_boxplot()+stat_compare_means(method ="wilcox.test")
p1
VlnPlot(immune.combined,features = 'gene_upautp1',pt.size = 0,sort = TRUE)
VlnPlot(immune.combined,features = 'gene_dnautp1',pt.size = 0,sort = TRUE)
p2<-VlnPlot(immune.combined,features = 'gene_dnautp1',pt.size = 0,group.by = "group")
p2 <- p2+geom_boxplot()+stat_compare_means(method ="wilcox.test")
p2
p.AddModuleScore_upautp <- ggviolin(immune.combined@meta.data, x = "group", y = "gene_upautp1",
                                 color = "group",add = 'mean_sd',fill = 'group',
                                 add.params = list(color = "black")) +
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") +
#scale_color_manual(values = my9color) + 
#scale_fill_manual(values = my9color) +
  theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)) +
  NoLegend() + labs(x = '')

gene_dnautp<-read.table("genelist_dnautophagy.txt",header = T)
gene<- as.list(gene_dnautp)
immune.combined <- AddModuleScore(object = immune.combined,features = gene,name = "gene_dnautp")
VlnPlot(immune.combined,features = 'gene_dnautp1',split.by = "group",pt.size = 0)
VlnPlot(immune.combined,features = 'gene_dnautp1',group.by = "group",pt.size = 0)
p.AddModuleScore_dnautp <- ggviolin(immune.combined@meta.data, x = "group", y = "gene_dnautp1",
                                    color = "group",add = 'mean_sd',fill = 'group',
                                    add.params = list(color = "black")) +
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") +
#scale_color_manual(values = my9color) + 
#scale_fill_manual(values = my9color) +
  theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)) +
  NoLegend() + labs(x = '')


VlnPlot(immune.combined,features = 'G2M.Score',split.by = "group",pt.size = 0)
VlnPlot(immune.combined,features = 'G2M.Score',group.by = "group",pt.size = 0)
VlnPlot(immune.combined,features = 'S.Score',group.by = "group",pt.size = 0)

p0<-VlnPlot(immune.combined,features = c('SPP1'),group.by = "group",pt.size = 0)
p0 <-p0+scale_fill_manual(values = c("#3B52FF", "#FF4319"))
p0
VlnPlot(immune.combined,features = c('SPP1'),pt.size = 0)
VlnPlot(immune.combined,features = c('FN1'),group.by = "group",pt.size = 0)
VlnPlot(immune.combined,features = c('DCN','SPARC','THBS1','PGK1','MEG3','RGCC','RHOA','KLF4','CTNNB1','FOXC1','SARS','GTF2I'),group.by = "group",pt.size = 0)
VlnPlot(immune.combined,features = 'vegfa1',pt.size = 0)
VlnPlot(immune.combined,features = immune.combined$Angio1,pt.size = 0)
VlnPlot(immune.combined,features = 'Angio1',pt.size = 0)
VlnPlot(immune.combined,features = c('PCSK9'),pt.size = 1,sort = T)
FeaturePlot(immune.combined,features = 'PCSK9',pt.size = 1)
VlnPlot(immune.combined,features = c('PCSK9'),pt.size = 1,sort = T,group.by = 'group')

FeaturePlot(immune.combined,features = 'SPP1',split.by = 'group',pt.size = 1)
DefaultAssay()
VlnPlot(immune.combined,features = c("HMOX1","NQO1","NFE2L2","KEAP1","SQSTM1"),group.by = 'group',pt.size = 0)
FeaturePlot(immune.combined,features = c("HMOX1","NQO1","NFE2L2","KEAP1","SQSTM1"),split.by = 'group',pt.size = 1)
FeaturePlot(immune.combined,features = c("SPP1","ACAN","COL2A1","MMP13","ADAMTS5","PTGS2","PDPN"),pt.size = 1)
FeaturePlot(immune.combined,'Angio1',split.by = 'group',pt.size = 1)
FeaturePlot(immune.combined,'Angio1',pt.size = 1)
p1<- VlnPlot(immune.combined,'senmayo1',group.by = "group",pt.size = 0)+geom_boxplot()
p1
mydata=FetchData(immune.combined,vars = c("UMAP_1","UMAP_2","group","seurat_clusters1","Angio1",'senmayo1'))
head(mydata)
class(mydata)
mydata_normal <- subset(mydata,group=="normal")
mydata_oa <- subset(mydata,group=="oa")


# =============================================================================
# 09. Module scores, expression subsets, and score visualization
# =============================================================================

# Configure the plotting palette.
library(viridis)

# Inspect palette colors.
library(scales)
show_col(viridis_pal(option = "B")(30))
# Visualize inflammation scores on the UMAP embedding.
p_normal<-ggplot(data = mydata_normal,aes(x=UMAP_1,y=UMAP_2,
                         colour=senmayo1))+geom_point(size=1)+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))+
  geom_point(size = 1.5)+theme_bw()+theme_bw()+theme(panel.grid.major =element_blank(),
                                                     panel.grid.minor = element_blank(),
                                                     panel.background = element_blank(),  # Remove the panel background.
                                                     panel.border = element_blank(),axis.line = element_line(color = 'black'))

p_normal

p_oa <- ggplot(data = mydata_oa,aes(x=UMAP_1,y=UMAP_2,
                                        colour=senmayo1))+geom_point(size=1)+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('blue','cyan','green','yellow','orange','red'))+
  geom_point(size = 1.5)+theme_bw()+theme(panel.grid.major =element_blank(),
                                        panel.grid.minor = element_blank(),
                                        panel.background = element_blank(),  # Remove the panel background.
                                        panel.border = element_blank(),axis.line = element_line(color = 'black'))
p_oa

p_normal|p_oa

# Explore expression-based subsets using the thresholds specified below.


subset_upspp1<-subset(immune.combined, subset = SPP1 > 1)
mydata_upspp1 <- FetchData(subset_upspp1,vars = c("UMAP_1","UMAP_2","group","seurat_clusters1","Angio1"))
p_upspp1 <-ggplot(data = mydata_upspp1,aes(x=UMAP_1,y=UMAP_2,
                                       colour=Angio1))+geom_point(size=1)+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#FCB216FF','#E25834FF','#992766FF','#450A69FF','#000004FF'))+
  geom_point(size = 1.5)+theme_bw()+theme(panel.grid.major =element_blank(),
                                          panel.grid.minor = element_blank(),
                                          panel.background = element_blank(),  # Remove the panel background.
                                          panel.border = element_blank(),axis.line = element_line(color = 'black'))
p_upspp1

subset_dnspp1<-subset(immune.combined, subset = SPP1 <= 1)
mydata_dnspp1 <- FetchData(subset_dnspp1,vars = c("UMAP_1","UMAP_2","group","seurat_clusters1","Angio1"))
p_dnspp1 <-ggplot(data = mydata_dnspp1,aes(x=UMAP_1,y=UMAP_2,
                                           colour=Angio1))+geom_point(size=1)+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#FCB216FF','#E25834FF','#992766FF','#450A69FF','#000004FF'))+
  geom_point(size = 1.5)+theme_bw()+theme(panel.grid.major =element_blank(),
                                          panel.grid.minor = element_blank(),
                                          panel.background = element_blank(),  # Remove the panel background.
                                          panel.border = element_blank(),axis.line = element_line(color = 'black'))
p_dnspp1

P1<-FeaturePlot(subset_upspp1,features = c('PDPN','COL2A1','OGN','MMP13','PTGS2','FN1'),pt.size = 1)
P2<-FeaturePlot(subset_dnspp1,features = c('PDPN','COL2A1','OGN','MMP13','PTGS2','FN1'),pt.size = 1)
P1|P2

P3 <- FeaturePlot(immune.combined,features = c('PDPN','SPP1'))

# Subset cartilage cells; verify the intended expression threshold before grouping.
subset_cartilage <- subset(immune.combined,idents = c("EC", "SPP1+ C_1","SPP1+ C_2",
                                                      "PreHTC_1", "PreHTC_2/FC",  "ProC", "HTC","RegC/HomoC", "Undecided"))
mydata_cc <- FetchData(subset_cartilage,vars = c("Phase","group","seurat_clusters1","Angio1","senmayo1"))
ccdata_normal <- subset(mydata_cc,group == 'normal')
table(ccdata_normal$seurat_clusters1,ccdata_normal$Phase)
ccdata_oa <- subset(mydata_cc,group == 'oa')
table(ccdata_oa$seurat_clusters1,ccdata_oa$Phase)

P3
senmayo_score <-as.data.frame(immune.combined$senmayo1)
angio_score <- as.data.frame(immune.combined$Angio1)
score_df <- cbind(senmayo_score,angio_score)

subset_cartilage1 <- subset(immune.combined,idents = c("EC/HTC/SPP1+C","PreHTC","FC"))
DefaultAssay(subset_cartilage1) <- "integrated"
subset_cartilage1 <- ScaleData(subset_cartilage1, verbose = FALSE)
subset_cartilage1 <- RunPCA(subset_cartilage1, npcs = 30, verbose = FALSE)
ElbowPlot(subset_cartilage1)
subset_cartilage1 <- RunUMAP(subset_cartilage1, reduction = "pca", dims = 1:10)
subset_cartilage1 <- FindNeighbors(subset_cartilage1, reduction = "pca", dims = 1:10)
subset_cartilage1 <- FindClusters(subset_cartilage1, resolution = 0.40)
DimPlot(subset_cartilage1,reduction = 'umap')
DefaultAssay(subset_cartilage1) <- "RNA"
FeaturePlot(subset_cartilage1,features = c("SPP1","COL10A1","COL1A2","CHRDL2","TGFBI","S100A4"))

# =============================================================================
# 10. Subset trajectories and expression-score plots
# =============================================================================

# Monocle 3 trajectory analysis using the RNA count matrix.

subset_cartilage13<- subset(subset_cartilage1,idents = c("1","2","3","5"))
DefaultAssay(subset_cartilage13) <- "integrated"
subset_cartilage13 <- ScaleData(subset_cartilage13, verbose = FALSE)
subset_cartilage13 <- RunPCA(subset_cartilage13, npcs = 30, verbose = FALSE)
ElbowPlot(subset_cartilage13)
subset_cartilage13 <- RunUMAP(subset_cartilage13, reduction = "pca", dims = 1:10)
subset_cartilage13 <- FindNeighbors(subset_cartilage13, reduction = "pca", dims = 1:10)
subset_cartilage13 <- FindClusters(subset_cartilage13, resolution = 0.40)
DimPlot(subset_cartilage13,reduction = 'umap')
DefaultAssay(subset_cartilage13) <- "RNA"
FeaturePlot(subset_cartilage13,features = c("SPP1","COL10A1","COL1A2","COL10A1","TGFBI","S100A4","ABI3BP"),order = T)
FeaturePlot(subset_cartilage13,features = c("SPP1","COL10A1"),pt.size =  0.5 , order = T)
new.cluster.ids <- c("PreHTC", "PreHTC", "PreHTC", "FC", "HTC/SPP1+C")
names(new.cluster.ids) <- levels(subset_cartilage13)
subset_cartilage13 <- RenameIdents(subset_cartilage13, new.cluster.ids)


immune.combined <- subset_cartilage13
immune.combined1 <- GetAssayData(immune.combined, assay = 'RNA',slot ='counts')
cell_metadata <- immune.combined@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(immune.combined1))
rownames(gene_annotation) <- rownames(immune.combined1)


# Create a Monocle cell_data_set object.
# Create a Monocle cell_data_set object.
cds <- new_cell_data_set(immune.combined1,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds@colData$replicate<- immune.combined$replicate
cds@colData$group<- immune.combined$group
# Run the preprocessing steps required by Monocle.

# Perform dimensionality reduction and clustering.
#NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)  # Preprocess the count matrix and compute principal components.
plot_pc_variance_explained(cds)  # Inspect the variance explained by principal components.
cds <- align_cds(cds, alignment_group = "replicate")

# Compute a UMAP embedding.
cds<- reduce_dimension(cds,preprocess_method = "PCA")  # Use PCA as the preprocessing method.


# Compute a t-SNE embedding.
#cds <- reduce_dimension(cds, reduction_method="tSNE")
#plot_cells(cds, reduction_method="tSNE", color_cells_by="seurat_clusters")

# Cluster cells.


cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(immune.combined, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed  # Reuse the Seurat UMAP coordinates after aligning cell names.
cds@colData$seurat_clusters <- immune.combined@active.ident

# Visualize the embedding.
plot_cells(cds, reduction_method="UMAP", color_cells_by="seurat_clusters")
cds <- cluster_cells(cds)
cds <- learn_graph(cds)  # Learn the principal graph for trajectory inference.

# Color cells using the selected metadata field.
head(colData(cds))
p_plot_cells<-plot_cells(cds,
                         color_cells_by = "seurat_clusters",
                         label_groups_by_cluster=FALSE,
                         label_leaves=FALSE,
                         label_branch_points=TRUE,
                         group_label_size=4,
                         cell_size=0.8)
#p_plot_cells <- p_plot_cells+scale_color_npg()
p_plot_cells
# Inspect the inferred cell-state trajectory.
# Black lines represent the principal graph.
# White circles label graph leaves or terminal states.
# Black circles label branch points.
# Control graph labels with label_leaves and label_branch_points.
# Order cells along the trajectory.
cds = order_cells(cds)  # Select a biologically justified root when ordering cells.
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1)


cds_3d <- reduce_dimension(cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
cds_3d <- order_cells(cds_3d, root_pr_nodes=get_earliest_principal_node(cds))
cds_3d <- order_cells(cds_3d)


position1 <-environment(cds_3d_plot_obj[["x"]][["visdat"]][["5f7455d1474"]])[["data"]]
data <- position1
library(plotly)
Component1 <- data[, 1]
Component2 <- data[, 2]
Component3 <- data[, 3]
p <- plot_ly(data, x = ~Component1, y = ~Component2, z = ~Component3, color = ~data$seurat_clusters,
             colors = c("#FF6DAE","#D4CA3A","#00BDFF"),
             marker = list(size = 2)) %>%
  add_markers()
p

cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="pseudotime")
cds_3d_plot_obj
plot_cells_3d( cds_3d,color_cells_by="seurat_clusters")
plot_cells_3d( cds_3d,color_cells_by="pseudotime")


# Test for genes whose expression varies along the principal graph.
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=6)
write.csv(Track_genes, "Trajectory_genes_all.csv", row.names = F)
Track_genes <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)

# Select the top ten genes for visualization.
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()

# Plot gene expression trends along pseudotime.
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="seurat_clusters",
                         min_expr=0.5, ncol = 2)
p3<-plot_genes_in_pseudotime(cds[c('HMOX1','KEAP1','NFE2L2','NQO1','SQSTM1'),], color_cells_by="seurat_clusters",
                             min_expr=0.5, ncol = 1)

p3


# Add regression line and confidence interval
# Add correlation coefficient: stat_cor()
ggscatter(score_df, x = "immune.combined$senmayo1", y = "immune.combined$Angio1",size = 0.5,
          add = "reg.line", conf.int = TRUE,
          add.params = list(color = "red",  # Set the color of the regression line and confidence interval.
                            fill = "#FF9900")

)+
  stat_cor(method = "pearson",
           label.x = 0.2, label.y = 0.4)

p <- ggplot(data=score_df,mapping = aes(x = "immune.combined$senmayo1",y = "immune.combined$Angio1")) +
  geom_point(shape=19 ,alpha = .01,size=0.5)+
  stat_binhex()+
  scale_fill_gradient(low = "lightblue", high = "red",limits =c(0,6000)) +  # Set the fill gradient and legend limits.
  theme_bw() +
  theme(
    legend.position = "left",  # Set the legend position.
    axis.title = element_text( family = "sans",size = 17,face = "bold",color = "gray25"),  # Set the axis-title font, size, and color.
    axis.text = element_text( family = "sans",size = 15,face = "bold"),  # Set the axis-label font and size.
    title = element_text( family = " sans ",size = 20,face = "bold"),  # Set the plot-title font and size.
    plot.title = element_text(margin = margin(t = 10, b = 10)),  # Set spacing around the plot title.
    axis.title.x = element_text(margin = margin(t = 15),family = "sans"),
    axis.title.y = element_text(margin = margin(r = 15) ,family = "sans"),  # Adjust spacing between the axis title and plotting area.
    legend.title = element_blank(),  # Remove the legend title.
    panel.border = element_rect( linetype = "solid", size = 2,color = "gray50")  # Set the panel-border color, line type, and width.
  ) +
  labs(title="Add marginal to scatter",hjust=0.5) + xlab("Senmayo") + ylab("Angiogenesis")  # Add the plot title and axis labels.
p

p <- ggplot(data=score_df,mapping = aes(x = immune.combined$senmayo1,y = immune.combined$Angio1)) +
  geom_point(shape=19 ,alpha = .01,size=0.5)+
  stat_binhex()
p


####


G2M_score <-as.data.frame(immune.combined$G2M.Score)
S1_score <- as.data.frame(immune.combined$S.Score)
score_df <- cbind(cbind(senmayo_score,angio_score),cbind(G2M_score,S1_score))

# Add regression line and confidence interval
# Add correlation coefficient: stat_cor()
ggscatter(score_df, x = "immune.combined$senmayo1", y = "immune.combined$Angio1",size = 0.5,
          add = "reg.line", conf.int = TRUE,
          add.params = list(color = "red",  # Set the color of the regression line and confidence interval.
                            fill = "#FF9900")

)+
  stat_cor(method = "pearson",
           label.x = 0.6, label.y = 0.4)

ggscatter(score_df, x = "immune.combined$Angio1", y = "immune.combined$S.Score",size = 0.5,
          add = "reg.line", conf.int = TRUE,
          add.params = list(color = "red",  # Set the color of the regression line and confidence interval.
                            fill = "#FF9900")

)+
  stat_cor(method = "pearson",
           label.x = 0.6, label.y = 0.4)
# set cluster to idents
Idents(immune.combined) <- immune.combined@active.ident


####
immune.combined <- irGSEA.score(object = immune.combined, assay = "RNA",
                             slot = "data", seeds = 123, ncores = 1,
                             min.cells = 3, min.feature = 0,
                             custom = F, geneset = NULL, msigdb = T,
                             species = "Homo sapiens", category = "H",
                             subcategory = NULL, geneid = "symbol",
                             method = c("AUCell", "UCell", "singscore",
                                        "ssgsea"),
                             aucell.MaxRank = NULL, ucell.MaxRank = NULL,
                             kcdf = 'Gaussian')
Seurat::Assays(immune.final)
result.dge <- irGSEA.integrate(object = immune.combined,

                               metadata = NULL, col.name = NULL,
                               method = c("AUCell", "UCell", "singscore",
                                          "ssgsea"))
save.image(file = "immune_result.RData")
irGSEA.heatmap.plot <- irGSEA.heatmap(object = result.dge,
                                      method = "RRA",
                                      top = 50,
                                      show.geneset = NULL)

irGSEA.heatmap.plot

irGSEA.bubble.plot <- irGSEA.bubble(object = result.dge,
                                    method = "RRA",
                                    top = 50)
irGSEA.bubble.plot


irGSEA.upset.plot <- irGSEA.upset(object = result.dge,
                                  method = "RRA")

#Warning in if (as.character(ta_call[[1]]) == "upset_top_annotation") {: the
#condition has length > 1 and only the first element will be used
irGSEA.upset.plot

irGSEA.barplot.plot <- irGSEA.barplot(object = result.dge,
                                      method = c("AUCell", "UCell", "singscore",
                                                 "ssgsea"))
irGSEA.barplot.plot

scatterplot <- irGSEA.density.scatterplot(object = immune.combined,
                                          method = "UCell",
                                          show.geneset = "HALLMARK-INFLAMMATORY-RESPONSE",
                                          reduction = "umap")

scatterplot

scatterplot <- irGSEA.density.scatterplot(object = immune.combined,
                                          method = "UCell",
                                          show.geneset = "HALLMARK-APOPTOSIS",
                                          reduction = "umap")

scatterplot

scatterplot <- irGSEA.density.scatterplot(object = immune.combined,
                                          method = "UCell",
                                          show.geneset = "HALLMARK-ANGIOGENESIS",
                                          reduction = "umap",size=2)
# Interactive help: ?irGSEA.density.scatterplot()

scatterplot

result.dge$ssgsea$Name

scatterplot <- irGSEA.density.scatterplot(object = immune.combined,
                                          method = "UCell",
                                          show.geneset = "HALLMARK-MTORC1-SIGNALING",
                                          reduction = "umap",size=2)
scatterplot

densityheatmap <-irGSEA.densityheatmap(object = immune.combined,
                                       method = "UCell",
                                       show.geneset = "HALLMARK-MTORC1-SIGNALING")

densityheatmap

scatterplot <- irGSEA.density.scatterplot(object = immune.combined,
                                          method = "UCell",
                                          show.geneset = "HALLMARK-TGF-BETA-SIGNALING",
                                          reduction = "umap",size=2)

scatterplot
class(irGSEA.barplot.plot)

densityheatmap <-irGSEA.densityheatmap(object = immune.combined,
                                       method = "UCell",
                                       show.geneset = "HALLMARK-TGF-BETA-SIGNALING")
densityheatmap

ridgeplot <- irGSEA.ridgeplot(object = immune.combined,
                              method = "UCell",
                              show.geneset = "HALLMARK-ANGIOGENESIS")
ridgeplot


scatterplot <- irGSEA.density.scatterplot(object = immune.combined,
                                          method = "UCell",
                                          show.geneset = "HALLMARK-ANDROGEN-RESPONSE",
                                          reduction = "umap",size=2)

scatterplot

densityheatmap <-irGSEA.densityheatmap(object = immune.combined,
                                       method = "UCell",
                                       show.geneset = "HALLMARK-ANDROGEN-RESPONSE")
densityheatmap

scatterplot <- irGSEA.density.scatterplot(object = immune.combined,
                                          method = "UCell",
                                          show.geneset = "HALLMARK-CHOLESTEROL-HOMEOSTASIS",
                                          reduction = "umap",size=2)

scatterplot

densityheatmap <-irGSEA.densityheatmap(object = immune.combined,
                                       method = "UCell",
                                       show.geneset = "HALLMARK-CHOLESTEROL-HOMEOSTASIS")
densityheatmap

scatterplot <- irGSEA.density.scatterplot(object = immune.combined,
                                          method = "UCell",
                                          show.geneset = "HALLMARK-OXIDATIVE-PHOSPHORYLATION",
                                          reduction = "umap",size=2)

scatterplot

densityheatmap <-irGSEA.densityheatmap(object = immune.combined,
                                       method = "UCell",
                                       show.geneset = "HALLMARK-OXIDATIVE-PHOSPHORYLATION")
densityheatmap

scatterplot <- irGSEA.density.scatterplot(object = immune.combined,
                                          method = "UCell",
                                          show.geneset = "HALLMARK-REACTIVE-OXYGEN-SPECIES-PATHWAY",
                                          reduction = "umap",size=2)

scatterplot

densityheatmap <-irGSEA.densityheatmap(object = immune.combined,
                                       method = "UCell",
                                       show.geneset = "HALLMARK-REACTIVE-OXYGEN-SPECIES-PATHWAY")
densityheatmap


scatterplot <- irGSEA.density.scatterplot(object = immune.combined,
                                          method = "UCell",
                                          show.geneset = "HALLMARK-PI3K-AKT-MTOR-SIGNALING",
                                          reduction = "umap",size=2)

scatterplot

densityheatmap <-irGSEA.densityheatmap(object = immune.combined,
                                       method = "UCell",
                                       show.geneset = "HALLMARK-PI3K-AKT-MTOR-SIGNALING")
densityheatmap


scatterplot <- irGSEA.density.scatterplot(object = immune.combined,
                                          method = "UCell",
                                          show.geneset = "HALLMARK-TNFA-SIGNALING-VIA-NFKB",
                                          reduction = "umap",size=2)

scatterplot

densityheatmap <-irGSEA.densityheatmap(object = immune.combined,
                                       method = "UCell",
                                       show.geneset = "HALLMARK-TNFA-SIGNALING-VIA-NFKB")
densityheatmap

scatterplot <- irGSEA.density.scatterplot(object = immune.combined,
                                          method = "UCell",
                                          show.geneset = "HALLMARK-APOPTOSIS",
                                          reduction = "umap",size=2)

scatterplot

densityheatmap <-irGSEA.densityheatmap(object = immune.combined,
                                       method = "UCell",
                                       show.geneset = "HALLMARK-APOPTOSIS")
densityheatmap
#> Picking joint bandwidth of 0.00533
densityheatmap <-irGSEA.densityheatmap(object = immune.combined,
                                       method = "UCell",
                                       show.geneset = "HALLMARK-ANGIOGENESIS",
                                      )

densityheatmap

densityheatmap <-irGSEA.densityheatmap(object = immune.combined,
                                       method = "UCell",
                                       show.geneset = "HALLMARK-TGF-BETA-SIGNALING",
)

densityheatmap


####draw enrichment results
DF <- read.table("nallenrichment.txt",header=T,encoding = "utf-8", sep = "\t")
library(forcats)
library(ggplot2)
DF$Term <- as.factor(DF$Term)
DF$Term <- fct_inorder(DF$Term)
DF$Generatio <- as.numeric(DF$Generatio)

P_result2<-ggplot(DF, aes(Group, Term)) +
  geom_point(aes(color=Value, size=Generatio))+theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5))+
  scale_color_gradient(low='blue',high='red')+  #(low='#6699CC',high='#CC3333')+
  labs(x=NULL,y=NULL)+guides(size=guide_legend(order=1))+
theme(legend.direction = "horizontal", legend.position = "top")+  # Place the legend horizontally above the plot.
  scale_y_discrete(position = "right")+scale_color_gradient(low = "blue", high = "red")+ theme(axis.text.y= element_text(angle =-30, hjust = 1))  # Place the y-axis labels on the right.

dev.off()
P_result2


###draw spp1
df<- read.table(file="spp1scrna.txt",sep="\t",header=T,check.names=FALSE)
# Alternative visualization: ridge plot.
p1 <- ggplot(df, aes(x = Value, y = Stage, fill = Stage)) +  # Specify the plotting data and aesthetic mappings.
  geom_density_ridges(stat = "density_ridges") +  # Draw density ridges.
  theme_prism(palette = "candy_bright",  # Select a plot theme.
              base_fontface = "plain",  # Set the font face.
              base_family = "serif",  # Set the font family.
              base_size = 16,  # Set the font size.
              base_line_size = 0.8)+  # Set the axis-line width.
  scale_fill_prism(palette = "candy_bright")+  # Set the fill colors.
  theme(legend.position = "none",  # Hide the legend.
        axis.line.y = element_blank(),  # Hide the y-axis line.
        axis.ticks.y = element_blank())+  # Hide y-axis ticks.
  scale_y_discrete(expand = c(0.05, 0))  # Adjust expansion around the discrete axis.
p1


# =============================================================================
# 11. Normal samples: integration, cell proportions, and exploration
# =============================================================================

# Trajectory analysis of normal samples.

normal1 <-Read10X(data.dir = paths$path_08)
normal1 <- CreateSeuratObject(normal1,min.cells = 3, min.features = 200)
normal2 <-Read10X(data.dir = paths$path_09)
normal2 <- CreateSeuratObject(normal2,min.cells = 3, min.features = 200)
normal3 <-Read10X(data.dir = paths$path_10)
normal3 <- CreateSeuratObject(normal3,min.cells = 3, min.features = 200)
normal1$replicate <- sample(c("normal1"), size = ncol(normal1), replace = TRUE)
normal2$replicate <- sample(c("normal2"), size = ncol(normal2), replace = TRUE)
normal3$replicate <- sample(c("normal3"), size = ncol(normal3), replace = TRUE)


list <-c(normal1,normal2,normal3)
# normalize and identify variable features for each dataset independently
list <- lapply(X = list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = list)
immune.anchors <- FindIntegrationAnchors(object.list = list, anchor.features = features)
# this command creates an 'integrated' data assay
DefaultAssay(immune.combined) <- "RNA"
immune.combined <- IntegrateData(anchorset = immune.anchors)
saveRDS(immune.combined,file ="normalizednormal.rds")
immune.combined[["percent.mt"]] <- PercentageFeatureSet(immune.combined, pattern ="^MT-")
VlnPlot(immune.combined, features = c("nFeatures_RNA","nCount_RNA","percent.mt"),ncol = 3)
immune.combined <- subset(immune.combined, subset = nCount_RNA<100000&percent.mt <5)
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
ElbowPlot(immune.combined)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:15)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:15)
immune.combined <- FindClusters(immune.combined, resolution = 0.5)
DimPlot(immune.combined,reduction = 'umap')

table(immune.combined$orig.ident,immune.combined$replicate)
#sample_normal <- sample(c("normal"),size = 9042,replace = TRUE)
#sample_oa <- sample(c("oa"),size = 57198,replace = TRUE)                           
#immune.combined$group <- c(sample_normal,sample_oa)
saveRDS(immune.combined,file="normalizednormal.rds")

immune.combined1 <-readRDS('normalizednormal.rds')
immune.combined <-readRDS('nallcombined1.rds')

table(scedata$orig.ident)  # Count cells in each cluster.
prop.table(table(Idents(scedata)))
table(Idents(scedata), scedata$orig.ident)  # Tabulate cell counts by cluster and sample.
Cellratio <- prop.table(table(Idents(scedata), scedata$orig.ident), margin = 2)  # Calculate cluster proportions within each group or sample.
Cellratio <- data.frame(Cellratio)
library(reshape2)
cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")  # Reshape cell proportions from long to wide format.
rownames(cellper) <- cellper[,1]
cellper <- cellper[,-1]


# Add sample and group metadata.
sample <- c("BM1","BM2","BM3","GM1","GM2","GM3")
group <- c("BM","BM","BM","GM","GM","GM")
samples <- data.frame(sample, group)  # Create a metadata data frame.

rownames(samples)=samples$sample
cellper$sample <- samples[rownames(cellper),'sample']  # Add a metadata column by matching row names.
cellper$group <- samples[rownames(cellper),'group']  # Add a metadata column by matching row names.

# Visualize the results.
pplist = list()
sce_groups = c('Endothelial','Fibroblast','Immune','Epithelial','Other')
library(ggplot2)
library(dplyr)
install.packages("ggpubr")
library(ggpubr)
for(group_ in sce_groups){
  cellper_  = cellper %>% select(one_of(c('sample','group',group_)))  # Select the sample, group, and current cell-population columns.
  colnames(cellper_) = c('sample','group','percent')  # Rename the selected columns.
  cellper_$percent = as.numeric(cellper_$percent)  # Convert proportions to numeric values.
  cellper_ <- cellper_ %>% group_by(group) %>% mutate(upper =  quantile(percent, 0.75),
                                                      lower = quantile(percent, 0.25),
                                                      mean = mean(percent),
                                                      median = median(percent))  # Summarize the distribution within each group.
  print(group_)
  print(cellper_$median)

  pp1 = ggplot(cellper_,aes(x=group,y=percent)) +  # Plot cell proportions by group.
    geom_jitter(shape = 21,aes(fill=group),width = 0.25) +
    stat_summary(fun=mean, geom="point", color="grey60") +
    theme_cowplot() +
    theme(axis.text = element_text(size = 10),axis.title = element_text(size = 10),legend.text = element_text(size = 10),
          legend.title = element_text(size = 10),plot.title = element_text(size = 10,face = 'plain'),legend.position = 'none') +
    labs(title = group_,y='Percentage') +
    geom_errorbar(aes(ymin = lower, ymax = upper),col = "grey60",width =  1)

# Compare groups using the statistical test specified below.
  labely = max(cellper_$percent)
  compare_means(percent ~ group,  data = cellper_)
  my_comparisons <- list( c("GM", "BM") )
  pp1 = pp1 + stat_compare_means(comparisons = my_comparisons,size = 3,method = "t.test")
  pplist[[group_]] = pp1
}

library(cowplot)
plot_grid(pplist[['Endothelial']],
          pplist[['Fibroblast']],
          pplist[['Immune']],
          pplist[['Epithelial']],
          pplist[['Other']])

#————————————————
# Attribution: adapted material credited to CSDN author TS de Mei Meng (TS's Dream), licensed under CC BY-SA 4.0. Retain attribution and the source link when redistributing adaptations.
# Source: https://blog.csdn.net/qq_42090739/article/details/123392130
##cell cycle
cc.genes.updated.2019
immune.combined <- CellCycleScoring(
  object = immune.combined,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes)
VlnPlot(immune.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB","G2M.Score","S.Score"), ncol = 6)  #+scale_color_npg()
####cell cycle
umapem<-immune.combined@reductions$umap@cell.embeddings
umapem<-as.data.frame(umapem)
metada= immune.combined@meta.data
dim(umapem);dim(metada)

metada$bar<-rownames(metada)
umapem$bar<-rownames(umapem)
ccdata<-merge(umapem,metada,by="bar")
head(ccdata)
library(ggplot2)
plot<-ggplot(ccdata, aes(UMAP_1, UMAP_2,label=Phase))+geom_point(aes(colour = factor(Phase)))+
#plot<-plot+scale_colour_manual(values=c("#CC33FF","Peru","#660000","#660099","#990033","black","red", "#666600", "green","#6699CC","#339900","#0000FF","#FFFF00","#808080"))+
  labs("@yunlai",x = "", y="")
plot=plot+scale_color_aaas()  +
  theme_bw()+theme(panel.grid=element_blank(),legend.title=element_blank(),legend.text = element_text(color="black", size = 10, face = "bold"))
plot<-plot+guides(colour = guide_legend(override.aes = list(size=5))) +theme(plot.title = element_text(hjust = 0.5))
plot
write.table(ccdata,file ="ccdata.txt",sep = "\t", quote = F)
table(immune.combined$group,immune.combined$Phase)
table(immune.combined$integrated_snn_res.0.4,immune.combined$Phase)
table(immune.combined$group,immune.combined$integrated_snn_res.0.4)


DefaultAssay(subset_cartilage) <- 'RNA'
##
DimPlot(subset_cartilage,pt.size = 1)
FeaturePlot(subset_cartilage,features = c("TF","CHRDL2"),pt.size = 1)
FeaturePlot(subset_cartilage,features = c("CKS2","TUBB4B"),pt.size = 1)
FeaturePlot(subset_cartilage,features = c("P3H2","DDX21","ATP13A3"),pt.size = 1,min.cutoff =0)  #0
FeaturePlot(subset_cartilage,features = c("TGFBI","S100A4","TNC","ABI3BP"),pt.size = 1)
FeaturePlot(subset_cartilage,features = c('JUN','FOS'),pt.size = 1,min.cutoff = 1)  #4
##EC
VlnPlot(immune.combined, features = c("TF","CHRDL2","TPM4","RARRES2","DSC2","CYSTM1","KIF9"),pt.size =0)  ##1,5
FeaturePlot(immune.combined,features = c("TF","CHRDL2"),pt.size = 1)
###RegC
VlnPlot(immune.combined, features = c("KIF9","CKS2","TUBB4B","HMOX1","TMEM176A"),pt.size =0)  # 2
FeaturePlot(immune.combined,features = c("KIF9","CKS2","TUBB4B","HMOX1","TMEM176A"),pt.size = 1)
FeaturePlot(immune.combined,features = c("TMEM176A"),pt.size = 2)
###ProC
VlnPlot(immune.combined, features = c("P3H2","DDX21","GLIS3","ATP13A3","NOP16","MTHFD2","TXNRD1","ACTN1","ASNS","ERCC1","FAM101B","BMP6"),pt.size =0)
FeaturePlot(immune.combined,features = c("P3H2","DDX21","ATP13A3"),pt.size = 1,min.cutoff =0)  #0
FeaturePlot(immune.combined,features = c("DDX21"),pt.size = 1,min.cutoff =1)  #0
##PreHTC
VlnPlot(immune.combined, features = c("TGFBI","S100A4","TNC","ABI3BP","CRTAC1","TNFAIP6","CRLF1"),pt.size =0)  #6,3
FeaturePlot(immune.combined,features = c("TGFBI","S100A4","TNC","ABI3BP","ABI3BP","CRTAC1","TNFAIP6","CRLF1"),pt.size = 1)
FeaturePlot(immune.combined,features = c("MMP1","CHI3L2","TNFAIP6"),pt.size = 1)
###FC
VlnPlot(immune.combined, features = c('COL1A1',	'COL1A2',	'TMSB4X',	'MTRNR2L2',	'S100A4',	'MGST3',	'MTRNR2L8',	'SERPINF1',	'ID3',	'IFITM2',	'MMP2',	'THY1',	'LGALS3BP',	'IFI16',	'PSME1',	'IGFBP7',	'IFI27',	'COL5A1',	'MTRNR2L1'),pt.size =0)  #6
FeaturePlot(immune.combined,features = c('COL1A1',	'COL1A2',	'TMSB4X',	'MTRNR2L2','OGN'),pt.size = 1)  #6

##HTC
VlnPlot(immune.combined, features = c('RNU6-2',	'WWP2',	'COL10A1',	'BHLHE41',	'IBSP',	'RMRP',	'PLCG2',	'SCIN',	'S100A1',	'MEF2C',	'CRISPLD1',	'VTRNA1-1'),pt.size =0)  #6
FeaturePlot(immune.combined, features = c('WWP2',	'COL10A1','VTRNA1-1','BHLHE41',	'IBSP',	'RMRP',	'PLCG2'),pt.size =1)  #4,6
##HomC
VlnPlot(immune.combined, features = c('JUN',	'PPP1R10',	'BRD2',	'ID3',	'ID1',	'RGS16',	'HES1',	'NR4A1',	'RSRP1',	'TRA2B',	'RSRC2',	'NR4A2',	'CLK1',	'RASD1',	'CCNL1',	'IER2',	'NFKBIA',	'SNAI1',	'HEXIM1',	'KLF2'),pt.size =0)  #3
FeaturePlot(immune.combined,features = c('JUN','PPP1R10',	'BRD2',	'ID3',	'ID1'),pt.size = 1,min.cutoff = 1)  #4
##HBB HBA1+ C
VlnPlot(immune.combined, features = c('HBB','HBA1'),pt.size =0,y.max =10,sort = TRUE)  #Red Blood Cell 8
FeaturePlot(immune.combined,features = c("FGFBP2"),min.cutoff =0,pt.size = 1)
FeaturePlot(immune.combined,features = c("MKI67","PCNA"),min.cutoff =0,pt.size = 1)
FeaturePlot(immune.combined,features = c("RUNX2"),min.cutoff =0,pt.size = 1)
FeaturePlot(immune.combined,features = c("RUNX2",'OGN','ALP','COL2A1','ACAN','SOX9'),min.cutoff =0,pt.size = 1)

FeaturePlot(immune.combined,features = c("SPP1","COL10A1" ),blend = T)
###
immune.combined <- RenameIdents(immune.combined, `0` = "ProC", `1` = "EC_1", `2` = "RegC",
                                `3` = "PreHTC/HomoC", `4` = "HTC", `5` = "EC_2", `6` = "PreHTC/HTC/FC")
DimPlot(immune.combined1,pt.size = 1,label = T,label.size = 5)
DimPlot(immune.combined1,pt.size = 1,label = T,label.size = 5)

immune.combined <- immune.combined1
DefaultAssay(immune.combined1) <- "RNA"
##EC
VlnPlot(immune.combined, features = c("TF","CHRDL2","TPM4","RARRES2","DSC2","CYSTM1","KIF9"),pt.size =0)  ##1,5
FeaturePlot(immune.combined,features = c("TF","CHRDL2"),pt.size = 1,min.cutoff =0)
FeaturePlot(immune.combined,features = c("MMP13"),pt.size = 1,min.cutoff =0)
FeaturePlot(immune.combined,features = c("TF","CHRDL2"),pt.size = 1,min.cutoff =0)
###RegC
VlnPlot(immune.combined, features = c("KIF9","CKS2","TUBB4B","HMOX1","TMEM176A"),pt.size =0)  # 2
FeaturePlot(immune.combined,features = c("KIF9","CKS2","TUBB4B","HMOX1","TMEM176A"),pt.size = 1)
FeaturePlot(immune.combined,features = c("KIF9","CKS2","TUBB4B","TMEM176A"),pt.size = 1,min.cutoff =0)
###ProC
VlnPlot(immune.combined, features = c("P3H2","DDX21","GLIS3","ATP13A3","NOP16","MTHFD2","TXNRD1","ACTN1","ASNS","ERCC1","FAM101B","BMP6"),pt.size =0)
FeaturePlot(immune.combined,features = c("P3H2"),pt.size = 1,min.cutoff =0)  #0
##HomoC
FeaturePlot(immune.combined,features = c('JUN','FOS'),pt.size = 1,min.cutoff = 1)  #4
##PreHTC
VlnPlot(immune.combined, features = c("TGFBI","S100A4","TNC","ABI3BP","CRTAC1","TNFAIP6","CRLF1"),pt.size =0)  #6,3
FeaturePlot(immune.combined,features = c("TGFBI","S100A4","TNC","ABI3BP","ABI3BP","CRTAC1","TNFAIP6","CRLF1"),pt.size = 1)
FeaturePlot(immune.combined,features = c("TGFBI","S100A4","TNC","ABI3BP"),pt.size = 1)
###FC
VlnPlot(immune.combined, features = c('COL1A1',	'COL1A2',	'TMSB4X',	'MTRNR2L2',	'S100A4',	'MGST3',	'MTRNR2L8',	'SERPINF1',	'ID3',	'IFITM2',	'MMP2',	'THY1',	'LGALS3BP',	'IFI16',	'PSME1',	'IGFBP7',	'IFI27',	'COL5A1',	'MTRNR2L1'),pt.size =0)  #6
FeaturePlot(immune.combined,features = c('COL1A1',	'COL1A2'),pt.size = 1)  #6

##HTC
VlnPlot(immune.combined, features = c('RNU6-2',	'WWP2',	'COL10A1',	'BHLHE41',	'IBSP',	'RMRP',	'PLCG2',	'SCIN',	'S100A1',	'MEF2C',	'CRISPLD1',	'VTRNA1-1'),pt.size =0)  #6
FeaturePlot(immune.combined, features = c(	'COL10A1','IBSP'),pt.size =1)  #4,6
VlnPlot(immune.combined, features = c(	'IBSP','SPP1'),pt.size =0,sort = T)  #4,6


FeaturePlot(immune.combined, features = c(	'MT1G','MT1M'),pt.size =1)  #4,6

DefaultAssay(immune.combined1) <- "integrated"
immune.combined1 <- RunUMAP(immune.combine, reduction = "pca", dims = 1:15)
immune.combined1 <- FindNeighbors(immune.combined1, reduction = "pca", dims = 1:15)
immune.combined1 <- FindClusters(immune.combined1, resolution = 0.40)
immune.combined1 <- FindClusters(immune.combined1, resolution = 0.80)
DimPlot(immune.combined1,reduction = "umap", pt.size = 1)

DefaultAssay(immune.combined1) <- "RNA"
##EC
VlnPlot(immune.combined, features = c("TF","CHRDL2","TPM4","RARRES2","DSC2","CYSTM1","KIF9"),pt.size =0)  ##1,5
FeaturePlot(immune.combined,features = c("TF","CHRDL2"),pt.size = 1,min.cutoff =0)
FeaturePlot(immune.combined1,features = c("MMP13"),pt.size = 1,min.cutoff =0)
FeaturePlot(immune.combined1,features = c("TF","CHRDL2"),pt.size = 1,min.cutoff =0)
FeaturePlot(immune.combined1,features = c("CHRDL2"),pt.size = 1,min.cutoff =0)
###RegC
VlnPlot(immune.combined, features = c("KIF9","CKS2","TUBB4B","HMOX1","TMEM176A"),pt.size =0)  # 2
FeaturePlot(immune.combined,features = c("KIF9","CKS2","TUBB4B","HMOX1","TMEM176A"),pt.size = 1)
FeaturePlot(immune.combined,features = c("KIF9","CKS2","TUBB4B","TMEM176A"),pt.size = 1,min.cutoff =0)
FeaturePlot(immune.combined1,features = c("CKS2"),pt.size = 1,min.cutoff =0)
###ProC
VlnPlot(immune.combined, features = c("P3H2","DDX21","GLIS3","ATP13A3","NOP16","MTHFD2","TXNRD1","ACTN1","ASNS","ERCC1","FAM101B","BMP6"),pt.size =0)
FeaturePlot(immune.combined1,features = c("P3H2"),pt.size = 1,min.cutoff =0)  #0
FeaturePlot(immune.combined1,features = c("ATP13A3"),pt.size = 1,min.cutoff =0)  #0
##HomoC
FeaturePlot(immune.combined1,features = c('JUN','FOS'),pt.size = 1,min.cutoff = 1)  #4
##PreHTC
VlnPlot(immune.combined, features = c("TGFBI","S100A4","TNC","ABI3BP","CRTAC1","TNFAIP6","CRLF1"),pt.size =0)  #6,3
FeaturePlot(immune.combined,features = c("TGFBI","S100A4","TNC","ABI3BP","ABI3BP","CRTAC1","TNFAIP6","CRLF1"),pt.size = 1)
FeaturePlot(immune.combined1,features = c("TGFBI","S100A4","TNC","ABI3BP"),pt.size = 1)
FeaturePlot(immune.combined1,features = c("TGFBI"),pt.size = 1,min.cutoff = 0.5)
FeaturePlot(immune.combined1,features = c("ABI3BP"),pt.size = 1,min.cutoff = 0.5)
###FC
VlnPlot(immune.combined, features = c('COL1A1',	'COL1A2',	'TMSB4X',	'MTRNR2L2',	'S100A4',	'MGST3',	'MTRNR2L8',	'SERPINF1',	'ID3',	'IFITM2',	'MMP2',	'THY1',	'LGALS3BP',	'IFI16',	'PSME1',	'IGFBP7',	'IFI27',	'COL5A1',	'MTRNR2L1'),pt.size =0)  #6
FeaturePlot(immune.combined1,features = c('COL1A1',	'COL1A2'),pt.size = 1)  #6

##HTC
VlnPlot(immune.combined, features = c('RNU6-2',	'WWP2',	'COL10A1',	'BHLHE41',	'IBSP',	'RMRP',	'PLCG2',	'SCIN',	'S100A1',	'MEF2C',	'CRISPLD1',	'VTRNA1-1'),pt.size =0)  #6
FeaturePlot(immune.combined1, features = c(	'COL10A1','IBSP'),pt.size =1)  #4,6

##SPP1
VlnPlot(immune.combined, features = c('RNU6-2',	'WWP2',	'COL10A1',	'BHLHE41',	'IBSP',	'RMRP',	'PLCG2',	'SCIN',	'S100A1',	'MEF2C',	'CRISPLD1',	'VTRNA1-1'),pt.size =0)  #6
FeaturePlot(immune.combined1, features = c(	'HBB'),pt.size =1)  #4,6
markers <- FindAllMarkers(immune.combined1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20

DefaultAssay(immune.combined1) <- "RNA"
DoHeatmap(immune.combined1, features = top10$gene)
immune.combined3 <- immune.combined1
immune.combined1 <- subset(x=immune.combined1,idents=c("0","1","2","3","4","5","6","7","8","9","11")  )

FeaturePlot(immune.combined1, features = c(	"TF","CHRDL2","CKS2","TUBB4B",'JUN','FOS',"DDX21","ATP13A3","TGFBI","S100A4","ABI3BP",'COL10A1','IBSP','COL1A1',	'COL1A2',"S100A4"),pt.size =1)  #4,6
FeaturePlot(immune.combined1, features = c(	"CHRDL2"),pt.size =1)  #4,6
FeaturePlot(immune.combined1, features = c(	"CKS2"),pt.size =1)  #4,6
FeaturePlot(immune.combined1, features = c(	"TUBB4B"),pt.size =1,min.cutoff = 2)  #4,6
FeaturePlot(immune.combined1, features = c(	"CKS2","TUBB4B"),pt.size =1)  #4,6
FeaturePlot(immune.combined1, features = c(	'JUN','FOS'),pt.size =1,min.cutoff = 2.5)  #4,6
FeaturePlot(immune.combined1, features = c(	"DDX21"),pt.size =1,min.cutoff = 1)  #4,6
FeaturePlot(immune.combined1, features = c(	"ATP13A3"),pt.size =1,min.cutoff = 0)  #4,6
FeaturePlot(immune.combined1, features = c(	"S100A4","ABI3BP"),pt.size =1,min.cutoff = 0.5)  #4,6
FeaturePlot(immune.combined1, features = c(	"TGFBI"),pt.size =1,min.cutoff = 0)  #4,6
FeaturePlot(immune.combined1, features = c(	'COL10A1','IBSP'),pt.size =1,min.cutoff = 0)  #4,6
FeaturePlot(immune.combined1, features = c(	'COL1A2'),pt.size =1,min.cutoff = 0)  #4,6
FeaturePlot(immune.combined, features = c(	'SPP1','MDM2','OMD'),pt.size =1,min.cutoff = 0)  #4,6
VlnPlot(immune.combined, features = c(	'SPP1','MDM2','OMD'),pt.size =0,sort =T)  #4,6
VlnPlot(immune.combined, features = c(	'DKK3',	'HIST1H1E',	'C1orf56',	'NRP2',	'CTSB',	'COL5A2',	'CITED4',	'CFH',	'SMOC2',	'OAF',
                                       'RFLNB',	'ST6GALNAC4',	'DNER',	'UBE2S',	'TFPI2',	'OMD',	'MDM2',	'IBSP',	'LRRC75A'),pt.size =0,sort = T)  #4,6
FeaturePlot(immune.combined1, features = c(	'S100A4'),pt.size =1,min.cutoff = 1)  #4,6
DimPlot(immune.combined1,reduction = "umap", pt.size = 1,label = F)

new.cluster.ids <- c("ProC","ProC", "RegC", "ProC","EC/HTC/SPP1+C", "EC", "EC",
                     "PreHTC", "ProC", "HomoC/RegC","FC")
names(new.cluster.ids) <- levels(immune.combined1)
immune.combined1 <- RenameIdents(immune.combined1, new.cluster.ids)

DimPlot(immune.combined1,reduction = "umap")
saveRDS(immune.combined1,file ="nnormalclusters.rds")
immune.combined <- readRDS(file ="nnormalclusters.rds")

as.data.frame(prop.table(table(Idents(immune.combined), immune.combined@meta.data[,"group"]), margin = 2))-> pdf -> td
library(tidyverse)
allcolour=c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00","#808000","#FF00FF","#FA8072","#7B68EE","#9400D3","#800080","#A0522D","#D2B48C","#D2691E","#87CEEB","#40E0D0","#5F9EA0","#FF1493","#0000CD","#008B8B","#FFE4B5","#8A2BE2","#228B22","#E9967A","#4682B4","#32CD32","#F0E68C","#FFFFE0","#EE82EE","#FF6347","#6A5ACD","#9932CC","#8B008B","#8B4513","#DEB887")
allcolour -> colour1
plt<- ggplot(td,aes(x=td[,2],y=td[,3],fill=td[,1]))+
  geom_bar(position = 'stack',stat="identity")+
  labs(x="replicate",y="Cells Ratio")+
  theme(panel.background=element_rect(fill='transparent', color='black'),
        legend.key=element_rect(fill='transparent', color='transparent'),axis.text = element_text(color="black"))+
  scale_y_continuous(expand=c(0.001,0.001))+
  scale_fill_manual(values=colour1)+
  guides(fill = guide_legend(keywidth = 1, keyheight = 1,ncol=1,title = 'Cluster'))

plt
#############
DefaultAssay(immune.combined) <- "integrated"
p3 <- DimPlot(immune.combined, reduction = "umap", group.by = "replicate",pt.size= 1)
p4 <- DimPlot(immune.combined, reduction = "umap", label = TRUE, repel = TRUE,pt.size= 1,cols = c("0" ="#DC143C","1"= "#0000FF","2"= "#20B2AA","3"= "#FFA500","4"= "#9370DB","5"= "#98FB98","6"= "#F08080","7"= "#1E90FF","8"= "#7CFC00","9"= "#FFFF00"))
p4 <- DimPlot(immune.combined, reduction = "umap", split.by = "replicate",pt.size= 1,cols = c("0" ="#DC143C","1"= "#0000FF","2"= "#20B2AA","3"= "#FFA500","4"= "#9370DB","5"= "#98FB98","6"= "#F08080","7"= "#1E90FF","8"= "#7CFC00","9"= "#FFFF00"))
p3|p4
p3
p4
table(immune.combined$integrated_snn_res.0.4,immune.combined$replicate)
####
# Interactive help: ?FindAllMarkers()
markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25,test.use ="wilcox")
markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.table(markers,file = "allmarkersres0.4.txt",sep ='\t',quote = F)
markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

DefaultAssay(immune.combined) <- "integrated"
DoHeatmap(immune.combined, features = top10$gene)
DoHeatmap(immune.combined, features = top10$gene,group.colors = c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00"))

##conserved cell type
DefaultAssay(immune.combined) <- "RNA"
conserved_markers <- FindConservedMarkers(immune.combined,ident.1= 0,grouping.var = "group")
conserved_markers0 <-conserved_markers

###

# =============================================================================
# 12. Normal-sample and selected-subpopulation trajectories
# =============================================================================

# Monocle 3 trajectory analysis using the RNA count matrix.
immune.combined2 <- immune.combined1
immune.combined <- immune.combined1
immune.combined1 <- GetAssayData(immune.combined, assay = 'RNA',slot ='counts')
cell_metadata <- immune.combined@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(immune.combined1))
rownames(gene_annotation) <- rownames(immune.combined1)

# Create a Monocle cell_data_set object.
# Create a Monocle cell_data_set object.
cds <- new_cell_data_set(immune.combined1,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds@colData$replicate<- immune.combined$replicate
cds@colData$group<- immune.combined$group
# Run the preprocessing steps required by Monocle.

# Perform dimensionality reduction and clustering.
#NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)  # Preprocess the count matrix and compute principal components.
plot_pc_variance_explained(cds)  # Inspect the variance explained by principal components.
cds <- align_cds(cds, alignment_group = "replicate")

# Compute a UMAP embedding.
cds<- reduce_dimension(cds,preprocess_method = "PCA")  # Use PCA as the preprocessing method.


# Compute a t-SNE embedding.
#cds <- reduce_dimension(cds, reduction_method="tSNE")
#plot_cells(cds, reduction_method="tSNE", color_cells_by="seurat_clusters")

# Cluster cells.


cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(immune.combined, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed  # Reuse the Seurat UMAP coordinates after aligning cell names.
cds@colData$seurat_clusters <- immune.combined@active.ident

# Visualize the embedding.
plot_cells(cds, reduction_method="UMAP", color_cells_by="seurat_clusters")
cds <- cluster_cells(cds)
cds <- learn_graph(cds)  # Learn the principal graph for trajectory inference.

# Color cells using the selected metadata field.
head(colData(cds))
p_plot_cells<-plot_cells(cds,
                         color_cells_by = "seurat_clusters",
                         label_groups_by_cluster=FALSE,
                         label_leaves=FALSE,
                         label_branch_points=TRUE,
                         group_label_size=4,
                         cell_size=0.8)
#p_plot_cells <- p_plot_cells+scale_color_npg()
p_plot_cells
# Inspect the inferred cell-state trajectory.
# Black lines represent the principal graph.
# White circles label graph leaves or terminal states.
# Black circles label branch points.
# Control graph labels with label_leaves and label_branch_points.
# Order cells along the trajectory.
cds = order_cells(cds)  # Select a biologically justified root when ordering cells.
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1)

# Test for genes whose expression varies along the principal graph.
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=6)
write.csv(Track_genes, "Trajectory_genes_all.csv", row.names = F)
Track_genes <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)

# Select the top ten genes for visualization.
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()

# Plot gene expression trends along pseudotime.
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="seurat_clusters",
                         min_expr=0.5, ncol = 2)
p3<-plot_genes_in_pseudotime(cds[c('HMOX1','KEAP1','NFE2L2','NQO1','SQSTM1'),], color_cells_by="seurat_clusters",
                             min_expr=0.5, ncol = 1)

immune.combined$pseudotime <- cds@principal_graph_aux@listData[["UMAP"]][["pseudotime"]]
immune.combined$seurat_clusters <- immune.combined@active.ident
gene_angio<-read.table("genelist_angio.txt",header = T)
gene<- as.list(gene_angio)
immune.combined <- AddModuleScore(object = immune.combined,features = gene,name = "Angio")

gene_senmayo<-read.table("genelist_senmayo.txt",header = T)
gene<- as.list(gene_senmayo)
immune.combined <- AddModuleScore(object = immune.combined,features = gene,name = "senmayo")
df1=FetchData(immune.combined, vars=c("pseudotime", "seurat_clusters","Angio1","senmayo1"))


p1 <- ggplot(data = df1)+
#geom_point(aes(x=pseudotime,y=Angio1,color=N,group=1))+
#scale_color_gradientn(colours=c("blue","green","yellow","red"))
  geom_smooth (aes(x=pseudotime,y=Angio1),method = lm,formula = y ~ x+I(x^2),linetype=1,se=T,span=1)
p1
p2 <- ggplot(data = df1)+
#geom_point(aes(x=pseudotime,y=Angio1,color=N,group=1))+
#scale_color_gradientn(colours=c("blue","green","yellow","red"))
  geom_smooth (aes(x=pseudotime,y=senmayo1),method = lm,formula = y ~ x+I(x^2),linetype=1,se=T,span=1)
p2

p1|p2
VlnPlot(immune.combined$pseudotime)
p3
# Interactive help: ?plot_genes_in_pseudotime()
p3 <-p3+scale_color_npg()
p3
p3 <- p3+scale_color_manual(values = c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00"))
p3
# Visualize selected genes on the cell embedding.
plot_cells(cds, genes=c("SQSTM1","COL2A1","NQO1","ACAN","MMP3","ADAMTS5","MMP13","COL10A1"), show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)
# Interpret trajectory-associated genes alongside genes of biological interest.
# Identify modules of coexpressed genes.
genelist <- pull(Track_genes, gene_short_name) %>% as.character()
gene_module <- find_gene_modules(cds[genelist,], resolution=1e-2, cores = 10)
cell_group <- tibble::tibble(cell=row.names(colData(cds)),
                             cell_group=colData(cds)$integrated_snn_res.0.4)
agg_mat <- aggregate_gene_expression(cds, gene_module, cell_group)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2",color = colorRampPalette(colors = c("blue","white","red"))(100))
write.csv(gene_module,"all_pseuodolgenemodule.csv",row.names = F )
###show module genes
plot_cells(cds_subset2,
           genes=gene_module %>% filter(module %in% c(102, 46)),
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

# Subset selected cell populations for trajectory analysis.

immune.combined <-readRDS('oa10clusters.rds')
DimPlot(immune.combined,reduction = "umap",label = T)
FeaturePlot(immune.combined, features = c(	"CHRDL2","TF","MMP13"),pt.size =1)  #4,6
FeaturePlot(immune.combined, features = c(	"CKS2"),pt.size =1)  #4,6
FeaturePlot(immune.combined, features = c(	"TUBB4B"),pt.size =1,min.cutoff = 2)  #4,6
FeaturePlot(immune.combined, features = c(	"CKS2","TUBB4B","MMP13"),pt.size =1)  #4,6
FeaturePlot(immune.combined, features = c(	'JUN','FOS'),pt.size =1,min.cutoff = 2.5)  #4,6
FeaturePlot(immune.combined, features = c(	"DDX21","P3H2"),pt.size =1,min.cutoff = 1)  #4,6
FeaturePlot(immune.combined, features = c(	"ATP13A3"),pt.size =1,min.cutoff = 0)  #4,6
FeaturePlot(immune.combined, features = c("TGFBI",	"S100A4","ABI3BP"),pt.size =1,min.cutoff = 0.5)  #4,6
FeaturePlot(immune.combined, features = c(	"TGFBI","ABI3BP"),pt.size =1,min.cutoff = 0)  #4,6
FeaturePlot(immune.combined, features = c(	'COL10A1','IBSP'),pt.size =1,min.cutoff = 0)  #4,6
FeaturePlot(immune.combined, features = c(	'COL1A2','COL1A1'),pt.size =1,min.cutoff = 0)  #4,6
FeaturePlot(immune.combined, features = c(	'SPP1'),pt.size =1,min.cutoff = 0)  #4,6
FeaturePlot(immune.combined, features = c(	'S100A4'),pt.size =1,min.cutoff = 1)  #4,6

new.cluster.ids <- c("EC", "SPP1+C", "PreHTC", "PreHTC", "RegC", "ProC",
                     "HTC", "FC", "HomoC","undecidedC","undecidedC")
names(new.cluster.ids) <- levels(immune.combined)
immune.combined <- RenameIdents(immune.combined, new.cluster.ids)

DimPlot(immune.combined,reduction = "umap",pt.size =1)

immune.combined$seurat_clusters <- immune.combined$integrated_snn_res.0.4
immune.combined_pro <- immune.combined[,immune.combined@meta.data$seurat_clusters %in% c(0,1,2,3,5)]
sce=immune.combined_pro
#sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 1e4) 
#sce <- FindVariableFeatures(sce, selection.method = 'vst', nfeatures = 2000)
#DefaultAssay(immune.combined) <- "integrated"
#immune.combined <- ScaleData(immune.combined, verbose = FALSE)
#immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
#ElbowPlot(immune.combined)
#immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:15)
#immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:15)
#immune.combined <- FindClusters(immune.combined, resolution = 0.40)
DefaultAssay(sce) <- "integrated"
sce <- ScaleData(sce, verbose = FALSE)
sce <- RunPCA(sce, npcs =30, verbose = FALSE)
ElbowPlot(sce)
sce <- FindNeighbors(sce, dims = 1:15)
sce <- FindClusters(sce, resolution = 0.4 )
# Look at cluster IDs of the first 5 cells
DefaultAssay(sce) <- "RNA"
head(Idents(sce), 5)
table(sce$seurat_clusters)
sce <- RunUMAP(sce, dims = 1:15)
DimPlot(sce, reduction = 'umap')
FeaturePlot(sce,features=c("SOD2","NQO1"))  # 0 ProC
FeaturePlot(sce,features=c("MMP3","KCNMA1"))  #  4 EC_1
FeaturePlot(sce,features=c("HSPA6","HSPA1A"))  #2 RegC 
FeaturePlot(sce,features=c("CYTOR"))  #3 EC_2 
FeaturePlot(sce,features = c("TF","CHRDL2","TPM4","RARRES2","DSC2","CYSTM1"),pt.size = 1)  #3 EC_2 
FeaturePlot(sce,features = c("FOS","MYOC","APOD"),pt.size = 1)  #1 "PreHTC/HomoC"
FeaturePlot(sce,features = c('JUN','PPP1R10',	'BRD2',	'ID3',	'ID1'),pt.size = 1,min.cutoff = 1)  #homo
FeaturePlot(sce,features = c("MMP1","CHI3L2","TNFAIP6"),pt.size = 1)  #prehtc
sce <- RenameIdents(sce, `0` = "ProC", `1` = "PreHTC/HomoC", `2` = "RegC",
                                `3` = "EC_2", `4` = "EC_1")
FeaturePlot(sce,features=c("HMOX1","NQO1","NFE2L2","KEAP1"))  # 0 ProC
###


immune.combined1 <- subset(x=immune.combined,idents=c("PreHTC","HTC","SPP1+C","FC")  )
DefaultAssay(immune.combined1) <- "integrated"
immune.combined1 <- ScaleData(immune.combined1, verbose = FALSE)
immune.combined1 <- RunPCA(immune.combined1, npcs = 30, verbose = FALSE)
ElbowPlot(immune.combined1)
immune.combined1 <- RunUMAP(immune.combined1, reduction = "pca", dims = 1:15)
immune.combined1 <- FindNeighbors(immune.combined1, reduction = "pca", dims = 1:15)
immune.combined1 <- FindClusters(immune.combined1, resolution = 0.20)
DimPlot(immune.combined1,reduction = 'umap')
DefaultAssay(immune.combined1) <- "RNA"
FeaturePlot(immune.combined1,features = c("SPP1","COL10A1","COL1A2","TGFBI","S100A4","ABI3BP","IBSP"))
new.cluster.ids <- c("SPP1+C", "PreHTC", "SPP1+C", "FC", "HTC")
names(new.cluster.ids) <- levels(immune.combined1)
immune.combined1 <- RenameIdents(immune.combined1, new.cluster.ids)

# Monocle 3 trajectory analysis using the RNA count matrix.

immune.combined2 <- immune.combined
immune.combined <- immune.combined1
immune.combined1 <- GetAssayData(immune.combined, assay = 'RNA',slot ='counts')
cell_metadata <- immune.combined@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(immune.combined1))
rownames(gene_annotation) <- rownames(immune.combined1)


# Create a Monocle cell_data_set object.
# Create a Monocle cell_data_set object.
cds <- new_cell_data_set(immune.combined1,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds@colData$replicate<- immune.combined$replicate
cds@colData$group<- immune.combined$group


#cds@colData$group<- immune.combined$group
# Run the preprocessing steps required by Monocle.

# Perform dimensionality reduction and clustering.
#NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)  # Preprocess the count matrix and compute principal components.
plot_pc_variance_explained(cds)  # Inspect the variance explained by principal components.
cds <- align_cds(cds, alignment_group = "replicate")

# Compute a UMAP embedding.
cds<- reduce_dimension(cds,preprocess_method = "PCA")  # Use PCA as the preprocessing method.

# Cluster cells.


cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(immune.combined, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed  # Reuse the Seurat UMAP coordinates after aligning cell names.
cds@colData$seurat_clusters <- immune.combined@active.ident

# Visualize the embedding.
plot_cells(cds, reduction_method="UMAP", color_cells_by="seurat_clusters")
cds <- cluster_cells(cds)
cds <- learn_graph(cds)  # Learn the principal graph for trajectory inference.

# Color cells using the selected metadata field.
head(colData(cds))
p_plot_cells<-plot_cells(cds,
                         color_cells_by = "seurat_clusters",
                         label_groups_by_cluster=FALSE,
                         label_leaves=FALSE,
                         label_branch_points=TRUE,
                         group_label_size=4,
                         cell_size=0.8)
#p_plot_cells <- p_plot_cells+scale_color_npg()
p_plot_cells
# Inspect the inferred cell-state trajectory.
# Black lines represent the principal graph.
# White circles label graph leaves or terminal states.
# Black circles label branch points.
# Control graph labels with label_leaves and label_branch_points.
# Order cells along the trajectory.
cds = order_cells(cds)  # Select a biologically justified root when ordering cells.
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1)


# Compute a t-SNE embedding.
#cds <- reduce_dimension(cds, reduction_method="tSNE")
#plot_cells(cds, reduction_method="tSNE", color_cells_by="seurat_clusters")

# Cluster cells.


cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(sce, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed  # Reuse the Seurat UMAP coordinates after aligning cell names.
cds@colData$integrated_snn_res.0.4 <- sce$integrated_snn_res.0.4

# Visualize the embedding.
plot_cells(cds, reduction_method="UMAP", color_cells_by="integrated_snn_res.0.4")
cds <- cluster_cells(cds)
cds <- learn_graph(cds)  # Learn the principal graph for trajectory inference.

# Color cells using the selected metadata field.
head(colData(cds))
p_plot_cells<-plot_cells(cds,
                         color_cells_by = "integrated_snn_res.0.4",
                         label_groups_by_cluster=FALSE,
                         label_leaves=FALSE,
                         label_branch_points=TRUE,
                         group_label_size=4,
                         cell_size=0.8)
#p_plot_cells <- p_plot_cells+scale_color_npg()
p_plot_cells
# Inspect the inferred cell-state trajectory.
# Black lines represent the principal graph.
# White circles label graph leaves or terminal states.
# Black circles label branch points.
# Control graph labels with label_leaves and label_branch_points.
# Order cells along the trajectory.
cds = order_cells(cds)  # Select a biologically justified root when ordering cells.
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5,
           group_label_size=4,cell_size=1)

# Test for genes whose expression varies along the principal graph.
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=6)
write.csv(Track_genes, "normal_5_clustes_Trajectory_genes_all.csv", row.names = F)
Track_genes <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)

# Select the top ten genes for visualization.
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()

# Plot gene expression trends along pseudotime.
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="integrated_snn_res.0.4",
                         min_expr=0.5, ncol = 2)
p3<-plot_genes_in_pseudotime(cds[c("HMOX1","NQO1","NFE2L2","KEAP1"),], color_cells_by="integrated_snn_res.0.4",
                             min_expr=0, ncol = 1)
# Interactive help: ?plot_genes_in_pseudotime()
p3 <-p3+scale_color_npg()
p3
p3 <- p3+scale_color_manual(values = c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00"))
p3
# Visualize selected genes on the cell embedding.
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="integrated_snn_res.0.4",
                         min_expr=0.5, ncol = 2)
plot_cells(cds, genes=c("HMOX1","NQO1","NFE2L2","KEAP1"), show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)
# Interpret trajectory-associated genes alongside genes of biological interest.
# Identify modules of coexpressed genes.
genelist <- pull(Track_genes, gene_short_name) %>% as.character()
gene_module <- find_gene_modules(cds[genelist,], resolution=1e-2, cores = 10)
cell_group <- tibble::tibble(cell=row.names(colData(cds)),
                             cell_group=colData(cds)$integrated_snn_res.0.4)
agg_mat <- aggregate_gene_expression(cds, gene_module, cell_group)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2",color = colorRampPalette(colors = c("blue","white","red"))(100))
write.csv(gene_module,"5normalclusters_all_pseuodolgenemodule.csv",row.names = F )
###
pseudotime <-pseudotime(cds,reduction_method = 'UMAP')
pseudotime <-pseudotime[rownames(sce@meta.data)]
sce$pseudotime <- pseudotime
FeaturePlot(sce,reduction = "umap",features = "pseudotime")

save(cds,Track_genes, file = "normal_5clusterd_pseudo.RData")
cds1 <-cds
load("normal_5clusterd_pseudo.RData")
saveRDS(sce,file ="normal5clusters.rds")
rm(list = ls())
##########normal 3 clusters


# Subset selected cell populations for trajectory analysis.
immune.combined$seurat_clusters <- immune.combined$integrated_snn_res.0.4
immune.combined_pro <- immune.combined[,immune.combined@meta.data$seurat_clusters %in% c(4,6)]
sce=immune.combined_pro
#sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 1e4) 
#sce <- FindVariableFeatures(sce, selection.method = 'vst', nfeatures = 2000)
#DefaultAssay(immune.combined) <- "integrated"
#immune.combined <- ScaleData(immune.combined, verbose = FALSE)
#immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
#ElbowPlot(immune.combined)
#immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:15)
#immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:15)
#immune.combined <- FindClusters(immune.combined, resolution = 0.40)
DefaultAssay(sce) <- "integrated"
sce <- ScaleData(sce, verbose = FALSE)
sce <- RunPCA(sce, npcs =30, verbose = FALSE)
ElbowPlot(sce)
sce <- FindNeighbors(sce, dims = 1:15)
sce <- FindClusters(sce, resolution = 0.15 )
# Look at cluster IDs of the first 5 cells
DefaultAssay(sce) <- "RNA"
head(Idents(sce), 5)
table(sce$seurat_clusters)
sce <- RunUMAP(sce, dims = 1:15)
DimPlot(sce, reduction = 'umap')
FeaturePlot(sce,features = c('COL1A1',	'COL1A2',	'TMSB4X',	'MTRNR2L2','OGN'),pt.size = 1)  #1 FC
FeaturePlot(sce,features = c("TGFBI","S100A4","TNC","ABI3BP","ABI3BP","CRTAC1","TNFAIP6","CRLF1"),pt.size = 1)  #0 prehtc/homoc
FeaturePlot(sce,features = c("MMP1","CHI3L2","TNFAIP6"),pt.size = 1)  #0 prehtc
FeaturePlot(sce, features = c('WWP2',	'COL10A1','VTRNA1-1','BHLHE41',	'IBSP',	'RMRP',	'PLCG2'),pt.size =1)  #3htc
FeaturePlot(sce,features = c('JUN','PPP1R10',	'BRD2',	'ID3',	'ID1'),pt.size = 1,min.cutoff = 1)  #3homoc
FeaturePlot(sce,features=c("SOD2","NQO1"))  # 0 ProC
FeaturePlot(sce,features=c("MMP3","KCNMA1"))  #  4 EC_1
FeaturePlot(sce,features=c("HSPA6","HSPA1A"))  #2 RegC 
FeaturePlot(sce,features=c("CYTOR"))  #3 EC_2 
FeaturePlot(sce,features = c("TF","CHRDL2","TPM4","RARRES2","DSC2","CYSTM1"),pt.size = 1)  #3 EC_2 
FeaturePlot(sce,features = c("FOS","MYOC","APOD"),pt.size = 1)  #1 "PreHTC/HomoC"
FeaturePlot(sce,features = c('JUN','PPP1R10',	'BRD2',	'ID3',	'ID1'),pt.size = 1,min.cutoff = 1)  #homo
FeaturePlot(sce,features = c("MMP1","CHI3L2","TNFAIP6"),pt.size = 1)  #prehtc
sce <- RenameIdents(sce, `0` = "PreHTC/HomoC", `1` = "FC", `2` = "Undecided",
                    `3` = "HTC")
markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
DefaultAssay(sce) <- "integrated"

markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(sce, features = top10$gene)


FeaturePlot(sce,features=c("HMOX1","NQO1","NFE2L2","KEAP1"))  # 0 ProC

sce <- subset(x=sce,idents=c("PreHTC/HomoC","FC","HTC")  )
DimPlot(sce, reduction = 'umap')
###
# Monocle 3 trajectory analysis using the RNA count matrix.
immune.combined1 <- GetAssayData(sce, assay = 'RNA',slot ='counts')
cell_metadata <- sce@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(immune.combined1))
rownames(gene_annotation) <- rownames(immune.combined1)

# Create a Monocle cell_data_set object.
# Create a Monocle cell_data_set object.
cds <- new_cell_data_set(immune.combined1,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds@colData$replicate<- sce$replicate
#cds@colData$group<- immune.combined$group
# Run the preprocessing steps required by Monocle.

# Perform dimensionality reduction and clustering.
#NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)  # Preprocess the count matrix and compute principal components.
plot_pc_variance_explained(cds)  # Inspect the variance explained by principal components.
cds <- align_cds(cds, alignment_group = "replicate")

# Compute a UMAP embedding.
cds<- reduce_dimension(cds,preprocess_method = "PCA")  # Use PCA as the preprocessing method.


# Compute a t-SNE embedding.
#cds <- reduce_dimension(cds, reduction_method="tSNE")
#plot_cells(cds, reduction_method="tSNE", color_cells_by="seurat_clusters")

# Cluster cells.


cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(sce, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed  # Reuse the Seurat UMAP coordinates after aligning cell names.
cds@colData$integrated_snn_res.0.15 <- sce$integrated_snn_res.0.15

# Visualize the embedding.
plot_cells(cds, reduction_method="UMAP", color_cells_by="integrated_snn_res.0.15",cell_size = 1)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)  # Learn the principal graph for trajectory inference.

# Color cells using the selected metadata field.
head(colData(cds))
p_plot_cells<-plot_cells(cds,
                         color_cells_by = "integrated_snn_res.0.15",
                         label_groups_by_cluster=FALSE,
                         label_leaves=FALSE,
                         label_branch_points=TRUE,
                         group_label_size=4,
                         cell_size=0.8)
#p_plot_cells <- p_plot_cells+scale_color_npg()
p_plot_cells
# Inspect the inferred cell-state trajectory.
# Black lines represent the principal graph.
# White circles label graph leaves or terminal states.
# Black circles label branch points.
# Control graph labels with label_leaves and label_branch_points.
# Order cells along the trajectory.
cds = order_cells(cds)  # Select a biologically justified root when ordering cells.
plot_cells(cds,
           color_cells_by = "pseudotime",
           show_trajectory_graph = FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=0,
           group_label_size=0,cell_size=1)
# Interactive help: ?plot_cells()

# Test for genes whose expression varies along the principal graph.
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=6)
write.csv(Track_genes, "normal_5_clustes_Trajectory_genes_all.csv", row.names = F)
Track_genes <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)

# Select the top ten genes for visualization.
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()

# Plot gene expression trends along pseudotime.
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="integrated_snn_res.0.15",
                         min_expr=0.5, ncol = 2)
p3<-plot_genes_in_pseudotime(cds[c("SPP1"),], color_cells_by="integrated_snn_res.0.15",
                             min_expr=0, ncol = 1)
# Interactive help: ?plot_genes_in_pseudotime()
p3 <-p3+scale_color_npg()
p3
p3 <- p3+scale_color_manual(values = c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00"))
p3
# Visualize selected genes on the cell embedding.
plot_cells(cds, genes=c("HMOX1","NQO1","NFE2L2","KEAP1"), show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)
# Interpret trajectory-associated genes alongside genes of biological interest.
# Identify modules of coexpressed genes.
genelist <- pull(Track_genes, gene_short_name) %>% as.character()
gene_module <- find_gene_modules(cds[genelist,], resolution=1e-2, cores = 10)
cell_group <- tibble::tibble(cell=row.names(colData(cds)),
                             cell_group=colData(cds)$integrated_snn_res.0.4)
agg_mat <- aggregate_gene_expression(cds, gene_module, cell_group)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2",color = colorRampPalette(colors = c("blue","white","red"))(100))
write.csv(gene_module,"5normalclusters_all_pseuodolgenemodule.csv",row.names = F )
###
pseudotime <-pseudotime(cds,reduction_method = 'UMAP')
pseudotime <-pseudotime[rownames(sce@meta.data)]
sce$pseudotime <- pseudotime
FeaturePlot(sce,reduction = "umap",features = "pseudotime")

save(cds,Track_genes, file = "normal_3clusterd_pseudo.RData")
cds1 <-cds
load("normal_3clusterd_pseudo.RData")
cds_3d <- reduce_dimension(cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d)

# =============================================================================
# 13. CellChat communication analysis and subsequent exploration
# =============================================================================

cds_3d <- learn_graph(cds_3d)
cds_3d <- order_cells(cds_3d)
colData(cds)$assigned_cell_type <- dplyr::recode(colData(cds)$integrated_snn_res.0.15,
                                                 "0"="PreHTC",
                                                 "1"="HTC",
                                                 "3"="FC")
plot_cells(cds,
           color_cells_by = "integrated_snn_res.0.15",
           show_trajectory_graph = FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=0,
           group_label_size=0,cell_size=1)

cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="pseudotime")
cds_3d_plot_obj1 <- plot_cells_3d(cds_3d, color_cells_by="integrated_snn_res.0.15")
plot_cells_3d( cds_3d,color_cells_by="integrated_snn_res.0.15")
cds_3d_plot_obj
plot2 <- (cds@colData@listData)
plot2$integrated_snn_res.0.15
plot2 <-plot2$integrated_snn_res.0.15
plot2 <- as.data.frame(plot2)
plot3 <- subset(plot2, plot2 =="0"|plot2 =="1"|plot2 =="3")
position1 <-environment(cds_3d_plot_obj[["x"]][["visdat"]][["5f7417504302"]])[["data"]]
data <- position1
library(plotly)
Component1 <- data[, 1]
Component2 <- data[, 2]
Component3 <- data[, 3]
p <- plot_ly(data, x = ~Component1, y = ~Component2, z = ~Component3, color = ~data$seurat_clusters,
             colors = c("#DC143C","#0000FF","#20B2AA","#FFA500"),
             marker = list(size = 2)) %>%
  add_markers()
# Interactive help: ?plot_ly
p
cds_3d_plot_obj
class(cds_3d_plot_obj1)

# Interactive help: ?plot_cells_3d()
#########load and read
load("normal_5clusterd_pseudo.RData")
saveRDS(sce,file ="normal5clusters.rds")
sce<-readRDS(file ="normal5clusters.rds")
p3<-plot_genes_in_pseudotime(cds[c("HMOX1","NQO1","NFE2L2","KEAP1","SQSTM1"),], color_cells_by="integrated_snn_res.0.4",
                             min_expr=0, ncol = 1)
p3
DefaultAssay(sce) <- "RNA"
FeaturePlot(sce,reduction = "umap",features = c("HMOX1","NQO1","NFE2L2","KEAP1","SQSTM1"),pt.size=1)


rm(list = ls())

##cellchat
data.input <- GetAssayData(immune.combined, assay = "RNA", slot = "data")  # normalized data matrix
labels <- Idents(immune.combined)
head(labels)
identity <- data.frame(group = labels, row.names = names(labels))  # create a dataframe of the cell labels
cellchat <- createCellChat( data.input)
# Interactive help: ?createCellChat()
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels")  # set "labels" as default cell identity
levels(cellchat@idents)  # show factor levels of the cell labels
CellChatDB <- CellChatDB.human  # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
#> Rows: 1,939
#> Columns: 11
#> $ interaction_name   <chr> "TGFB1_TGFBR1_TGFBR2", "TGFB2_TGFBR1_TGFBR2", "TGF…
#> $ pathway_name       <chr> "TGFb", "TGFb", "TGFb", "TGFb", "TGFb", "TGFb", "T…
#> $ ligand             <chr> "TGFB1", "TGFB2", "TGFB3", "TGFB1", "TGFB1", "TGFB…
#> $ receptor           <chr> "TGFbR1_R2", "TGFbR1_R2", "TGFbR1_R2", "ACVR1B_TGF…
#> $ agonist            <chr> "TGFb agonist", "TGFb agonist", "TGFb agonist", "T…
#> $ antagonist         <chr> "TGFb antagonist", "TGFb antagonist", "TGFb antago…
#> $ co_A_receptor      <chr> "", "", "", "", "", "", "", "", "", "", "", "", ""…
#> $ co_I_receptor      <chr> "TGFb inhibition receptor", "TGFb inhibition recep…
#> $ evidence           <chr> "KEGG: hsa04350", "KEGG: hsa04350", "KEGG: hsa0435…
#> $ annotation         <chr> "Secreted Signaling", "Secreted Signaling", "Secre…
#> $ interaction_name_2 <chr> "TGFB1 - (TGFBR1+TGFBR2)", "TGFB2 - (TGFBR1+TGFBR2…

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")  # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)  # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4)  # do parallel
#> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
#> in future (>= 1.13.0) when running R from RStudio, because it is
#> considered unstable. Because of this, plan("multicore") will fall
#> back to plan("sequential"), and plan("multiprocess") will fall back to
#> plan("multisession") - not plan("multicore") as in the past. For more details,
#> how to control forked processing or not, and how to silence this warning in
#> future R sessions, see ?future::supportsMulticore
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat@netP$pathways
head(cellchat@LR$LRsig)
cellchat@netP$pathways
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
save(cellchat,file="normalcellchat.Rdata")
# Visualize outgoing interactions from each cell population.
mat <- cellchat@net$count
par(mfrow =c(3,3),xpd=T)

for (i in 1:nrow(mat)){
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat),title.name = rownames(mat)[i])
}

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = "TGFb", color.heatmap = "Reds")
#> Do heatmap based on a single object

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")  # the slot 'netP' means the inferred intercellular communication network of signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

# Bubble plot of significant ligand-receptor interactions.

levels(cellchat@idents)
# show all the significant interactions (L-R pairs)

# Specify the source and target cell populations.
p=netVisual_bubble(cellchat, sources.use = c(4),
                   targets.use = c(1,2,3,5,6,7),remove.isolate = FALSE)
p
##
# Circ diagram
pairLR.TGFb <- extractEnrichedLR(cellchat, signaling = 'TGFb', geneLR.return = FALSE)
LR.show <- pairLR.TGFb[1,]  # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,2,4,6)  # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
#> [[1]]
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
##TGFb
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = 'TGFb', layout = "circle")

##FGF
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = 'FGF', layout = "circle")

##VEGF
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = 'VEGF', layout = "circle")

##ANGPTL
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling =  "ANGPTL", layout = "circle")
##SPP1
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling =  "SPP1", layout = "circle")

##BMP
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling =  "BMP", layout = "circle")

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")  # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups TGFb
netAnalysis_signalingRole_network(cellchat, signaling = 'TGFb', width = 8, height = 2.5, font.size = 10)

# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups FGF
netAnalysis_signalingRole_network(cellchat, signaling = 'FGF', width = 8, height = 2.5, font.size = 10)

# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups VEGF
netAnalysis_signalingRole_network(cellchat, signaling = 'VEGF', width = 8, height = 2.5, font.size = 10)

# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups ANGPTL
netAnalysis_signalingRole_network(cellchat, signaling = 'ANGPTL', width = 8, height = 2.5, font.size = 10)

# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups BMP
netAnalysis_signalingRole_network(cellchat, signaling = 'BMP', width = 8, height = 2.5, font.size = 10)

##############
# this command creates an 'integrated' data assay
DefaultAssay(immune.combined) <- "RNA"
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- IntegrateData(anchorset = immune.anchors)
saveRDS(immune.combined,file ="all_normalizednormal.rds")

# =============================================================================
# 14. Reload combined data and repeat selected analyses
# =============================================================================

immune.combined <-readRDS(file =paths$path_11)
immune.combined[["percent.mt"]] <- PercentageFeatureSet(immune.combined, pattern ="^MT-")
immune.combined <- subset(immune.combined, subset = nCount_RNA<80000&percent.mt <5)
VlnPlot(immune.combined, features = c("nFeatures_RNA","nCount_RNA","percent.mt"),ncol = 3)
DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
ElbowPlot(immune.combined)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:15)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:15)
immune.combined <- FindClusters(immune.combined, resolution = 0.40)
table(immune.combined$orig.ident,immune.combined$replicate)
sample_normal <- sample(c("normal"),size = 9042,replace = TRUE)
sample_oa <- sample(c("oa"),size = 57198,replace = TRUE)
immune.combined$group <- c(sample_normal,sample_oa)

DimPlot(immune.combined, label = TRUE,pt.size = 1)

##
DimPlot(immune.combined,pt.size = 1,label = TRUE)
##EC
VlnPlot(immune.combined, features = c("TF","CHRDL2","TPM4","RARRES2","DSC2","CYSTM1","KIF9","FST","SNORC",'S100A2'),pt.size =0)  ##0
FeaturePlot(immune.combined,features = c("TF","CHRDL2","TPM4","RARRES2","DSC2","CYSTM1"),pt.size = 1)
FeaturePlot(immune.combined,features = c("FST","SNORC",'S100A2'),pt.size = 1)
###RegC
VlnPlot(immune.combined, features = c("KIF9","CKS2","TUBB4B","HMOX1","TMEM176A"),pt.size =0)  # 9
FeaturePlot(immune.combined,features = c("KIF9","CKS2","TUBB4B","HMOX1","TMEM176A"),pt.size = 1)
FeaturePlot(immune.combined,features = c("TMEM176A"),pt.size = 2)
###ProC
VlnPlot(immune.combined, features = c("P3H2","DDX21","GLIS3","ATP13A3","NOP16","MTHFD2","TXNRD1","ACTN1","ASNS","ERCC1","FAM101B","BMP6"),pt.size =0)
FeaturePlot(immune.combined,features = c("P3H2","DDX21","GLIS3","ATP13A3"),pt.size = 1,min.cutoff =0)  #5
##PreHTC
VlnPlot(immune.combined, features = c("TGFBI","S100A4","TNC","ABI3BP","CRTAC1","TNFAIP6","CRLF1"),pt.size =0)  #1,2
FeaturePlot(immune.combined,features = c("TGFBI","S100A4","TNC","ABI3BP","ABI3BP","CRTAC1","TNFAIP6","CRLF1"),pt.size = 1)
FeaturePlot(immune.combined,features = c("MMP1","CHI3L2","TNFAIP6"),pt.size = 1)
FeaturePlot(immune.combined,features = c("COL1A2"),pt.size = 1)
###FC
VlnPlot(immune.combined, features = c('COL1A1',	'COL1A2',	'TMSB4X',	'MTRNR2L2',	'S100A4',	'MGST3',	'MTRNR2L8',	'SERPINF1',	'ID3',	'IFITM2',	'MMP2',	'THY1',	'LGALS3BP',	'IFI16',	'PSME1',	'IGFBP7',	'IFI27',	'COL5A1',	'MTRNR2L1'),pt.size =0)  #6
FeaturePlot(immune.combined,features = c('COL1A1',	'COL1A2',	'TMSB4X',	'MTRNR2L2','OGN','S100A4'),pt.size = 1)  #6,3

##HTC
VlnPlot(immune.combined, features = c('RNU6-2',	'WWP2',	'COL10A1',	'BHLHE41',	'IBSP',	'RMRP',	'PLCG2',	'SCIN',	'S100A1',	'MEF2C',	'CRISPLD1',	'VTRNA1-1'),pt.size =0)  #6
FeaturePlot(immune.combined, features = c('WWP2',	'COL10A1','VTRNA1-1','BHLHE41',	'IBSP',	'RMRP',	'PLCG2'),pt.size =1)  #4
##HomC
VlnPlot(immune.combined, features = c('JUN',	'PPP1R10',	'BRD2',	'ID3',	'ID1',	'RGS16',	'HES1',	'NR4A1',	'RSRP1',	'TRA2B',	'RSRC2',	'NR4A2',	'CLK1',	'RASD1',	'CCNL1',	'IER2',	'NFKBIA',	'SNAI1',	'HEXIM1',	'KLF2'),pt.size =0)  #3
FeaturePlot(immune.combined,features = c('JUN','PPP1R10',	'BRD2',	'ID3',	'ID1'),pt.size = 1,min.cutoff = 1)  #8
##HBB HBA1+ C
VlnPlot(immune.combined, features = c('HBB','HBA1'),pt.size =0,y.max =10,sort = TRUE)  #Red Blood Cell 7
FeaturePlot(immune.combined,features = c('HBB','HBA1'),min.cutoff =0,pt.size = 1)
##spp1
FeaturePlot(immune.combined,features = c('SPP1'),min.cutoff =0,pt.size = 1)  # 4
FeaturePlot(immune.combined,features = c('ANGPTL4'),min.cutoff =0,pt.size = 1)  # 4
markers <- FindAllMarkers(immune.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)
write.table(markers,file = "oaallmarkersres0.4.txt",sep ='\t',quote = F)
markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10

DefaultAssay(immune.combined) <- "integrated"
DoHeatmap(immune.combined, features = top10$gene)
immune.combined@active.ident<- immune.combined$integrated_snn_res.0.4
immune.combined <- RenameIdents(immune.combined, `0` = "EC", `1` = "SPP1_C", `2` = "PreHTC",
                                `3` = "FC_2", `4` = "SPP1_C/HTC", `5` = "ProC", `6` = "FC_1", `7` = "RBC", `8` = "HomoC", `9` = "RegC")
immune.combined <- RenameIdents(immune.combined, "PreHTC_2" = "PreHTC_2/SPP1_C",
                                 "HTC" = "SPP1_C/HTC", "ProC" = "ProC", "FC_1" = "FC_1", "RBC" = "RBC", "HomoC" = "HomoC", "RegC" = "RegC")
saveRDS(immune.combined,file ="oa10clusters.rds")
immune.combined <- readRDS(file='oa10clusters.rds')

# Monocle 3 trajectory analysis using the RNA count matrix.
immune.combined1 <- GetAssayData(immune.combined, assay = 'RNA',slot ='counts')
cell_metadata <- immune.combined@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(immune.combined1))
rownames(gene_annotation) <- rownames(immune.combined1)

# Create a Monocle cell_data_set object.
# Create a Monocle cell_data_set object.
cds <- new_cell_data_set(immune.combined1,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds@colData$replicate<- immune.combined$replicate
#cds@colData$group<- immune.combined$group
# Run the preprocessing steps required by Monocle.

# Perform dimensionality reduction and clustering.
#NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)  # Preprocess the count matrix and compute principal components.
plot_pc_variance_explained(cds)  # Inspect the variance explained by principal components.
cds <- align_cds(cds, alignment_group = "replicate")

# Compute a UMAP embedding.
cds<- reduce_dimension(cds,preprocess_method = "PCA")  # Use PCA as the preprocessing method.


# Compute a t-SNE embedding.
#cds <- reduce_dimension(cds, reduction_method="tSNE")
#plot_cells(cds, reduction_method="tSNE", color_cells_by="seurat_clusters")

# Cluster cells.


cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(immune.combined, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed  # Reuse the Seurat UMAP coordinates after aligning cell names.
cds@colData$integrated_snn_res.0.4 <- immune.combined$integrated_snn_res.0.4

# Visualize the embedding.
plot_cells(cds, reduction_method="UMAP", color_cells_by="integrated_snn_res.0.4")
cds <- cluster_cells(cds)
cds <- learn_graph(cds)  # Learn the principal graph for trajectory inference.

# Color cells using the selected metadata field.
head(colData(cds))
p_plot_cells<-plot_cells(cds,
                         color_cells_by = "integrated_snn_res.0.15",
                         label_groups_by_cluster=FALSE,
                         label_leaves=FALSE,
                         label_branch_points=TRUE,
                         group_label_size=4,
                         cell_size=0.8)
#p_plot_cells <- p_plot_cells+scale_color_npg()
p_plot_cells
# Inspect the inferred cell-state trajectory.
# Black lines represent the principal graph.
# White circles label graph leaves or terminal states.
# Black circles label branch points.
# Control graph labels with label_leaves and label_branch_points.
# Order cells along the trajectory.
cds = order_cells(cds)  # Select a biologically justified root when ordering cells.
plot_cells(cds,
           color_cells_by = "pseudotime",
           show_trajectory_graph = FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=0,
           group_label_size=0,cell_size=1)
# Interactive help: ?plot_cells()

# Test for genes whose expression varies along the principal graph.
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=6)
write.csv(Track_genes, "normal_5_clustes_Trajectory_genes_all.csv", row.names = F)
Track_genes <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)

# Select the top ten genes for visualization.
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()

# Plot gene expression trends along pseudotime.
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="integrated_snn_res.0.15",
                         min_expr=0.5, ncol = 2)
p3<-plot_genes_in_pseudotime(cds[c("SPP1"),], color_cells_by="integrated_snn_res.0.15",
                             min_expr=0, ncol = 1)
# Interactive help: ?plot_genes_in_pseudotime()
p3 <-p3+scale_color_npg()
p3
p3 <- p3+scale_color_manual(values = c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00"))
p3
# Visualize selected genes on the cell embedding.
plot_cells(cds, genes=c("HMOX1","NQO1","NFE2L2","KEAP1"), show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)
# Interpret trajectory-associated genes alongside genes of biological interest.
# Identify modules of coexpressed genes.
genelist <- pull(Track_genes, gene_short_name) %>% as.character()
gene_module <- find_gene_modules(cds[genelist,], resolution=1e-2, cores = 10)
cell_group <- tibble::tibble(cell=row.names(colData(cds)),
                             cell_group=colData(cds)$integrated_snn_res.0.4)
agg_mat <- aggregate_gene_expression(cds, gene_module, cell_group)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2",color = colorRampPalette(colors = c("blue","white","red"))(100))
write.csv(gene_module,"5normalclusters_all_pseuodolgenemodule.csv",row.names = F )
###
pseudotime <-pseudotime(cds,reduction_method = 'UMAP')
pseudotime <-pseudotime[rownames(sce@meta.data)]
sce$pseudotime <- pseudotime
FeaturePlot(sce,reduction = "umap",features = "pseudotime")

save(cds,Track_genes, file = "normal_3clusterd_pseudo.RData")
cds1 <-cds
load("normal_5clusterd_pseudo.RData")
saveRDS(sce,file ="normal5clusters.rds")
rm(list = ls())


library(CellChat)

##cellchat oa
##cellchat
data.input <- GetAssayData(immune.combined, assay = "RNA", slot = "data")  # normalized data matrix
labels <- Idents(immune.combined)
head(labels)
identity <- data.frame(group = labels, row.names = names(labels))  # create a dataframe of the cell labels
cellchat <- createCellChat( data.input)
# Interactive help: ?createCellChat()
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels")  # set "labels" as default cell identity
levels(cellchat@idents)  # show factor levels of the cell labels
CellChatDB <- CellChatDB.human  # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
#> Rows: 1,939
#> Columns: 11
#> $ interaction_name   <chr> "TGFB1_TGFBR1_TGFBR2", "TGFB2_TGFBR1_TGFBR2", "TGF…
#> $ pathway_name       <chr> "TGFb", "TGFb", "TGFb", "TGFb", "TGFb", "TGFb", "T…
#> $ ligand             <chr> "TGFB1", "TGFB2", "TGFB3", "TGFB1", "TGFB1", "TGFB…
#> $ receptor           <chr> "TGFbR1_R2", "TGFbR1_R2", "TGFbR1_R2", "ACVR1B_TGF…
#> $ agonist            <chr> "TGFb agonist", "TGFb agonist", "TGFb agonist", "T…
#> $ antagonist         <chr> "TGFb antagonist", "TGFb antagonist", "TGFb antago…
#> $ co_A_receptor      <chr> "", "", "", "", "", "", "", "", "", "", "", "", ""…
#> $ co_I_receptor      <chr> "TGFb inhibition receptor", "TGFb inhibition recep…
#> $ evidence           <chr> "KEGG: hsa04350", "KEGG: hsa04350", "KEGG: hsa0435…
#> $ annotation         <chr> "Secreted Signaling", "Secreted Signaling", "Secre…
#> $ interaction_name_2 <chr> "TGFB1 - (TGFBR1+TGFBR2)", "TGFB2 - (TGFBR1+TGFBR2…

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")  # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)  # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4)  # do parallel
#> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
#> in future (>= 1.13.0) when running R from RStudio, because it is
#> considered unstable. Because of this, plan("multicore") will fall
#> back to plan("sequential"), and plan("multiprocess") will fall back to
#> plan("multisession") - not plan("multicore") as in the past. For more details,
#> how to control forked processing or not, and how to silence this warning in
#> future R sessions, see ?future::supportsMulticore
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat@netP$pathways
head(cellchat@LR$LRsig)
cellchat@netP$pathways
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Visualize outgoing interactions from each cell population.
mat <- cellchat@net$count
par(mfrow =c(3,4),xpd=T)
for (i in 1:nrow(mat)){
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat),title.name = rownames(mat)[i])
}

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = "TGFb", color.heatmap = "Reds")
#> Do heatmap based on a single object

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")  # the slot 'netP' means the inferred intercellular communication network of signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

# Analyze SPP1 ligand-receptor interactions.
netAnalysis_contribution(cellchat, signaling = "SPP1")

# Bubble plot of significant ligand-receptor interactions.

levels(cellchat@idents)
# show all the significant interactions (L-R pairs)

# Specify the source and target cell populations.
p=netVisual_bubble(cellchat, sources.use = c(2),
                   targets.use = c(1,3,4,5,6,7,8,9),remove.isolate = FALSE)
p
##
# Circ diagram
pairLR.TGFb <- extractEnrichedLR(cellchat, signaling = 'TGFb', geneLR.return = FALSE)
LR.show <- pairLR.TGFb[1,]  # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,2,4,6)  # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
#> [[1]]
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
#contribution
netAnalysis_contribution(cellchat, signaling = 'SPP1')
##TGFb
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = 'TGFb', layout = "circle")

##FGF
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = 'FGF', layout = "circle")

##VEGF
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = 'VEGF', layout = "circle")

##ANGPTL
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling =  "ANGPTL", layout = "circle")

##BMP
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling =  "BMP", layout = "circle")

##SPP1
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling =  "SPP1", layout = "circle")

# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")  # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups TGFb
netAnalysis_signalingRole_network(cellchat, signaling = 'TGFb', width = 8, height = 2.5, font.size = 10)

# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups FGF
netAnalysis_signalingRole_network(cellchat, signaling = 'FGF', width = 8, height = 2.5, font.size = 10)

# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups VEGF
netAnalysis_signalingRole_network(cellchat, signaling = 'SPP1', width = 8, height = 2.5, font.size = 10)

# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups ANGPTL
netAnalysis_signalingRole_network(cellchat, signaling = 'ANGPTL', width = 8, height = 2.5, font.size = 10)

# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups BMP
netAnalysis_signalingRole_network(cellchat, signaling = 'BMP', width = 8, height = 2.5, font.size = 10)

# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups BMP
netAnalysis_signalingRole_network(cellchat, signaling = 'SPP1', width = 8, height = 2.5, font.size = 10)
plotGeneExpression(cellchat, signaling = "SPP1")
plotGeneExpression(cellchat, signaling = "FGF")
plotGeneExpression(cellchat, signaling = "ANGPTL")
plotGeneExpression(cellchat, signaling = "TGFb")

###outgoing and incoming
library(NMF)
library(ggalluvial)

selectK(cellchat, pattern = "outgoing")
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
# river plot
netAnalysis_river(cellchat, pattern = "outgoing")
#> Please make sure you have load `library(ggalluvial)` when running this function

selectK(cellchat, pattern = "incoming")
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
netAnalysis_river(cellchat, pattern = "incoming")

###
cellchat <- computeNetSimilarity(cellchat, type = "structural")
cellchat <- netEmbedding(cellchat, type = "structural")
#> Manifold learning of the signaling networks for a single dataset
cellchat <- netClustering(cellchat, type = "structural")
#> Classification learning of the signaling networks for a single dataset
# Visualization in 2D-space
netVisual_embedding(cellchat, type = "structural", label.size = 3.5)
netVisual_embeddingZoomIn(cellchat, type = "structural", nCol = 2)


# =============================================================================
# 15. Red blood cell filtering and OA cartilage trajectories
# =============================================================================

# Remove residual red blood cells using HBB expression after the initial filtering step.
#sce <- subset(x=sce,idents=c( "EC", "SPP1_C", "PreHTC","FC_2","SPP1_C/HTC","ProC","FC_1","HomoC","RegC")  )
sce <- subset(sce,subset = HBB<6)
DimPlot(sce, reduction = 'umap')
DefaultAssay(sce) <- "integrated"
sce <- ScaleData(sce, verbose = FALSE)
sce <- RunPCA(sce, npcs =30, verbose = FALSE)
ElbowPlot(sce)
sce <- FindNeighbors(sce, dims = 1:15)
sce <- FindClusters(sce, resolution = 0.4 )

# Look at cluster IDs of the first 5 cells
DefaultAssay(sce) <- "RNA"
head(Idents(sce), 5)
table(sce$seurat_clusters)
sce <- RunUMAP(sce, dims = 1:15)
DimPlot(sce, reduction = 'umap')
##
DimPlot(sce,pt.size = 1,label = TRUE)
##EC
VlnPlot(sce, features = c("TF","CHRDL2","TPM4","RARRES2","DSC2","CYSTM1","KIF9","FST","SNORC",'S100A2'),pt.size =0)  ##0
FeaturePlot(sce,features = c("TF","CHRDL2","TPM4","RARRES2","DSC2","CYSTM1"),pt.size = 1)
FeaturePlot(sce,features = c("FST","SNORC",'S100A2'),pt.size = 1)
###RegC
VlnPlot(sce, features = c("KIF9","CKS2","TUBB4B","HMOX1","TMEM176A"),pt.size =0)  # 10
FeaturePlot(sce,features = c("KIF9","CKS2","TUBB4B","HMOX1","TMEM176A"),pt.size = 1)
FeaturePlot(sce,features = c("TMEM176A"),pt.size = 2)
###ProC
VlnPlot(sce, features = c("P3H2","DDX21","GLIS3","ATP13A3","NOP16","MTHFD2","TXNRD1","ACTN1","ASNS","ERCC1","FAM101B","BMP6"),pt.size =0)
FeaturePlot(sce,features = c("P3H2","DDX21","GLIS3","ATP13A3"),pt.size = 1,min.cutoff =0)  #5,4
##PreHTC
VlnPlot(sce, features = c("TGFBI","S100A4","TNC","ABI3BP","CRTAC1","TNFAIP6","CRLF1"),pt.size =0)  #2,3,
FeaturePlot(sce,features = c("TGFBI","S100A4","TNC","ABI3BP","ABI3BP","CRTAC1","TNFAIP6","CRLF1"),pt.size = 1)
FeaturePlot(sce,features = c("MMP1","CHI3L2","TNFAIP6"),pt.size = 1)
FeaturePlot(sce,features = c("COL1A2"),pt.size = 1)
###FC
VlnPlot(sce, features = c('COL1A1',	'COL1A2',	'TMSB4X',	'MTRNR2L2',	'S100A4',	'MGST3',	'MTRNR2L8',	'SERPINF1',	'ID3'),pt.size =0)  #
FeaturePlot(sce,features = c('COL1A1',	'COL1A2',	'TMSB4X',	'MTRNR2L2','OGN','S100A4'),pt.size = 1)  #7,4,2

##HTC
VlnPlot(sce, features = c('RNU6-2',	'WWP2',	'COL10A1',	'BHLHE41',	'IBSP',	'RMRP'),pt.size =0)  #6
FeaturePlot(sce, features = c('WWP2',	'COL10A1','VTRNA1-1','BHLHE41',	'IBSP',	'RMRP',	'PLCG2'),pt.size =1)  #6
##HomC
VlnPlot(sce, features = c('JUN',	'PPP1R10',	'BRD2',	'ID3',	'ID1',	'RGS16',	'HES1',	'NR4A1',	'RSRP1',	'TRA2B',	'RSRC2',	'NR4A2',	'CLK1',	'RASD1',	'CCNL1',	'IER2',	'NFKBIA',	'SNAI1',	'HEXIM1',	'KLF2'),pt.size =0)  #8
FeaturePlot(sce,features = c('JUN','PPP1R10',	'BRD2',	'ID3',	'ID1'),pt.size = 1,min.cutoff = 1)  #8
##HBB HBA1+ C
VlnPlot(immune.combined, features = c('HBB','HBA1'),pt.size =0,y.max =10,sort = TRUE)  #Red Blood Cell 7
FeaturePlot(sce,features = c('HBB','HBA1'),min.cutoff =0,pt.size = 1)
##spp1
FeaturePlot(sce,features = c('SPP1'),min.cutoff =0,pt.size = 1)  # 1,6
sce@active.ident <- sce$integrated_snn_res.0.4
sce <- RenameIdents(sce, `0` = "EC", `1` = "SPP1_C", `2` = "PreHTC/FC",
                    `3` = "PreHTC", `4` = "FC/ProC", `5` = "ProC", `6` = "HTC/SPP1_C", `7` = "FC_2", `8` = "HomoC", `9` = "undecided",`10` = "RegC")
FeaturePlot(sce,features=c("HMOX1","NQO1","NFE2L2","KEAP1"))  # 0 ProC

FeaturePlot(immune.combined,features=c("CD29", "CD44", "CD73", "CD90", "CD105"),pt.size = 1)  # 0 ProC
FeaturePlot(immune.combined,features=c("CD24"),pt.size = 1)  # 0 ProC
saveRDS(sce,file ="removedRBC_oa10clusters.rds")

immune.combined <- readRDS(file ="removedRBC_oa10clusters.rds")
DimPlot(immune.combined,reduction = "umap")
# Run Monocle 3 trajectory analysis on OA cartilage cells.

immune.combined1 <- GetAssayData(sce, assay = 'RNA',slot ='counts')
cell_metadata <- sce@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(immune.combined1))
rownames(gene_annotation) <- rownames(immune.combined1)

# Create a Monocle cell_data_set object.
# Create a Monocle cell_data_set object.
cds <- new_cell_data_set(immune.combined1,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds@colData$replicate<- sce$replicate
#cds@colData$group<- immune.combined$group
# Run the preprocessing steps required by Monocle.

# Perform dimensionality reduction and clustering.
#NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)  # Preprocess the count matrix and compute principal components.
plot_pc_variance_explained(cds)  # Inspect the variance explained by principal components.
cds <- align_cds(cds, alignment_group = "replicate")

# Compute a UMAP embedding.
cds<- reduce_dimension(cds,preprocess_method = "PCA")  # Use PCA as the preprocessing method.


# Compute a t-SNE embedding.
#cds <- reduce_dimension(cds, reduction_method="tSNE")
#plot_cells(cds, reduction_method="tSNE", color_cells_by="seurat_clusters")

# Cluster cells.


cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(sce, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed  # Reuse the Seurat UMAP coordinates after aligning cell names.
cds@colData$integrated_snn_res.0.4 <- sce$integrated_snn_res.0.4

# Visualize the embedding.
plot_cells(cds, reduction_method="UMAP", color_cells_by="integrated_snn_res.0.4",cell_size=0.8)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)  # Learn the principal graph for trajectory inference.

# Color cells using the selected metadata field.
head(colData(cds))
p_plot_cells<-plot_cells(cds,
                         color_cells_by = "integrated_snn_res.0.4",
                         label_groups_by_cluster=FALSE,
                         label_leaves=FALSE,
                         label_branch_points=TRUE,
                         group_label_size=4,
                         cell_size=0.8)
#p_plot_cells <- p_plot_cells+scale_color_npg()
p_plot_cells
# Inspect the inferred cell-state trajectory.
# Black lines represent the principal graph.
# White circles label graph leaves or terminal states.
# Black circles label branch points.
# Control graph labels with label_leaves and label_branch_points.
# Order cells along the trajectory.
cds = order_cells(cds)  # Select a biologically justified root when ordering cells.
plot_cells(cds,
           color_cells_by = "pseudotime",
           show_trajectory_graph = TRUE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           graph_label_size=0,
           group_label_size=0,cell_size=1)
# Interactive help: ?plot_cells()

# Test for genes whose expression varies along the principal graph.
write.csv(Track_genes, "oa_removedrbc_10clustes_Trajectory_genes_all.csv", row.names = F)
Track_genes <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)

# Select the top ten genes for visualization.
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()

# Plot gene expression trends along pseudotime.
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="integrated_snn_res.0.4",
                         min_expr=0.5, ncol = 2)
p3<-plot_genes_in_pseudotime(cds[c("FTH1","C2orf40","COL2A1","TNFRSF11B","CLU","SNORC","SPARC","ACAN","LUM","SPP1" ),], color_cells_by="seurat_clusters",
                           min_expr=0, ncol = 2)
Track_genes_sig
p3
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=6)

plot_genes_in_pseudotime(cds[c("FTH1","C2orf40"),], color_cells_by="integrated_snn_res.0.4",
                         min_expr=0, ncol = 2)
# Interactive help: ?plot_genes_in_pseudotime()
p3 <-p3+scale_color_npg()
p3
p3 <- p3+scale_color_manual(values = c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00"))
p3
# Visualize selected genes on the cell embedding.
plot_cells(cds, genes=c("HMOX1","NQO1","NFE2L2","KEAP1"), show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)
# Interpret trajectory-associated genes alongside genes of biological interest.
# Identify modules of coexpressed genes.
genelist <- pull(Track_genes, gene_short_name) %>% as.character()
gene_module <- find_gene_modules(cds[genelist,], resolution=1e-2, cores = 10)
cell_group <- tibble::tibble(cell=row.names(colData(cds)),
                             cell_group=colData(cds)$seurat_clusters)
agg_mat <- aggregate_gene_expression(cds, gene_module, cell_group)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2",color = colorRampPalette(colors = c("blue","white","red"))(100))
plot_cells(cds,
           genes=gene_module %>% filter(module %in% c(61)),
           group_cells_by="partition",
           color_cells_by="partition",
           show_trajectory_graph=FALSE)
write.csv(gene_module,"OA10clustersremovedRBC_pseuodolgenemodulen.csv",row.names = F )
save(cds,Track_genes, file = "oal_removedRBC10clusterd_pseudo.RData")
load("oal_removedRBC10clusterd_pseudo.RData")

# =============================================================================
# 16. Selected cartilage populations and additional trajectories
# =============================================================================

# Subset pre-HTC, FC, and SPP1 populations after red blood cell removal.
sce <- subset(x=sce,idents=c("SPP1_C", "PreHTC/FC","PreHTC","HTC/SPP1_C","FC_2")  )

DimPlot(sce, reduction = 'umap')
DefaultAssay(sce) <- "integrated"
sce <- ScaleData(sce, verbose = FALSE)
sce <- RunPCA(sce, npcs =30, verbose = FALSE)
ElbowPlot(sce)
sce <- FindNeighbors(sce, dims = 1:15)
sce <- FindClusters(sce, resolution = 0.4 )

# Look at cluster IDs of the first 5 cells
DefaultAssay(sce) <- "RNA"
head(Idents(sce), 5)
table(sce$seurat_clusters)
sce <- RunUMAP(sce, dims = 1:15)
DimPlot(sce, reduction = 'umap')
##
#DimPlot(sce,pt.size = 1,label = TRUE)
DefaultAssay(sce) <- "integrated"
markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)


markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(sce, features = top10$gene)
##EC
VlnPlot(sce, features = c("TF","CHRDL2","TPM4","RARRES2","DSC2","CYSTM1","KIF9","FST","SNORC",'S100A2'),pt.size =0)  ##0
FeaturePlot(sce,features = c("TF","CHRDL2","TPM4","RARRES2","DSC2","CYSTM1"),pt.size = 1)
FeaturePlot(sce,features = c("FST","SNORC",'S100A2'),pt.size = 1)
###RegC
VlnPlot(sce, features = c("KIF9","CKS2","TUBB4B","HMOX1","TMEM176A"),pt.size =0)  # 10
FeaturePlot(sce,features = c("KIF9","CKS2","TUBB4B","HMOX1","TMEM176A"),pt.size = 1)
FeaturePlot(sce,features = c("TMEM176A"),pt.size = 2)
###ProC
VlnPlot(sce, features = c("P3H2","DDX21","GLIS3","ATP13A3","NOP16","MTHFD2","TXNRD1","ACTN1","ASNS","ERCC1","FAM101B","BMP6"),pt.size =0)
FeaturePlot(sce,features = c("P3H2","DDX21","GLIS3","ATP13A3"),pt.size = 1,min.cutoff =0)  #5,4
##PreHTC
VlnPlot(sce, features = c("TGFBI","S100A4","TNC","ABI3BP","CRTAC1","TNFAIP6","CRLF1"),pt.size =0)  #0.3
FeaturePlot(sce,features = c("TGFBI","S100A4","TNC","ABI3BP","ABI3BP","CRTAC1","TNFAIP6","CRLF1"),pt.size = 1)
FeaturePlot(sce,features = c("MMP1","CHI3L2","TNFAIP6","VCAM1","ACKR1"),pt.size = 1)
FeaturePlot(sce,features = c("VCAM1","BMP4"),pt.size = 1)
FeaturePlot(sce,features = c("COL1A2"),pt.size = 1)
###FC
VlnPlot(sce, features = c('COL1A1',	'COL1A2',	'TMSB4X',	'MTRNR2L2',	'S100A4',	'MGST3',	'MTRNR2L8',	'SERPINF1',	'ID3'),pt.size =0)  #
FeaturePlot(sce,features = c('COL1A1',	'COL1A2',	'TMSB4X',	'MTRNR2L2','OGN','S100A4'),pt.size = 1)  #5

##HTC
VlnPlot(sce, features = c('RNU6-2',	'WWP2',	'COL10A1',	'BHLHE41',	'IBSP',	'RMRP'),pt.size =0)  #6
FeaturePlot(sce, features = c('WWP2',	'COL10A1','VTRNA1-1','BHLHE41',	'IBSP',	'RMRP',	'PLCG2'),pt.size =1)  #6,
##HomC
VlnPlot(sce, features = c('JUN',	'PPP1R10',	'BRD2',	'ID3',	'ID1',	'RGS16',	'HES1',	'NR4A1',	'RSRP1',	'TRA2B',	'RSRC2',	'NR4A2',	'CLK1',	'RASD1',	'CCNL1',	'IER2',	'NFKBIA',	'SNAI1',	'HEXIM1',	'KLF2'),pt.size =0)  #8
FeaturePlot(sce,features = c('JUN','PPP1R10',	'BRD2',	'ID3',	'ID1'),pt.size = 1,min.cutoff = 1)  #8
##HBB HBA1+ C
VlnPlot(immune.combined, features = c('HBB','HBA1'),pt.size =0,y.max =10,sort = TRUE)  #Red Blood Cell 7
FeaturePlot(sce,features = c('HBB','HBA1'),min.cutoff =0,pt.size = 1)
##spp1
FeaturePlot(sce,features = c('SPP1'),min.cutoff =0,pt.size = 1)  # 1,2,7,6
VlnPlot(sce, features = c('SPP1'),pt.size =0)  #6
sce@active.ident <- sce$integrated_snn_res.0.4
sce <- RenameIdents(sce, `0` = "PreHTC_1", `1` = "SPP1_C1", `2` = "SPP1_C2",
                    `3` = "PreHTC_2", `4` = "PreHTC/FC", `5` = "FC", `6` = "HTC/SPP1_C", `7` = "PreHTC/SPP1_C")
DimPlot(sce, reduction = 'umap',label = TRUE,pt.size = 1)
FeaturePlot(sce,features=c("HMOX1","NQO1","NFE2L2","KEAP1"))  # 0 ProC
saveRDS(sce,file ="removedRBCpre_htc_spp1_fc_oa8clusters.rds")

# Run Monocle 3 trajectory analysis on OA cartilage cells.

immune.combined1 <- GetAssayData(sce, assay = 'RNA',slot ='counts')
cell_metadata <- sce@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(immune.combined1))
rownames(gene_annotation) <- rownames(immune.combined1)

# Create a Monocle cell_data_set object.
# Create a Monocle cell_data_set object.
load(file = "oal_removedRBCprehtchtcfc8clusterd_pseudo.RData")
cds <- new_cell_data_set(immune.combined1,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds@colData$replicate<- sce$replicate
#cds@colData$group<- immune.combined$group
# Run the preprocessing steps required by Monocle.

# Perform dimensionality reduction and clustering.
#NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)  # Preprocess the count matrix and compute principal components.
plot_pc_variance_explained(cds)  # Inspect the variance explained by principal components.
cds <- align_cds(cds, alignment_group = "replicate")

# Compute a UMAP embedding.
cds<- reduce_dimension(cds,preprocess_method = "PCA")  # Use PCA as the preprocessing method.


# Compute a t-SNE embedding.
#cds <- reduce_dimension(cds, reduction_method="tSNE")
#plot_cells(cds, reduction_method="tSNE", color_cells_by="seurat_clusters")

# Cluster cells.


cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(sce, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed  # Reuse the Seurat UMAP coordinates after aligning cell names.
cds@colData$integrated_snn_res.0.4 <- sce$integrated_snn_res.0.4

# Visualize the embedding.
plot_cells(cds, reduction_method="UMAP", color_cells_by="integrated_snn_res.0.4")
cds <- cluster_cells(cds)
cds <- learn_graph(cds)  # Learn the principal graph for trajectory inference.

# Color cells using the selected metadata field.
head(colData(cds))
p_plot_cells<-plot_cells(cds,
                         color_cells_by = "integrated_snn_res.0.4",
                         label_groups_by_cluster=FALSE,
                         label_leaves=FALSE,
                         label_branch_points=TRUE,
                         group_label_size=4,
                         cell_size=0.8)
#p_plot_cells <- p_plot_cells+scale_color_npg()
p_plot_cells
###3d

cds_3d <- reduce_dimension(cds, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
cds_3d <- order_cells(cds_3d)

cds_3d_plot_obj <- plot_cells_3d(cds_3d, color_cells_by="pseudotime")
cds_3d_plot_obj1 <- plot_cells_3d(cds_3d, color_cells_by="integrated_snn_res.0.4")
cds_3d_plot_obj1 <- plot_cells_3d(cds_3d)
cds_3d_plot_obj1
cds_3d_plot_obj
plot2 <- (cds@colData@listData)
plot2$integrated_snn_res.0.4
plot2 <-plot2$integrated_snn_res.4
plot2 <- as.data.frame(plot2)
position1 <-   environment(cds_3d_plot_obj[["x"]][["visdat"]][["5f74116a1f28"]])[["data"]]
position2 <-   environment(cds_3d_plot_obj[["x"]][["visdat"]][["5f74116a1f28"]])[["data"]]
data <- position1
library(plotly)
Component1 <- data[, 1]
Component2 <- data[, 2]
Component3 <- data[, 3]
p <- plot_ly(data, x = ~Component1, y = ~Component2, z = ~Component3, color = ~data$seurat_clusters,
             colors = c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF"),
             marker = list(size = 2)) %>%
  add_markers()
# Interactive help: ?plot_ly
p
cds_3d_plot_obj1
class(cds_3d_plot_obj1)
# Inspect the inferred cell-state trajectory.
# Black lines represent the principal graph.
# White circles label graph leaves or terminal states.
# Black circles label branch points.
# Control graph labels with label_leaves and label_branch_points.
# Order cells along the trajectory.
cds = order_cells(cds)  # Select a biologically justified root when ordering cells.
plot_cells(cds,
           color_cells_by = "pseudotime",
           show_trajectory_graph = FALSE,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           graph_label_size=0,
           group_label_size=0,cell_size=1)
# Interactive help: ?plot_cells()

# Test for genes whose expression varies along the principal graph.
Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=6)
write.csv(Track_genes, "oa_removedrbc_prehtchtcspp17clustes_Trajectory_genes_all.csv", row.names = F)
Track_genes <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)

# Select the top ten genes for visualization.
Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()

# Plot gene expression trends along pseudotime.
plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="integrated_snn_res.0.4",
                         min_expr=0.5, ncol = 2)
Track_genes_sig
p3<-plot_genes_in_pseudotime(cds[c( "SPP1","ACAN","ADAMTS5")], color_cells_by="integrated_snn_res.0.4",
                             min_expr=0, ncol = 1)
# Interactive help: ?plot_genes_in_pseudotime()
p3 <-p3+scale_color_npg()
p3
p3 <- p3+scale_color_manual(values = c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00"))
p3
# Visualize selected genes on the cell embedding.
plot_cells(cds, genes=c("HMOX1","NQO1","NFE2L2","KEAP1"), show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)
# Interpret trajectory-associated genes alongside genes of biological interest.
# Identify modules of coexpressed genes.
genelist <- pull(Track_genes, gene_short_name) %>% as.character()
gene_module <- find_gene_modules(cds[genelist,], resolution=1e-2, cores = 10)
cell_group <- tibble::tibble(cell=row.names(colData(cds)),
                             cell_group=colData(cds)$integrated_snn_res.0.4)
agg_mat <- aggregate_gene_expression(cds, gene_module, cell_group)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat, scale="column", clustering_method="ward.D2",color = colorRampPalette(colors = c("blue","white","red"))(100))
write.csv(gene_module,"OA10clustersremovedRBC_pseuodolgenemodule.csv",row.names = F )
save(cds,Track_genes, file = "oal_removedRBCprehtchtcfc8clusterd_pseudo.RData")
###
cds1 <- cds

# =============================================================================
# 17. DDRTree exploration and gene-expression correlations
# =============================================================================

# Reduce dimensionality using the selected genes and the DDRTree method.
cds1 <- reduceDimension(cds1, max_components = 2, method = 'DDRTree')
FeaturePlot(immune.combined,features = "CD44",split.by = "group")
p1<-VlnPlot(immune.combined,features =c("CD44") ,split.by = "group",pt.size = 0)
p1 <- p1 +
  stat_compare_means(comparisons = my_comparisons,label = "p.signif") +
#scale_color_manual(values = my9color) + 
#scale_fill_manual(values = my9color) +
  theme(axis.text.x.bottom = element_text(angle = 90,vjust = 0.5,hjust = 1)) +
  NoLegend() + labs(x = '')
p1 + stat_compare_means(aes(group = group),method = "wilcox.test",label = "p.signif")
p1

p1<-p1 + stat_compare_means(method = "wilcox.test",label = "p.signif")+geom_boxplot()+NoLegend()
p1
p1<-VlnPlot(immune.combined,features=c("HMOX1","NQO1","NFE2L2","KEAP1","SQSTM1"),group.by = "group",pt.size = 0)

p1<-VlnPlot(immune.combined,features=c("HMOX1"),group.by = "group",pt.size = 0)
p1<-p1 + stat_compare_means(method = "wilcox.test",label = "p.signif")+geom_boxplot()+NoLegend()
p1

p2<-VlnPlot(immune.combined,features=c("NQO1"),group.by = "group",pt.size = 0)
p2<-p2 + stat_compare_means(method = "wilcox.test",label = "p.signif")+geom_boxplot()+NoLegend()
p2

p3<-VlnPlot(immune.combined,features=c("NFE2L2"),group.by = "group",pt.size = 0)
p3<-p3 + stat_compare_means(method = "wilcox.test",label = "p.signif")+geom_boxplot()+NoLegend()
p3

p4<-VlnPlot(immune.combined,features=c("KEAP1"),group.by = "group",pt.size = 0)
p4<-p4 + stat_compare_means(method = "wilcox.test",label = "p.signif")+geom_boxplot()+NoLegend()
p4
.libPaths()
# Interactive help: ?Vlnplot()

p5<-VlnPlot(immune.combined,features=c("COL2A1"),pt.size = 0)
p5
p6<-VlnPlot(immune.combined,features=c("SPP1"),pt.size = 0)
p5|p6
gene <- c("SPP1","COL10A1",	"COL13A1",	"COL14A1",	"COL16A1",	"COL17A1",	"COL19A1",	"COL1A1",	"COL1A2",	"COL21A1",	"COL23A1",	"COL25A1",	"COL2A1",	"COL3A1",	"COL4A1",	"COL4A3",	"COL4A4",	"COL4A5",	"COL4A6",	"COL5A1",	"COL5A2",	"COL5A3",	"COL6A1",	"COL7A1",	"COL9A1",	"COL9A2",	"COL9A3")
exprs <- data.frame(FetchData(object = immune.combined, vars = gene))
cor_exprs <- round(cor(exprs), 3)
pexprs <- cor_pmat(exprs)
ggcorrplot(cor_exprs,hc.order = T,  # Cluster cells.
           ggtheme = ggplot2::theme_void(base_size = 15),  # Customize the plot theme.
           colors = c("CornflowerBlue","white","Salmon"),  # Set a custom color palette.
           lab = T,lab_size = 5,  # Set the font size.
           tl.cex = 15,  # Set the font size.
           p.mat = pexprs,  # Provide the matrix of significance values.
           sig.level = 0.01,  # Set the significance threshold.
           pch = 4,  # Select the plotting symbol for nonsignificant correlations.
           pch.cex = 10,
           insig = "blank")  # Leave nonsignificant correlations blank.

exprs <-data.frame(GetAssayData(immune.combined, assay = 'RNA',slot ='counts'))
exprs
df_n2 <- data.frame(FetchData(immune.combined,vars = c("SPP1","CD44","ANGPTL4","FGF2")))
library(IDPmisc)
exprs <- data.frame(FetchData(object = immune.combined, vars = c("SPP1")))
df <-cbind(immune.combined$senmayo1,immune.combined$Angio1)
df_n <-cbind(exprs,df)
p <- ggplot(data=score_df,mapping = aes(x = "immune.combined$senmayo1",y = "immune.combined$Angio1")) +
  geom_point(shape=19 ,alpha = .01,size=0.5)+
  stat_binhex()+
  scale_fill_gradient(low = "lightblue", high = "red",limits =c(0,6000)) +  # Set the fill gradient and legend limits.
  theme_bw() +
  theme(
    legend.position = "left",  # Set the legend position.
    axis.title = element_text( family = "sans",size = 17,face = "bold",color = "gray25"),  # Set the axis-title font, size, and color.
    axis.text = element_text( family = "sans",size = 15,face = "bold"),  # Set the axis-label font and size.
    title = element_text( family = " sans ",size = 20,face = "bold"),  # Set the plot-title font and size.
    plot.title = element_text(margin = margin(t = 10, b = 10)),  # Set spacing around the plot title.
    axis.title.x = element_text(margin = margin(t = 15),family = "sans"),
    axis.title.y = element_text(margin = margin(r = 15) ,family = "sans"),  # Adjust spacing between the axis title and plotting area.
    legend.title = element_blank(),  # Remove the legend title.
    panel.border = element_rect( linetype = "solid", size = 2,color = "gray50")  # Set the panel-border color, line type, and width.
  ) +
  labs(title="Add marginal to scatter",hjust=0.5) + xlab("Senmayo") + ylab("Angiogenesis")  # Add the plot title and axis labels.
p

df_n2[df_n2==0] <- NA
df_n2<-na.omit(df_n2)
IDPmisc::iplot(x =score_df$`immune.combined$senmayo1`,y= score_df$`immune.combined$Angio1`,pix = 0.4,
           colramp = colorRampPalette(colors = c('black','blue','yellow','red')),xlab="Senmayo Score", ylab="Angiogenesis Score",main ='correlation plot with density')
IDPmisc::iplot(x =df_n$`SPP1`,y= df_n$`1`,pix = 0.4,
               colramp = colorRampPalette(colors = c('black','blue','yellow','red')),xlab="SPP1 Expression", ylab="Senmayo Score",main ='correlation plot with density')
IDPmisc::iplot(x =df_n2$`SPP1`,y= df_n2$`CD44`,pix = 0.4,
               colramp = colorRampPalette(colors = c('black','blue','yellow','red')),xlab="SPP1 Expression", ylab="Senmayo Score",main ='correlation plot with density')
###homogene
library(homologene)
a <- c("SPP1","AKT")
class(a)
class(gene_angio)
gene_angio <- as.character(gene_angio$gene)
gene_senmayo <- as.character(gene_senmayo$Gene)
genelistangio<-homologene(gene_senmayo, inTax = 9606, outTax = 10090)
write.table(genelistangio,file ="mmu_genelist_senmayo.txt",quote =F,sep = "\t")

##sc_metabolism
install.packages(c("devtools", "data.table", "wesanderson", "Seurat", "devtools", "AUCell", "GSEABase", "GSVA", "ggplot2","rsvd"))
devtools::install_github("YosefLab/VISION")
devtools::install_github("wu-yc/scMetabolism")
library(scMetabolism)
library(ggplot2)
library(rsvd)
load(file = "pbmc_demo.rda")
count <- GetAssayData(immune.combined, assay = 'RNA',slot ='counts')

# =============================================================================
# 18. Metabolic scores and group comparisons
# =============================================================================

immune.combined <-readRDS(file =paths$path_12)
countexp.Seurat<-sc.metabolism.Seurat(obj = immune.combined, method = "VISION", imputation = F, ncores = 2, metabolism.type = "KEGG")
saveRDS(countexp.Seurat,file ="countexp.Seurat.rds")
DotPlot.metabolism(obj = countexp.Seurat,
                   pathway = input.pathway,
                   phenotype = "ident",
                   norm = "y")
saveRDS(countexp.Seurat,file ="scmetabilsmn.rds")
countexp.Seurat <- readRDS(file ="scmetabilsmn.rds")
new.cluster.ids <- c("EC", "PreHTC", "SPP1+C", "PreHTC", "FC", "ProC",
                     "HTC", "RegC/HomoC", "RBC","undecidedC")
names(new.cluster.ids) <- levels(countexp.Seurat)
countexp.Seurat <- RenameIdents(countexp.Seurat, new.cluster.ids)
metabolism <- countexp.Seurat@assays[["METABOLISM"]]
metabolism_score<-metabolism$score
inputpathway <- row.names(metabolism_score)
countexp.Seurat $seurat_clusters <- countexp.Seurat@active.ident
BoxPlot.metabolism(countexp.Seurat,pathway = c ("Linoleic acid metabolism",	"alpha-Linolenic acid metabolism",
"Steroid hormone biosynthesis",	"Pentose and glucuronate interconversions",	"Sulfur metabolism",	"Steroid biosynthesis",
"Pantothenate and CoA biosynthesis",	"Galactose metabolism",	"Terpenoid backbone biosynthesis",	"Arachidonic acid metabolism",
"Nicotinate and nicotinamide metabolism",	"Metabolism of xenobiotics by cytochrome P450",	"Oxidative phosphorylation",
"Other types of O-glycan biosynthesis",	"Ether lipid metabolism",	"Riboflavin metabolism",	"Phenylalanine metabolism",
"Glycosaminoglycan biosynthesis - keratan sulfate",	"Porphyrin and chlorophyll metabolism",	"Starch and sucrose metabolism"
), phenotype = "seurat_clusters", ncol = 4)
p1<-BoxPlot.metabolism(countexp.Seurat,pathway = c ("Linoleic acid metabolism"), phenotype = "group", ncol = 1)
c<-row.names(test12)
class(c)
p1<-  p1
BoxPlot.metabolism(countexp.Seurat,pathway = c (c), phenotype = "group", ncol = 4)
p1<-BoxPlot.metabolism(countexp.Seurat,pathway = c ("Vitamin B6 metabolism"), phenotype = "group", ncol = 1)
p1 <- p1+geom_signif(
  comparisons = list(c("normal", "oa")),
  y_position = 1.0,
  map_signif_level = F,
  textsize = 7,)
p1
VlnPlot(countexp.Seurat)
DimPlot.metabolism(obj = countexp.Seurat, pathway = "Linoleic acid metabolism", dimention.reduction.type = "umap", dimention.reduction.run = F, size = 1)
p1<-BoxPlot.metabolism(countexp.Seurat,pathway = c(c),phenotype = "integrated_snn_res.0.4", ncol = 4)
p1
input.pathway<-c(inputpathway)
BoxPlot.metabolism(countexp.Seurat,pathway = input.pathway,phenotype = "group", ncol = 8)
p1<-BoxPlot.metabolism(countexp.Seurat,pathway = input.pathway,phenotype = "integrated_snn_res.0.4", ncol = 1)
DotPlot.metabolism(obj = countexp.Seurat,
                   pathway = input.pathway,
                   phenotype = "group")
DotPlot.metabolism(obj = countexp.Seurat,
                   pathway = input.pathway,
                   phenotype = "integrated_snn_res.0.4")
group1 <-as.data.frame(immune.combined$integrated_snn_res.0.4)
group2 <- as.data.frame(immune.combined$group)
group3 <- cbind(group1,group2)
metaboism_score2 <- t(metabolism_score)
metabolism_score3 <- cbind(group3,metaboism_score2)


group1 <-as.data.frame(immune.combined$seurat_clusters1)
group2 <- as.data.frame(immune.combined$group)
group3 <- cbind(group1,group2)
metaboism_score2 <- t(metabolism_score)
metabolism_score3 <- cbind(group3,metaboism_score2)
aggregate(metabolism_score3$`Glycolysis / Gluconeogenesis`, by=list(type=metabolism_score3$`immune.combined$group`),mean)
###
aggregate(metabolism_score3$`Glycolysis / Gluconeogenesis`, by=list(type=metabolism_score3$`immune.combined$group`),mean)
aggregate(metabolism_score3$`Citrate cycle (TCA cycle)`, by=list(type=metabolism_score3$`immune.combined$group`),mean)
aggregate(metabolism_score3$`Pentose phosphate pathway`, by=list(type=metabolism_score3$`immune.combined$group`),mean)
aggregate(metabolism_score3$`Pentose and glucuronate interconversions`, by=list(type=metabolism_score3$`immune.combined$group`),mean)
aggregate(metabolism_score3$`Glycolysis / Gluconeogenesis`, by=list(type=metabolism_score3$`immune.combined$group`),mean)
class(metabolism_score3$`Glycolysis / Gluconeogenesis`)
class(metabolism_score3[,1])
class(aggregate(as.numeric(metabolism_score3[,3]), by=list(type=metabolism_score3$`immune.combined$group`),mean))
a<-as.data.frame( metabolism_score3[,1])

df<-metabolism_score3
output <- list()  # Initialize a list for the summarized results.
for (i in seq_along(df)) {  # Iterate over the selected columns.
  output[[i]] <- aggregate(as.numeric(metabolism_score3[,i]), by=list(type=df$`immune.combined$group`),mean)  # Calculate the mean score within each group.
}
output  # Inspect the summarized results.


output2 <- list()  # Initialize a list for the summarized results.
for (i in seq_along(df)) {  # Iterate over the selected columns.
  output2[[i]] <- aggregate(as.numeric(metabolism_score3[,i]), by=list(type=df$`immune.combined$seurat_clusters1`),mean)  # Calculate the mean score within each group.
}
output2  # Inspect the summarized results.
output2[[80]]
output2[[80]][order(output2[[80]][,2]),]

df2<- as.data.frame(output)
write.table(df2,file = "metabolism_oa_ctrl.txt",sep = "\t",quote = F, row.names = F)
df1 <- as.data.frame(output2)
write.table(df1,file = "metabolism_10clusters.txt",sep = "\t",quote = F, row.names = F)
write.table(df1,file = "metabolism_n9clusters.txt",sep = "\t",quote = F, row.names = F)
colnames(metabolism_score3)
write.table(colnames(metabolism_score3),file = "metabolism_colname.txt",sep = "\t",quote = F, row.names = F)
test<- read.table(file = "TOP20_metabolism_normal_oa.txt",header= T,sep = "\t", row.names = 1)
test10<- read.table(file = "TOP10_metabolism_normal_oa.txt",header= T,sep = "\t", row.names = 1)
test11<-read.table(file = "TOP10_metabolism_10clusters.txt",header= T,sep = "\t", row.names = 1)

test$Term <- row.names(test)
test1$Term <- row.names(test1)
test1 <- read.table(file = "new_metabolism_clusters.txt",header= T,sep = "\t", row.names = 1)
test12<-read.table(file = "TOP20_metabolism_10clusters.txt",header= T,sep = "\t", row.names = 1)
test$Term <- row.names(test)
test1$Term <- row.names(test1)
test_merge <- merge(test,test1)
test_merge <- test_merge[,-c(2:3)]
rownames(test_merge) <- test_merge$Term
test_merge <- test_merge[,-1]
write.table(test_merge,file = "metabolism_clusters_unnormalto1.txt",sep = "\t",quote = F, row.names = T)
test<-read.table(file = "metabolism_clustersto1.txt",header= T,sep = "\t", row.names = 1)
test <- test_merge
library(pheatmap)
pheatmap(test,color = colorRampPalette(colors = c("blue","white","red"))(100),display_numbers = TRUE,  # Display numeric values in heatmap cells.
         number_format = "%.2f")  # Format numeric labels to two decimal places.
pheatmap(test10,color = colorRampPalette(colors = c("blue","white","red"))(100))
pheatmap(test11,color = colorRampPalette(colors = c("blue","white","red"))(100))
# Interactive help: ?pheatmap()
pheatmap(test12,color = colorRampPalette(colors = c("blue","white","red"))(100),display_numbers = TRUE,  # Display numeric values in heatmap cells.
         number_format = "%.2f")  # Format numeric labels to two decimal places.

df_cor <-cbind(immune.combined$Angio1,metabolism_score3$`Arachidonic acid metabolism`)
df_cor<- cbind(immune.combined$senmayo1,df_cor)
df_cor <- as.data.frame(df_cor)

IDPmisc::iplot(x =df_cor$'V1',y= df_cor$'V3',pix = 0.4,
               colramp = colorRampPalette(colors = c('black','blue','yellow','red')),xlab="Senmayo Score", ylab="Arachidonic acid metabolism Score",main ='correlation plot with density')

IDPmisc::iplot(x =df_cor$'V2',y= df_cor$'V3',pix = 0.4,
               colramp = colorRampPalette(colors = c('black','blue','yellow','red')),xlab="Angiogenesis Score", ylab="Arachidonic acid metabolism Score",main ='correlation plot with density')


p <- ggplot(data=df_cor,mapping = aes(x = "V1",y = "V3")) +
  geom_point(shape=19 ,alpha = .01,size=0.5)+
  stat_binhex()+
  scale_fill_gradient(low = "lightblue", high = "red",limits =c(0,6000)) +  # Set the fill gradient and legend limits.
  theme_bw() +
  theme(
    legend.position = "left",  # Set the legend position.
    axis.title = element_text( family = "sans",size = 17,face = "bold",color = "gray25"),  # Set the axis-title font, size, and color.
    axis.text = element_text( family = "sans",size = 15,face = "bold"),  # Set the axis-label font and size.
    title = element_text( family = " sans ",size = 20,face = "bold"),  # Set the plot-title font and size.
    plot.title = element_text(margin = margin(t = 10, b = 10)),  # Set spacing around the plot title.
    axis.title.x = element_text(margin = margin(t = 15),family = "sans"),
    axis.title.y = element_text(margin = margin(r = 15) ,family = "sans"),  # Adjust spacing between the axis title and plotting area.
    legend.title = element_blank(),  # Remove the legend title.
    panel.border = element_rect( linetype = "solid", size = 2,color = "gray50")  # Set the panel-border color, line type, and width.
  ) +
  labs(title="Add marginal to scatter",hjust=0.5) + xlab("Senmayo Score") + ylab("Arachidonic acid metabolism Score")  # Add the plot title and axis labels.
p

ggscatter(df_cor, x = "V1", y = "V3",size = 0.5,
          add = "reg.line", conf.int = TRUE,
          add.params = list(color = "red",  # Set the color of the regression line and confidence interval.
                            fill = "#FF9900")

)+
  stat_cor(method = "pearson",
           label.x = 0.2, label.y = 1.0)

ggscatter(df_cor, x = "V2", y = "V3",size = 0.5,
          add = "reg.line", conf.int = TRUE,
          add.params = list(color = "red",  # Set the color of the regression line and confidence interval.
                            fill = "#FF9900")

)+
  stat_cor(method = "spearman",
           label.x = 0.2, label.y = 1.0)

ggscatter(df_cor, x = "V1", y = "V2",size = 0.5,
          add = "reg.line", conf.int = TRUE,
          add.params = list(color = "red",  # Set the color of the regression line and confidence interval.
                            fill = "#FF9900")

)+
  stat_cor(method = "spearman",
           label.x = 0.2, label.y = 1.0)
FeaturePlot(immune.combined,features = c("CD44"),split.by = "group")
FeaturePlot(immune.combined,features = c("FKBP11"),split.by = "group")
VlnPlot(immune.combined,features = c("FKBP11"))
VlnPlot(immune.combined,features = c("FKBP11"),group.by = "group",pt.size = 0)
VlnPlot(immune.combined,features = c("KLF4"),group.by = "group",pt.size = 0)
VlnPlot(immune.combined,features = c("FABP4"),group.by = "group",pt.size = 0)
VlnPlot(immune.combined,features = c("KLF2"),group.by = "group",pt.size = 0)
VlnPlot(immune.combined,features = c("MMP3"),group.by = "group",pt.size = 0)
VlnPlot(immune.combined,features = c("CHI3L1"),group.by = "group",pt.size = 0)
VlnPlot(immune.combined,features = c("VEGFA"),group.by = "group",pt.size = 0)
VlnPlot(immune.combined,features = c("VEGFA"),pt.size = 0)
VlnPlot(immune.combined,features = c("PTGS2"),group.by = "group",,pt.size = 0)
p1<-VlnPlot(immune.combined,features = c("PTGS1","PTGS2","LOX"),group.by = "group",pt.size = 0)
p1<-VlnPlot(immune.combined,features = c("PTGS1","PTGS2","LOX"),pt.size = 0)
p1
p1<-VlnPlot(immune.combined,features = c("PTGS1"),group.by = "group",pt.size = 0.5)
p1<-VlnPlot(immune.combined,features = c("PTGS2"),group.by = "group",pt.size = 0.5)
p1<-VlnPlot(immune.combined,features = c("LOX"),group.by = "group",pt.size = 0.5)
p1 <-p1+geom_signif(
  comparisons = list(c("normal", "oa")),
  y_position = 4,
  map_signif_level = F,
  textsize = 5,)
p1
FeaturePlot(immune.combined,features =  c("PTGS1","PTGS2","LOX"),split.by = "group",pt.size = 1)
FeaturePlot(immune.combined,features = c("CXCL8","IL11"))
FeaturePlot(immune.combined,features =  c("PIEZO1"),split.by = "group",pt.size = 1)
VlnPlot(immune.combined,features = c("PIEZO1"),group.by = "group",pt.size = 0)
VlnPlot(immune.combined,features = c("PECAM1","EMCN"),pt.size = 0)
FeaturePlot(immune.combined,features =  c("PECAM1","EMCN"),split.by = "group",pt.size = 2)
FeaturePlot(immune.combined,features =  c("COL2A1","ACAN","SPP1","SOX9"),pt.size = 1)
FeaturePlot(immune.combined,features =  c("CYP7A1","CYP7B1","CH25H"),pt.size = 1)
VlnPlot(immune.combined,features = c("CYP7A1","CYP7B1","CH25H"),pt.size = 1,sort = T)
###
BiocManager::install('Biobase')
library(Biobase)
library(NMF)
library(Seurat)
vm <- immune.combined@assays$integrated@scale.data
res <- nmf(vm,10,method = "snmf/r")
rm(immune.combined)
immune.combined_nmf <- NormalizeData(immune.combined) %>% FindVariableFeatures() %>% ScaleData(do.center = FALSE)
# Interactive help: ?ScaleData()
vm <- immune.combined_nmf@assays$RNA@scale.data
res <- nmf(vm,11,method = "snmf/r")
runtime(res)
res_10 <- res
res_11 <- res
fs <- extractFeatures(res,30L)
fs <- lapply(fs,function(x)rownames(res)[x])
fs <- do.call("rbind",fs)
rownames(fs) <- paste0("cluster",1:11)


immune.combined_nmf@reductions$nmf <- immune.combined_nmf$pca
immune.combined_nmf@reductions$nmf@cell.embeddings <- t(coef(res))
immune.combined_nmf@reductions$nmf@features.loading <- basis(res)

s.f <- 1:11
cell1 <- colnames(immune.combined_nmf)
cell2 <- colnames(coef(res))
cells <- intersect(cell1,cell2)

immune.combined_nmf <- immune.combined_nmf[,cells]
immune.combined_nmf <- RunPCA(immune.combined_nmf,verbose = F)
immune.combined_nmf@reductions$nmf <- immune.combined_nmf@reductions$pca
immune.combined_nmf@reductions$nmf@cell.embeddings <- t(coef(res)[,cells])
immune.combined_nmf@reductions$nmf@feature.loadings <- basis(res)
immune.combined_nmf <- RunUMAP(immune.combined_nmf,reduction = 'nmf',dims = s.f)

immune.combined_nmf <- FindNeighbors(immune.combined_nmf,dims = s.f)%>% FindClusters()
p <- DimPlot(immune.combined_nmf)
p
set.seed
undeciedC = immune.combined[, Idents(immune.combined) %in% c( "undecidedC" )]
DimPlot(undeciedC,reduction = "umap", label = TRUE, repel = TRUE)
P1<-DimPlot(undeciedC,reduction = "umap", label = TRUE, repel = TRUE,cols = c("undecidedC"= "#7CFC00"))
VlnPlot(immune.combined,features = c("TF"),pt.size = 0)
P1|P2

P1 <-VlnPlot(immune.combined,features = 'Angio1',pt.size = 0,sort = T)
P1 <- P1 +geom_boxplot()
P1

##cellchat_all
immune.combined <- readRDS(file = "nallcombined1.rds")
immune.combined_normal <- subset(immune.combined, group=="normal")
immune.combined_oa <- subset(immune.combined, group=="oa")
cellchat1


# =============================================================================
# 19. Additional cell communication analysis
# =============================================================================


##cellchat oa
##cellchat
data.input <- GetAssayData(immune.combined_oa, assay = "RNA", slot = "data")  # normalized data matrix
labels <- Idents(immune.combined_oa)
head(labels)
identity <- data.frame(group = labels, row.names = names(labels))  # create a dataframe of the cell labels
cellchat <- createCellChat( data.input)
# Interactive help: ?createCellChat()
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels")  # set "labels" as default cell identity
levels(cellchat@idents)  # show factor levels of the cell labels
CellChatDB <- CellChatDB.human  # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)
#> Rows: 1,939
#> Columns: 11
#> $ interaction_name   <chr> "TGFB1_TGFBR1_TGFBR2", "TGFB2_TGFBR1_TGFBR2", "TGF…
#> $ pathway_name       <chr> "TGFb", "TGFb", "TGFb", "TGFb", "TGFb", "TGFb", "T…
#> $ ligand             <chr> "TGFB1", "TGFB2", "TGFB3", "TGFB1", "TGFB1", "TGFB…
#> $ receptor           <chr> "TGFbR1_R2", "TGFbR1_R2", "TGFbR1_R2", "ACVR1B_TGF…
#> $ agonist            <chr> "TGFb agonist", "TGFb agonist", "TGFb agonist", "T…
#> $ antagonist         <chr> "TGFb antagonist", "TGFb antagonist", "TGFb antago…
#> $ co_A_receptor      <chr> "", "", "", "", "", "", "", "", "", "", "", "", ""…
#> $ co_I_receptor      <chr> "TGFb inhibition receptor", "TGFb inhibition recep…
#> $ evidence           <chr> "KEGG: hsa04350", "KEGG: hsa04350", "KEGG: hsa0435…
#> $ annotation         <chr> "Secreted Signaling", "Secreted Signaling", "Secre…
#> $ interaction_name_2 <chr> "TGFB1 - (TGFBR1+TGFBR2)", "TGFB2 - (TGFBR1+TGFBR2…

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")  # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)  # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4)  # do parallel
#> Warning: [ONE-TIME WARNING] Forked processing ('multicore') is disabled
#> in future (>= 1.13.0) when running R from RStudio, because it is
#> considered unstable. Because of this, plan("multicore") will fall
#> back to plan("sequential"), and plan("multiprocess") will fall back to
#> plan("multisession") - not plan("multicore") as in the past. For more details,
#> how to control forked processing or not, and how to silence this warning in
#> future R sessions, see ?future::supportsMulticore
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

cellchat_oa <- cellchat
cellchat_normal <- cellchat

object.list <- list(normal = cellchat_normal, oa = cellchat_oa)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
# Interactive help: ?mergeCellChat()
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2

saveRDS(cellchat, file = "cellchat_comparisonAnalysis.rds")

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
gg1 + gg2

cellchat@netP$pathways
head(cellchat@LR$LRsig)
cellchat@netP$pathways
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Visualize outgoing interactions from each cell population.
mat <- cellchat@net$count
par(mfrow =c(3,4),xpd=T)
for (i in 1:nrow(mat)){
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat),title.name = rownames(mat)[i])
}

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = "TGFb", color.heatmap = "Reds")
#> Do heatmap based on a single object

# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
# Compute the network centrality scores
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")  # the slot 'netP' means the inferred intercellular communication network of signaling pathways
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming")
ht1 + ht2

# Analyze SPP1 ligand-receptor interactions.
netAnalysis_contribution(cellchat, signaling = "SPP1")

# Bubble plot of significant ligand-receptor interactions.

levels(cellchat@idents)
# show all the significant interactions (L-R pairs)

# Specify the source and target cell populations.
p=netVisual_bubble(cellchat, sources.use = c(2),
                   targets.use = c(1,3,4,5,6,7,8,9),remove.isolate = FALSE)
p
##
# Circ diagram
pairLR.TGFb <- extractEnrichedLR(cellchat, signaling = 'TGFb', geneLR.return = FALSE)
LR.show <- pairLR.TGFb[1,]  # show one ligand-receptor pair
# Hierarchy plot
vertex.receiver = seq(1,2,4,6)  # a numeric vector
netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
#> [[1]]
# Circle plot
netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
#contribution
netAnalysis_contribution(cellchat, signaling = 'SPP1')
##TGFb
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = 'TGFb', layout = "circle")

##FGF
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = 'FGF', layout = "circle")

##VEGF
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = 'VEGF', layout = "circle")

##ANGPTL
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling =  "ANGPTL", layout = "circle")

##BMP
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling =  "BMP", layout = "circle")

##SPP1
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling =  "SPP1", layout = "circle")

remove.packages('rlang')
install.packages('rlang')

# =============================================================================
# 20. Scissor demonstration and OA phenotype exploration
# =============================================================================

# Install and load Scissor.
remove.packages('vctrs')
install.packages("vctrs")

# Install the package using devtools.

# Load the package.
remotes::install_local(paths$path_13)
library("Scissor")
library(Scissor)

# Load the LUAD demonstration single-cell dataset.

# Load single-cell data.
location <- "https://xialab.s3-us-west-2.amazonaws.com/Duanchen/Scissor_data/"
load(url(paste0(location, 'scRNA-seq.RData')))
load('scRNA-seq.RData')
sc1 <- sc_dataset
## bulk-seq data
load(url(paste0(location, 'TCGA_LUAD_exp1.RData')))
load( 'TCGA_LUAD_exp1.RData')
# Load phenotype data.
load(url(paste0(location, 'TCGA_LUAD_survival.RData')))

# Load single-cell data.
sc_dataset <- Seurat_preprocessing(sc_dataset, verbose = F)

# Prepare phenotype data.
phenotype <- bulk_survival[,2:3]
colnames(phenotype) <- c("time", "status")

# Run Scissor to associate cells with bulk-sample phenotypes.
infos1 <- Scissor(bulk_dataset, sc_dataset, phenotype, alpha = 0.05,
                  family = "cox", Save_file = 'Scissor_LUAD_survival.RData')


Scissor_select <- rep(0, ncol(sc_dataset))
names(Scissor_select) <- colnames(sc_dataset)
# Scissor_pos contains cells with positive regression coefficients (beta > 0).
# Scissor_neg contains cells with negative regression coefficients (beta < 0).
Scissor_select[infos1$Scissor_pos] <- 1
Scissor_select[infos1$Scissor_neg] <- 2
sc_dataset <- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = "scissor")
DimPlot(sc_dataset, reduction = 'umap', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))
bulk <- bulk_dataset
phenotype2 <- phenotype
###scissor SPP1
sc2 <- readRDS(file = "nallcombined1.rds")
sc2 <- subset(x = sc2, downsample = 100)
# TODO: The original draft contains an incomplete assignment. Define sc3 before using it.
# sc3 <- NULL
bulk1 <- read.table(file = paths$path_14,header =T,sep ='\t')
bulk2 <- read.table(file = paths$path_15,header =T,sep ='\t')
bulk1 <- bulk1[,1:19]
bulk2 <- bulk2[,1:21]
bulk_merge <- read.table(file = paths$path_16,header =T,sep ='\t')
bulk_merge <- as.matrix(bulk_merge)
row.names(bulk_merge) <- bulk_merge[,1]
bulk_merge <- bulk_merge[,-1]
boxplot(bulk_merge)
class(bulk_merge)
summary(bulk_merge)
f <- as.data.frame( c(rep("1", 20), rep("0", 18)))
phenotype <-f
colnames(f)[1] <-"disease"
class(phenotype)


# Load single-cell data.
sc_dataset <- sc2

sc_dataset <- as.matrix( GetAssayData(object = sc_dataset, slot = "counts"))
head(sc_dataset)
dim(sc_dataset)
sc_dataset <- Seurat_preprocessing(sc_dataset, verbose = F)

sc_dataset1 <- sc_dataset

# Prepare phenotype data.
phenotype <- bulk_survival[,2:3]
colnames(phenotype) <- c("time", "status")

# Run Scissor to associate cells with bulk-sample phenotypes.
bulk_dataset2 <- bulk_dataset
bulk_dataset <- bulk_merge
# Interactive help: ?Scissor()

tag <- c('oa', 'normal')
infos4 <- Scissor(bulk_dataset, sc_dataset, phenotype, tag = tag, alpha = 0.5,
                  family = "binomial", Save_file = "Scissor_LUAD_TP53_mutation.RData")
infos1 <- Scissor(bulk_dataset, sc_dataset, phenotype,tag = tag, alpha = 0.05,
                  family = "binomial")


Scissor_select <- rep(0, ncol(sc_dataset))
names(Scissor_select) <- colnames(sc_dataset)
# Scissor_pos contains cells with positive regression coefficients (beta > 0).
# Scissor_neg contains cells with negative regression coefficients (beta < 0).
Scissor_select[infos1$Scissor_pos] <- 1
Scissor_select[infos1$Scissor_neg] <- 2
sc_dataset <- AddMetaData(sc_dataset, metadata = Scissor_select, col.name = "scissor")
DimPlot(sc_dataset, reduction = 'umap', group.by = 'scissor', cols = c('grey','indianred1','royalblue'), pt.size = 1.2, order = c(2,1))
bulk <- bulk_dataset
phenotype2 <- phenotype

newLimit <- 1200000
memory.limit()
memory.limit(newLimit)

####

# =============================================================================
# 21. ESR1 expression subsets and cartilage-cell comparisons
# =============================================================================

# Explore ESR1 expression-based subsets; inspect the thresholds below.
VlnPlot(immune.combined,features = "ESR1",group.by = "group")

subset_ESR1<-subset(immune.combined, subset = ESR1 > 1)
mydata_upspp1 <- FetchData(subset_upspp1,vars = c("UMAP_1","UMAP_2","group","seurat_clusters1","Angio1"))
p_upspp1 <-ggplot(data = mydata_upspp1,aes(x=UMAP_1,y=UMAP_2,
                                           colour=Angio1))+geom_point(size=1)+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#FCB216FF','#E25834FF','#992766FF','#450A69FF','#000004FF'))+
  geom_point(size = 1.5)+theme_bw()+theme(panel.grid.major =element_blank(),
                                          panel.grid.minor = element_blank(),
                                          panel.background = element_blank(),  # Remove the panel background.
                                          panel.border = element_blank(),axis.line = element_line(color = 'black'))
p_upspp1

subset_dnspp1<-subset(immune.combined, subset = SPP1 <= 1)
mydata_dnspp1 <- FetchData(subset_dnspp1,vars = c("UMAP_1","UMAP_2","group","seurat_clusters1","Angio1"))
p_dnspp1 <-ggplot(data = mydata_dnspp1,aes(x=UMAP_1,y=UMAP_2,
                                           colour=Angio1))+geom_point(size=1)+
  scale_color_gradientn(values = seq(0,1,0.2),colours = c('#FCB216FF','#E25834FF','#992766FF','#450A69FF','#000004FF'))+
  geom_point(size = 1.5)+theme_bw()+theme(panel.grid.major =element_blank(),
                                          panel.grid.minor = element_blank(),
                                          panel.background = element_blank(),  # Remove the panel background.
                                          panel.border = element_blank(),axis.line = element_line(color = 'black'))
p_dnspp1

P1<-FeaturePlot(subset_upspp1,features = c('PDPN','COL2A1','OGN','MMP13','PTGS2','FN1'),pt.size = 1)
P2<-FeaturePlot(subset_dnspp1,features = c('PDPN','COL2A1','OGN','MMP13','PTGS2','FN1'),pt.size = 1)
P1|P2

P3 <- FeaturePlot(immune.combined,features = c('PDPN','SPP1'))


###immune.combined chondrocyte
FeaturePlot(immune.combined,features = c('COL2A1','ACAN'),pt.size = 1)
new.cluster.ids <- c("Chondrocyte", "Chondrocyte", "Chondrocyte", "Chondrocyte", "Chondrocyte", "Chondrocyte",
                     "Chondrocyte", "Chondrocyte", "RBC","Chondrocyte")
names(new.cluster.ids) <- levels(immune.combined)
immune.combined <- RenameIdents(immune.combined, new.cluster.ids)

DotPlot(immune.combined, features = c('COL2A1','ACAN','HBB','HBA1'))+RotatedAxis()+
  scale_x_discrete("")+scale_y_discrete("")
VlnPlot(immune.combined,features = c('COL2A1','ACAN'),pt.size = 0)
DimPlot(immune.combined, reduction = "umap", label = FALSE, repel = TRUE)
# Subset cartilage cells; verify the intended expression threshold before grouping.
subset_cartilage <- subset(immune.combined,idents = c("EC", "SPP1+ C_1","SPP1+ C_2",
                                                      "PreHTC_1", "PreHTC_2/FC",  "ProC", "HTC","RegC/HomoC", "Undecided"))
mydata_cc <- FetchData(subset_cartilage,vars = c("Phase","group","seurat_clusters1","Angio1","senmayo1"))
ccdata_normal <- subset(mydata_cc,group == 'normal')
table(ccdata_normal$seurat_clusters1,ccdata_normal$Phase)
ccdata_oa <- subset(mydata_cc,group == 'oa')
table(ccdata_oa$seurat_clusters1,ccdata_oa$Phase)

P3

# Analyze ESR1 expression and associated cell states.
immune.combined@active.ident
subset_cartilage<-subset(immune.combined, idents = c("EC", "SPP1+ C_1","SPP1+ C_2",
                                                     "PreHTC_1", "PreHTC_2/FC",  "ProC", "HTC","RegC/HomoC", "Undecided"))
FeaturePlot(subset_cartilage,features = c('ESR1'),pt.size = 1,order = T)
subset_esr1<-subset(subset_cartilage, subset = ESR1 > 0)
VlnPlot(immune.combined,features = c("ESR1"))
DimPlot(subset_esr1)
subset_esr1$esr1 <- sample(c("ESR1+"), size = ncol(subset_esr1), replace = TRUE)
meta1 <- subset_esr1$esr1
meta1
subset_esr1_0<-subset(subset_cartilage, subset = ESR1 <= 0)
subset_esr1_0$esr1 <- sample(c("ESR1-"), size = ncol(subset_esr1_0), replace = TRUE)
meta2 <- subset_esr1_0$esr1
meta2 <- c(meta1,meta2)
subset_cartilage$esr1 <- meta2
library("ggpubr")
library("ggplot2")
DefaultAssay(subset_cartilage) <- "integrated"
subset_cartilage <- ScaleData(subset_cartilage, verbose = FALSE)
subset_cartilage <- RunPCA(subset_cartilage, npcs = 30, verbose = FALSE)
ElbowPlot(subset_cartilage)
subset_cartilage <- RunUMAP(subset_cartilage, reduction = "pca", dims = 1:20)
subset_cartilage <- FindNeighbors(subset_cartilage, reduction = "pca", dims = 1:20)
subset_cartilage <- FindClusters(subset_cartilage, resolution = 0.2)
DimPlot(subset_cartilage,reduction = 'umap',pt.size = 1)
DefaultAssay(subset_cartilage) <- "RNA"
FeaturePlot(subset_cartilage,features="ESR1",order = T,pt.size = 1)
# NOTE: This process can take a long time for big datasets, comment out for expediency. More
# approximate techniques such as those implemented in ElbowPlot() can be used to reduce
# computation time

p1<-VlnPlot(subset_cartilage,features = c("PTGS2","MMP13","ADAMTS5","COL2A1","ACAN","SOX9"),group.by = "esr1",pt.size = 0)
p1
p1<-VlnPlot(subset_cartilage,features = c("PTGS2"),group.by = "esr1",pt.size = 0)
p1<-p1+geom_boxplot()
p1
p2<-VlnPlot(subset_cartilage,features = c("PTGS2"),group.by = "esr1",pt.size = 0)
p2<-p1+geom_boxplot()
p2
p3<-VlnPlot(subset_cartilage,features = c("senmayo1"),group.by = "esr1",pt.size = 0)
p3<-p3+geom_boxplot()+stat_compare_means(method ="wilcox.test")
p3
p2<-VlnPlot(subset_cartilage,features = c("COL2A1"),group.by = "esr1",pt.size = 0)
p2<-p2+geom_boxplot()
p2
p3<-VlnPlot(immune.combined,features = c("COL2A1"),group.by = "esr1")
p3<-p3+geom_boxplot()
p3
p4<-VlnPlot(subset_cartilage,features = c("ACAN"),group.by = "esr1",pt.size = 0)
p4<-p4+geom_boxplot()
p4
p5<-VlnPlot(immune.combined,features = c("SOX9"),group.by = "esr1")
p5<-p5+geom_boxplot()
p5
p6<-VlnPlot(subset_cartilage,features = c("MMP3"),group.by = "esr1",pt.size = 0)
p6<-p6+geom_boxplot()
p6
p7<-VlnPlot(immune.combined,features = c("ADAMTS5"),group.by = "esr1")
p7<-p7+geom_boxplot()
p7


gene_bp<-read.table("Go_BP_C.txt",header = T)
gene<- as.list(gene_bp)
subset_cartilage <- AddModuleScore(object = subset_cartilage,features = gene,name = "GO_BP_C")
p13<-VlnPlot(subset_cartilage,features = 'GO_BP_C1',group.by = "group",pt.size = 0)
p13
p14<-VlnPlot(subset_cartilage,features = 'GO_BP_C1',group.by = "esr1",pt.size = 0)
p14 <- p14+geom_boxplot()+stat_compare_means(method ="wilcox.test")
p14


##cell cycle
cc.genes.updated.2019
subset_cartilage <- CellCycleScoring(
  object = subset_cartilage,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes)
VlnPlot(subset_cartilage, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB","G2M.Score","S.Score"), ncol = 6)  #+scale_color_npg()
####cell cycle
umapem<-immune.combined@reductions$umap@cell.embeddings
umapem<-as.data.frame(umapem)
metada= immune.combined@meta.data
dim(umapem);dim(metada)

metada$bar<-rownames(metada)
umapem$bar<-rownames(umapem)
ccdata<-merge(umapem,metada,by="bar")
head(ccdata)
library(ggplot2)
plot<-ggplot(ccdata, aes(UMAP_1, UMAP_2,label=Phase))+geom_point(aes(colour = factor(Phase)))+
#plot<-plot+scale_colour_manual(values=c("#CC33FF","Peru","#660000","#660099","#990033","black","red", "#666600", "green","#6699CC","#339900","#0000FF","#FFFF00","#808080"))+
  labs("@yunlai",x = "", y="")
plot=plot+scale_color_aaas()  +
  theme_bw()+theme(panel.grid=element_blank(),legend.title=element_blank(),legend.text = element_text(color="black", size = 10, face = "bold"))
plot<-plot+guides(colour = guide_legend(override.aes = list(size=5))) +theme(plot.title = element_text(hjust = 0.5))
plot


subset_cartilage@active.ident<- subset_cartilage$esr1

saveRDS(subset_cartilage,file = "subset_cartilage.rds")

