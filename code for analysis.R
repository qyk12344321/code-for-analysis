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
library(ggsignif)
library(Seurat)
library(SeuratData)
install.packages('dpylr')
library(patchwork)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c('dpylr'))
library(dplyr)
library(magrittr)
library(limma)
library(monocle3)
library("ggsci")
install.packages("psych")
library(psych)
library(pheatmap)
library(corrplot)
BiocManager::install('multtest')
install.packages('metap')
library(multtest)
library(metap)
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
#安装R包
install.packages("ggplot2")
install.packages("ggridges")
install.packages("reshape")
install.packages("ggprism")
#加载R包
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
library(ggcorrplot) 

normal1 <-Read10X(data.dir = "D:/beifen/QYK/scrna/oa/gse169454/normal/filtered/filtered/normal/normal1")
normal1 <- CreateSeuratObject(normal1,min.cells = 3, min.features = 200)
normal2 <-Read10X(data.dir = "D:/beifen/QYK/scrna/oa/gse169454/normal/filtered/filtered/normal/normal2")
normal2 <- CreateSeuratObject(normal2,min.cells = 3, min.features = 200)
normal3 <-Read10X(data.dir = "D:/beifen/QYK/scrna/oa/gse169454/normal/filtered/filtered/normal/normal3")
normal3 <- CreateSeuratObject(normal3,min.cells = 3, min.features = 200)
normal1$replicate <- sample(c("normal1"), size = ncol(normal1), replace = TRUE)
normal2$replicate <- sample(c("normal2"), size = ncol(normal2), replace = TRUE)
normal3$replicate <- sample(c("normal3"), size = ncol(normal3), replace = TRUE)
OA1 <-Read10X(data.dir = "D:/beifen/QYK/scrna/oa/gse169454/OA/filter/OA1")
OA1 <- CreateSeuratObject(OA1,min.cells = 3, min.features = 200)
OA2 <-Read10X(data.dir = "D:/beifen/QYK/scrna/oa/gse169454/OA/filter/OA2")
OA2 <- CreateSeuratObject(OA2,min.cells = 3, min.features = 200)
OA3 <-Read10X(data.dir = "D:/beifen/QYK/scrna/oa/gse169454/OA/filter/OA3")
OA3 <- CreateSeuratObject(OA3,min.cells = 3, min.features = 200)
OA4 <-Read10X(data.dir = "D:/beifen/QYK/scrna/oa/gse169454/OA/filter/OA4")
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


new.cluster.ids <- c("EC", "PreHTC", "SPP1+C", "PreHTC", "FC", "ProC",
                     "HTC", "RegC/HomoC", "RBC","undecidedC")
names(new.cluster.ids) <- levels(immune.combined)
immune.combined <- RenameIdents(immune.combined, new.cluster.ids)

DefaultAssay(immune.combined) <- "integrated"
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 30, verbose = FALSE)
ElbowPlot(immune.combined)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:15)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:15)
immune.combined <- FindClusters(immune.combined, resolution = 0.40)
table(immune.combined$orig.ident,immune.combined$replicate)

cc.genes.updated.2019
immune.combined <- CellCycleScoring(
  object = immune.combined,
  g2m.features = cc.genes$g2m.genes,
  s.features = cc.genes$s.genes)
VlnPlot(immune.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB","G2M.Score","S.Score"), ncol = 6)#+scale_color_npg()
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
 
plot=plot+scale_color_aaas()  +
  theme_bw()+theme(panel.grid=element_blank(),legend.title=element_blank(),legend.text = element_text(color="black", size = 10, face = "bold"))


immune.combined1 <- GetAssayData(immune.combined, assay = 'RNA',slot ='counts')
cell_metadata <- immune.combined@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(immune.combined1))
rownames(gene_annotation) <- rownames(immune.combined1)


cds <- new_cell_data_set(immune.combined1,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds@colData$replicate<- immune.combined$replicate
cds@colData$group<- immune.combined$group


#NormalizeData+ScaleData+RunPCA
cds <- preprocess_cds(cds, num_dim = 50)     #preprocess_cds函数相当于seurat中NormalizeData+ScaleData+RunPCA
plot_pc_variance_explained(cds)    #像seurat一样展示pc数
cds <- align_cds(cds, alignment_group = "replicate")


cds<- reduce_dimension(cds,preprocess_method = "PCA") #preprocess_method默认是PCA



#cds <- reduce_dimension(cds, reduction_method="tSNE")
#plot_cells(cds, reduction_method="tSNE", color_cells_by="seurat_clusters")




cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(immune.combined, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed  
cds@colData$integrated_snn_res.0.4 <- immune.combined$integrated_snn_res.0.4


plot_cells(cds, reduction_method="UMAP", color_cells_by="integrated_snn_res.0.4") 
cds <- cluster_cells(cds) 
cds <- learn_graph(cds)   


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

cds = order_cells(cds)  #
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=1.5,   
           group_label_size=4,cell_size=1.5)


Track_genes <- graph_test(cds, neighbor_graph="principal_graph", cores=6)
write.csv(Track_genes, "Trajectory_genes_all.csv", row.names = F)
Track_genes <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value < 1e-3)


Track_genes_sig <- Track_genes %>% top_n(n=10, morans_I) %>%
  pull(gene_short_name) %>% as.character()

plot_genes_in_pseudotime(cds[Track_genes_sig,], color_cells_by="integrated_snn_res.0.4", 
                         min_expr=0.5, ncol = 2)
p3<-plot_genes_in_pseudotime(cds[c("SQSTM1","COL2A1","NQO1","ACAN","MMP3","ADAMTS5","MMP13","COL10A1"),], color_cells_by="integrated_snn_res.0.4", 
                             min_expr=0.5, ncol = 2)
?plot_genes_in_pseudotime()
p3 <-p3+scale_color_npg()
p3
p3 <- p3+scale_color_manual(values = c("#DC143C","#0000FF","#20B2AA","#FFA500","#9370DB","#98FB98","#F08080","#1E90FF","#7CFC00","#FFFF00"))
p3

plot_cells(cds, genes=c("SQSTM1","COL2A1","NQO1","ACAN","MMP3","ADAMTS5","MMP13","COL10A1"), show_trajectory_graph=FALSE,
           label_cell_groups=FALSE,  label_leaves=FALSE)

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
?irGSEA.density.scatterplot()

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


##cellchat 
##cellchat
data.input <- GetAssayData(immune.combined_oa, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(immune.combined_oa)
head(labels)
identity <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels
cellchat <- createCellChat( data.input)
?createCellChat()
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
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
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat) # subset the expression data of signaling genes for saving computation cost
future::plan("multiprocess", workers = 4) # do parallel
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
?mergeCellChat()
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
