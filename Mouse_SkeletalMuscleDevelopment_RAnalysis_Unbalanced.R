.libPaths("/home/sannald/projects/thesis/babel/.guix-profile/site-library")

library(Seurat)
library(dplyr)
library(ggplot2)
library(SeuratData)
library(SeuratDisk)
library(scater)

##### SO - Reference Data

so <- readRDS("/fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_mouse/muscle_labels/so.rds")
DefaultAssay(so) <- "integratedRNA"
so <- RunUMAP(so, reduction = "pca", reduction.name = "umap.RNA", dims = 1:30, return.model = TRUE)

svg("/fast/AG_Haghverdi/Siddharth_Annaldasula/Section2_Myofiber/Myofiber_GroupedCelltypes.svg", width = 7, height = 7)
DimPlot(so, reduction = "umap.RNA", group.by = "celltype", label = TRUE, label.size = 5, repel = TRUE) + 
  NoLegend() + ggtitle("Skeletal Muscle Development, by Celltypes") + xlim(-15, 15) + ylim(-14, 12) +
  scale_colour_discrete(drop=TRUE, limits = levels(so$celltype))
dev.off()

svg("/fast/AG_Haghverdi/Siddharth_Annaldasula/Section2_Myofiber/Myofiber_GroupedTimepoints.svg", width = 7, height = 7)
DimPlot(so, reduction = "umap.RNA", group.by = "orig.ident", label = TRUE, label.size = 5, repel = TRUE) + 
  ggtitle("Skeletal Muscle Development, by Timepoints") + xlim(-15, 15) + ylim(-14, 12) 
dev.off()

##### E14
Convert("/fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_mouse/muscle_e14/muscle_e14_unbalanced_integrated_celltypes/atac_rna_test_preds.h5ad", dest = "h5seurat", overwrite = TRUE)
musclee14_test <- LoadH5Seurat("/fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_mouse/muscle_e14/muscle_e14_unbalanced_integrated_celltypes/atac_rna_test_preds.h5seurat", assays = "RNA")

musclee14_test_cells <- read.table(
  file = "/fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_mouse/muscle_e14/muscle_e14_unbalanced_integrated_celltypes/atac_rna_test_cells.txt",
  col.names = c("cells")
)
musclee14_test <- RenameCells(musclee14_test, new.names = musclee14_test_cells$cells)

DefaultAssay(musclee14_test) <- 'RNA'
musclee14_test@meta.data["nCount_RNA"] = colSums(x = musclee14_test, slot = "counts") 
musclee14_test@meta.data["nFeature_RNA"] = colSums(x = GetAssayData(object = musclee14_test, slot = "counts") > 0)
musclee14_test[["percent.mt"]] <- PercentageFeatureSet(musclee14_test, pattern = "^mt-")
VlnPlot(musclee14_test, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(musclee14_test, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(musclee14_test, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
musclee14_test <- subset(musclee14_test,subset = nCount_RNA < 25000 & nCount_RNA > 400 & percent.mt < 20)

musclee14_test <- SCTransform(object = musclee14_test, method = 'glmGamPoi', variable.features.n = 3000,verbose = FALSE)
musclee14_test <- RunPCA(musclee14_test, npcs = 30, verbose = FALSE)
musclee14_test <- RunUMAP(musclee14_test, reduction = "pca", reduction.name = "umap.RNA",dims = 1:30)

musclee14_test.anchors <- FindTransferAnchors(reference = so, query = musclee14_test, normalization.method = "SCT", reference.reduction = "pca")
#predictionse14 <- TransferData(anchorset = musclee14_test.anchors, refdata = so_e14_test_cells$celltype, dims = 1:30)
#musclee14_test <- AddMetaData(musclee14_test, metadata = predictionse14)
musclee14_test <- MapQuery(anchorset = musclee14_test.anchors, reference = so, query = musclee14_test,
                           refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap.RNA")

musclee14_test_cells$cells <- paste("E14", musclee14_test_cells$cells, sep="_")

so_e14 <- subset(so, cells = names(so$orig.ident[so$orig.ident == 'E14']))
so_e14_test_cells <- subset(so, cells = musclee14_test_cells$cells)

svg("/fast/AG_Haghverdi/Siddharth_Annaldasula/Section2_Myofiber_Unbalanced/Mouse_E14_Original_Truth.svg", width = 7, height = 7)
DimPlot(so_e14, reduction = "umap.RNA", group.by = "celltype", label = TRUE, label.size = 5,repel = TRUE) + 
  NoLegend() + ggtitle("E14 All Data, Truth Labels")  + xlim(-15, 15) + ylim(-14, 12) +
  scale_colour_discrete(drop=TRUE, limits = levels(so$celltype))
dev.off()

svg("/fast/AG_Haghverdi/Siddharth_Annaldasula/Section2_Myofiber_Unbalanced/Mouse_E14_Test_Truth.svg", width = 7, height = 7)
DimPlot(so_e14_test_cells, reduction = "umap.RNA", group.by = "celltype", label = TRUE, label.size = 5,repel = TRUE) + 
  NoLegend() + ggtitle("E14 Test Data, Truth Labels")  + xlim(-15, 15) + ylim(-14, 12) +
  scale_colour_discrete(drop=TRUE, limits = levels(so$celltype))
dev.off()

svg("/fast/AG_Haghverdi/Siddharth_Annaldasula/Section2_Myofiber_Unbalanced/Mouse_E14_Inferred_Transferred.svg", width = 7, height = 7)
DimPlot(musclee14_test, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, label.size = 5, repel = TRUE) + 
  NoLegend() + ggtitle("E14 Inferred Test Data, Transferred Labels") + xlim(-15, 15) + ylim(-14, 12) +
  scale_colour_discrete(drop=TRUE, limits = levels(so$celltype))
dev.off() 

musclee14_test_celltypes <- data.frame(actual = so_e14_test_cells$celltype, predicted = musclee14_test$predicted.celltype)
musclee14_test_celltypes$same <- musclee14_test_celltypes$actual == musclee14_test_celltypes$predicted
length(which(musclee14_test_celltypes$same == "TRUE"))/length(so_e14_test_cells$celltype)
common_gene_names <- intersect(row.names(so_e14_test_cells@assays$RNA@counts),row.names(musclee14_test@assays$RNA@counts))
e14_ref_test_count_mtx <- so_e14_test_cells@assays$RNA@counts[common_gene_names,]
e14_test_count_mtx <- musclee14_test@assays$RNA@counts[common_gene_names,]
musclee14_test_celltypes$cor <- diag(cor(log1p(as.matrix(e14_ref_test_count_mtx)), log1p(as.matrix(e14_test_count_mtx)), method = 'pearson'))

test <- musclee14_test_celltypes %>%
  group_by(actual) %>%
  summarise(truth = sum(same), cells=n(), r = mean(cor)) 
print(test, n = 30)


##### E18
Convert("/fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_mouse/muscle_e18/muscle_e18_unbalanced_integratedcelltypes/atac_rna_test_preds.h5ad", dest = "h5seurat", overwrite = TRUE)
muscleE18_test <- LoadH5Seurat("/fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_mouse/muscle_e18/muscle_e18_unbalanced_integratedcelltypes/atac_rna_test_preds.h5seurat", assays = "RNA")

muscleE18_test_cells <- read.table(
  file = "/fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_mouse/muscle_e18/muscle_e18_unbalanced_integratedcelltypes/atac_rna_test_cells.txt",
  col.names = c("cells")
)
muscleE18_test <- RenameCells(muscleE18_test, new.names = muscleE18_test_cells$cells)

DefaultAssay(muscleE18_test) <- 'RNA'
muscleE18_test@meta.data["nCount_RNA"] = colSums(x = muscleE18_test, slot = "counts") 
muscleE18_test@meta.data["nFeature_RNA"] = colSums(x = GetAssayData(object = muscleE18_test, slot = "counts") > 0)
muscleE18_test[["percent.mt"]] <- PercentageFeatureSet(muscleE18_test, pattern = "^mt-")
muscleE18_test <- subset(muscleE18_test,subset = nCount_RNA < 25000 & nCount_RNA > 400 & percent.mt < 20)

muscleE18_test <- SCTransform(object = muscleE18_test, method = 'glmGamPoi', variable.features.n = 3000, vars.to.regress = "percent.mt", verbose = FALSE)
muscleE18_test <- RunPCA(muscleE18_test, npcs = 30, verbose = FALSE)
muscleE18_test <- RunUMAP(muscleE18_test, reduction = "pca", reduction.name = "umap.RNA",dims = 1:30)

DefaultAssay(muscleE18_test) <- "SCT"

muscleE18_test.anchors <- FindTransferAnchors(reference = so, query = muscleE18_test, normalization.method = "SCT", reference.reduction = "pca")
#predictionsE18 <- TransferData(anchorset = muscleE18_test.anchors, refdata = so$celltype)
#muscleE18_test <- AddMetaData(muscleE18_test, metadata = predictionsE18)
muscleE18_test <- MapQuery(anchorset = muscleE18_test.anchors, reference = so, query = muscleE18_test,
                           refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap.RNA")

muscleE18_test_cells$cells <- paste("E18", muscleE18_test_cells$cells, sep="_")
so_E18_test_cells <- subset(so, cells = muscleE18_test_cells$cells)


so_E18 <- subset(so, cells = names(so$orig.ident[so$orig.ident == 'E18']))

svg("/fast/AG_Haghverdi/Siddharth_Annaldasula/Section2_Myofiber_Unbalanced/Mouse_E18_Original_Truth.svg", width = 7, height = 7)
DimPlot(so_E18, reduction = "umap.RNA", group.by = "celltype", label = TRUE, label.size = 5,repel = TRUE) + 
  NoLegend() + ggtitle("E18 All Data, Truth Labels")  + xlim(-15, 15) + ylim(-14, 12) +
  scale_colour_discrete(drop=TRUE, limits = levels(so$celltype))
dev.off()

svg("/fast/AG_Haghverdi/Siddharth_Annaldasula/Section2_Myofiber_Unbalanced/Mouse_E18_Test_Truth.svg", width = 7, height = 7)
DimPlot(so_E18_test_cells, reduction = "umap.RNA", group.by = "celltype", label = TRUE, label.size = 5,repel = TRUE) + 
  NoLegend() + ggtitle("E18 Test Data, Truth Labels")  + xlim(-15, 15) + ylim(-14, 12) +
  scale_colour_discrete(drop=TRUE, limits = levels(so$celltype))
dev.off()

svg("/fast/AG_Haghverdi/Siddharth_Annaldasula/Section2_Myofiber_Unbalanced/Mouse_E18_Inferred_Transferred.svg", width = 7, height = 7)
DimPlot(muscleE18_test, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, label.size = 5, repel = TRUE) + 
  NoLegend() + ggtitle("E18 Inferred Test Data, Transferred Labels") + xlim(-15, 15) + ylim(-14, 12) +
  scale_colour_discrete(drop=TRUE, limits = levels(so$celltype))
dev.off() 

muscleE18_test_celltypes <- data.frame(actual = so_E18_test_cells$celltype, predicted = muscleE18_test$predicted.celltype)
muscleE18_test_celltypes$same <- muscleE18_test_celltypes$actual == muscleE18_test_celltypes$predicted
common_gene_names <- intersect(row.names(so_E18_test_cells@assays$RNA@counts),row.names(muscleE18_test@assays$RNA@counts))
E18_ref_test_count_mtx <- so_E18_test_cells@assays$RNA@counts[common_gene_names,]
E18_test_count_mtx <- muscleE18_test@assays$RNA@counts[common_gene_names,]
muscleE18_test_celltypes$cor <- diag(cor(log1p(as.matrix(E18_ref_test_count_mtx)), log1p(as.matrix(E18_test_count_mtx)), method = 'pearson'))

test <- muscleE18_test_celltypes %>%
  group_by(actual) %>%
  summarise(truth = sum(same), cells=n(), r = mean(cor)) 
print(test, n = 30)

##### P5
Convert("/fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_mouse/muscle_p5/muscle_p5_unbalanced_integratedcelltypes/atac_rna_test_preds.h5ad", dest = "h5seurat", overwrite = TRUE)
muscleP5_test <- LoadH5Seurat("/fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_mouse/muscle_p5/muscle_p5_unbalanced_integratedcelltypes/atac_rna_test_preds.h5seurat", assays = "RNA")

muscleP5_test_cells <- read.table(
  file = "/fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_mouse/muscle_p5/muscle_p5_unbalanced_integratedcelltypes/atac_rna_test_cells.txt",
  col.names = c("cells")
)
muscleP5_test <- RenameCells(muscleP5_test, new.names = muscleP5_test_cells$cells)

DefaultAssay(muscleP5_test) <- 'RNA'
muscleP5_test@meta.data["nCount_RNA"] = colSums(x = muscleP5_test, slot = "counts") 
muscleP5_test@meta.data["nFeature_RNA"] = colSums(x = GetAssayData(object = muscleP5_test, slot = "counts") > 0)
muscleP5_test[["percent.mt"]] <- PercentageFeatureSet(muscleP5_test, pattern = "^mt-")
muscleP5_test <- subset(muscleP5_test,subset = nCount_RNA < 25000 & nCount_RNA > 400 & percent.mt < 20)

muscleP5_test <- SCTransform(object = muscleP5_test, method = 'glmGamPoi', variable.features.n = 3000, verbose = FALSE)
muscleP5_test <- RunPCA(muscleP5_test, npcs = 30, verbose = FALSE)
muscleP5_test <- RunUMAP(muscleP5_test, reduction = "pca", reduction.name = "umap.RNA",dims = 1:30)

DefaultAssay(muscleP5_test) <- "SCT"

muscleP5_test.anchors <- FindTransferAnchors(reference = so, query = muscleP5_test, normalization.method = "SCT", reference.reduction = "pca")
#predictionsP5 <- TransferData(anchorset = muscleP5_test.anchors, refdata = so$celltype)
#muscleP5_test <- AddMetaData(muscleP5_test, metadata = predictionsP5)
muscleP5_test <- MapQuery(anchorset = muscleP5_test.anchors, reference = so, query = muscleP5_test,
                          refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap.RNA")

#pbmc10k_ref_test_cells <- subset(pbmc10k_ref, cells = pbmc10k_test_cells$cells)

muscleP5_test_cells$cells <- paste("P5", muscleP5_test_cells$cells, sep="_")
so_P5_test_cells <- subset(so, cells = muscleP5_test_cells$cells)

so_P5 <- subset(so, cells = names(so$orig.ident[so$orig.ident == 'P5']))

svg("/fast/AG_Haghverdi/Siddharth_Annaldasula/Section2_Myofiber_Unbalanced/Mouse_P5_Original_Truth.svg", width = 7, height = 7)
DimPlot(so_P5, reduction = "umap.RNA", group.by = "celltype", label = TRUE, label.size = 5,repel = TRUE) + 
  NoLegend() + ggtitle("P5 All Data, Truth Labels")  + xlim(-15, 15) + ylim(-14, 12) +
  scale_colour_discrete(drop=TRUE, limits = levels(so$celltype))
dev.off()

svg("/fast/AG_Haghverdi/Siddharth_Annaldasula/Section2_Myofiber_Unbalanced/Mouse_P5_Test_Truth.svg", width = 7, height = 7)
DimPlot(so_P5_test_cells, reduction = "umap.RNA", group.by = "celltype", label = TRUE, label.size = 5,repel = TRUE) + 
  NoLegend() + ggtitle("P5 Test Data, Truth Labels")  + xlim(-15, 15) + ylim(-14, 12) +
  scale_colour_discrete(drop=TRUE, limits = levels(so$celltype))
dev.off()

svg("/fast/AG_Haghverdi/Siddharth_Annaldasula/Section2_Myofiber_Unbalanced/Mouse_P5_Inferred_Transferred.svg", width = 7, height = 7)
DimPlot(muscleP5_test, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, label.size = 5, repel = TRUE) + 
  NoLegend() + ggtitle("P5 Inferred Test Data, Transferred Labels") + xlim(-15, 15) + ylim(-14, 12) +
  scale_colour_discrete(drop=TRUE, limits = levels(so$celltype))
dev.off() 

muscleP5_test_celltypes <- data.frame(actual = so_P5_test_cells$celltype, predicted = muscleP5_test$predicted.celltype)
muscleP5_test_celltypes$same <- muscleP5_test_celltypes$actual == muscleP5_test_celltypes$predicted

common_gene_names <- intersect(row.names(so_P5_test_cells@assays$RNA@counts),row.names(muscleP5_test@assays$RNA@counts))
P5_ref_test_count_mtx <- so_P5_test_cells@assays$RNA@counts[common_gene_names,]
P5_test_count_mtx <- muscleP5_test@assays$RNA@counts[common_gene_names,]
muscleP5_test_celltypes$cor <- diag(cor(log1p(as.matrix(P5_ref_test_count_mtx)), log1p(as.matrix(P5_test_count_mtx)), method = 'pearson'))

test <- muscleP5_test_celltypes %>%
  group_by(actual) %>%
  summarise(truth = sum(same), cells=n(), r = mean(cor)) 
print(test, n = 30)


##### Adult
Convert("/fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_mouse/muscle_adult/muscle_adult_unbalanced_integratedcelltypes/atac_rna_test_preds.h5ad", dest = "h5seurat", overwrite = TRUE)
muscleAdult_test <- LoadH5Seurat("/fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_mouse/muscle_adult/muscle_adult_unbalanced_integratedcelltypes/atac_rna_test_preds.h5seurat", assays = "RNA")

muscleAdult_test_cells <- read.table(
  file = "/fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_mouse/muscle_adult/muscle_adult_unbalanced_integratedcelltypes/atac_rna_test_cells.txt",
  col.names = c("cells")
)
muscleAdult_test <- RenameCells(muscleAdult_test, new.names = muscleAdult_test_cells$cells)

DefaultAssay(muscleAdult_test) <- 'RNA'
muscleAdult_test@meta.data["nCount_RNA"] = colSums(x = muscleAdult_test, slot = "counts") 
muscleAdult_test@meta.data["nFeature_RNA"] = colSums(x = GetAssayData(object = muscleAdult_test, slot = "counts") > 0)
muscleAdult_test[["percent.mt"]] <- PercentageFeatureSet(muscleAdult_test, pattern = "^mt-")
muscleAdult_test <- subset(muscleAdult_test,subset = nCount_RNA < 25000 & nCount_RNA > 400 & percent.mt < 20)

muscleAdult_test <- SCTransform(object = muscleAdult_test, method = 'glmGamPoi', variable.features.n = 3000, verbose = FALSE)
muscleAdult_test <- RunPCA(muscleAdult_test, npcs = 30, verbose = FALSE)
muscleAdult_test <- RunUMAP(muscleAdult_test, reduction = "pca", reduction.name = "umap.RNA",dims = 1:30)

DefaultAssay(muscleAdult_test) <- "SCT"

muscleAdult_test.anchors <- FindTransferAnchors(reference = so, query = muscleAdult_test, normalization.method = "SCT", reference.reduction = "pca")
#predictionsAdult <- TransferData(anchorset = muscleAdult_test.anchors, refdata = so$celltype)
#muscleAdult_test <- AddMetaData(muscleAdult_test, metadata = predictions)
muscleAdult_test <- MapQuery(anchorset = muscleAdult_test.anchors, reference = so, query = muscleAdult_test,
                             refdata = list(celltype = "celltype"), reference.reduction = "pca", reduction.model = "umap.RNA")

pbmc10k_ref_test_cells <- subset(pbmc10k_ref, cells = pbmc10k_test_cells$cells)

muscleAdult_test_cells$cells <- paste("Adult", muscleAdult_test_cells$cells, sep="_")
so_Adult_test_cells <- subset(so, cells = muscleAdult_test_cells$cells)

so_Adult <- subset(so, cells = names(so$orig.ident[so$orig.ident == 'Adult']))

svg("/fast/AG_Haghverdi/Siddharth_Annaldasula/Section2_Myofiber_Unbalanced/Mouse_Adult_Original_Truth.svg", width = 7, height = 7)
DimPlot(so_Adult, reduction = "umap.RNA", group.by = "celltype", label = TRUE, label.size = 5,repel = TRUE) + 
  NoLegend() + ggtitle("Adult All Data, Truth Labels")  + xlim(-15, 15) + ylim(-14, 12) +
  scale_colour_discrete(drop=TRUE, limits = levels(so$celltype))
dev.off()

svg("/fast/AG_Haghverdi/Siddharth_Annaldasula/Section2_Myofiber_Unbalanced/Mouse_Adult_Test_Truth.svg", width = 7, height = 7)
DimPlot(so_Adult_test_cells, reduction = "umap.RNA", group.by = "celltype", label = TRUE, label.size = 5,repel = TRUE) + 
  NoLegend() + ggtitle("Adult Test Data, Truth Labels")  + xlim(-15, 15) + ylim(-14, 12) +
  scale_colour_discrete(drop=TRUE, limits = levels(so$celltype))
dev.off()

svg("/fast/AG_Haghverdi/Siddharth_Annaldasula/Section2_Myofiber_Unbalanced/Mouse_Adult_Inferred_Transferred.svg", width = 7, height = 7)
DimPlot(muscleAdult_test, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, label.size = 5, repel = TRUE) + 
  NoLegend() + ggtitle("Adult Inferred Test Data, Transferred Labels") + xlim(-15, 15) + ylim(-14, 12) +
  scale_colour_discrete(drop=TRUE, limits = levels(so$celltype))
dev.off() 

muscleAdult_test_celltypes <- data.frame(actual = so_Adult_test_cells$celltype, predicted = muscleAdult_test$predicted.celltype)
muscleAdult_test_celltypes$same <- muscleAdult_test_celltypes$actual == muscleAdult_test_celltypes$predicted
length(which(muscleAdult_test_celltypes$same == "TRUE"))/length(so_Adult_test_cells$celltype)

common_gene_names <- intersect(row.names(so_Adult_test_cells@assays$RNA@counts),row.names(muscleAdult_test@assays$RNA@counts))
Adult_ref_test_count_mtx <- so_Adult_test_cells@assays$RNA@counts[common_gene_names,]
Adult_test_count_mtx <- muscleAdult_test@assays$RNA@counts[common_gene_names,]
muscleAdult_test_celltypes$cor <- diag(cor(log1p(as.matrix(Adult_ref_test_count_mtx)), log1p(as.matrix(Adult_test_count_mtx)), method = 'pearson'))

test <- muscleAdult_test_celltypes %>%
  group_by(actual) %>%
  summarise(truth = sum(same), cells=n(), r = mean(cor)) 
print(test, n = 30)
