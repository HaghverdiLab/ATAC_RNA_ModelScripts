.libPaths("/home/sannald/projects/thesis/babel/.guix-profile/site-library")

library(Seurat)
library(dplyr)
library(ggplot2)
library(SeuratData)
library(SeuratDisk)
library(scater)
library(RColorBrewer)
library(qlcMatrix)

set.seed(42)

# Reference
reference <- LoadH5Seurat("/fast/AG_Haghverdi/Siddharth_Annaldasula/data/PBMC_reference/pbmc_multimodal.h5seurat")

# Original RNA
counts <- Read10X_h5("/fast/AG_Haghverdi/Siddharth_Annaldasula/data/PBMC_Healthy_10X_10k/pbmc_granulocyte_sorted_10k_filtered_feature_bc_matrix.h5")
pbmc10k_ref <- CreateSeuratObject(counts = counts$`Gene Expression`, project = "pbmc10k10XRNA")

pbmc10k_ref[["percent.mt"]] <- PercentageFeatureSet(pbmc10k_ref, pattern = "^MT-")
VlnPlot(pbmc10k_ref, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(pbmc10k_ref, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc10k_ref, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

pbmc10k_ref <- subset(pbmc10k_ref, subset = nCount_RNA > 1000 & nCount_RNA < 25000)

pbmc10k_ref <- SCTransform(pbmc10k_ref, verbose = FALSE)
pbmc10k_ref <- RunPCA(pbmc10k_ref)
ElbowPlot(pbmc10k)


DefaultAssay(pbmc10k_ref) <- "SCT"
pbmc.anchors <- FindTransferAnchors(reference = reference, query = pbmc10k_ref, dims = 1:50, recompute.residuals = FALSE, normalization.method = "SCT", reference.reduction = "spca")
#predictions <- TransferData(anchorset = pbmc.anchors, refdata = reference$celltype.l2, weight.reduction = pbmc10k_ref[['pca']], dims = 1:50)
#pbmc10k_ref <- AddMetaData(pbmc10k_ref, metadata = predictions)

reference <- RunUMAP(reference, reduction = "pca", dims = 1:50, return.model = TRUE)
pbmc10k_ref <- MapQuery(anchorset = pbmc.anchors, reference = reference, query = pbmc10k_ref,
                    refdata = list(celltype = "celltype.l2"), reference.reduction = "pca", reduction.model = "umap")

Idents(pbmc10k_ref) <- "predicted.celltype"
levels(pbmc10k_ref) <- c("CD4 Naive", "CD4 TCM", "CD4 CTL", "CD4 TEM", "CD4 Proliferating",
                  "CD8 Naive", "dnT",
                  "CD8 TEM", "CD8 TCM", "CD8 Proliferating", "MAIT", "NK", "NK_CD56bright",
                  "NK Proliferating", "gdT",
                  "Treg", "B naive", "B intermediate", "B memory", "Plasmablast",
                  "CD14 Mono", "CD16 Mono",
                  "cDC1", "cDC2", "pDC", "HSPC", "Eryth", "ASDC", "ILC", "Platelet")

pbmc10k_ref <- RunUMAP(pbmc10k_ref, dims = 1:50)
DimPlot(pbmc10k_ref, label = TRUE, repel = TRUE, reduction = "umap") + NoLegend()

p1 <- DimPlot(reference, reduction = "umap", group.by = "celltype.l2", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations") 
p2 <- DimPlot(pbmc10k_ref, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Query transferred labels")
p1 + p2

celltypes <- pbmc10k_ref$predicted.celltype
write.csv(celltypes,"/fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_human/pbmc10k_celltypes_new.csv", quote=FALSE)
celltypes <- read.csv("/fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_human/pbmc10k_celltypes_new.csv")
pbmc10k_ref$predicted.celltype <- celltypes$x
remove_celltypes <- c("Platelet","CD4 Proliferating","CD8 Proliferating","ILC")
pbmc10k_ref <- subset(pbmc10k_ref, cells = names(pbmc10k_ref$predicted.celltype[!pbmc10k_ref$predicted.celltype %in% remove_celltypes]))

#### Distribution of Data Celltypes before and after Balancing

pbmc10k_ref_celltypesl1_counts <- data.frame(name = names(table(celltypel2_l1[pbmc10k_ref$predicted.celltype])), count = as.numeric(table(celltypel2_l1[pbmc10k_ref$predicted.celltype])))
ggplot(pbmc10k_ref_celltypesl1_counts, aes(x = name, y = count)) + 
  ggtitle("Distribution of Celltypes, PBMC Data") +
  xlab("PBMC Celltypes") +
  ylab("Cell Count") + 
  geom_col(aes(reorder(name, -count))) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5), legend.position="none") 


pbmc10k_ref_celltypes_counts <- data.frame(name = names(table(pbmc10k_ref$predicted.celltype)), count = as.numeric(table(pbmc10k_ref$predicted.celltype)))
svg("/fast/AG_Haghverdi/Siddharth_Annaldasula/Section2_PBMC/PBMC_Data_Celltypes.svg", width = 7, height = 7)
ggplot(pbmc10k_ref_celltypes_counts, aes(x = name, y = count)) + 
  ggtitle("Distribution of Celltypes, PBMC Data") +
  xlab(NULL) +
  ylab("Cell Count") + 
  geom_col(aes(reorder(name, -count), fill = levels(as.factor(pbmc10k_ref$predicted.celltype))),color = "black") + 
  theme(text = element_text(size=18), axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5), legend.position="none") 
dev.off()

pbmc10k_ref_test_celltypes_counts <- pbmc10k_ref_celltypes_counts
pbmc10k_ref_test_celltypes_counts$count <- sapply(pbmc10k_ref_celltypes_counts$count, function(x) ifelse(x > 200, 200, x))
svg("/fast/AG_Haghverdi/Siddharth_Annaldasula/Section2_PBMC/PBMC_Data_Balanced_Celltypes.svg", width = 7, height = 7)
ggplot(pbmc10k_ref_test_celltypes_counts, aes(x = name, y = count)) + 
  ggtitle("Distribution of Celltypes after Balancing, PBMC Data") + 
  xlab(NULL) +
  ylab("Cell Count") + 
  geom_col(aes(reorder(name, -count), fill = levels(as.factor(pbmc10k_ref$predicted.celltype))),color = "black") +
  theme(text = element_text(size=18),axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position="none")
dev.off()

####
Convert("/fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_human/pbmc10k_balanced200new_filter/atac_rna_test_preds.h5ad", dest = "h5seurat", overwrite = TRUE)
pbmc10k_test <- LoadH5Seurat("/fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_human/pbmc10k_balanced200new_filter/atac_rna_test_preds.h5seurat", assays = "RNA")
pbmc10k_test_cells <- read.table(
  file = "/fast/AG_Haghverdi/Siddharth_Annaldasula/projects/atac_rna_enc_human/pbmc10k_balanced200new_filter/atac_rna_test_cells.txt",
  col.names = c("cells")
)
pbmc10k_test <- RenameCells(pbmc10k_test, new.names = pbmc10k_test_cells$cells)

DefaultAssay(pbmc10k_test) <- 'RNA'
nCount = colSums(x = pbmc10k_test, slot = "counts")  # nCount_RNA
nFeature = colSums(x = GetAssayData(object = pbmc10k_test, slot = "counts") > 0)  # nFeatureRNA
pbmc10k_test@meta.data["nCount_RNA"] = colSums(x = pbmc10k_test, slot = "counts") 
pbmc10k_test@meta.data["nFeature_RNA"] = colSums(x = GetAssayData(object = pbmc10k_test, slot = "counts") > 0)
pbmc10k_test[["percent.mt"]] <- PercentageFeatureSet(pbmc10k_test, pattern = "^MT-")
pbmc10k_test <- SCTransform(pbmc10k_test, verbose = FALSE)
pbmc10k_test <- RunPCA(pbmc10k_test)

pbmc10kref_10ktest.anchors <- FindTransferAnchors(reference = pbmc10k_ref, query = pbmc10k_test, dims = 1:50, normalization.method = "SCT",reference.reduction = "pca")
#predictions_10ktest <- TransferData(anchorset = pbmc10kref_10ktest.anchors, refdata = pbmc10k_ref$predicted.id, dims = 1:50)
#pbmc10k_test <- AddMetaData(pbmc10k_test, metadata = predictions_10ktest)
pbmc10k_ref <- RunUMAP(pbmc10k_ref, reduction = "pca", dims = 1:50, return.model = TRUE)
pbmc10k_test <- MapQuery(anchorset = pbmc10kref_10ktest.anchors, reference = pbmc10k_ref, query = pbmc10k_test,
                         refdata = list(celltype = "predicted.celltype"), reference.reduction = "pca", reduction.model = "umap")
pbmc10k_ref_test_cells <- subset(pbmc10k_ref, cells = pbmc10k_test_cells$cells)

DimPlot(so_Adult_test_cells, reduction = "umap.RNA", group.by = "celltype", label = TRUE, label.size = 5,repel = TRUE) + 
  NoLegend() + ggtitle("Adult Test Data, Truth Labels")  + xlim(-9, 11) + ylim(-13, 12) +
  scale_colour_discrete(drop=TRUE, limits = levels(so$celltype))

svg("/fast/AG_Haghverdi/Siddharth_Annaldasula/Section2_PBMC/Human_PBMC_Original_Truth.svg", width = 7, height = 7)
DimPlot(pbmc10k_ref, reduction = "umap", group.by = "predicted.celltype", label = TRUE, label.size = 4, repel = TRUE) + 
  NoLegend() + ggtitle("PBMC, Truth Labels") + xlim(-9, 11) + ylim(-12, 13) +
  scale_colour_discrete(drop=TRUE, limits = levels(as.factor(pbmc10k_ref$predicted.celltype)))
dev.off()

svg("/fast/AG_Haghverdi/Siddharth_Annaldasula/Section2_PBMC/Human_PBMC_Test_Truth.svg", width = 7, height = 7)
DimPlot(pbmc10k_ref_test_cells, reduction = "umap", group.by = "predicted.celltype", label = TRUE, label.size = 4, repel = TRUE) + 
  NoLegend() + ggtitle("PBMC Test Data, Truth Labels") + xlim(-9, 11) + ylim(-12, 13) +
  scale_colour_discrete(drop=TRUE, limits = levels(as.factor(pbmc10k_ref$predicted.celltype)))
dev.off()

svg("/fast/AG_Haghverdi/Siddharth_Annaldasula/Section2_PBMC/Human_PBMC_Inferred_Transferred.svg", width = 7, height = 7)
DimPlot(pbmc10k_test, reduction = "ref.umap", group.by = "predicted.celltype", label = TRUE, label.size = 4, repel = TRUE) + 
  NoLegend() + ggtitle("PBMC Inferred Test Data, Transferred Labels")  + xlim(-9, 11) + ylim(-12, 13)
dev.off()
p1 + p2 + p3




common_cells_pbmc10k <- intersect(names(pbmc10k_ref_test_cells$predicted.celltype),names(pbmc10k_test$predicted.celltype))
pbmc10k_ref_test_cells <- subset(pbmc10k_ref_test_cells, cells = common_cells_pbmc10k)
pbmc10k_test <- subset(pbmc10k_test, cells = common_cells_pbmc10k)
pbmc10k_test_celltypes <- data.frame(actual = pbmc10k_ref_test_cells$predicted.celltype, predicted = pbmc10k_test$predicted.celltype)

common_gene_names <- intersect(row.names(pbmc10k_ref_test_cells@assays$RNA@counts),row.names(pbmc10k_test@assays$RNA@counts))
pbmc10k_ref_test_count_mtx <- pbmc10k_ref_test_cells@assays$RNA@counts[common_gene_names,]
pbmc1k0_test_count_mtx <- pbmc10k_test@assays$RNA@counts[common_gene_names,]
test <- corSparse(pbmc10k_ref_test_count_mtx,pbmc1k0_test_count_mtx)
test <- cor(log1p(as.matrix(pbmc10k_ref_test_count_mtx)), log1p(as.matrix(pbmc1k0_test_count_mtx)), method = 'pearson')
pbmc10k_test_celltypes$cor <- diag(test)

pbmc10k_test_celltypes$same <- pbmc10k_test_celltypes$actual == pbmc10k_test_celltypes$predicted
print(pbmc10k_test_celltypes %>% group_by(actual) %>% summarise(truth = sum(same), cells=n(), r = mean(cor)),n=30)  

length(which(pbmc10k_test_celltypes$same == "TRUE"))/length(common_cells_pbmc10k)

test <- data.frame(l2 = reference$celltype.l2, l1 = reference$celltype.l1)
unique(test)

celltypel2_l1 <- c(
  'B naive' = 'B',
  'B intermediate' = 'B',
  'B memory' = 'B',
  'Plasmablast' = 'B',
  'CD14 Mono' = 'Mono',
  'CD16 Mono' = 'Mono',
  'CD4 Proliferating' = 'CD4 T',
  'CD4 CTL' = 'CD4 T',
  'CD4 Naive' = 'CD4 T',
  'CD4 TCM' = 'CD4 T',
  'CD4 TEM' = 'CD4 T',
  'Treg' = 'CD4 T',
  'CD8 Proliferating' = 'CD8 T',
  'CD8 Naive' = 'CD8 T',
  'CD8 TCM' = 'CD8 T',
  'CD8 TEM' = 'CD8 T',
  'ASDC' = 'DC',
  'cDC2' = 'DC',
  'cDC1' = 'DC',
  'pDC' = 'DC',
  'NK' = 'NK',
  'NK_CD56bright' = 'NK',
  'NK Proliferating' = 'NK',
  'dnT' = 'other T',
  'gdT' = 'other T',
  'MAIT' = 'other T',
  'Eryth' = 'other',
  'HSPC' = 'other',
  'Platelet' = 'other',
  'ILC' = 'other',
  'Doublet' = 'other'
)

pbmc10k_test_celltypes$actual_l1cat <- celltypel2_l1[pbmc10k_test_celltypes$actual]
pbmc10k_test_celltypes$predicted_l1cat <- celltypel2_l1[pbmc10k_test_celltypes$predicted]
pbmc10k_test_celltypes$same_l1cat <- pbmc10k_test_celltypes$actual_l1cat == pbmc10k_test_celltypes$predicted_l1cat
length(which(pbmc10k_test_celltypes$same_l1cat == "TRUE"))/length(common_cells_pbmc10k)

print(pbmc10k_test_celltypes[order(pbmc10k_test_celltypes$actual_l1cat),] %>% group_by(actual) %>% summarise(truth = sum(same), cells=n(),r = mean(cor)),n=30)
print(pbmc10k_test_celltypes %>% group_by(actual_l1cat) %>% summarise(truth = sum(same_l1cat), cells=n(),r = mean(cor)),n=30)  


test <- pbmc10k_test_celltypes %>%
  group_by(predicted) %>%
  summarise(cells=n()) 
print(test, n = 30)





pbmc_euc_dist <- sqrt(rowSums((pbmc10k_ref_test_cells@reductions$umap@cell.embeddings - pbmc10k_test@reductions$ref.umap@cell.embeddings)^2))
pbmc_euc_dist_df <- data.frame(eucdist = pbmc_euc_dist)
sum(sqrt(rowSums((pbmc10k_ref_test_cells@reductions$umap@cell.embeddings - pbmc10k_test@reductions$ref.umap@cell.embeddings)^2)))
sum(sqrt(rowSums((pbmc10k_ref_test_cells@reductions$umap@cell.embeddings - pbmc10k_test@reductions$ref.umap@cell.embeddings)^2)))/length(pbmc10k_test$predicted.celltype)
svg("/fast/AG_Haghverdi/Siddharth_Annaldasula/Section2_PBMC/Human_PBMC_EucDist.svg", width = 7, height = 7)
ggplot(pbmc_euc_dist_df, aes(eucdist)) + geom_histogram(binwidth = 0.5,center=0.25,fill = 'sky blue',color = 'black') +
  ggtitle("PBMC, Ref. UMAP Euclidean Distance between Original and Inferred") +
  xlab("Euclidean Distance") + ylab("Count")
dev.off() 
pbmc_euc_dist_df$celltype <- pbmc10k_ref_test_cells$predicted.id
print(pbmc_euc_dist_df %>% group_by(celltype) %>% summarise(euc = sum(eucdist)/n(), cells=n()),n=30) 
