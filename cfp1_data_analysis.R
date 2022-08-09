library(Seurat)
library(harmony)


# CFP1, B6 and Old mouse single cell RNA-seq data analysis

#read merged gene expression data

seurat.obj <- readRDS("seurat.obj.add.mt.percent.unnormalized.sample.Rds")
seurat.obj$Sample_Names <- paste0(seurat.obj$batch, "_", seurat.obj$Samples)

##------------------------------------------------------------------##
##  select used samples
#------------------------------------------------------------------##
cfp1 <- subset(seurat.obj, Sample_Names %in% c("20200618_B6",  "20200618_CFP1", "20200618_Old", "20200415_CFP1"))
cfp1 <- subset(cfp1, subset = nFeature_RNA > 200 & percent.mt < 20)

cfp1$condition <- cfp1$Samples
cfp1$condition[cfp1$condition %in% c("B6", "WT2")] <- "control"
##-----------------------------------------------------------------##
# Data normalization
cfp1 <- NormalizeData(cfp1)
cfp1 <- FindVariableFeatures(cfp1)
cfp1 <- ScaleData(vars.to.regress = c("percent.mt"), verbose = FALSE) 

cfp1 <- RunHarmony(cfp1, "Sample_Names", plot_convergence = FALSE)


##---------------------- Clustering ------------------------------------##

cfp1 <- RunTSNE(cfp1, reduction = "harmony", dims = 1:20) 
cfp1 <- FindNeighbors(cfp1, reduction = "harmony", dims = 1:20)
cfp1 <- FindClusters(cfp1, resolution = 0.8)

DimPlot(cfp1, label = T, reduction = "tsne")

FeaturePlot(cfp1, features = c("Rorc", "Cd3d", "Tcf7","Cd4", "Klf4", "Ncr1", "Ccr6"),  reduction = "tsne") 
DimPlot(cfp1, label = T, reduction = "tsne") & NoLegend()

cfp1_filter <- subset(cfp1, seurat_clusters %in% c(1,4,5,10,12,13,15,18,19,20), invert=T)



#------------------ downstream analysis------------------------##

pdf("fg_a_dimplot.pdf")
DimPlot(cfp1_filter, label = T) & NoLegend()
DimPlot(cfp1_filter, label = T, reduction = "tsne") & NoLegend()
DimPlot(cfp1_filter, label = T, group.by = "celltype", reduction = "tsne") & NoLegend()
dev.off()

genes <- c("Rorc", "Cd3d", "Tcf7", "Klf4", "Ncr1", "Ccr6", "Cd4", "Gata3", "Tbx21", "Maf", "Il17a", "Il22", "Ifng", "Tnf")


pdf("fg_b_maker_gene_dimplot_tsne.pdf")
FeaturePlot(cfp1_filter, features = genes, cols = c('grey80', '#FF6347'), reduction = "tsne")  & NoAxes() & NoLegend()
dev.off()

genes2 <- readxl::read_xlsx("senescence genelist.xlsx")$Senescence
seurat.obj <- AddModuleScore(seurat.obj, features = list(genes2), name="senescence")

pdf("vlnplot.pdf")
VlnPlot(cfp1_filter, features = senescence1, group.by = "clusters")
dev.off()


pdf("fg_e_umap_tsne_group_conditon.pdf")
DimPlot(cfp1_filter,  group.by = "condition", reduction = "tsne")
dev.off()


df2 <- table(seurat.obj$condition, seurat.obj$celltype)
df2 <- reshape2::melt(df2)
names(df2)[1:2] <- c('sample','cluster')


ggplot(df2, aes(sample, value, fill=cluster)) + 
	geom_bar(position="fill", stat="identity"5) +
	labs(x="", y="Fraction of total cells") + theme_classic() + 
ggsave("fg_c_barplot.pdf")


##-----------------------------------------------------------##
## pseudotime trajectory  ##
##-----------------------------------------------------------##
library(monocle)
library(Seurat)

ilc3_obj <- readRDS("cluster_data.rdata")
#pheno data
data6 <- GetAssayData(ilc3_obj, slot = "counts")
cell.ids <- colnames(data6)
cell.info <- ilc3_obj@meta.data
pd6 <- new("AnnotatedDataFrame", data= data.frame(cell.info))
fData <- data.frame(gene_short_name = row.names(data6), row.names = row.names(data6))


mc_6_repeat <- newCellDataSet(as.matrix(data6), 
                              phenoData = pd6, 
                              featureData = fd6)

mc_6_repeat <- estimateSizeFactors(mc_6_repeat)
mc_6_repeat <- estimateDispersions(mc_6_repeat)


diff_test_res <- differentialGeneTest(mc_6_repeat,fullModelFormulaStr = "~seurat_clusters")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
mc_6_repeat <- setOrderingFilter(mc_6_repeat, ordering_genes)


mc_6_repeat <- orderCells(mc_6_repeat, root_state = 4)



pdf("ilc3_all_obj_pseudotime1.pdf")
plot_cell_trajectory(mc_6_repeat, color_by="celltype", show_branch_points = F)
plot_cell_trajectory(mc_6_repeat, color_by="Pseudotime", show_branch_points = F)
plot_cell_trajectory(mc_6_repeat, color_by="seurat_clusters2", show_branch_points = F)
dev.off()

pdf("ilc3_all_obj_pseudotime2.pdf", w=5.5, h=9)
plot_cell_trajectory(mc_6_repeat, color_by="celltype") +
  facet_wrap(~celltype, nrow = 3)
dev.off()


index <- c("Ccr6", "Cd4", "Il17a", "Il22", "Ifng", "Klf4", "Tcf7", "Rorc", "Cxxc1",
           "Ahr", "Gata3", "Runx3", "Tbx21", "Maf",  "Washc1" ,"Washc2", 
           "Washc3", "Washc4", "Washc5")

index <- c("Il17f", "Tbx21", "Setd2", "Csf2", "Setd1a", "Setd1b", "Ncr1", "GM-CSF")

my_genes <- row.names(subset(fData(mc_6_repeat), gene_short_name %in% index))
cds_subset <- mc_6_repeat[my_genes, ]


pdf("ilc3_obj_pseudotime_all.pdf", w=10, h=8)
plot_cell_trajectory(mc_6_repeat, color_by="seurat_clusters") +
  facet_wrap(~condition, nrow = 1)
dev.off()


pdf("genes_plot_trajectory2.pdf")
plot_trajectory(cds_subset, ncol = 3)
dev.off()




##---------------------------------------------------------------##
##---------------------------------------------------------------##
# B6 and Old mouse single cell RNA-seq data analysis


seurat.obj <- readRDS("seurat.obj.add.mt.percent.unnormalized.sample.Rds")
seurat.obj$Sample_Names <- paste0(seurat.obj$batch, "_", seurat.obj$Samples)
seurat.obj <- subset(seurat.obj, Sample_Names %in% c("20211102_WT2", "20200702_B6", "20200918_Old", "20200618_Old"))
seurat.obj <- subset(seurat.obj, subset = nFeature_RNA > 200 & percent.mt < 20)

seurat.obj <- NormalizeData(seurat.obj, verbose = FALSE)

seurat.obj <- FindVariableFeatures(seurat.obj) 
seurat.obj <- ScaleData(seurat.obj, vars.to.regress = c("percent.mt"), verbose = FALSE) %>%
seurat.obj <- RunPCA(seurat.obj, verbose = FALSE)

seurat.obj <- RunHarmony(seurat.obj, "Sample_Names", plot_convergence = FALSE)


seurat.obj <- RunTSNE(seurat.obj, reduction = "harmony", dims = 1:20) 
seurat.obj <- FindNeighbors(seurat.obj, reduction = "harmony", dims = 1:20) 
seurat.obj <- FindClusters(seurat.obj, resolution = 1.2) 


DimPlot(seurat.obj, group.by = "Sample_Names")
DimPlot(seurat.obj, label = T)
DimPlot(seurat.obj, label = T, reduction = "tsne")


old_mouse2 <- subset(seurat.obj, seurat_clusters %in% c(13, 14, 18,9, 16, 19, 10, 0,17), invert=T)

old_mouse2$condition <- old_mouse2$Samples
old_mouse2$condition[grep("B6|WT", old_mouse2$condition)] <- "control"
VlnPlot(old_mouse2, features = c("Cxxc1", "Klf4"), group.by = "condition")

#----- downsample -------#
df2 <- as.data.frame(table(downsample_dat$clusters, downsample_dat$condition))

ggplot(df2, aes(cluster, value, fill=sample)) + 
  geom_bar(position="fill",width = 0.85 , stat="identity", alpha=0.85) +
  labs(x="", y="Fraction of total cells") +  
  geom_hline(yintercept=0.5, linetype="dashed", color = "grey30") +
  theme_classic() + scale_fill_manual(values = c("#5F4B8BFF", "#E69A8DFF")) +
  theme(axis.text = element_text(color = "black", size=10),


ggsave("old_mouse_cell_fraction.pdf", w=8,h=5)


downsample_dat <- subset(old_mouse2, downsample = 5016)
table(downsample_dat$condition)

downsample_dat$clusters <- as.character(downsample_dat$clusters)
downsample_dat$clusters <- as.numeric(downsample_dat$clusters)

pdf("cluster_dimplot.pdf")
DimPlot(downsample_dat, label = T, reduction = "tsne", group.by = "clusters") 
dev.off()

pdf("cluster_dimplot_splitby_condition.pdf", h=5)
DimPlot(downsample_dat, label = T, reduction = "tsne", group.by = "clusters", split.by = "condition") 
dev.off()


pdf("cluster_dimplot_splitby_condition2.pdf")
FeaturePlot(downsample_dat, reduction = "tsne", features = c("Rorc", "Tcf7", "Cd4", "Klf4", "Ccr6", "Ncr1"), 
            cols = c("grey80", "red"), pt.size = 0.1, order = T) & NoLegend() & NoAxes()
dev.off()

pdf("cluster_dimplot_splitby_condition_split_condition.pdf", w=6, h=9)
FeaturePlot(downsample_dat, reduction = "tsne", features = c("Tcf7","Klf4" ,"Cd4", "Cdkn2a", "Cdkn1a"), 
            cols = c("grey80", "red"), pt.size = 0.05, split.by = "condition", order = F) & NoLegend() & NoAxes()
dev.off()










