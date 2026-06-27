
## ==============================================================================
## Project: Hanna A. et al. Circulation 2026; Cardiac Fibroblast Heterogeneity & Smad1 Integration
## ==============================================================================

## %% 0. LOAD PACKAGES
# Suppress package startup messages for cleaner execution
suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(SCP)       
  library(orthogene) 
  library(SCpubr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(patchwork)
})

# set seed
set.seed(42)

## %% 1. LOAD RAW UMI COUNT MATRICES & CREATE SEURAT OBJECTS
# Dataset 1: GSE132146 (Col1a1 GFP+ Cardiac Fibroblasts)
counts_7dpi  <- Read10X_h5("7dMI_counts.h5")
counts_14dpi <- Read10X_h5("14dMI_counts.h5")
counts_30dpi <- Read10X_h5("30dMI_counts.h5")

seurat_7dpi  <- CreateSeuratObject(counts = counts_7dpi, project = "GSE132146_7d", min.cells = 3, min.features = 200)
seurat_14dpi <- CreateSeuratObject(counts = counts_14dpi, project = "GSE132146_14d", min.cells = 3, min.features = 200)
seurat_30dpi <- CreateSeuratObject(counts = counts_30dpi, project = "GSE132146_30d", min.cells = 3, min.features = 200)
str(seurat_7dpi)
seurat_7dpi$condition  <- "7d-PO"
seurat_14dpi$condition <- "14d-PO"
seurat_30dpi$condition <- "30d-PO"
seurat_7dpi$study      <- "GSE132146"
seurat_14dpi$study     <- "GSE132146"
seurat_30dpi$study     <- "GSE132146"

## [Dataset 2: E-MTAB-7365] Pdgfra+/Sca1+/Cd31- Infarct Fibroblasts

# Load raw count matrix 
counts_7365 <- read.table("EMTAB7365_days3_7_counts.txt", header = TRUE, row.names = 1, sep = "\t")

# Create Seurat Object
seurat_7365 <- CreateSeuratObject(counts = counts_7365, project = "EMTAB7365", min.cells = 3, min.features = 200)

# Load the cell metadata table to identify 3d-PO and 7d-PO
meta_7365 <- read.table("Timepoint_table.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Seurat's AddMetaData requires the rownames of the metadata to match the cell barcodes (colnames of counts).
# Assuming the first column contains the cell barcodes if they aren't already set as rownames:
if(colnames(meta_7365)[1] != "row.names" && !identical(rownames(meta_7365), colnames(seurat_7365))) {
  rownames(meta_7365) <- meta_7365[, 1]
}

# Add the imported metadata to the Seurat object
seurat_7365 <- AddMetaData(seurat_7365, metadata = meta_7365)

# Standardize the study column for downstream batch correction
seurat_7365$study <- "EMTAB7365"

# 1. Rename metadata column "experiment" to "condition"
seurat_7365$condition <- seurat_7365$experiment
seurat_7365$experiment <- NULL # Removes the old column to keep metadata clean

View(seurat_7365@meta.data)
# 2. Rename condition "MI-Day7" to "7d-PO"
seurat_7365$condition[seurat_7365$condition == "MI-day 7"] <- "7d-PO"

# 3. Subset Seurat object to keep only cells where condition is "7d-PO"
# Using the subset function directly on the metadata column
seurat_7365_7d <- subset(seurat_7365, subset = condition == "7d-PO")

## %% 2. MERGE, QC & UNIFORM PREPROCESSING
fb <- merge(x = seurat_7dpi, 
            y = c(seurat_14dpi, seurat_30dpi, seurat_7365_7d),
            add.cell.ids = c("GSE_7d", "GSE_14d", "GSE_30d", "EMTAB_7d"))


# Join layers if using Seurat V5
if(inherits(fb[["RNA"]], "Assay5")) {
  fb[["RNA"]] <- JoinLayers(fb[["RNA"]])
}

# Quality Control (Mitochondrial genes)
fb[["percent.mt"]] <- PercentageFeatureSet(fb, pattern = "^mt-")

# Standard scRNA-seq filtering threshold
fb <- subset(fb, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 10)

# Normalization & Scaling
fb <- NormalizeData(fb, normalization.method = "LogNormalize", scale.factor = 10000)
fb <- FindVariableFeatures(fb, selection.method = "vst", nfeatures = 2000)
fb <- ScaleData(fb, features = rownames(fb))
fb <- RunPCA(object = fb, npcs = 100)
## %% 3.1 HARMONY BATCH CORRECTION & CLUSTERING
ElbowPlot(fb, ndims = 100)

fb <- RunHarmony(fb, group.by.vars = "study", dims.use = 1:50, max.iter = 100)

# Unsupervised Clustering
fb <- FindNeighbors(fb, reduction = "harmony", dims = 1:50)
fb <- FindClusters(fb, resolution = c(0.5, 0.7, 1)) 
fb <- RunUMAP(fb, reduction = "harmony", dims = 1:50, reduction.name = "umap.harmony")
DimPlot(fb, reduction = "umap.harmony", label = TRUE, group.by = "study")

## %% 3.2 Clustering resolutions
# Visually inspect clustering resolutions on UMAP
resolutions <- c(0.5, 0.7, 1.0)
plots <- list()
for (res in resolutions) {
  fb <- FindClusters(fb, resolution = res) # 1 is default Louvain
  p <- DimPlot(fb, label = TRUE, group.by = "seurat_clusters") + ggtitle(paste("Resolution =", res))
  plots[[as.character(res)]] <- p
}

# Combine UMAP plots by resolution
wrap_plots(plots)

## %% 3.3 DIFFERENTIAL EXPRESSION & ENRICHMENT BY RESOLUTION
de_results_list <- list()
enrichment_plots <- list()

resolutions <- c(0.5, 0.7, 1.0)
assay_name <- DefaultAssay(fb) 

for (res in resolutions) {
  # Set the active identity to the current resolution
  cluster_col <- paste0(assay_name, "_snn_res.", res)
  Idents(fb) <- fb@meta.data[[cluster_col]]
  
  # Find all markers for this resolution
  markers <- FindAllMarkers(fb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  de_results_list[[as.character(res)]] <- markers
  
  # Extract Top 50 significant genes per cluster for enrichment
  top_genes <- markers %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    top_n(n = 50, wt = avg_log2FC)
  
  # Convert to a named list of gene vectors 
  cluster_gene_list <- split(top_genes$gene, top_genes$cluster)
  
  # Run Comparative GO Enrichment Analysis (Biological Process)
  comp_enrich <- compareCluster(
    geneClusters = cluster_gene_list,
    fun = "enrichGO",
    OrgDb = org.Mm.eg.db,
    keyType = "SYMBOL",
    ont = "BP",       
    pvalueCutoff = 0.05
  )
  
  # 5. Generate and store the DotPlot
  p_enrich <- dotplot(comp_enrich, showCategory = 4) +
    ggtitle(paste("GO Biological Process - Resolution", res)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  enrichment_plots[[as.character(res)]] <- p_enrich
}

print(enrichment_plots[["0.5"]])
print(enrichment_plots[["0.7"]])
print(enrichment_plots[["1"]])

write.csv(de_results_list[["0.5"]], "DEGs_Resolution_0.5.csv")
write.csv(de_results_list[["0.7"]], "DEGs_Resolution_0.7.csv")
write.csv(de_results_list[["1"]], "DEGs_Resolution_1.csv")

## The resolution that gives the best separation of cell states is 0.7
fb <- FindClusters(fb, resolution = 0.7)
# Add the resolution 0.7 to the metadata
fb$seurat_clusters <- fb@meta.data$seurat_clusters_0.7
View(fb@meta.data)
## %% 4. IDENTIFYING MYOFIBROBLASTS & Fibroblasts
FeaturePlot(fb, features = c("Postn", "Cthrc1", "Acta2"), pt.size = 0.1)
ggsave("Canonical_Myofibroblast_Markers.jpeg", width = 15, height = 10, dpi = 300)

# Find all markers to assign clusters (fib1-fib9, FMF, RMF)
all_markers <- FindAllMarkers(fb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(all_markers, "All_Clusters_Markers.csv")

cluster_mapping <- c(`3` = "RMF", `6` = "FMF", `0`="fib1", `1`="fib2", `2`="fib3", 
                     `4`="fib4", `5`="fib5", `7`="fib6", `8`="fib7", `9`="fib8", `10`="fib9")
fb <- RenameIdents(fb, cluster_mapping)
fb$clusters <- Idents(fb)

# Figure S16: All Clusters Heatmap
top5_all <- all_markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DoHeatmap(fb, features = top5_all$gene, group.by = "clusters", size = 3, angle = 45) +
  scale_fill_viridis_c(option = "plasma")
ggsave("Fig_S16_All_Clusters_Heatmap.jpeg", width = 14, height = 12, dpi = 300)

## %% 5. RMF vs FMF 
# Subset purely the myofibroblasts

myo_pool <- subset(fb, idents = c("RMF", "FMF"))
myo_pool <- ScaleData(myo_pool, features = rownames(myo_pool))

# Differential Expression
rmf_vs_fmf_de <- FindMarkers(myo_pool, ident.1 = "RMF", ident.2 = "FMF", logfc.threshold = 0.1, min.pct = 0.1)
write.csv(rmf_vs_fmf_de, "RMF_vs_FMF_Differential_Expression.csv")

# Density Plots for key FMF and RMF markers
key_markers <- c("Sfrp2", "Wisp2", "Id2", "Acvrl1", "Col15a1", "Tgfb1", "Cilp", "Crlf1")

do_NebulosaPlot(sample = myo_pool, 
                             features = key_markers, 
                             joint = TRUE, font.size = 22)
ggsave("myo_marker_density.jpeg", width = 10, height = 15, units = "in", dpi = 600)

## %% 6. TEMPORAL ANALYSIS OF CLUSTERS 

CellDimPlot(
  srt = fb, group.by = "clusters",
  label = TRUE, label_insitu = TRUE,
  stat_plot_type = "ring", stat_plot_label = TRUE, stat_plot_size = 0.15
)
ggsave("fb_dimplot_temporal.jpeg", width = 4, height = 4)

CellStatPlot(
  srt = fb,
  group.by = "condition",
  stat.by = "clusters",
  label = TRUE,
  plot_type = "trend"
)
ggsave("fb_statplot_temporal.jpeg", width = 4, height = 4)

## %% 7. GSEA: RMF VS FMF
# Un-thresholded DE for proper GSEA ranking
gsea_de <- FindMarkers(myo_pool, ident.1 = "FMF", ident.2 = "RMF", logfc.threshold = 0, min.pct = 0.1)
gsea_de <- gsea_de[order(-gsea_de$avg_log2FC), ] # Rank by LogFC descending
gsea_de$group1 <- ifelse(gsea_de$avg_log2FC > 0, "FMF", "RMF")

# Format for SCP RunGSEA
geneID <- row.names(gsea_de)
geneScore <- gsea_de$avg_log2FC

# Run Reactome GSEA (Ensure SCP Assay conversion if required by package version)

gsea_reactome <- RunGSEA(
  geneID = geneID, 
  geneScore = geneScore, 
  geneID_groups = gsea_de$group1, 
  db = "Reactome", 
  species = "Mus_musculus"
)

GSEAPlot(res = gsea_reactome, db = "Reactome", plot_type = "dot", direction = "both", topTerm = 20)
ggsave("Reactome_GSEA_DotPlot.jpeg", width = 10, height = 8, dpi = 300)

## %% 8. SMAD1 INHIBITION MODULE
# Load curated intersection of upregulated genes in Smad1 loss (in vivo & in vitro)
vitro_degs <- read.csv("vitro_degs.csv")
vivo_degs  <- read.csv("vivo_degs.csv")

up_vitro <- subset(vitro_degs, log2fc > 0)
up_vivo  <- subset(vivo_degs, log2fc > 0)

# intersect by mgi_symbol
s1module <- intersect(up_vitro$mgi_symbol, up_vivo$mgi_symbol)

# Calculate Module Scores
myo_pool <- AddModuleScore(myo_pool, features = list(s1module), 
                           name = "Smad1_Upregulated_Module", 
                           weight.by.variance = TRUE)

# Ensure precise column name reference (Seurat appends '1' to the name)
VlnPlot(myo_pool, features = "Smad1_Upregulated_Module1", group.by = "clusters", pt.size = 0.05) +
  stat_summary(fun.data = "mean_sdl", geom = "crossbar", width = 0.2) +
  ggtitle("Smad1 Inhibition Module") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("Smad1_Module_Score_Violin.jpeg", width = 5, height = 6, dpi = 300)

# Statistical validation of module score (Wilcoxon Rank Sum Test)
mf1_scores <- myo_pool@meta.data[myo_pool@meta.data$clusters == "RMF", "Smad1_Upregulated_Module1"]
mf2_scores <- myo_pool@meta.data[myo_pool@meta.data$clusters == "FMF", "Smad1_Upregulated_Module1"]
wilcox.test(mf2_scores, mf1_scores)

# Heatmap of Smad1 Module genes present in the object
genes_in_obj <- intersect(s1module, rownames(myo_pool))
DoHeatmap(myo_pool, features = genes_in_obj, group.by = "clusters", label = TRUE) + 
  scale_fill_viridis_c(option = "plasma")
ggsave("Smad1_Module_Heatmap.jpeg", width = 12, height = 10, dpi = 300)

# Custom DB GSEA specific to the Smad1 Inhibition Module
custom_TERM2GENE <- data.frame(
  term = "Smad1KO_Upregulated_Module",
  gene = s1module
)

# Prepare Database via SCP function
SCP_custom_db <- PrepareDB(
  species = "Mus_musculus", 
  db = "SMAD1MOD",
  custom_TERM2GENE = custom_TERM2GENE,
  custom_species = "Mus_musculus",
  custom_IDtype = "symbol",
  custom_version = "v1",
  convert_species = FALSE,
  Ensembl_version = "current_release"
)


# Run GSEA using the custom database
gsea_smad1 <- RunGSEA(
  geneID = geneID, 
  geneScore = geneScore, 
  db = "SMAD1MOD", 
  species = "Mus_musculus"
)

GSEAPlot(res = gsea_smad1, db = "SMAD1MOD", id_use = "Smad1KO_Upregulated_Module", plot_type = "line")
ggsave("Smad1_Module_GSEA.jpeg", width = 6, height = 5, dpi = 300)

# Save Processed Object
saveRDS(fb, file = "Final_Integrated_Fibroblast_Seurat.rds")
