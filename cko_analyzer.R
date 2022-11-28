#AUTHOR: Brikha R Shrestha,PhD. Harvard Medical School, Boston, MA, USA. brikha_shrestha@meei.harvard.edu
#LAST UPDATED: 27 November 2022

library(Seurat)
library(ggplot2)
library(dplyr)
library(DoubletFinder)
library(wesanderson)
library(cowplot)

# PREPARE DATA ---------------------------------------------------------------
# Import raw read counts after downloading from GEO. Assumes data is in the working directory; add filepath if necessary
countMatrix <- read.csv("rawCounts_CKO.csv", row.names = 1)

# Remove anything that is not an SGN (this removes VGNs specifically)
sgnGeneExpression <- apply(countMatrix[c("Cd24a", "C1ql3", "Mafb"),], 2, sum)
sgn <- names(which(sgnGeneExpression>3))
countMatrixSGN <- countMatrix[,sgn]

# FILTERING, NORMALIZATION & SCALING  ----------------------------------------------------------
# Creating a Seurat object
runx1CKO <- CreateSeuratObject(counts=data.matrix(countMatrixSGN), project = "Runx1_LOF", min.cells=5, min.features = 4000)
runx1CKO <-  NormalizeData(runx1CKO, normalization.method = "LogNormalize", scale.factor = 1e5)
runx1CKO <- ScaleData(runx1CKO)
runx1CKO <- FindVariableFeatures (runx1CKO, nfeatures=235)

# CURATE METADATA ----------------------------------------------------------------
runx1CKO <- AddMetaData(runx1CKO, read.csv("metaData_CKO_v2.csv")$GROUP, col.name = "group")
metaData <- runx1CKO@meta.data

# REMOVE GLIAL GENES -------------------------------------------------------------
# Find joint glia+neuron profiles based on expression of known glia-enriched markers 
dataMat <- runx1CKO@assays$RNA@data
glialCells <- intersect(which(dataMat["Mpz",] > 0.9), intersect(which(dataMat["Mbp",] > 0.9), which(dataMat["Plp1",] > 0.9)))
glialCells <- colnames(dataMat)[glialCells]
neuronOnly <- colnames(dataMat)[!colnames(dataMat) %in% glialCells]

# Identify genes expressed at higher levels in glia than in neurons. This is preferred over a presence/absence call due to the sparse nature of scRNA-seq data.
glialGenes <- apply(dataMat, 1, function(x) mean(x[glialCells]) > (2 * mean(x[neuronOnly])))
glialGenes <- names(glialGenes)[glialGenes]

# filter out glia-enriched genes from the list of variable features
runx1CKO@assays$RNA@var.features <- runx1CKO@assays$RNA@var.features[!runx1CKO@assays$RNA@var.features %in% glialGenes]


# PCA & CLUSTERING -----------------------------------------------------
runx1CKO <- RunPCA(runx1CKO, npcs = 20, ndims.print = 1:20, nfeatures.print = 20)
runx1CKO <- FindNeighbors(runx1CKO, reduction = "pca", dims = 1:20, nn.eps = 0.5)
runx1CKO <- FindClusters(runx1CKO, resolution = 0.15, n.start = 10)


# UMAP & FEATURE PLOTS --------------------------------------------------------------------
runx1CKO <- RunUMAP(runx1CKO, dims = 1:10, min.dist = 0.2)

# See how it looks...Note that since UMAP is a stochastic algorithm, the plot below may look slightly different from that in Shrestha et al.
DimPlot(runx1CKO, reduction = "umap", pt.size = 0.6, group.by="seurat_clusters")

# Check expression of genes of interest...
FeaturePlot(runx1CKO, features = "Lypd1", cols = c("grey70","skyblue", "red"), reduction = "umap", pt.size = 1.5) + 
    scale_color_gradientn(colours = c("gray90", "dodgerblue1", "pink1", "red2"), values = c(0,0.5,0.8,3))


# SUBSETS OF INTEREST --------------------------------------------------------------------
conOnly <- rownames(runx1CKO@meta.data)[grep("con", runx1CKO@meta.data$group)]
ckoOnly <- rownames(runx1CKO@meta.data)[grep("cko", runx1CKO@meta.data$group)]

umap1 <- runx1CKO@reductions$umap@cell.embeddings[,"UMAP_1"]
allT1_only <- names(umap1[umap1>-5])

# DIFF. EXP GENES ----------------------------------------------------------
# What is differentially expressed among CON Type I SGNs?
conDiffGenes <- FindAllMarkers(subset(runx1CKO, cells = conOnly), test.use="bimod", return.thresh = 0.001, min.pct=0.3)
conDiffGenes[,"foldDiff"] <- abs(conDiffGenes$avg_log2FC)

# Make a subset of the most differentially expressed
top_conDiffGenes <- conDiffGenes[tail(order(conDiffGenes$foldDiff), 144),]
top_conDiffGeneNames <- (top_conDiffGenes %>% arrange(cluster, avg_log2FC, pct.1)) [,"gene"]
top_conDiffGeneNames <- unique(top_conDiffGeneNames)

# Most enriched among Type I subtypes
aGenes <- (conDiffGenes %>% filter(cluster == "0") %>% arrange(-avg_log2FC)) [1:20,"gene"]
bGenes <- (conDiffGenes %>% filter(cluster == "1") %>% arrange(-avg_log2FC)) [1:20,"gene"]
cGenes <- (conDiffGenes %>% filter(cluster == "2") %>% arrange(-avg_log2FC)) [1:20,"gene"]
combinedDiffGenes <- c(aGenes, bGenes, cGenes)

# What is different between CON and CKO?
mutMarkers <- FindMarkers(runx1CKO, ident.1=conOnly, ident.2=ckoOnly, test.use="bimod")
mutMarkers[,"gene"] <- rownames(mutMarkers)

# DOUBLET DETECTION -------------------------------------------------------
# Determine the optimal pK value; Look for a maxima in BCMVN.Based on McGinnis et al (2019).
sweep.res.list_SGNs <- paramSweep_v3(runx1CKO, PCs = 1:10, sct = FALSE)
sweep.stats_SGNs <- summarizeSweep(sweep.res.list_SGNs, GT = FALSE)

# generate plot for pK sweep
bcmvn_SGNs <- find.pK(sweep.stats_SGNs)
bcmvn_SGNs$pK  <-  as.numeric(as.vector(bcmvn_SGNs$pK))
pal2 <- wes_palette("Darjeeling1", 2, type = "discrete")
pK_sweep <- ggplot(bcmvn_SGNs, aes(x=pK, y=BCmetric)) + 
              geom_point(color=pal2[2]) + 
              geom_line(color=pal2[2]) + 
              theme_classic()
print(pK_sweep)

# Homotypic Doublet Proportion Estimate 
annotations <- runx1CKO@meta.data$seurat_clusters
homotypic.prop <- modelHomotypic(annotations)           ## ex: annotations <- runx1CKO@meta.data$ClusteringResults
nExp_poi <- round(0.075*length(dim(runx1CKO@assays$RNA@counts)[2]))  ## Assuming 7.5% doublet formation rate
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

# Run DoubletFinder 
runx1CKO <- doubletFinder_v3(runx1CKO, PCs = 1:10, pN = 0.25, pK = 0.01, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE)
DimPlot(runx1CKO, reduction = "umap", pt.size = 0.6, group.by="DF.classifications_0.25_0.01_0")
FeaturePlot(runx1CKO, features = "pANN_0.25_0.01_0", reduction = "umap", pt.size = 1.5) +
  scale_color_gradientn(colours = c("gray90", "dodgerblue1", "pink1", "red"), values = c(0,0.3,0.5,1)) +  ylim(-20,12.5) + xlim(-4,15)


# MAKE PLOTS ------------------------------------
# organize data for easy access
metaData <- runx1CKO@meta.data[allT1_only,c("group", "seurat_clusters")]
expMatrix <- runx1CKO@assays$RNA@data[,allT1_only]
coordinates <- runx1CKO@reductions$umap@cell.embeddings[allT1_only,]
genename <- "Calb2"
expData <- cbind(coordinates, data.frame(expMatrix[genename,]), data.frame(metaData))
colnames(expData) <- c("UMAP1", "UMAP2", "counts", "group", "cluster")

# offset CKO cells so they can be seen  easily, offset values can  be adjusted as needed. Resultant UMAP coordinates should only be used for feature visualization.
expData[which(expData$group=="cko"), "UMAP1"] <- expData[which(expData$group=="cko"), "UMAP1"] + 6.6
expData[which(expData$group=="cko"), "UMAP2"] <- expData[which(expData$group=="cko"), "UMAP2"] - 1

# show cluster identities
p1 <- ggplot(expData, aes(x=UMAP1, y=UMAP2)) +
  geom_point(aes(color=cluster),shape=16, alpha=0.6, size=2) +
  scale_color_manual(values = c("springgreen4", "magenta2","dodgerblue1")) +
  theme_void() + theme(legend.position = "right", panel.border = element_rect(colour = "gray90", fill=NA, size=0.6)) 

print(p1) 

#show genotype
p2 <- ggplot(expData, aes(x=UMAP1, y=UMAP2)) +
  geom_point(aes(color=group),shape=16, alpha=0.6, size=2) +
  scale_color_manual(values = c("limegreen", "royalblue3"))  +
  theme_void() + theme(legend.position = "right", panel.border = element_rect(colour = "gray90", fill=NA, size=0.6)) 

print(p2) 

#show gene expression
p3 <- ggplot(expData, aes(x=UMAP1, y=UMAP2)) +
    geom_point(aes(color=counts),shape=16, alpha=0.6, size=2) +
    scale_color_gradientn(colours = c("gray90", "dodgerblue1", "pink1", "red2"), values = c(0,0.2,0.3,1)) +
    theme_void() + theme(legend.position = "right", panel.border = element_rect(colour = "gray90", fill=NA, size=0.6)) 
  
print(p3) 

# HEATMAPS ----------------------------------------------------------------
# Heatmaps showing Top 20 markers in Type I SGN clusters
hMapCon <- DoHeatmap(subset(runx1CKO, cells = intersect(names(runx1CKO@active.ident[which(runx1CKO@active.ident!=3)]), conOnly)), features = combinedDiffGenes, size = 3) + 
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")) )

hMapCko <- DoHeatmap(subset(runx1CKO, cells = intersect(names(runx1CKO@active.ident[which(runx1CKO@active.ident!=3)]), ckoOnly)), features = combinedDiffGenes, size = 3) +
  scale_fill_gradientn(colors = rev(RColorBrewer::brewer.pal(n = 10, name = "RdBu")))

plot_grid(hMapCon,hMapCko)

# SUBTYPE PROPORTIONS --------------
subGroups <- runx1CKO@active.ident
conSubGroup <- subGroups[rownames(metaData)[which(metaData$group=="con")]]
ckoSubGroup <- subGroups[rownames(metaData)[which(metaData$group=="cko")]]

# subtype proportions in CON; dividing by total for Type I SGNs
table(conSubGroup)[1:3] / sum(table(conSubGroup)[1:3]) 

# subtype proportions in CKO
table(ckoSubGroup)[1:3] / sum(table(ckoSubGroup)[1:3]) 


