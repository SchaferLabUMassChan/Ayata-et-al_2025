####Code for MERFISH Analysis #########
### Analysis run using R v4.2.2 and Seurat v5.0.3 ######

#load programs
library(Seurat)
library(dplyr)
library(Matrix)
library(reticulate)
library(scCustomize)
library(ggplot2)
library(sp)

options(scipen = 10000000)

### load in data tables and create Seurat object #######

#read in data table
FADTV411 <- read.table(file='counts_and_metadata_FADTV_411.csv', sep=",", header=TRUE,row.names=1)
# transpose to rows are genes and columns are cells
FADTV411 <- t(FADTV411)
# Remove "Blank" genes
FADTV411.Genes <- FADTV411[1:398,]
# For each cell, calculate the total transcript count
FADTV411[472,] <- colSums(FADTV411.Genes)
# Remove cells with total transcripts <= 40
Real.Cells <- which(FADTV411[472,] > 40)
FADTV411 <- FADTV411[,Real.Cells]
# Remove cells with volume <= 100 um^3
Big.Cells <- which(FADTV411[464,]>100)
FADTV411 <- FADTV411[,Big.Cells]
# Store volume, x position, y position of each cell
volume <- FADTV411[464,]
center.x <- FADTV411[465,]
center.y <- FADTV411[466,]
# Calculate the mean transcript count per cell across the entire sample
mean.RNA <-mean(FADTV411[472,])
# Remove metadata and "blank" genes
FADTV411 <- FADTV411[1:398,]
# Normalize data by the mean transcript count per cell across the entire sample
FADTV411 <- FADTV411/mean.RNA
# Normalize data by each cell's volume
FADTV411 <- FADTV411/volume
# create Seurat Object with volume, x position, and y position as metadata
FADTV411_s <-CreateSeuratObject(counts=FADTV411)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
FADTV411_s@meta.data <- cbind(FADTV411_s@meta.data,Volumes)
FADTV411_s@meta.data <- cbind(FADTV411_s@meta.data,Center.x)
FADTV411_s@meta.data <- cbind(FADTV411_s@meta.data,Center.y)

## repeat for samples 2-8
FADTV348 <- read.table(file='counts_and_metadata_FADTV_348.csv', sep=",", header=TRUE,row.names=1)
FADTV348 <- t(FADTV348)
FADTV348.Genes <- FADTV348[1:398,]
FADTV348[472,] <- colSums(FADTV348.Genes)
Real.Cells <- which(FADTV348[472,] > 40)
FADTV348 <- FADTV348[,Real.Cells]
Big.Cells <- which(FADTV348[464,]>100)
FADTV348 <- FADTV348[,Big.Cells]
volume <- FADTV348[464,]
center.x <- FADTV348[465,]
center.y <- FADTV348[466,]
mean.RNA <-mean(FADTV348[472,])
FADTV348 <- FADTV348[1:398,]
FADTV348 <- FADTV348/mean.RNA
FADTV348 <- FADTV348/volume
FADTV348_s <-CreateSeuratObject(counts=FADTV348)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
FADTV348_s@meta.data <- cbind(FADTV348_s@meta.data,Volumes)
FADTV348_s@meta.data <- cbind(FADTV348_s@meta.data,Center.x)
FADTV348_s@meta.data <- cbind(FADTV348_s@meta.data,Center.y)

FADPU538 <- read.table(file='counts_and_metadata_FADPU_538.csv', sep=",", header=TRUE,row.names=1)
FADPU538 <- t(FADPU538)
FADPU538.Genes <- FADPU538[1:398,]
FADPU538[472,] <- colSums(FADPU538.Genes)
Real.Cells <- which(FADPU538[472,] > 40)
FADPU538 <- FADPU538[,Real.Cells]
Big.Cells <- which(FADPU538[464,]>100)
FADPU538 <- FADPU538[,Big.Cells]
volume <- FADPU538[464,]
center.x <- FADPU538[465,]
center.y <- FADPU538[466,]
mean.RNA <-mean(FADPU538[472,])
FADPU538 <- FADPU538[1:398,]
FADPU538 <- FADPU538/mean.RNA
FADPU538 <- FADPU538/volume
FADPU538_s <-CreateSeuratObject(counts=FADPU538)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
FADPU538_s@meta.data <- cbind(FADPU538_s@meta.data,Volumes)
FADPU538_s@meta.data <- cbind(FADPU538_s@meta.data,Center.x)
FADPU538_s@meta.data <- cbind(FADPU538_s@meta.data,Center.y)

FADTV635 <- read.table(file='counts_and_metadata_FADTV_635.csv', sep=",", header=TRUE,row.names=1)
FADTV635 <- t(FADTV635)
FADTV635.Genes <- FADTV635[1:398,]
FADTV635[472,] <- colSums(FADTV635.Genes)
Real.Cells <- which(FADTV635[472,] > 40)
FADTV635 <- FADTV635[,Real.Cells]
Big.Cells <- which(FADTV635[464,]>100)
FADTV635 <- FADTV635[,Big.Cells]
volume <- FADTV635[464,]
center.x <- FADTV635[465,]
center.y <- FADTV635[466,]
mean.RNA <-mean(FADTV635[472,])
FADTV635 <- FADTV635[1:398,]
FADTV635 <- FADTV635/mean.RNA
FADTV635 <- FADTV635/volume
FADTV635_s <-CreateSeuratObject(counts=FADTV635)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
FADTV635_s@meta.data <- cbind(FADTV635_s@meta.data,Volumes)
FADTV635_s@meta.data <- cbind(FADTV635_s@meta.data,Center.x)
FADTV635_s@meta.data <- cbind(FADTV635_s@meta.data,Center.y)

FADPU521 <- read.table(file='counts_and_metadata_FADPU_521.csv', sep=",", header=TRUE,row.names=1)
FADPU521 <- t(FADPU521)
FADPU521.Genes <- FADPU521[1:398,]
FADPU521[472,] <- colSums(FADPU521.Genes)
Real.Cells <- which(FADPU521[472,] > 40)
FADPU521 <- FADPU521[,Real.Cells]
Big.Cells <- which(FADPU521[464,]>100)
FADPU521 <- FADPU521[,Big.Cells]
volume <- FADPU521[464,]
center.x <- FADPU521[465,]
center.y <- FADPU521[466,]
mean.RNA <-mean(FADPU521[472,])
FADPU521 <- FADPU521[1:398,]
FADPU521 <- FADPU521/mean.RNA
FADPU521 <- FADPU521/volume
FADPU521_s <-CreateSeuratObject(counts=FADPU521)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
FADPU521_s@meta.data <- cbind(FADPU521_s@meta.data,Volumes)
FADPU521_s@meta.data <- cbind(FADPU521_s@meta.data,Center.x)
FADPU521_s@meta.data <- cbind(FADPU521_s@meta.data,Center.y)

FADPU663 <- read.table(file='counts_and_metadata_FADPU_663.csv', sep=",", header=TRUE,row.names=1)
FADPU663 <- t(FADPU663)
FADPU663.Genes <- FADPU663[1:398,]
FADPU663[472,] <- colSums(FADPU663.Genes)
Real.Cells <- which(FADPU663[472,] > 40)
FADPU663 <- FADPU663[,Real.Cells]
Big.Cells <- which(FADPU663[464,]>100)
FADPU663 <- FADPU663[,Big.Cells]
volume <- FADPU663[464,]
center.x <- FADPU663[465,]
center.y <- FADPU663[466,]
mean.RNA <-mean(FADPU663[472,])
FADPU663 <- FADPU663[1:398,]
FADPU663 <- FADPU663/mean.RNA
FADPU663 <- FADPU663/volume
FADPU663_s <-CreateSeuratObject(counts=FADPU663)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
FADPU663_s@meta.data <- cbind(FADPU663_s@meta.data,Volumes)
FADPU663_s@meta.data <- cbind(FADPU663_s@meta.data,Center.x)
FADPU663_s@meta.data <- cbind(FADPU663_s@meta.data,Center.y)

FADPU662 <- read.table(file='counts_and_metadata_FADPU_662.csv', sep=",", header=TRUE,row.names=1)
FADPU662 <- t(FADPU662)
FADPU662.Genes <- FADPU662[1:398,]
FADPU662[472,] <- colSums(FADPU662.Genes)
Real.Cells <- which(FADPU662[472,] > 40)
FADPU662 <- FADPU662[,Real.Cells]
Big.Cells <- which(FADPU662[464,]>100)
FADPU662 <- FADPU662[,Big.Cells]
volume <- FADPU662[464,]
center.x <- FADPU662[465,]
center.y <- FADPU662[466,]
mean.RNA <-mean(FADPU662[472,])
FADPU662 <- FADPU662[1:398,]
FADPU662 <- FADPU662/mean.RNA
FADPU662 <- FADPU662/volume
FADPU662_s <-CreateSeuratObject(counts=FADPU662)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
FADPU662_s@meta.data <- cbind(FADPU662_s@meta.data,Volumes)
FADPU662_s@meta.data <- cbind(FADPU662_s@meta.data,Center.x)
FADPU662_s@meta.data <- cbind(FADPU662_s@meta.data,Center.y)

FADPU668 <- read.table(file='counts_and_metadata_FADPU_668.csv', sep=",", header=TRUE,row.names=1)
FADPU668 <- t(FADPU668)
FADPU668.Genes <- FADPU668[1:398,]
FADPU668[472,] <- colSums(FADPU668.Genes)
Real.Cells <- which(FADPU668[472,] > 40)
FADPU668 <- FADPU668[,Real.Cells]
Big.Cells <- which(FADPU668[464,]>100)
FADPU668 <- FADPU668[,Big.Cells]
volume <- FADPU668[464,]
center.x <- FADPU668[465,]
center.y <- FADPU668[466,]
mean.RNA <-mean(FADPU668[472,])
FADPU668 <- FADPU668[1:398,]
FADPU668 <- FADPU668/mean.RNA
FADPU668 <- FADPU668/volume
FADPU668_s <-CreateSeuratObject(counts=FADPU668)
Volumes <- data.frame(volume)
Center.x <- data.frame(center.x)
Center.y <- data.frame(center.y)
FADPU668_s@meta.data <- cbind(FADPU668_s@meta.data,Volumes)
FADPU668_s@meta.data <- cbind(FADPU668_s@meta.data,Center.x)
FADPU668_s@meta.data <- cbind(FADPU668_s@meta.data,Center.y)

#create metadata for conditions
FADTV411_s$sample <- "FADTV411"
FADTV411_s$genotype <- "5xFAD_PUhigh"
FADTV348_s$sample <- "FADTV348"
FADTV348_s$genotype <- "5xFAD"
FADPU538_s$sample <- "FADPU538"
FADPU538_s$genotype <- "5xFAD"
FADTV635_s$sample <- "FADTV635"
FADTV635_s$genotype <- "Ctrl"
FADPU521_s$sample <- "FADPU521"
FADPU521_s$genotype <- "5xFAD_PUlow"
FADPU663_s$sample <- "FADPU663"
FADPU663_s$genotype <- "Ctrl"
FADPU662_s$sample <- "FADPU662"
FADPU662_s$genotype <- "5xFAD"
FADPU668_s$sample <- "FADPU668"
FADPU668_s$genotype <- "5xFAD_PUlow"


### Clustering and annotation of all cell types ######

#merge all samples into a single object 
all.combined <- merge(FADTV411_s, y = c(FADTV348_s, FADPU538_s, FADTV635_s, FADPU521_s,FADPU663_s ,FADPU662_s,FADPU668_s ),merge.data=TRUE, add.cell.ids = c("1", "2", "3", "4", "5","6","7","8"))
all.combined <- JoinLayers(all.combined)

#remove cells with nFeautures_RNA <= 10
all.combined <- subset(all.combined,nFeature_RNA > 10)
#normalize data and find variable features (all transcripts included as variable features since only 398 genes are on the panel)
all.combined <- NormalizeData(all.combined)
all.combined <- FindVariableFeatures(all.combined, selection.method = "vst", nfeatures = 398)
#scale data and run PCA
all.genes <- rownames(all.combined)
all.combined <- ScaleData(all.combined, features = all.genes)
all.combined <- RunPCA(all.combined, features = VariableFeatures(object = all.combined),npcs = 50)
# test of number of significant Principal Components. 39 (p-value < 0.05) are used for further analysis
all.combined <- JackStraw(all.combined, num.replicate = 100,dims=40)
all.combined <- ScoreJackStraw(all.combined, dims = 1:40)
JackStrawPlot(all.combined, dims = 1:40)
#clustering and linear dimensionality reduction
all.combined <- FindNeighbors(all.combined, dims = 1:39)
all.combined <- FindClusters(all.combined, resolution = 1.8)
all.combined <- RunUMAP(all.combined, dims = 1:39)

### cell-type annotation of clusters
Idents(all.combined) <- all.combined@meta.data$seurat_clusters
new.cluster.ids <- c("Excitatory Neurons","Inhibitory Neurons","Endothelial Cells","Oligodendrocytes","Astrocytes",
                     "Excitatory Neurons","Excitatory Neurons","Excitatory Neurons","Oligodendrocytes","Excitatory Neurons",
                     "Inhibitory Neurons","Oligodendrocytes","Inhibitory Neurons","Excitatory Neurons","Astrocytes",
                     "Microglia","Inhibitory Neurons","OPCs","Oligodendrocyte/Neuron Hybrids","Microglia",
                     "Astrocyte/Neuron Hybrids","Meningeal Cells","Excitatory Neurons","Pericytes","Oligodendrocyte/Neuron Hybrids",
                     "Inhibitory Neurons","Inhibitory Neurons","Endothelial Cells","Oligodendrocytes","Oligodendrocyte/Astrocyte Hybrids",
                     "Oligodendrocytes","Border-Associated Macrophages","Excitatory Neurons","Microglia/Oligodendrocyte Hybrids","Endothelial Cell/Oligodendrocyte Hybrids",
                     "Excitatory Neurons","Oligodendrocytes","Meningeal Cells","Choroid Plexus","Oligodendrocytes",
                     "Inhibitory Neurons","Inhibitory Neurons","Inhibitory Neurons","Oligodendrocyte/Neuron Hybrids","Ependymal Cells",
                     "Astrocyte/Neuron Hybrids","Astrocyte/Neuron Hybrids","Astrocyte/Neuron Hybrids","Excitatory Neurons","Microglia/Astrocyte Hybrids",
                     "Microglia/Endothelial Cell Hybrids","Astrocyte/Endothelial Cell Hybrids","FADPU521 Artifacts","Inhibitory Neurons","Oligodendrocytes",
                     "Lymphoid Cells","OPCs","OPC/Endothelial Cell Hybrids","Microglia","Microglia/OPC Hybrids",
                     "Ifng Plaques","OPC/Astrocyte/Neuron Hybrids","Excitatory Neurons","Inhibitory Neurons","Oligodendrocytes",
                     "Inhibitory Neurons","Oligodendrocyte/Astrocyte Hybrids","Oligodendrocytes","Excitatory Neurons","Neuron/Oligodendocyte Hybrids",
                     "Neuron/OPC Hybrids","Broad Hybrids","Broad Hybrids","Broad Hybrids","Excitatory Neurons",
                     "Neuron/Oligodendrocyte Hybrids","Oligodendrocytes","Inhibitory Neurons"
)
names(new.cluster.ids) <- levels(all.combined)
all.combined <- RenameIdents(all.combined, new.cluster.ids)
all.combined$celltype <- Idents(all.combined)

## data output
saveRDS(all.combined, "All celltypes_annotated.rds")

### Iterative reclustering and annotation after removal of clusters containing artifacts or hybrids####
all.combined <- readRDS("250122_all_cells_annotated.rds")

#subset to remove clusters containing artifacts or hybrids - round 1
Idents(all.combined) <- all.combined$celltype
all.combined.clean <- subset(all.combined, idents = c("Excitatory Neurons","Inhibitory Neurons","Endothelial Cells","Oligodendrocytes","Astrocytes",
                                                       "Microglia","OPCs","Meningeal Cells","Pericytes","Border-Associated Macrophages","Choroid Plexus","Ependymal Cells",
                                                       "Lymphoid Cells","Ifng Plaques"))

#reclustering - round 1
all.combined.clean <- NormalizeData(all.combined.clean)
all.combined.clean <- FindVariableFeatures(all.combined.clean, selection.method = "vst", nfeatures = 398)
all.genes <- rownames(all.combined.clean)
all.combined.clean <- ScaleData(all.combined.clean, features = all.genes)
all.combined.clean <- RunPCA(all.combined.clean, features = VariableFeatures(object = all.combined.clean),npcs = 50)
all.combined.clean <- JackStraw(all.combined.clean, num.replicate = 100,dims=40)
all.combined.clean <- ScoreJackStraw(all.combined.clean, dims = 1:40)
JackStrawPlot(all.combined.clean, dims = 1:40)
all.combined.clean <- FindNeighbors(all.combined.clean, dims = 1:39)
all.combined.clean <- FindClusters(all.combined.clean, resolution = 1.8)
all.combined.clean <- RunUMAP(all.combined.clean, dims = 1:39)

# annotate after clean round 1 ####
Idents(all.combined.clean) <- all.combined.clean@meta.data$seurat_clusters
new.cluster.ids <- c("Inhibitory Neurons","Excitatory Neurons","Oligodendrocytes","Inhibitory Neurons","Astrocytes",
                     "Excitatory Neurons","Excitatory Neurons","Excitatory Neurons","Oligodendrocytes","Excitatory Neurons",
                     "Excitatory Neurons","Oligodendrocytes","Microglia","Endothelial Cells","Astrocytes",
                     "Inhibitory Neurons","Endothelial Cells","OPCs","Excitatory Neurons","Pericytes",
                     "Excitatory Neurons","Pericyte/Endothelial/Meninges Hybrids","Microglia","Meningeal Cells","Inhibitory Neurons",
                     "Inhibitory Neurons","Oligodendrocytes","Inhibitory Neurons","Excitatory Neurons","Border-Associated Macrophages",
                     "Oligodendrocytes","Excitatory Neurons","Oligodendrocytes","Oligodendrocytes","Excitatory Neurons",
                     "Choroid Plexus","Inhibitory Neurons","Meningeal Cells","Inhibitory Neurons","Inhibitory Neurons",
                     "Ependymal Cells","Excitatory Neurons","Antigen-Presenting Cells","Inhibitory Neurons","Oligodendrocytes",
                     "OPCs","T Cells","Astrocytes","Microglia","Oligodendrocytes",
                     "Endothelial Cells","Ifng Plaques","Excitatory Neurons","Inhibitory Neurons","Excitatory Neurons",
                     "Inhibitory Neurons","Hybrid","Hybrid","Hybrid","Hybrid",
                     "Hybrid","Hybrid","Hybrid","Hybrid","Hybrid",
                     "Hybrid","Excitatory Neurons","Hybrid"
)
names(new.cluster.ids) <- levels(all.combined.clean)
all.combined.clean <- RenameIdents(all.combined.clean, new.cluster.ids)
all.combined.clean$celltype_clean1 <- Idents(all.combined.clean)

#subset to remove clusters containing artifacts or hybrids - Round 2 #
Idents(all.combined.clean) <- all.combined.clean$celltype_clean1
all.combined.clean2 <- subset(all.combined.clean, idents = c("Excitatory Neurons","Inhibitory Neurons","Endothelial Cells","Oligodendrocytes","Astrocytes",
                                                      "Microglia","OPCs","Meningeal Cells","Pericytes","Border-Associated Macrophages","Choroid Plexus","Ependymal Cells",
                                                      "Antigen-Presenting Cells","Ifng Plaques","T Cells"))

#Reclustering - Round 2
all.combined.clean2 <- NormalizeData(all.combined.clean2)
all.combined.clean2 <- FindVariableFeatures(all.combined.clean2, selection.method = "vst", nfeatures = 398)
all.genes <- rownames(all.combined.clean2)
all.combined.clean2 <- ScaleData(all.combined.clean2, features = all.genes)
all.combined.clean2 <- RunPCA(all.combined.clean2, features = VariableFeatures(object = all.combined.clean2),npcs = 50)
all.combined.clean2 <- JackStraw(all.combined.clean2, num.replicate = 100,dims=40)
all.combined.clean2 <- ScoreJackStraw(all.combined.clean2, dims = 1:40)
JackStrawPlot(all.combined.clean2, dims = 1:40)
all.combined.clean2 <- FindNeighbors(all.combined.clean2, dims = 1:39)
all.combined.clean2 <- FindClusters(all.combined.clean2, resolution = 1.8)
all.combined.clean2 <- RunUMAP(all.combined.clean2, dims = 1:39)

# annotation after round 2
Idents(all.combined.clean2) <- all.combined.clean2@meta.data$seurat_clusters
new.cluster.ids <- c("Excitatory Neurons","Excitatory Neurons","Oligodendrocytes","Excitatory Neurons","Astrocytes",
                     "Inhibitory Neurons","Inhibitory Neurons","Excitatory Neurons","Oligodendrocytes","Oligodendrocytes",
                     "Astrocytes","Excitatory Neurons","Excitatory Neurons","Endothelial Cells","Endothelial Cells",
                     "Inhibitory Neurons","Microglia","Microglia","Inhibitory Neurons","OPCs",
                     "Inhibitory Neurons","Inhibitory Neurons","Meningeal Cells","Pericytes","Inhibitory Neurons",
                     "Inhibitory Neurons","Oligodendrocytes","Excitatory Neurons","Border-Associated Macrophages","Oligodendrocytes",
                     "Excitatory Neurons","Oligodendrocytes","Inhibitory Neurons","Choroid Plexus","Oligodendrocytes",
                     "Inhibitory Neurons","Inhibitory Neurons","Ependymal Cells","Excitatory Neurons","Antigen-Presenting Cells",
                     "Inhibitory Neurons","OPC/Neurons","T Cells","Astrocytes","Inhibitory Neurons",
                     "Microglia","Oligodendrocytes","Ifng Plaques","Inhibitory Neurons","Excitatory Neurons",
                     "Inhibitory Neurons","Hybrids","Excitatory Neurons","Excitatory Neurons","Excitatory Neurons",
                     "Inhibitory Neurons","Oligodendrocytes","Excitatory Neurons","Excitatory Neurons","Oligodendrocytes",
                     "Hybrids","Oligodendrocytes","Excitatory Neurons","Excitatory Neurons","Excitatory Neurons"
)
names(new.cluster.ids) <- levels(all.combined.clean2)
all.combined.clean2 <- RenameIdents(all.combined.clean2, new.cluster.ids)
all.combined.clean2$celltype_clean1 <- Idents(all.combined.clean2)

#subset to remove clusters containing artifacts or hybrids - Round 3
Idents(all.combined.clean2) <- all.combined.clean2$celltype_clean1
all.combined.clean3 <- subset(all.combined.clean2, idents = c("Excitatory Neurons","Inhibitory Neurons","Endothelial Cells","Oligodendrocytes","Astrocytes",
                                                             "Microglia","OPCs","Meningeal Cells","Pericytes","Border-Associated Macrophages","Choroid Plexus","Ependymal Cells",
                                                             "Antigen-Presenting Cells","Ifng Plaques","T Cells"))

#Reclustering - Round 3
all.combined.clean3 <- NormalizeData(all.combined.clean3)
all.combined.clean3 <- FindVariableFeatures(all.combined.clean3, selection.method = "vst", nfeatures = 398)
all.genes <- rownames(all.combined.clean3)
all.combined.clean3 <- ScaleData(all.combined.clean3, features = all.genes)
all.combined.clean3 <- RunPCA(all.combined.clean3, features = VariableFeatures(object = all.combined.clean3),npcs = 50)
all.combined.clean3 <- JackStraw(all.combined.clean3, num.replicate = 100,dims=40)
all.combined.clean3 <- ScoreJackStraw(all.combined.clean3, dims = 1:40)
JackStrawPlot(all.combined.clean3, dims = 1:40)
all.combined.clean3 <- FindNeighbors(all.combined.clean3, dims = 1:36)
all.combined.clean3 <- FindClusters(all.combined.clean3, resolution = 1.8)
all.combined.clean3 <- RunUMAP(all.combined.clean3, dims = 1:36)

# annotate after Round 3 ####
Idents(all.combined.clean3) <- all.combined.clean3@meta.data$seurat_clusters
new.cluster.ids <- c("Inhibitory Neurons","Excitatory Neurons","Oligodendrocytes","Astrocytes","Excitatory Neurons",
                     "Oligodendrocytes","Inhibitory Neurons","Excitatory Neurons","Oligodendrocytes","Excitatory Neurons",
                     "Excitatory Neurons","Excitatory Neurons","Astrocytes","Microglia","Endothelial Cells",
                     "Inhibitory Neurons","Endothelial Cells","Inhibitory Neurons","OPCs","Excitatory Neurons",
                     "Excitatory Neurons","Microglia","Pericytes","Meningeal Fibroblasts","Inhibitory Neurons",
                     "Inhibitory Neurons","Excitatory Neurons","Oligodendrocytes","Excitatory Neurons","Border-Associated Macrophages",
                     "Oligodendrocytes","Antigen-Presenting Cells","Choroid Plexus Cells","Oligodendrocytes","Inhibitory Neurons",
                     "Inhibitory Neurons","Inhibitory Neurons","Ependymal Cells","Excitatory Neuron - Artifacts","Inhibitory Neurons",
                     "Astrocyte/Neuron Hybrids","Oligodendrocytes","T Cells","Neuron/BAM Hybrids","Ifng Plaques",
                     "Hybrids","Hybrids","Hybrids","Hybrids","Inhibitory Neurons",
                     "Hybrids","Hybrids"
)
names(new.cluster.ids) <- levels(all.combined.clean3)
all.combined.clean3 <- RenameIdents(all.combined.clean3, new.cluster.ids)
all.combined.clean3$celltype_clean3 <- Idents(all.combined.clean3)
#subset to remove clusters containing artifacts or hybrids - Round 4 #
Idents(all.combined.clean3) <- all.combined.clean3$celltype_clean3
all.combined.clean4 <- subset(all.combined.clean3, idents = c("Excitatory Neurons","Inhibitory Neurons","Endothelial Cells","Oligodendrocytes","Astrocytes",
                                                              "Microglia","OPCs","Meningeal Fibroblasts","Pericytes","Border-Associated Macrophages","Choroid Plexus Cells","Ependymal Cells",
                                                              "Antigen-Presenting Cells","Ifng Plaques","T Cells"))
all.combined.clean4 <- NormalizeData(all.combined.clean4)
all.combined.clean4 <- FindVariableFeatures(all.combined.clean4, selection.method = "vst", nfeatures = 398)
all.genes <- rownames(all.combined.clean4)
all.combined.clean4 <- ScaleData(all.combined.clean4, features = all.genes)
all.combined.clean4 <- RunPCA(all.combined.clean4, features = VariableFeatures(object = all.combined.clean4),npcs = 50)
all.combined.clean4 <- JackStraw(all.combined.clean4, num.replicate = 100,dims=45)
all.combined.clean4 <- ScoreJackStraw(all.combined.clean4, dims = 1:45)
JackStrawPlot(all.combined.clean4, dims = 1:45)
all.combined.clean4 <- FindNeighbors(all.combined.clean4, dims = 1:38)
all.combined.clean4 <- FindClusters(all.combined.clean4, resolution = 2.4)
all.combined.clean4 <- RunUMAP(all.combined.clean4, dims = 1:38)

#### annotate after Round 4
Idents(all.combined.clean4) <- all.combined.clean4@meta.data$seurat_clusters
new.cluster.ids <- c("Oligodendrocytes","Excitatory Neurons","Inhibitory Neurons","Excitatory Neurons","Excitatory Neurons",
                     "Inhibitory Neurons","Excitatory Neurons","Excitatory Neurons","Oligodendrocytes","Endothelial Cells",
                     "Excitatory Neurons","Inhibitory Neurons","Astrocytes","Excitatory Neurons","Endothelial Cells",
                     "Astrocytes","Oligodendrocytes","Oligodendrocytes","Inhibitory Neurons","Astrocytes",
                     "OPCs","Excitatory Neurons","Pericytes","Microglia","Inhibitory Neurons",
                     "Microglia","Excitatory Neurons","Inhibitory Neurons","Meningeal Fibroblasts","Microglia",
                     "Astrocytes","Excitatory Neurons","Inhibitory Neurons","Oligodendrocytes","Inhibitory Neurons",
                     "Border-Associated Macrophages","Excitatory Neurons","Oligodendrocytes","Inhibitory Neurons","Oligodendrocytes",
                     "Excitatory Neurons","Choroid Plexus Cells","Oligodendrocytes","Inhibitory Neurons","Microglia/Neuron Hybrids",
                     "Inhibitory Neurons","Meningeal Fibroblasts","Ependymal Cells","Excitatory Neurons","Excitatory Neurons",
                     "Oligodendrocytes","Antigen-Presenting Cells","Inhibitory Neurons","Oligodendrocytes","Excitatory Neurons",
                     "T Cells","Inhibitory Neurons","Excitatory/BAMs","OL/T Cells","Oligodendrocytes",
                     "Plaques","Hybrids","Hybrids","Inhibitory Neurons","64",
                     "Hybrids", "Hybrids","Excitatory Neurons","Hybrids","Excitatory Neurons",
                     "Hybrids")
names(new.cluster.ids) <- levels(all.combined.clean4)
all.combined.clean4 <- RenameIdents(all.combined.clean4, new.cluster.ids)
all.combined.clean4$celltype_clean4 <- Idents(all.combined.clean4)

#### remove portions of FADPU668 that show artifact near edge of sample####
#subset out FADPU668
Idents(object=all.combined.clean4) <- all.combined.clean4@meta.data$sample
FADPU668 <- subset(x = all.combined.clean4, idents = c("FADPU668"),invert = FALSE)
Other_Samples <- subset(x = all.combined.clean4, idents = c("FADPU668"),invert = TRUE)

#define X,Y coordinates of polygon outlining the region containing the artifact
FADPU668.x <- FADPU668@meta.data[,5]
FADPU668.y <- FADPU668@meta.data[,6]
FADPU668_exclude.x <- c(4432,4388,4201,4185,3993,3985,3797,3789,3594,3586,2975,2971,2779,2763,2388,2364,1353,1353)
FADPU668_exclude.y <- c(6328,6184,6186,6008,5992,5793,5789,5585,5585,5178,5170,4870,4966,4774,4782,4379,4371,6603)

#identify cells contained within polygon region
FADPU668_exclude <- point.in.polygon(FADPU668.x, FADPU668.y, FADPU668_exclude.x, FADPU668_exclude.y)
FADPU668.info <- cbind.data.frame(FADPU668_exclude)
FADPU668_exclude.info <- rowSums(FADPU668.info)
FADPU668.info <- cbind.data.frame(FADPU668_exclude.info )
FADPU668.info <- FADPU668.info %>% mutate(Region = if_else(FADPU668_exclude.info > 0,"Exclude","Include"))
Region_FADPU668 <-FADPU668.info[,2]

#Remove cells contained within polygon region
FADPU668 <- AddMetaData(FADPU668, metadata = Region_FADPU668, col.name = "Region")
FADPU668.include <- subset(FADPU668, Region == c("Include"))

##merge FADPU668 back with other samples
all.combined.clean4_no_artifact <- merge(FADPU668.include,y=c(Other_Samples))
all.combined.clean4_no_artifact <- JoinLayers(all.combined.clean4_no_artifact)

#####subset to remove clusters containing artifacts or hybrids - Round 5
Idents(all.combined.clean4_no_artifact) <- all.combined.clean4_no_artifact$celltype_clean4
all.combined.clean5 <- subset(all.combined.clean4_no_artifact, idents = c("Excitatory Neurons","Inhibitory Neurons","Endothelial Cells","Oligodendrocytes","Astrocytes",
                                                              "Microglia","OPCs","Meningeal Fibroblasts","Pericytes","Border-Associated Macrophages","Choroid Plexus Cells","Ependymal Cells",
                                                              "Antigen-Presenting Cells","T Cells"))
#Reclustering - round 5
all.combined.clean5 <- NormalizeData(all.combined.clean5)
all.combined.clean5 <- FindVariableFeatures(all.combined.clean5, selection.method = "vst", nfeatures = 398)
all.genes <- rownames(all.combined.clean5)
all.combined.clean5 <- ScaleData(all.combined.clean5, features = all.genes)
all.combined.clean5 <- RunPCA(all.combined.clean5, features = VariableFeatures(object = all.combined.clean5),npcs = 50)
all.combined.clean5 <- JackStraw(all.combined.clean5, num.replicate = 100,dims=45)
all.combined.clean5 <- ScoreJackStraw(all.combined.clean5, dims = 1:45)
JackStrawPlot(all.combined.clean5, dims = 1:45)
all.combined.clean5 <- FindNeighbors(all.combined.clean5, dims = 1:38)
all.combined.clean5 <- FindClusters(all.combined.clean5, resolution = 2.4)
all.combined.clean5 <- RunUMAP(all.combined.clean5, dims = 1:38)

#### annotate after round 5
Idents(all.combined.clean6) <- all.combined.clean6@meta.data$seurat_clusters
new.cluster.ids <- c("Oligodendrocytes","Inhibitory Neurons","Excitatory Neurons","Excitatory Neurons","Inhibitory Neurons",
                     "Excitatory Neurons","Excitatory Neurons","Excitatory Neurons","Oligodendrocytes","Astrocytes",
                     "Endothelial Cells","Inhibitory Neurons","Endothelial Cells","Microglia","Excitatory Neurons",
                     "Excitatory Neurons","Astrocytes","Oligodendrocytes","Oligodendrocytes","Excitatory Neurons",
                     "Excitatory Neurons","OPCs","Astrocytes","Pericytes","Excitatory Neurons",
                     "Inhibitory Neurons","Inhibitory Neurons","Astrocytes","Inhibitory Neurons","Meningeal Fibroblasts",
                     "Excitatory Neurons","Microglia","Inhibitory Neurons","Oligodendrocytes","Microglia",
                     "Inhibitory Neurons","Border-Associated Macrophages","Excitatory Neurons","Oligodendrocytes","Oligodendrocytes",
                     "Choroid Plexus Cells","Inhibitory Neurons","Newly-Formed Oligodendrocytes","Inhibitory Neurons","Meningeal Fibroblasts",
                     "Inhibitory Neurons","Ependymal Cells","Excitatory Neurons","Excitatory Neurons","Inhibitory Neurons",
                     "Antigen-Presenting Cells","Oligodendrocytes","Oligodendrocytes","T Cells","Oligodendrocytes",
                     "Oligodendrocytes","Excitatory Neurons","Inhibitory Neurons","Excitatory Neurons","Excitatory Neurons",
                     "Oligodendrocytes","Excitatory Neurons","Inhibitory Neurons","Inhibitory Neurons","Inhibitory Neurons",
                     "Inhibitory Neurons","Inhibitory Neurons","Oligodendrocytes","Excitatory Neurons","Inhibitory Neurons",
                     "Excitatory Neurons","Excitatory Neurons","Excitatory Neurons"
                                      )
names(new.cluster.ids) <- levels(all.combined.clean6)
all.combined.clean6 <- RenameIdents(all.combined.clean6, new.cluster.ids)
all.combined.clean6$celltype_clean <- Idents(all.combined.clean6)

## Data Output
saveRDS(all.combined.clean6, "All celltypes_no hybrids_annotated.rds")

#### Subclustering of microglia, starting from Seurat object containing all cells including hybrids #####
all.combined <- readRDS("250122_all_cells_annotated.rds")

## Subset to include clusters containing microglia, microglia hybrids, and border-associated macrophages 
microglia.8sample.unclean <- subset(all.combined, idents = c("Microglia","Microglia/Oligodendrocyte Hybrids","Microglia/Endothelial Cell Hybrids","Microglia/Astrocyte Hybrids","Microglia/OPC Hybrids","Border-Associated Macrophages"))

#### remove portions of FADPU668 that show artifact near edge of sample ####
Idents(object=microglia.8sample.unclean) <- microglia.8sample.unclean@meta.data$sample
FADPU668 <- subset(x = microglia.8sample.unclean, idents = c("FADPU668"),invert = FALSE)
Other_Samples <- subset(x = microglia.8sample.unclean, idents = c("FADPU668"),invert = TRUE)
FADPU668.x <- FADPU668@meta.data[,5]
FADPU668.y <- FADPU668@meta.data[,6]

# X, Y coordinates of polygon containing region to be removed
FADPU668_exclude.x <- c(4432,4388,4201,4185,3993,3985,3797,3789,3594,3586,2975,2971,2779,2763,2388,2364,1353,1353)
FADPU668_exclude.y <- c(6328,6184,6186,6008,5992,5793,5789,5585,5585,5178,5170,4870,4966,4774,4782,4379,4371,6603)

# Identify microglia contained within polygon
FADPU668_exclude <- point.in.polygon(FADPU668.x, FADPU668.y, FADPU668_exclude.x, FADPU668_exclude.y)
FADPU668.info <- cbind.data.frame(FADPU668_exclude)
FADPU668_exclude.info <- rowSums(FADPU668.info)
FADPU668.info <- cbind.data.frame(FADPU668_exclude.info )
FADPU668.info <- FADPU668.info %>% mutate(Region = if_else(FADPU668_exclude.info > 0,"Exclude","Include"))
Region_FADPU668 <-FADPU668.info[,2]

# remove microglia contained within polygon
FADPU668 <- AddMetaData(FADPU668, metadata = Region_FADPU668, col.name = "Region")
FADPU668.include <- subset(FADPU668, Region == c("Include"))

# merge microglia from FADPU668 back with other samples
microglia.8sample.unclean <- merge(FADPU668.include,y=c(Other_Samples))
microglia.8sample.unclean <- JoinLayers(microglia.8sample.unclean)

#### subclustering of microglia - round 1
microglia.8sample.unclean <-NormalizeData(microglia.8sample.unclean)
microglia.8sample.unclean <- FindVariableFeatures(microglia.8sample.unclean, selection.method = "vst", nfeatures = 398)
all.genes <- rownames(microglia.8sample.unclean)
microglia.8sample.unclean <- ScaleData(microglia.8sample.unclean, features = all.genes)
microglia.8sample.unclean <- RunPCA(microglia.8sample.unclean, features = VariableFeatures(object = microglia.8sample.unclean),npcs = 50)
microglia.8sample.unclean <- JackStraw(microglia.8sample.unclean, num.replicate = 100,dims=40)
microglia.8sample.unclean <- ScoreJackStraw(microglia.8sample.unclean, dims = 1:40)
JackStrawPlot(microglia.8sample.unclean, dims = 1:40)
microglia.8sample.unclean <- FindNeighbors(microglia.8sample.unclean, dims = 1:29)
microglia.8sample.unclean <- FindClusters(microglia.8sample.unclean, resolution = 2.8)
microglia.8sample.unclean <- RunUMAP(microglia.8sample.unclean, dims = 1:29)

### Removal of clusters containing border-associated macrophages or hybrids - round 1
microglia.8sample.unclean2 <- subset(x = microglia.8sample.unclean, idents =  c("0","1","2","3","4","7","9","10","12","14","15","19","23"))

#subclustering of microglia - round 2
microglia.8sample.unclean2 <-NormalizeData(microglia.8sample.unclean2)
microglia.8sample.unclean2 <- FindVariableFeatures(microglia.8sample.unclean2, selection.method = "vst", nfeatures = 398)
all.genes <- rownames(microglia.8sample.unclean2)
microglia.8sample.unclean2 <- ScaleData(microglia.8sample.unclean2, features = all.genes)
microglia.8sample.unclean2 <- RunPCA(microglia.8sample.unclean2, features = VariableFeatures(object = microglia.8sample.unclean2),npcs = 50)
microglia.8sample.unclean2 <- JackStraw(microglia.8sample.unclean2, num.replicate = 100,dims=30)
microglia.8sample.unclean2 <- ScoreJackStraw(microglia.8sample.unclean2, dims = 1:30)
JackStrawPlot(microglia.8sample.unclean2, dims = 1:30)
microglia.8sample.unclean2 <- FindNeighbors(microglia.8sample.unclean2, dims = 1:18)
microglia.8sample.unclean2 <- FindClusters(microglia.8sample.unclean2, resolution = 2.8)
microglia.8sample.unclean2 <- RunUMAP(microglia.8sample.unclean2, dims = 1:18)

##Removal of clusters containing border-associated macrophages or hybrids - round 2
microglia.8sample.clean <- subset(x = microglia.8sample.unclean2, idents =  c("0","1","2","3","4","5","6","7","8","9","10","11","12","14","16","17","19","20","21","22","23","24","25","26","27","28","29","30"))
microglia.8sample.clean <-NormalizeData(microglia.8sample.clean)
microglia.8sample.clean <- FindVariableFeatures(microglia.8sample.clean, selection.method = "vst", nfeatures = 398)
all.genes <- rownames(microglia.8sample.clean)
microglia.8sample.clean <- ScaleData(microglia.8sample.clean, features = all.genes)
microglia.8sample.clean <- RunPCA(microglia.8sample.clean, features = VariableFeatures(object = microglia.8sample.clean),npcs = 50)
microglia.8sample.clean <- JackStraw(microglia.8sample.clean, num.replicate = 100,dims=25)
microglia.8sample.clean <- ScoreJackStraw(microglia.8sample.clean, dims = 1:25)
JackStrawPlot(microglia.8sample.clean, dims = 1:25)
microglia.8sample.clean <- FindNeighbors(microglia.8sample.clean, dims = 1:17)
microglia.8sample.clean <- FindClusters(microglia.8sample.clean, resolution = 0.7)
microglia.8sample.clean <- RunUMAP(microglia.8sample.clean, dims = 1:17)

##data output
saveRDS(microglia.8sample.clean,"Microglia_all samples_no hybrids.rds")

###### Analysis of microglia-plaque nearest neighbor distance (NND) ######
microglia.8sample.clean <- readRDS("250216_8sample_microglia_subclusters_clean.rds")

###analysis of FADPU538
FADPU538 <- subset(microglia.8sample.clean, subset = (sample == "FADPU538"))

#define microglia contained within the cortex
FADPU538.x <- FADPU538@meta.data[,5]
FADPU538.y <- FADPU538@meta.data[,6]
FADPU538_cx_x <- c(3573, 3893, 3770, 4161, 4493, 4814, 5080, 5195, 5219, 5174, 5084, 5039, 5035, 5068, 5019, 4959, 4757, 4611, 4439, 4016, 3700, 3758, 3936, 4065, 4446, 6898, 7171, 4861, 3609, 3230, 3089, 2980, 3054, 3437)
FADPU538_cx_y <- c(6953, 7315, 7350, 7694, 7239, 6875, 6416, 5952, 5737, 5669, 4865, 4684, 4037, 3114, 2344, 2052, 1744, 1599, 1593, 1667, 1659, 1138, 870, 563, 183, 152, 6813, 9488, 8418, 7844, 7618, 7275, 7022, 7020)
FADPU538_pol_cx <- point.in.polygon(FADPU538.x, FADPU538.y, FADPU538_cx_x, FADPU538_cx_y)
FADPU538.info <- cbind.data.frame(FADPU538_pol_cx)
FADPU538.info <- data.frame(FADPU538.info)
FADPU538@meta.data <- cbind(FADPU538@meta.data,FADPU538.info)
FADPU538_cx <- subset(FADPU538,subset = FADPU538_pol_cx == 1)
FADPU538_rest <- subset(FADPU538,subset = FADPU538_pol_cx == 0)
FADPU538_cx$region <- "Cx"
FADPU538_rest$region <- "rest"
Idents(FADPU538_cx)<-FADPU538_cx$region
levels(FADPU538_cx)
FADPU538_all <- merge(FADPU538_cx, y = c( FADPU538_rest))
FADPU538_all <- JoinLayers(FADPU538_all)

## output csv files containing X and Y coordinates of microglia contained within the cortex for NND analysis in python
Idents(FADPU538_cx) <- FADPU538_cx$sample
table(Idents(FADPU538_cx))
output_dataframe=data.frame(FADPU538_cx@active.ident)
output_dataframe=cbind(output_dataframe,FADPU538_cx@meta.data$center.x)
output_dataframe=cbind(output_dataframe,FADPU538_cx@meta.data$center.y)
output_dataframe <- cbind(rownames(output_dataframe), data.frame(output_dataframe, row.names=NULL))
colnames(output_dataframe) <- c('cell', 'cluster', 'x', 'y')
write.csv(output_dataframe, 'all_cells.csv', row.names=FALSE)

###run NND analysis using "python find_nearest_neighbors2.py"

###read output file of "python find_nearest_neighbors2.py" back into Rstudio and attach as metadata to cortical microglia
NND<-read.csv("FADPU538_plaque_nearest_neighbors.csv",header=FALSE)
NND<- NND[,3]
NND <- data.frame(NND)
FADPU538_cx@meta.data <- cbind(FADPU538_cx@meta.data,NND)

##Annotate cortical microglia as plaque-associated (NND < 15) or distal (NND >= 15)
FADPU538_cx$JP <- FADPU538_cx$NND < 15
Idents(FADPU538_cx) <- "JP"

### Repeat for FADTV348
FADTV348 <- subset(x = microglia.8sample.clean, subset = (sample == "FADTV348"))
FADTV348.x <- FADTV348@meta.data[,5]
FADTV348.y <- FADTV348@meta.data[,6]
FADTV348_cx_x <- c(1802, 1866, 2059, 2254, 2370, 2452, 2415, 2248, 1969, 1077, 266, -81, 919, 1247, 2372, 3139, 3193, 3195, 3271, 3391, 3435, 2720, 2535, 2237, 2048, 1833, 1548, 1417, 1368, 1367, 1455, 1755)
FADTV348_cx_y <- c(4103, 4382, 4567, 4659, 4719, 4702, 4121, 3555, 2987, 3041, 4301, 7507, 10339, 11143, 11613, 11296, 11104, 10952, 10719, 10394, 9823, 9719, 9635, 9060, 8425, 7613, 6749, 6320, 5915, 5485, 5098, 4069)
FADTV348_pol_cx <- point.in.polygon(FADTV348.x, FADTV348.y, FADTV348_cx_x, FADTV348_cx_y)
FADTV348.info <- cbind.data.frame(FADTV348_pol_cx)
FADTV348.info <- data.frame(FADTV348.info)
FADTV348@meta.data <- cbind(FADTV348@meta.data,FADTV348.info)
FADTV348_cx <- subset(FADTV348,subset = FADTV348_pol_cx == 1)
FADTV348_rest <- subset(FADTV348,subset = FADTV348_pol_cx == 0)
FADTV348_cx$region <- "Cx"
FADTV348_rest$region <- "rest"
FADTV348_all <- merge(FADTV348_cx, y = c( FADTV348_rest))
FADTV348_all <- JoinLayers(FADTV348_all)
Idents(FADTV348_cx) <- FADTV348_cx$sample
table(Idents(FADTV348_cx))
output_dataframe=data.frame(FADTV348_cx@active.ident)
output_dataframe=cbind(output_dataframe,FADTV348_cx@meta.data$center.x)
output_dataframe=cbind(output_dataframe,FADTV348_cx@meta.data$center.y)
output_dataframe <- cbind(rownames(output_dataframe), data.frame(output_dataframe, row.names=NULL))
colnames(output_dataframe) <- c('cell', 'cluster', 'x', 'y')
write.csv(output_dataframe, 'all_cells.csv', row.names=FALSE)

NND<-read.csv("FADTV348_Plaques_nearest_neighbors.csv",header=FALSE)
NND<- NND[,3]
NND <- data.frame(NND)
FADTV348_cx@meta.data <- cbind(FADTV348_cx@meta.data,NND)
FADTV348_cx$JP <- FADTV348_cx$NND < 15
Idents(FADTV348_cx) <- "JP"

## merge cortical microglia from FADTV348 and FADPU538
Microglia_cx_5xFAD <- merge(FADTV348_cx, FADPU538_cx)
Microglia_cx_5xFAD <- JoinLayers(Microglia_cx_5xFAD)

##Data Output
saveRDS(Microglia_cx_5xFAD, "Microglia_5xFAD cortex_with plaque NND")

