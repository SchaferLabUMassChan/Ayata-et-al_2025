find_nearest_neighbors2.py takes as input a csv filed named "all_cells.csv" that contains a list of cells, their cluster identity, and their X and Y positions. The script analyzes all cluster pairings and calculates the nearest neighbor distance for each cell in a cluster, for all pairings of clusters contained within all_cells.csv. Data is output to a folder named "python_nearest_neighbor_outputs" as .csv files for each cluster-cluster pairing.

In Ayata_et_al_2025, "all_cells.csv" is generated from the metadata of a Seurat object using RStudio:

###example R code to generate "all.cells.csv"### 
Idents(FADPU538_cx) <- FADPU538_cx$celltype
table(Idents(FADPU538_cx))
output_dataframe=data.frame(FADPU538_cx@active.ident)
output_dataframe=cbind(output_dataframe,FADPU538_cx@meta.data$center.x)
output_dataframe=cbind(output_dataframe,FADPU538_cx@meta.data$center.y)
output_dataframe <- cbind(rownames(output_dataframe), data.frame(output_dataframe, row.names=NULL))
colnames(output_dataframe) <- c('cell', 'cluster', 'x', 'y')
write.csv(output_dataframe, 'all_cells.csv', row.names=FALSE)
