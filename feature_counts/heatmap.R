#!/usr/bin/Rscript
z_mat<- read.table("z_score_df.csv", header=T, row.names=1)
library(ComplexHeatmap)
dge<- read.table("../deseq2/differential_gene_exp_GSE213637_20.csv", header=T, sep=",", row.names=1)
z_mat[, 1:ncol(z_mat)][z_mat[, 1:ncol(z_mat)] == "NULL"] <- 0
is.na(z_mat)<- sapply(z_mat, is.infinite)
z_mat[is.na(z_mat)]<-0
z_mat<- merge(dge, z_mat, by="row.names")
z_mat<- na.omit(z_mat)
z_mat<- z_mat[order(z_mat$padj.BH),]
z_mat<- z_mat[z_mat$padj.BH < 0.05,]
z_mat<- z_mat[(z_mat$log2FoldChange >= 1)||(z_mat$log2FoldChange <= -1),]
row.names(z_mat)<- z_mat$Row.names
z_mat[,1]<- NULL
z_mat<- z_mat[,c(8:13)]
z_mat<- as.matrix(z_mat)
df <- data.frame(Condition = c(rep("Normal", 3), rep("Viral_infected", 3)))
ha = HeatmapAnnotation(df = df, col = list(Condition = c("Viral_infected" = "#DD292A", "Normal" = "#4FBFAD")), show_annotation_name = TRUE)

jpeg("heatmap_all_genes_0.01.jpeg", height=1200, width=1200, res=300)

Heatmap(z_mat, top_annotation = ha, show_row_names = FALSE, heatmap_legend_param = list(title="Z-score"))

dev.off()
