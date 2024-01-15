#!/usr/bin/Rscript

#if (!require("BiocManager", quietly = TRUE))
 #   install.packages("BiocManager")

#BiocManager::install("annotate")
#BiocManager::install("geneplotter")
#BiocManager::install("genefilter")
#BiocManager::install("DESeq2")

library(XML)
library(annotate)
library(geneplotter)
library(genefilter)
library(DESeq2)
library(dplyr)
library(ggplot2)
library(reshape2)
#library(ggfortify)
#library(ComplexHeatmap)

countdata<- read.csv('../feature_counts/raw_counts_feature_counts.csv', header = T, check.names=F)
countdata<- countdata[!duplicated(countdata$Geneid),] #remove duplicates
row.names(countdata)<- countdata[,1]
countdata<- countdata[,-1]
metadata<- read.table('../metadata.csv', header=T, row.names = 1, check.names=F, sep=',')
countdata_1<- countdata[,-1]
dds<- DESeqDataSetFromMatrix(countData=countdata_1,
                              colData=metadata, 
                              design=~condition, tidy = F)

dds<- dds[rowSums(counts(dds)) > 20,] ##removing low read count genes
dds<- DESeq(dds)
res<- results(dds, contrast=c("condition","viral_infected", "normal"))
res$padj.BH <- p.adjust(res$pvalue, method="BH")
#res<- res[, c(1:5,7)]
res <- res[order(res$padj.BH),]
res <- res[(res$padj.BH< 0.05),]
summary(res)

write.table(res, 'differential_gene_exp_GSE213637_20.csv', sep=",", row.names=T, col.names=NA)

## volcano plot ##

df<- read.table("differential_gene_exp_GSE213637_20.csv", header=T, sep=",")
names(df)[names(df)== "X"]<- "Geneid"
jpeg("volcano_plot_GSE213637.jpeg", width=2000, height=2000, res=300)
with(df, plot(log2FoldChange, -log10(padj.BH), pch=20, ylim=c(0,310), xlim=c(-7,10), cex=1, col="red", xlab="Log2 fold change", ylab="-log10 Padj value", cex.lab=1.5, cex.axis=1.5))
with(subset(df, padj.BH >= 0.05 | (res$log2FoldChange < 1 & res$log2FoldChange > -1)), points(log2FoldChange, -log10(padj.BH), pch=20, col="gray", cex=1))
gene1<- "RNF125"
gene2<- "CTGF"
gene3<- "LOC105603123"
gene4<- "AMOT"
with(subset(df, Geneid %in% gene1), points(log2FoldChange, -log10(padj.BH), pch=20, col="blue", cex=1.5), text(log2FoldChange, -log10(padj.BH), labels = Geneid, cex=1.2, pos=1))
with(subset(df, Geneid %in% gene2), points(log2FoldChange, -log10(padj.BH), pch=20, col="green", cex=1.5))
with(subset(df, Geneid %in% gene3), points(log2FoldChange, -log10(padj.BH), pch=20, col="purple", cex=1.5))
with(subset(df, Geneid %in% gene4), points(log2FoldChange, -log10(padj.BH), pch=20, col="black", cex=1.5))
legend("topright", legend=c("RNF125", "CTGF", "LOC105603123","AMOT","Significant DEGs"), pch=20, col=c("blue", "green", "purple", "black", "red"), cex=1)
abline(h=1.30103, col = "blue",lty=2)
abline(v=1, col = "blue",lty=2)
abline(v=-1, col = "blue",lty=2)
dev.off()


## PCA plot ##

rld <- rlog(dds, blind=TRUE)

jpeg("pca_plot.jpeg", res=300, width=1500, height=1500)
pcaData <- plotPCA(rld, intgroup="condition", returnData=TRUE)
# Create the PCA plot using ggplot2
ggplot(pcaData, aes(x=PC1, y=PC2, color=condition, shape=condition)) +
  geom_point(size=4) +
  labs(title="PCA Plot of Samples by Condition",
       x=paste0("PC1 (", round(attr(pcaData, "percentVar")[1]*100, 2), "% variance)"),
       y=paste0("PC2 (", round(attr(pcaData, "percentVar")[2]*100, 2), "% variance)")) +
  scale_color_manual(values=c("#00AFBB", "#E7B800", "#FC4E07")) +  # Customize colors
  theme_bw() +
  theme(plot.title=element_text(size=12, face="bold"),
        axis.title=element_text(size=10),
        legend.title=element_text(size=10),
        legend.text=element_text(size=8))
dev.off()





