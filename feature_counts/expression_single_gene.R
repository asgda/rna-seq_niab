#!/usr/bin/Rscript
library(reshape2)
library(dplyr)
library(ggplot2)
a<- read.table("tpm_counts.csv", header=T, check.names=F)
names(a)[1]<- "Geneid"
b<- read.table("../metadata.csv", check.names=F, header=T, sep=',')
gene_exp<- data.frame(a[a$Geneid=="TRIM4",])
rownames(gene_exp)<- gene_exp[,1]
gene_exp <- gene_exp[,-1]
gene_exp<- t(gene_exp)
colnames(gene_exp)[1]<- "TPM_values"
condition <- data.frame(Condition= c(rep("Normal",3), rep("Viral_infected",3)))

gene_exp<- cbind(gene_exp, condition)

jpeg("TRIM4_expression_TPM.jpeg", width=1200, height= 1200, res=300)
p <- ggplot(gene_exp, aes(Condition, TPM_values)) + geom_boxplot() + geom_jitter(size=2, alpha=1/1.5, width = 0.2, aes(colour=Condition)) + theme(text = element_text(size=20)) + theme(legend.position="none") + theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + theme(axis.title.x=element_blank()) + scale_x_discrete(limits=c("Normal", "Viral_infected")) + labs(y = "TRIM4 expression (TPM)") + ylim(10, 50)

p_value <- t.test(TPM_values ~ Condition, data = gene_exp)$p.value  # Replace with your actual p-value calculation
p <- p + annotate("text", x = 2, y = 50, label = paste("p = ", format(p_value, digits = 3)), size = 3)

print(p)
dev.off()
