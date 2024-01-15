#!/usr/bin/Rscript

## data preparation 

t1<- read.table("counts_restricted_niab.csv", sep='\t', header=T, check.names=F)
t1<- t1[,c(1,6:ncol(t1))]
col1<- colnames(t1)
new_col1<- gsub("../STAR_run/(.*)Aligned.sortedByCoord.out.bam", "\\1",col1)
names(t1)<- new_col1

## round the values

t1[,-1]<- round(t1[,-1])
write.table(t1, "raw_counts_feature_counts.csv", sep=',', row.names=F)

## RAW COUNTS TO TPM
a<- read.table("raw_counts_feature_counts.csv", header=T, check.names=F, sep=',')
len<- a$Length
a<- a[,-2]
rownames(a)<- a$Geneid
a<- a[,-1]
x <- a / len
tpm.mat <- t( t(x) * 1e6 / colSums(x) )
tpm.mat <- data.frame(tpm.mat)
write.table(tpm.mat, "tpm_counts.csv", row.names=T, col.names=NA)

## TPM TO Z_SCORE
log_df<- log(tpm.mat+1)
medians<- apply(log_df,1,median)
std_dev<- apply(log_df,1,sd)
z_score<- (log_df - medians)/(std_dev)
z_score<- na.omit(z_score)
write.table(z_score, "z_score_df.csv", row.names=T, col.names=NA)


