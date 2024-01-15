#!/usr/bin/Rscript

library("ggplot2")
table<- read.table("STAR_results.csv", header=T, sep=',', check.names=F)

jpeg("unique_map.jpeg", res=300, width=1500, height=1500)
ggplot(table, aes(x = Sample, y = Uniquely_mapped_reads)) +
  geom_bar(stat = "identity", fill = "steelblue") +  # Use bars for visualization
  labs(title = "STAR Aligner Results", x = "Samples", y = "Uniquely mapped reads %") +  # Add labels and title
  theme_bw() +  # Use a clean background theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_hline(yintercept = 75, linetype = "dashed", color = "red")
dev.off()

jpeg("unmapped.jpeg", res=300, width=1500, height=1500)
ggplot(table, aes(x = Sample, y = `%_of_unmapped_reads`)) +
  geom_bar(stat = "identity", fill = "steelblue") +  # Use bars for visualization
  labs(title = "STAR Aligner Results", x = "Samples", y = "% of unmapped reads") +  # Add labels and title
  theme_bw() +  # Use a clean background theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

jpeg("mismatch.jpeg", res=300, width=1500, height=1500)
ggplot(table, aes(x = Sample, y =Mismatch_rate )) +
  geom_bar(stat = "identity", fill = "steelblue") +  # Use bars for visualization
  labs(title = "STAR Aligner Results", x = "Samples", y = "Mismatch Rate (%)") +  # Add labels and title
  theme_bw() +  # Use a clean background theme
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()
