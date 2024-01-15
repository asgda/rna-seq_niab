#!/bin/bash

#SBATCH -N 2
#SBATCH --ntasks-per-node=20
#SBATCH --time=12:00:00
#SBATCH --job-name=feature_counts
#SBATCH --error=job.%J.err_node_20
#SBATCH --output=job.%J.out_node_20
#SBATCH --partition=standard

cores=32
bam_files=../STAR_run/*out.bam
counts_out=counts_restricted_niab.csv
gtf=/scratch/avik.bio.iith/gencode_gtf_files/sheep/hgdownload.soe.ucsc.edu/goldenPath/oviAri4/bigZips/genes/oviAri4.ncbiRefSeq.gtf

echo ==== Start ====
date
    # -M multimapped, --fraction count read fractions, -O overlapping reads, -p paired-end reads
    featureCounts \
    -t exon \
    -g gene_name \
    -a ${gtf} \
    -p \
    -o ${counts_out} \
    -M \
    -O \
    --fraction \
    -T ${cores} \
    ${bam_files}

##remove first column from the file

awk 'NR>1' ${counts_out} > temp.txt && mv temp.txt ${counts_out}
