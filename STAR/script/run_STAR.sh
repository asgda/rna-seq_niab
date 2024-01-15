#!/bin/bash
#SBATCH -N 4
#SBATCH --ntasks-per-node=20
#SBATCH --time=24:00:00
#SBATCH --job-name=star_sheep
#SBATCH --error=job.%J.err_node_20
#SBATCH --output=job.%J.out_node_20
#SBATCH --partition=standard

sample="../data/sample.txt"
while IFS= read -r line
do
	line=`echo ${line}`
	genomeDir="/scratch/avik.bio.iith/make_genome_STAR_sheep"
        fastq="../data"
        	echo $line[Start]`date`
                        	## STAR
                        	/scratch/avik.bio.iith/STAR-2.7.10a/source/STAR \
                        	--runThreadN 16 \
                        	--genomeDir $genomeDir \
                        	--readFilesIn ${fastq}/${line}_1.fastq ${fastq}/${line}_2.fastq \
                        	--readFilesCommand cat \
                        	--outSAMtype BAM SortedByCoordinate \
                        	--outFileNamePrefix ${line} \
                        	--outSAMunmapped Within \
                        	--outFilterType BySJout \
				--outFilterMultimapNmax 20 \
                        	--alignSJoverhangMin 8 \
                        	--alignSJDBoverhangMin 1 \
                        	--outFilterMismatchNmax 999 \
                        	--outFilterMismatchNoverLmax 0.04 \
                        	--twopassMode Basic \
                        	--outSAMattributes All \
                        	--outSAMstrandField intronMotif
        	echo $line[End]`date`
done <$sample
