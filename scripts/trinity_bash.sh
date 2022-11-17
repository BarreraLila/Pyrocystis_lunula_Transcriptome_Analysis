#!/bin/bash

#SBATCH --nodes=1
#SBATCH --time=08-00:00:00
#SBATCH --mem=500G
#SBATCH --ntasks=24
#SBATCH -J trinity
#SBATCH -o trinity.%A.out
#SBATCH -e trinity.%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lilalbar@ad.unc.edu

## The output of Trinity is a fasta file that contains contigs 
## (contiguous sequences of DNA/RNA) it was able to assemble.

##Change the below - figure out what files we are using
cd /pine/scr/l/i/lilalbar/p.lunula_fastq/data/trimmed_paired

module load trinity

#Easy way to condense filenames is to make variables that contain the path to the files
#That way you can organize better and still easily use them in scripts
sampA2="/pine/scr/l/i/lilalbar/p.lunula_fastq/data/trimmed_paired/A2"
sampA3="/pine/scr/l/i/lilalbar/p.lunula_fastq/data/trimmed_paired/A3"
sampA4="/pine/scr/l/i/lilalbar/p.lunula_fastq/data/trimmed_paired/A4"
sampA5="/pine/scr/l/i/lilalbar/p.lunula_fastq/data/trimmed_paired/A5"
sampB2="/pine/scr/l/i/lilalbar/p.lunula_fastq/data/trimmed_paired/B2"
sampB3="/pine/scr/l/i/lilalbar/p.lunula_fastq/data/trimmed_paired/B3"
sampB4="/pine/scr/l/i/lilalbar/p.lunula_fastq/data/trimmed_paired/B4"
sampB5="/pine/scr/l/i/lilalbar/p.lunula_fastq/data/trimmed_paired/B5"



Trinity --seqType fq \
	--max_memory 500G \
	--left "${sampA2}/R1/A2_R1.gz","${sampA3}/R1/A3_R1.gz","${sampA4}/R1/A4_R1.gz","${sampA5}/R1/A5_R1.gz","${sampB2}/R1/B2_R1.gz","${sampB3}/R1/B3_R1.gz","${sampB4}/R1/B4_R1.gz","${sampB5}/R1/B5_R1.gz" \
	--right "${sampA2}/R2/A2_R2.gz","${sampA3}/R2/A3_R2.gz","${sampA4}/R2/A4_R2.gz","${sampA5}/R2/A5_R2.gz","${sampB2}/R2/B2_R2.gz","${sampB3}/R2/B3_R2.gz","${sampB4}/R2/B4_R2.gz","${sampB5}/R2/B5_R2.gz" \
	--CPU 24 \
	--output /pine/scr/l/i/lilalbar/p.lunula_fastq/github/trinity_assembly
