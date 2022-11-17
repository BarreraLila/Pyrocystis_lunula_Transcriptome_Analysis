#!/bin/bash

#BATCH -p general
#SBATCH --nodes=1
#SBATCH --time=0-10:00:00
#SBATCH --mem=300G
#SBATCH --ntasks=24
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lilalbar@ad.unc.edu
#SBATCH -J salmon_quant
#SBATCH -o salmon_quant.%A.out
#SBATCH -e salmon_quant.%A.err

module load salmon

cd /pine/scr/l/i/lilalbar/p.lunula_fastq/github/trinity_assembly

sampA2="/pine/scr/l/i/lilalbar/p.lunula_fastq/data/trimmed_paired/A2"
sampA3="/pine/scr/l/i/lilalbar/p.lunula_fastq/data/trimmed_paired/A3"
sampA4="/pine/scr/l/i/lilalbar/p.lunula_fastq/data/trimmed_paired/A4"
sampA5="/pine/scr/l/i/lilalbar/p.lunula_fastq/data/trimmed_paired/A5"
sampB2="/pine/scr/l/i/lilalbar/p.lunula_fastq/data/trimmed_paired/B2"
sampB3="/pine/scr/l/i/lilalbar/p.lunula_fastq/data/trimmed_paired/B3"
sampB4="/pine/scr/l/i/lilalbar/p.lunula_fastq/data/trimmed_paired/B4"
sampB5="/pine/scr/l/i/lilalbar/p.lunula_fastq/data/trimmed_paired/B5"

salmon quant -l A -i ./Alignment_Index \
	-1 ${sampA2}/R1/A2_R1.gz \
	-2 ${sampA2}/R2/A2_R2.gz \
	-p 24 \
	-o A2_quants

salmon quant -l A -i ./Alignment_Index \
	-1 ${sampA3}/R1/A3_R1.gz \
	-2 ${sampA3}/R2/A3_R2.gz \
	-p 24 \
	-o A3_quants

salmon quant -l A -i ./Alignment_Index \
	-1 ${sampA4}/R1/A4_R1.gz \
	-2 ${sampA4}/R2/A4_R2.gz \
	-p 24 \
	-o A4_quants

salmon quant -l A -i ./Alignment_Index \
	-1 ${sampA5}/R1/A5_R1.gz \
	-2 ${sampA5}/R2/A5_R2.gz \
	-p 24 \
	-o A5_quants

salmon quant -l A -i ./Alignment_Index \
	-1 ${sampB2}/R1/B2_R1.gz \
	-2 ${sampB2}/R2/B2_R2.gz \
	-p 24 \
	-o B2_quants

salmon quant -l A -i ./Alignment_Index \
	-1 ${sampB3}/R1/B3_R1.gz \
	-2 ${sampB3}/R2/B3_R2.gz \
	-p 24 \
	-o B3_quants

salmon quant -l A -i ./Alignment_Index \
	-1 ${sampB4}/R1/B4_R1.gz \
	-2 ${sampB4}/R2/B4_R2.gz \
	-p 24 \
	-o B4_quants

salmon quant -l A -i ./Alignment_Index \
	-1 ${sampB5}/R1/B5_R1.gz \
	-2 ${sampB5}/R2/B5_R2.gz \
	-p 24 \
	-o B5_quants
