#!/bin/bash

#SBATCH -p general
#SBATCH --nodes=1

# only asking for 1 node 

#SBATCH --time=0-16:00:00
#SBATCH --mem=90G
#SBATCH --ntasks=16
#SBATCH -J clean_align
#SBATCH -o cleanalign.%A.out
#SBATCH -e cleanalign.%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lilalbar@unc.edu

# load your required tools for QC and trimming
module load trimmomatic
module load fastqc

# now cd into your dir with the files
cd /pine/scr/l/i/lilalbar/p.lunula_fastq/

# here we are creating an output fastq file on stats before we trim and saves - 
# visualize if not provided in the file

fastqc -t 8 --outdir ../data/before_trimming/

# here it looks like we are taking all forward (R1) reads and trimmming them

# 'for f in ' looks like a loop - "for files in" maybe 
# * means zero or more characters
# *R1* we're specifing for R1
# ls /pine/scr/l/i/lilalbar/p.lunula_fastq/data/*/*R1* ls the forward reads 
# of each sample file in the data dir on the scratch space

for f in `ls /pine/scr/l/i/lilalbar/p.lunula_fastq/data/*/*R1* | cut -f1,2,3 -d'_'`
do
        echo $f
        trimmomatic PE -threads 8 \
                ${f}_R1.fq ${f}_R2.fq \
                ${f}_fwd_paired_trimmed.fq.gz ${f}_fwd_unpaired_trimmed.fq.gz \
                ${f}_rev_paired_trimmed.fq.gz ${f}_rev_unpaired_trimmed.fq.gz \
                #ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads \
                LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:36

done

# this gives fastq output file after trimming and saves it
fastqc -t 8 --outdir ../data/after_trimming/