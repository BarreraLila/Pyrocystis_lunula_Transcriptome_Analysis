# Pyrocystis_lunula_Transcriptome_Analysis

## Introduction

Plan: Mini explanation of sections + script

## Cleaning Raw Data 

- stats?

## Transcriptome Assembly

## Clustered Assembly

- uses - KAAS, InterproScan and BLAST, phylo and kegg
Abundance Estimations

After completing our Trinity assembly it is a good idea to cluster similar contigs together. This is to help remove redundancy and will improve our alignment by removing the chance of multi-mapping due to poorly reconstructed transcripts. Clustering can easily be completed by using CD-HIT-EST from the CD-HIT suite. Our clustered assembly creation script is below.
<code>#!/bin/bash
#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --time=0-10:00:00
#SBATCH --mem=500G
#SBATCH --ntasks=12
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lilalbar@ad.unc.edu
#SBATCH -J cdhitremoveredundancy
#SBATCH -o cdhitremoveredundanc.%A.out
#SBATCH -e cdhitremoveredundancy.%A.err

module load cdhit

cd /pine/scr/l/i/lilalbar/p.lunula_fastq/github/trinity_assembly

cd-hit-est -i Trinity.fasta -o clustered_assembly.fa -c 0.98 -n 10 -d 100 -M 500000 -T 12
</code>
Here we are clustering contigs based on whether they have a 98% similarity with each other. We went with the standard option of choosing 10 for our word size. For our description we allowed a length of 100. 500Gb of memory was specified in the SLURM headers so we allowed the clustering to utilize all of that memory if needed. In a similar manner we specified 12 tasks so we allow our clustering to utilize up to 12 threads.

## Annotation
- Functional
- Taxonomic

## DESeq
## BLAST

[link]
(https://github.com/BarreraLila/Pyrocystis_lunula_Transcriptome_Analysis)

![KEGG pathway of photosynthesis comparing *P.lunula* in the Dark phase vs the Light phase. Image rendered using Pathview package in Rstudio](/images/0907_ko00195.phegg.pathview.png)