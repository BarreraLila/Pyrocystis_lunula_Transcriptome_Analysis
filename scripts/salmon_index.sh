#!/bin/bash

#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --time=0-2:00:00
#SBATCH --mem=200G
#SBATCH --ntasks=24
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lilalbar@ad.unc.edu
#SBATCH -J salmon_index
#SBATCH -o salmon_index.%A.out
#SBATCH -e salmon_index.%A.err

module load salmon

cd /pine/scr/l/i/lilalbar/p.lunula_fastq/github/trinity_assembly/Annotation

salmon index -i ../Alignment_Index -t clustered_assembly.fa
