#!/bin/bash

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

