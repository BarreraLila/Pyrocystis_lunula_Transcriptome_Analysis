#!/bin/bash
#SBATCH -p general
#SBATCH --nodes=1
#SBATCH --time=0-24:00:00
#SBATCH --mem=300G
#SBATCH --ntasks=12
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lilalbar@ad.unc.edu
#SBATCH -J diamond
#SBATCH -o diamond.%A.out
#SBATCH -e diamond.%A.err

module load python3

python3 phylodbannotation/fastannotation.py clustered_assembly.fa clustered_assembly_phylo_annotation.m8 phylo_annotated.fa output.tsv