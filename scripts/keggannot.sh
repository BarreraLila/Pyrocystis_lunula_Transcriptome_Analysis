#!/bin/bash

#BATCH -p general
#SBATCH --nodes=1
#SBATCH --time=0-5:00:00
#SBATCH --mem=300G
#SBATCH --ntasks=24
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=lilalbar@ad.unc.edu
#SBATCH -J keggannot_k18
#SBATCH -o keggannot_k18.%A.out
#SBATCH -e keggannot_k18.%A.err

cd /pine/scr/l/i/lilalbar/p.lunula_fastq/github/trinity_assembly/Annotation
#git clone https://github.com/ctberthiaume/keggannot.git
#git clone https://github.com/AlexanderLabWHOI/keggannot.git
cd keggannot/
module purge

#module load python/3.5.1
module load python/2.7.12 

#python2 setup.py install --user
#emailed marchetti github user and he suggest to change from /nas/longleaf/data/KEGG/KEGG to /nas/longleaf/data/KEGG/2018.11.15
./bin/keggannot_genes2ko -m /nas/longleaf/data/KEGG/2018.11.15 /pine/scr/l/i/lilalbar/p.lunula_fastq/github/trinity_assembly/Annotation/clustered_assembly_kegg_annotation.m8 > ../annot_retry_kegg2018.11.15/annotated_kegg2018_information.tsv
