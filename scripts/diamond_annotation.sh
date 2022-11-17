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

#https://github.com/bbuchfink/diamond/wiki/3.-Command-line-options
module load diamond

cd /pine/scr/l/i/lilalbar/p.lunula_fastq/github/trinity_assembly/Annotation

#Need to blast against KEGG

#--in stdin if parameter is omitted
#/pine/scr/l/i/lilalbar/p.lunula_fastq/github/trinity_assembly/Annotation
#blastx Align translated DNA query sequences against a protein reference database.

diamond makedb --in /nas/longleaf/data/KEGG/KEGG/genes/fasta/genes.pep.fasta -d keggdb
diamond blastx -d keggdb \
	-q clustered_assembly.fa \
	-o clustered_assembly_kegg_annotation.m8 \
	-p 12 -e 0.000001 -k 1
