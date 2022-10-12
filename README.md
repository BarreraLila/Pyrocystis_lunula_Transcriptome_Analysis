# Pyrocystis_lunula_Transcriptome_Analysis

## Introduction

Plan: Mini explanation of sections + script

### File Structure
```
P. Lunula
|- Github
|	|- Trinity_Assembly
|	|	|- Annontation
|	|	|- Assembly		
|	|	|- PhyloDB
|	|	|- salmon_quants
|	|	|	|- A2_quants
|	|	|	|- A3_quants
|	|	|	|- A4_quants
|	|	|	|- A5_quants
|	|	|	|- B2_quants
|	|	|	|- B3_quants
|	|	|	|- B4_quants
|	|	|	|- B5_quants
|	|	|- trinity_assembly_output	
|	| - Scripts
|	|	|- Job Fails
|- PhyloDB
|- PhyloDBAnnotation
|- Data
|	|- A2
|	|	|- R1
|	|	|- R2
|	|- A3
|	|	|- R1
|	|	|- R2
|	|- A4
|	|	|- R1
|	|	|- R2
|	|- A5
|	|	|- R1
|	|	|- R2
|	|- B2
|	|	|- R1
|	|	|- R2
|	|- B3
|	|	|- R1
|	|	|- R2
|	|- B4
|	|	|- R1
|	|	|- R2
|	|- B5
|	|	|- R1
|	|	|- R2
```

## Cleaning Raw Data 

- stats?

## Transcriptome Assembly
Now that the data has been cleaned up, we can go about reassembling our transcripts from our data using Trinity. Trinity is a modular software package that reconstructs plausible transcripts from raw RNASeq data. Trinity uses k-mer and de Bruijn graphs and sequences to reconstruct these transcripts. Trinity will output a FASTA file containing all full length transcripts as well as partial length isoforms. Our Trinity assembly script is below:

```
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

Trinity --seqType fq --max_memory 500G \
	--left "${sampA2}/R1/A2_R1.gz","${sampA3}/R1/A3_R1.gz",\
"${sampA4}/R1/A4_R1.gz","${sampA5}/R1/A5_R1.gz",\
"${sampB2}/R1/B2_R1.gz","${sampB3}/R1/B3_R1.gz",\
"${sampB4}/R1/B4_R1.gz","${sampB5}/R1/B5_R1.gz" \
	--right "${sampA2}/R2/A2_R2.gz","${sampA3}/R2/A3_R2.gz",\
"${sampA4}/R2/A4_R2.gz","${sampA5}/R2/A5_R2.gz",\
"${sampB2}/R2/B2_R2.gz","${sampB3}/R2/B3_R2.gz",\
"${sampB4}/R2/B4_R2.gz","${sampB5}/R2/B5_R2.gz" \
	--CPU 24 \
	--output /pine/scr/l/i/lilalbar/p.lunula_fastq/github/trinity_assembly
```
We start by created variable for each of our samples for simplicity and less confusion when writing the actual Trinity call. We then pass all our left read samples into the Trinity call using the <i>--left</i> flag. We then repeat the process for our right reads using the <i>--right</i> flag. We then specify an output folder for our Trinity FASTA file, this file by default is named <i>Trinity.fasta</i>. This file contains all our transcripts and uses a custom naming convention, this convention is as follows:
* TRINITY_DNxxxx_cx_gx_ix

In this naming convention the DNxxxx_cx are referring to the read cluster to which a contig belongs to. The cx is contig number of that read cluster. The gx is the “gene” of the contig, the use of gene here does not correspond to an actual gene, but instead is just how Trinity decides to separate similar transcripts. The ix represents the isoform of that gene. An example of a line in a Trinity output file is provided below.

```
>TRINITY_DN30_c0_g1_i1 len=2677 path=[2675:0-1921 2676:1922-2676]
ACCATGA
>TRINITY_DN30_c0_g2_i1 len=1940 path=[2673:0-1921 2674:1922-1939]
ACGATGA
```

The FASTA file is formatted in such that the first line is the FASTA header containing the name of the transcript, the length of the transcript, as well as the paths taken to create that transcript. The next line contains the sequence for the transcript listed in the line above it.

## Abundance Estimation
In order to be able to downstream analysis of our data so far, we will also need to quantify the counts of our transcripts across all our samples. There were many different packages we explored to complete this step but the two packages we attempted to use were Salmon and RSEM. This section will focus on our Salmon methodology as this was the package we used to quantify our transcripts moving forward. For Salmon we first must create an index with our clustered assembly, this only needs to be done once and can be reused for all samples. After this we must run a salmon quant of each of our samples, we did all our samples in one bash file for efficiency. Our bash file is provided below.

```
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

```
## Clustered Assembly

- uses - KAAS, InterproScan and BLAST, phylo and kegg
Abundance Estimations

After completing our Trinity assembly it is a good idea to cluster similar contigs together. This is to help remove redundancy and will improve our alignment by removing the chance of multi-mapping due to poorly reconstructed transcripts. Clustering can easily be completed by using CD-HIT-EST from the CD-HIT suite. Our clustered assembly creation script is below.
```
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
```
Here we are clustering contigs based on whether they have a 98% similarity with each other. We went with the standard option of choosing 10 for our word size. For our description we allowed a length of 100. 500Gb of memory was specified in the SLURM headers so we allowed the clustering to utilize all of that memory if needed. In a similar manner we specified 12 tasks so we allow our clustering to utilize up to 12 threads.

## Annotation
In annotating our clustered assembly, we ran into many challenges and issues. These will be discussed more in a later section. We ended up using two different databases to annotate our information as we wanted both functional and taxonomic annotations. The two databases that we used to annotate our data were Kegg and PhyloDB.
###Functional
To start our Kegg annotation we needed to blast our sequences against the database using diamond. We used Kegg files provided on the Longleaf server to build our database to annotate from. To get these annotations we used a modified version of the program <a href ="https://github.com/ctberthiaume/keggannot">KeggAnnot</a>. The modified version can be found <a href ="https://github.com/reismgadsden/keggannot">here</a>. The modifications simply let us skip over unresolved pathways and allow the program to run to completion. The output of this step was a annotated fasta file as well as a tsv file containing the annotation information for each of the reconstructed contigs. The script for constructing this database is provided below:

```
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
#blastx Align translated DNA query sequences against a protein reference database.
diamond makedb --in /nas/longleaf/data/KEGG/KEGG/genes/fasta/genes.pep.fasta -d keggdb
diamond blastx -d keggdb \
	-q clustered_assembly.fa \
	-o clustered_assembly_kegg_annotation.m8 \
	-p 12 -e 0.000001 -k 1

```

###Taxonomic
Our taxonomic annotation used PhyloDB, a marine life database used for metagenomics. Similar to our Kegg steps we used to diamond to construct a database using files provided to us. We then used the <a href ="https://github.com/Lswhiteh/phylodbannotation">PhyloDB Annotation</a> program that was designed by Logan Whitehouse in order to annotate our transcripts.

```
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
#Need to blast against PhyloDB
#blastx Align translated DNA query sequences against a protein reference database.
diamond makedb --in /path/to/Phylo/Database -d keggdb
diamond blastx -d keggdb \
	-q clustered_assembly.fa \
	-o clustered_assembly_kegg_annotation.m8 \
	-p 12 -e 0.000001 -k 1
```

## DESeq
## BLAST

[link](https://github.com/BarreraLila/Pyrocystis_lunula_Transcriptome_Analysis)

![KEGG pathway of photosynthesis comparing *P.lunula* in the Dark phase vs the Light phase. Image rendered using Pathview package in Rstudio](/images/0907_ko00195.phegg.pathview.png)