# Pyrocystis lunula Transcriptome Analysis
By: Lila L. Barrera & Reis M. Gadsden

## Table of Contents
1. [Introductions](#1)
2. [File Structure](#2)
3. [Cleaning Raw Data](#3)
4. [Transcriptome Assembly](#4)
5. [Abundance Estimation](#5)
6. [Clustered Assembly](#6)
7. [Annotation](#7)
	a. [Functional](#7a)
	b. [Taxonomic](#7b)
8. [DESeq2](#8)
9. [Pathview](#9)

## Introduction<a name="1"></a>

Plan: Mini explanation of sections + script

### File Structure<a name="2"></a>
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

## Cleaning Raw Data <a name="3"></a>

The raw data was cleaned on the Galaxy platform using the trimmomatic tool. Initially the data was cleaned using trimmomatic and FastQC on Longleaf using the following script:

```
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
```

However the data from the galaxy platform was choosen to be used in the following steps.

## Transcriptome Assembly<a name="4"></a>
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

## Abundance Estimation<a name="5"></a>
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
## Clustered Assembly<a name="6"></a>

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

## Annotation<a name="7"></a>
In annotating our clustered assembly, we ran into many challenges and issues. These will be discussed more in a later section. We ended up using two different databases to annotate our information as we wanted both functional and taxonomic annotations. The two databases that we used to annotate our data were Kegg and PhyloDB.
### Functional<a name="7a"></a>
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

KeggAnnot was then ran on the output files using the following script:

```
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

```

### Taxonomic<a name="7b"></a>
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
diamond makedb --in /path/to/Phylo/Database -d phylo
diamond blastx -d phylo \
	-q clustered_assembly.fa \
	-o clustered_assembly_phylo_annotation.m8 \
	-p 12 -e 0.000001 -k 1
```

<a href="https://github.com/Lswhiteh/phylodbannotation">Fastannotation.py</a> was then ran on the output files using the following script:

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

module load python3

python3 phylodbannotation/fastannotation.py clustered_assembly.fa clustered_assembly_phylo_annotation.m8 phylo_annotated.fa output.tsvs
```

## DESeq<a name="8"></a>

## KEGG Pathview<a name="9"></a>
During our analysis we used the bioconductor package in R, <a href="https://bioconductor.org/packages/release/bioc/html/pathview.html">Pathview</a>, in order to visualize the up- and down-regulation of proteins within certain biological pathways. Some of the significant pathways that we choose to model were:
* <a href ="https://www.genome.jp/pathway/map00195">Photosynthesis</a>
<p align="center">
	<img src ="/images/0907_ko00195.phegg.pathview.png" alt="KEGG pathway of photosynthesis comparing *P.lunula* in the Dark phase vs the Light phase. Image rendered using Pathview package in Rstudio" width=75%>
</p>

* <a href="https://www.genome.jp/entry/ko04710">Circadian Rhythm</a>
<p align="center">
	<img src ="/images/0922_ko04710.phegg.pathview.png" alt="KEGG pathway of fluid shear stress comparing *P.lunula* in the Dark phase vs the Light phase. Image rendered using Pathview package in Rstudio" width=75%>
</p>

* <a href="https://www.genome.jp/dbget-bin/www_bget?pathway:map05418">Fluid Shear Stress</a>
<p align="center">
	<img src ="/images/0922_ko05418.phegg.pathview.png" alt="KEGG pathway of fluid shear stress comparing *P.lunula* in the Dark phase vs the Light phase. Image rendered using Pathview package in Rstudio" width=75%>
</p>


The R script that was used in order to generate these visuals is provided below:
```
##################################################################
#                                                                #
#                       pathview_script.R                        #
#                                                                #
#                                                                #
# This script first cleans our input data then precedes to build #
# the molecule simulation data needed for the Pathview package.  #
#                                                                #
# author: Reis Gadsden                                           #
#                                                                #
##################################################################

# install pathview
BiocManager::install("pathview")
library(pathview) 

# YOU WILL NEED TO UPDATE THIS TO BE THE CORRECT PATH
path_to_csv <- "~/lunula/PheggdDESeq_KO.csv"

# load the csv
phegg <- read.csv(path_to_csv)

# cast our data frame to a matrix
phegg.matrix <- as.matrix(phegg)

# set the row names of our matrix to the KO number
rownames(phegg.matrix) <- substring(phegg.matrix[, 1], 2)

# create a list that unjankifies our rownames (yeah this is dumb but it doesn't work without)
phegg.matrix.data <- sim.mol.data(mol.type="gene.ko", nmol = 5000)

# generate pathview graph
phegg.matrix.pathview <- pathview(gene.data = phegg.matrix.data, pathway.id = "05418", species = "ko", out.suffix = "phegg.pathview", gene.idtype = "kegg", kegg.native = TRUE, low = list(gene="#FFC02A", cpd = "yellow"), high = list(gene = "#0C7BDC", cpd = "blue"))
```