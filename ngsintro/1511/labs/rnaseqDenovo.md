---
layout: default
title:  'RNAseq'
---

# Transcriptome De Novo Assembly

## Trinity

Trinity is one of several de novo transcriptome assemblers. By efficiently constructing and analyzing sets of de Bruijn graphs, Trinity reconstructs a large fraction of transcripts, including alternatively spliced isoforms and transcripts from recently duplicated genes.
This approach provides a unified solution for transcriptome reconstruction in any sample, especially in the absence of a reference genome.

Grabherr MG, Haas BW, Yassour M et al. (2011) Full-length transcriptome assembly from RNA-Seq data without a reference genome.
Nature Biotechnology.
2011 May 15;29(7):644-52.

## Getting started

Trinity combines three independent software modules: Inchworm, Chrysalis, and Butterfly, applied sequentially to process large volumes of RNA-Seq reads.
Trinity partitions the sequence data into many individual de Bruijn graphs, each representing the transcriptional complexity at at a given gene or locus, and then processes each graph independently to extract full-length splicing isoforms and to tease apart transcripts derived from paralogous genes.
Briefly, the process works like so:

Inchworm assembles the RNA-Seq data into the unique sequences of transcripts, often generating full-length transcripts for a dominant isoform, but then reports just the unique portions of alternatively spliced transcripts.

Chrysalis clusters the Inchworm contigs into clusters and constructs complete de Bruijn graphs for each cluster.
Each cluster represents the full transcriptonal complexity for a given gene (or sets of genes that share sequences in common).
Chrysalis then partitions the full read set among these disjoint graphs.

Butterfly then processes the individual graphs in parallel, tracing the paths that reads and pairs of reads take within the graph, ultimately reporting full-length transcripts for alternatively spliced isoforms, and teasing apart transcripts that corresponds to paralogous genes.

A basic recommendation is to have 1G of RAM per 1M pairs of Illumina reads in order to run the Inchworm and Chrysalis steps.
Simpler transcriptomes require less memory than complex transcriptomes.
Butterfly requires less memory and can also be spread across multiple processors.

The entire process can require ~1 hour per million pairs of reads in the current implementation.
There are various things that can be done to modify performance.
Please review the guidelines in the official Trinity documentation for more advice on this topic.
Typical Trinity usage is as follows:

```bash
$ Trinity --seqType (fq for fastq or fa for fast) --left ~/path/to/reads_1.fq --right ~/path/to/reads_2.fq (or --single for single reads) --CPU 4 --output ~/path/to/output_dir
```

## Exercise 1: Running Trinity

Get yourself familiar with Trinity by having a look at the manual: http://trinityrnaseq.sourceforge.net/

Have a look at the example data used in this exercise.
The data is obtained from mouse dendritic cells (mouse_left.fasta and mouse_right.fasta and) and a whitefly (whitefly_both.fasta), and the files are located in `/proj/g2015045/labs/transcriptome_assembly/`.
The mouse data is strand-specific (RF), the whitefly data is unstranded.
For strand-specific data, specify the library type.
There are four library types:

Paired reads:  
RF: first read (/1) of fragment pair is sequenced as anti-sense (reverse(R)), and second read (/2) is in the sense strand (forward(F)); typical of the dUTP/UDG sequencing method.
FR: first read (/1) of fragment pair is sequenced as sense (forward), and second read (/2) is in the antisense strand (reverse)  

Unpaired (single) reads:  
F: the single read is in the sense (forward) orientation  
R: the single read is in the antisense (reverse) orientation

By setting the --SS_lib_type parameter to one of the above, you are indicating that the reads are strand-specific.
By default, reads are treated as not strand-specific.

### Trinity on Uppmax example command line:
```bash
$ salloc -A g2015045 -t 04:00:00 -p core -n 8 --no-shell --reservation=g2015045_20151120 &
$ module load bioinfo-tools java/sun_jdk1.7.0_25 bowtie/1.1.0 samtools/0.1.19 trinity/2014-07-17
$ Trinity --seqType fa --left /proj/g2015045/labs/transcriptome_assembly/mouse_left.fasta --right /proj/g2015045/labs/transcriptome_assembly/mouse_right.fasta --SS_lib_type RF --CPU 8 --max_memory 16G --output trinity_out/
```

NB! -It is recommended to use fully specified paths for sequence files with Trinity.
    -Depending on version of Trinity used --max_memory is sometime given by the command --JM

Data sets: mouse (mouse_left.fasta, mouse_right.fasta), whitefly (whitefly_both.fasta)

## Exercise 2: Assess the data

Explore the Trinity output file Trinity.fasta located in the trinity_out_dir/output directory (or output directory you specify).
Transcripts are grouped as follows: * components: the set of all sequences that share at least one k-mer (including paralogs) * contigs: transcripts that share a number of k-mers (the set of isoforms of a gene) * sequences (isoforms and allelic variation)

2.1 Count the number of sequences in the Trinity.fasta file (hint: use 'grep' and 'wc')

2.2 Get basic information about the assembly with TrinityStats.
```bash
$ /sw/apps/bioinfo/trinity/2014-07-17/milou/util/TrinityStats.pl Trinity.fasta
```
- How many "genes" did Trinity assemble? 
- How many transcripts?
- How large is the assembly? (nr of bases)
- What is N50?

2.3 Filter out sequences shorter than 1000 nucleotides (hint: do a web search for appropriate tools. Someone else must have had the exact same problem.). Count the number of sequences again.

2.4 Align some sequences to a protein database and assess full-lengthness of a couple of sequences (hint: NCBI has an online blast version).

2.5 Find alternatively spliced genes (hint: see above) You can verify alternative splicing by using the UCSC genome browser (do a web search to find it): - Select BLAT from the menu at the top of the page and paste in a mouse transcript sequence from Trinity.fasta.
- Select the mouse/mm9 genome and click “submit”.
- Click on the top scoring hit.

Examine the alignments by clicking “details” on the resulting page.
- Your sequences will be displayed in the browser.
- Enable the mouse annotations (ENSEMBL gene build, UCSC genes, human proteins etc.).

Optional: Do a new transcriptome assembly of whitefly RNAseq data using above code as help.

