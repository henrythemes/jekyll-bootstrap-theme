---
layout: default
title:  'RNAseq'
---

# Transcriptome De Novo Assembly

## Trinity

By efficiently constructing and analyzing sets of de Bruijn graphs, Trinity fully reconstructs a large fraction of transcripts, including alternatively spliced isoforms and transcripts from recently duplicated genes.
Compared with other de novo transcriptome assemblers, Trinity recovers more full-length transcripts across a broad range of expression levels, with a sensitivity similar to methods that rely on genome alignments.
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
The data is obtained from mouse dendritic cells (mouse_left.fasta and mouse_right.fasta and) and a whitefly (whitefly_both.fasta), and the files are located in `/proj/g2015005/labs/transcriptome_assembly/`.
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
$ module load bioinfo-tools java/sun_jdk1.7.0_25 bowtie/1.1.0 samtools/0.1.19 trinity/2014-07-17
$ Trinity --seqType fa --left /proj/g2015005/labs/transcriptome_assembly/mouse_left.fasta --right /proj/g2015005/labs/transcriptome_assembly/mouse_right.fasta --SS_lib_type RF --CPU 8 --max_memory 16G --output trinity_out/
```

NB! -It is recommended to use fully specified paths for sequence files with Trinity.
    -Depending on version of Trinity used --max_memory is sometime given by the command --JM

Data sets: mouse (mouse_left.fasta, mouse_right.fasta), whitefly (whitefly_both.fasta)

## Exercise 2: Assess the data

Explore the Trinity output file Trinity.fasta located in the trinity_out_dir/output directory (or output directory you specify).
Transcripts are grouped as follows: * components: the set of all sequences that share at least one k-mer (including paralogs) * contigs: transcripts that share a number of k-mers (the set of isoforms of a gene) * sequences (isoforms and allelic variation)

2.1 Count the number of sequences in the Trinity.fasta file (hint: use 'grep' and 'wc')

2.2 Count the number of nucleotides in the Trinity.fasta file (hint: use 'grep' and 'wc')

2.3 Filter out sequences shorter than 1000 nucleotides (hint: do a web search for appropriate tools.
Someone else must have had the exact same problem.).
Count the number of sequences again.

2.4 Align some sequences to a protein database and assess full-lengthness of a couple of sequences (hint: NCBI has an online blast version).

2.5 Find alternatively spliced genes (hint: see above) You can verify alternative splicing by using the UCSC genome browser (do a web search to find it): - Select BLAT from the menu at the top of the page and paste in a mouse transcript sequence from Trinity.fasta.  
- Select the mouse/mm9 genome and click “submit”.  
- Click on the top scoring hit.

Optional: examine the alignments by clicking “details” on the resulting page.  
- Your sequences will be displayed in the browser.  
- Enable the mouse annotations (ENSEMBL gene build, UCSC genes, human proteins etc.).

Redo above exercises with the whitefly data.

## Final commments

Many tools that are available on most unix systems can be used to count featueres (like sequences etc) in text files.

```bash
$ grep '>' file.fasta | wc 
```

Will identify all fasta header rows in a fasta file and parse the output to the command wc which count the number of rows and characters and can hence be used to count the number of sequences in a file

```bash
$ module load Fastx/0.0.14
$ fasta-formatter -i file.fasta > filenew.fasta
$ grep '>' filenew.fasta | tr -d "len=" | awk '{sum+=$2} END {print sum}'
```

Will use the program fasta-formatter ot create a fasta file that instead of being 80 characters wide will have the header on 1 row and then the full sequence on the second row.
Header rows are then parse to the tr command, where the characters len= are deleted and the second column that now contains the length reported by trinity are summed up using awk.

```bash
$ module load ucsc-utlities/v287
$ faFilter -minSize=1000 file.fasta filenew.fasta
```

Filter the fasta file to only contain reads longer or equal to 1000bp long.

```bash
$  grep -E '.{1000}' filenew.fasta -B1 | greop '>' -c
```

NB! This needs a fasta file that have each single sequence on one line (eg. the output from fasta-formatter).
It identifies all rows which have 100 or more characters and grabs the row before that (-B1) pipes this to grep and count the number of header rows that are followed by a 1000bp or longer sequence.
