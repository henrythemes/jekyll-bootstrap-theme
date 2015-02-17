---
layout: default
title:  'RNAseq'
---

# Transcriptome Mapping First

<font color="red">**Please read everything carefully!**</font>

# Reference-genome based: Tophat and Cufflinks

A common problem in the analysis of RNA-seq data is to relate it to a known genome sequence and use that information to study the expression of genes - within a sample or across multiple conditions, such as tissues or time points. A popular pipeline to perform such an analysis is the Tuxedo protocol, which consist of set of programs that can be used to go from mapping of short reads to reference genomes all the way to detection of differential gene expression. The two main programs included in package are Tophat - a short read mapper and Cufflinks that performs analysis of the mapped reads. In this exercise you will learn how to use some of these tools to study gene expression differences between different human tissues.

Note: Do not simply copy the various unix commands quoted throughout this tutorial. Many include placeholders (folder names etc), so make sure you use whatever file names you have created.

[Illumina Bodymap2.0](http://www.ebi.ac.uk/gxa/experiments/E-MTAB-513) data consists of 16 human tissues that were sequenced using both single-end and pair-end technologies. The mapped reads can also be visualised at the [Ensembl genome browser](http://www.ensembl.info/blog/2011/05/24/human-bodymap-2-0-data-from-illumina/)   

In this tutorial we will focus on a limited set of tissues and only do the analysis for chromosome 1 of the human genome. The files are in addition just a small subset of the original files as the analysis of them would take too long to fit at a course lab.

The main goal goal with this tutorial is to detect differential gene expression between two different tissues (pick any two tissues that you want to compare) from human. For all included tissues there is one single-end library and one pair-end available. In order to test for significance in gene expression more than one replicate from each tissue is needed. In this lab we will use the two different library types as replicates in the detection of differential gene expression. Below is a summary of the data and tissues available.

* Single-end reads, 75bp
  * ERR030888: Female adipose
  * ERR030890: Female brain
  * ERR030892: Female colon
  * ERR030893: Female kidney
  * ERR030901: Female ovary
* Pair-end reads, 2 x 50bp
  * ERR030880: Female adipose
  * ERR030882: Female brain
  * ERR030884: Female colon
  * ERR030885: Female kidney
  * ERR030874: Female ovary
* human reference genome (or in this lab, for the sake of saving time only chromosome 1 named: rm.chr.1.fa)
* genome index for aligning reads with Bowtie2
* reference genome annotation based on the [EnsEMBL](http://www.ensembl.org/index.html) database named: Homo_sapiens.GRCh38_Chr1.77.gtf
* NB! All intermediate files are available at Uppmax so if any of the steps fails you can pick up the analysis at the next step pf the tutorial

## Tophat

Tophat is a script pipeline built on-top of the popular read-aligner Bowtie. It is used to align RNA-Seq reads to a genome in order to identify exon-exon splice junctions. More specifically, it produces data that we can use to not only study the expression of genes, but also their different isoforms. You will have a bit of waiting time during the exercises as the more complex analyses are running, so please check out some of the details of tophat here during those times: [here](http://ccb.jhu.edu/software/tophat/index.shtml)

## Cufflinks

Cufflinks is a collection of programs that perform different steps in the analysis of aligned RNA-seq reads ([Details](http://cole-trapnell-lab.github.io/cufflinks/cufflinks/index.html)). The output is usually a list of transcribed loci (primarily ‘genes’) and their expression levels within and/or between samples. For the analysis of multiple data sets, the general workflow in cufflinks consists of the following steps:
   * Cufflinks: Assemble the aligned reads of a given sample, identify transcribed loci and determine expression
   * Cuffmerge: Reconcile data on transcribed loci across multiple samples to produce a consensus annotation of loci
   * Cuffdiff: Compare read data across samples, guided by consensus annotation, and determine differential expression of loci, test for significance The main output we are interested in comes from the cuffdiff analysis and consists of differential expression estimates for a set of genes. In the following, we will be going through the necessary steps to accomplish this.
---++ Step-by-Step
   1 Prepare your data
   1 Load software
   1 Run Tophat on individual samples
   1 Run Cufflinks on individual samples
   1 Run Cuffmerge merge detected transcript over all samples
   1 Run Cuffdiff
## 0) Book a node

We have reserved half a node for each student during this course. By now, you are probably already familiar with the procedure:

<verbatim>$ salloc -A g2015005 -t 08:00:00 -p core -n 8 --no-shell --reservation=g2015005_thu &</verbatim>

Make sure you ony do this once, otherwise you will take away resources from the other course participants!
## 1) Prepare your data

*%RED% Note: It is completely up to your how you organize your data - what follows below is merely a suggestion:%ENDCOLOR%*

   * create a folder for your project <verbatim>$ cd ~/glob2
$ mkdir transcriptome
$ cd transcriptome</verbatim>

   * sym-link the required files and folders (this will ceate a symbolic link to the original folders/files and saves you the trouble of always typing the full path - BUT: Do not write into these linked folders, because that data is shared across everyone working with these folders...)
   <verbatim>ln -s /proj/g2015005/labs/transcriptome_map/reads/PE/</verbatim>
   <verbatim>ln -s /proj/g2015005/labs/transcriptome_map/reads/SE</verbatim>
   <verbatim>ln -s /proj/g2015005/labs/transcriptome_map/results</verbatim>
   <verbatim>ln -s /proj/g2015005/labs/transcriptome_map/reference</verbatim>

Your directory structure should look like this:

   * working directory (your choice)
      * PE/
      * SE/
      * results/
      * reference/

We are skipping a few steps here, namely obtaining a reference annotation and genome sequence and preparing the latter for use with the Bowtie2 aligner. We have taken care of that for you (located in the subfolder /reference). This is due to two main factors. First, it takes a lot of CPU hours to convert a genome sequence into a Bowtie index. Second, finding the latest release of a genome sequence free of unmapped fragments and haplotype data as well as a fully Tophat/Cufflinks-compatible annotation of that sequence is an exercise in frustration for beginners. One useful resource here is the FTP server of Illumina <a href="http://support.illumina.com/sequencing/sequencing_software/igenome.html" target="_top">here</a> . Finally, also note that we are providing you with the outputs of the different steps. This is to make sure that if you run into some trouble, like software crashing half-way through analysis, you can still continue with your exercises.
## 2) Load software

You have done this before, but here is a quick reminder:

<verbatim>$ module load bioinfo-tools samtools/0.1.19 bowtie2/2.2.3 tophat/2.0.12 cufflinks/2.2.1</verbatim>

If any of these packages does load as expected, you can check that module names are correct using the command

<verbatim>$ module avail</verbatim>

It may be that you need to load a different version.

## 3) Run Tophat

What goes in, what comes out:

In:   
   * One or several <a href="http://en.wikipedia.org/wiki/FASTQ_format" target="_top">FastQ</a> files (one for single-end reads, two for paired-end)
   * An indexed reference genome (and optional a reference annnotation (eg. a [[http://www.ensembl.org/info/website/upload/gff.html][GTF or GFF file]]))

Out:
   * A read alignment (BAM)<hr />

NOTE: The /reads folder contains a small subset of an actual FASTQ, limited to a specific chromosome. For the next step, please chose two tissues that you want to analyse and make sure to align both a pair-end and a single-end library from you tissues of choice. We do this since the time needed to align an actual read file with data from all human chromosomes can take several hours. With these sub-sampled data sets, it should be possible to align them with tophat in a reasonable time frame. However, if this should still take too long (&gt; 20mins), you may wish to abort this step. You already have the corresponding output the subfolder results/tophat/.

Tophat will take one or multiple FASTQ files and align the reads therein to a genomic reference. A common command may look like this:

<verbatim>$ tophat -o tophat_outputSE30888 --solexa-quals -p 8 --library-type=fr-unstranded reference/rm.chr.1 SE/ERR030888.fq.gz
$ tophat -o tophat_outputPE30880 --solexa-quals -p 8 -r 200 --mate-std-dev 90 --library-type=fr-unstranded reference/rm.chr.1 PE/ERR030880_1.fq.gz PE/ERR030880_2.fq.gz</verbatim>

We specify the output location (-o), the number of CPUs to use (-p), which type of sequencing library was used to produce the data (here ‘fr-unstranded’), in which format the quality information was stored (here ‘solexa’, pre 1.3), the location of the reference annotation, the location of the Bowtie2-formatted index file for the genome sequence and finally a FASTQ file. For pair-end data we have two FASTQ files and also define the expected size and variation in size of the fragments sequenced.

While this is running, you may want to head over to the <a href="http://cole-trapnell-lab.github.io/cufflinks/cufflinks/index.html" target="_top">tophat manual</a> and have a look at the available options and technical details.

The aligned reads are found in the output directory *you have chosen*, e.g. accepted_hits.bam. For convenience, you may want to sym-link this file into your main project folder:

<verbatim>$ ln -s tophat_outputSE30888/accepted_hits.bam SE30888.bam</verbatim>
## 4) Cufflinks: Assembly and transcript calling

What goes in, what comes out:

In: A read alignment in BAM format (SAM is also an option, but should not be used due to it being uncrompressed)

Out: A number of files, including a transcriptome annotation reconstructed from the read distribution
---

General command format:

<verbatim>$ cufflinks -o my_output_folder -p 8 -g reference/Homo_sapiens.GRCh38_Chr1.77.gtf my_infile.bam</verbatim>

Note: The Cufflinks step can take a while - so this is now a good time to get a coffee, or read through the available documentation on the Cufflinks website. While we have tried to cover the technical details of these tools to some degree in the lecture, there are a lot of details that will help you use these programs to their greatest effect.

Here we specify where to store the output, how many CPUs to use as well as where to find the reference files (genome sequence and annotation). Depending on the size of the BAM/SAM file, this step may require several hours to complete. To make this analysis feasible within the time limits of the course, we have created the chromosome-limited files you have been using.

The command line output will read something like:

&gt; Processed 48858 loci. [* * * * *] 100%

&gt; Map Properties:

&gt; Total Map Mass: 1893657.20

&gt; Fragment Length Distribution: Truncated Gaussian (default)

&gt; Default Mean: 200

&gt; Default Std Dev: 80

One important thing that can be noted here:<br /><br />
   * Processed loci - these are the transcribed regions, or 'genes'. How does this number compare to offical estimates of human gene content?
   * Total Map Mass - a measure for the size of your read library
   * Length distribution - this value measures the distance between mate-paired reads. It can either be specified (if known) or will be determined by Cufflinks. Since we are using single-end reads, this value should not matter. The output of this run can then be found under my_output_folder/ and includes a total of 4 files:

genes.fpkm_tracking

isoforms.fpkm_tracking

skipped.gtf

transcripts.gtf

The first two files contain basic information about expressed genes and transcripts, respectively - those known from the annotation file as well as novel loci identified by cufflinks -and the strength of their expression, given in FPKM. FPKM stands for ‘Fragments Per Kilobase of exon per Million fragments mapped’ and is a normalized measure of transcript abundance. That is the short explanation. The longer version for the more mathematically inclined among us can be found at <a href="http://www.nature.com/nmeth/journal/v5/n7/abs/nmeth.1226.html" target="_top">Mortazavi et al. 2008</a> .

These output files are tab-delimited and can e.g. be opened in e.g. Microsoft Excel (or similiar) to be analyzed and/or visualized.
## 5) Cuffmerge: Reconciling different transcript models

What goes in, what comes out:

In: An optional reference annotation and a list of transcript annotations to merge

Out: A consensus annotation, taking into account all input annotations
---

It is important to keep in mind that reference annotations are very likely incomplete. This is because some genes or individual exons may be expressed at very low levels or under specific conditions, thus having evaded prior detection. Moreover, many vertebrate genomes have only been annotated by reference to other genomes, which themselves may only be poorly characterized. Using the expression data obtained through cufflinks may hence allow us to improve existing annotations. Cuffmerge is a tool that takes cufflinks-derived annotation files (known & ‘novel’ loci) and reconciles them into a consensus annotation, discarding e.g. spuriously transcribed loci and merging overlapping loci into larger transcription units where possible.

Again, the commands below are<b> just examples</b>, your files and folder may be called differently.

<verbatim>$ cd ~/glob2
$ mkdir cuffmerge
$ cd cuffmerge
$ ln -s ../cufflinks.brainSE/transcripts.gtf brainSE.gtf
$ ln -s ../cufflinks.kidneyPE/transcripts.gtf kidneyPE.gtf</verbatim>

<i> *Note:* If this didn't work (check that the linked files actually exist and have content), then you probably chose a different way of organizing your folders and will have to figure out where your source files are yourself <img alt="wink" border="0" src="http://array.medsci.uu.se/twiki/pub/TWiki/SmiliesPlugin/wink.gif" title="wink" /></i>

Now that we have both transcript model files in one location, we can attempt to merge them. For this, we first have to create a text file that contains a list of GTF files to merge (quite inconvenient, I know). Use whichever tool you feel comfortable with and write the name of each gtf file in one line, then save it as transcripts.txt.

<verbatim>$ cuffmerge -o merged -g reference/Homo_sapiens.GRCh38_Chr1.77.gtf -p 8 -s reference/genome.fa transcripts.txt</verbatim>

This will save the reconciled annotation file as merged/merged.gtf. Symlink this file into your main project folder.

<verbatim>$ cd ~/glob2/transcription
$ ln -s cuffmerge/merged/merged.gtf</verbatim>

Now we are ready to check for differential expression in our read data from chromosome 1.
## 6) Cuffdiff: Differential expression analysis

What goes in, what comes out:

In: A consensus annotation (or just the reference annotation), read alignments for all samples that are to be compared and quantified

Out: A number of files, including tests for differential expression for all pairwise comparisons (gene_exp.diff)
---

Cuffdiff takes aligned reads from two or more samples, estimates comparable expression values and performs statistical analysis of the resulting data to determine which genes exhibit significantly different activity profiles between any two samples. The nitty-gritty details of the underlying mathematics can be found here: [[http://cole-trapnell-lab.github.io/cufflinks/cuffdiff/index.html][cuffdiff]] .

For running Cuffdiff, we type something like this (being in the main directory of our project):

<verbatim>$ cuffdiff -o cuffdiff.brain_vs_kidney -L brain,kidney -p 8 -u merged.gtf brainPE.bam,brainSE.bam kidneyPE.bam,kidneySE.bam</verbatim>

Adopt this to your data - if uncertain, run cuffdiff -h to learn more about the options. Note that the labels you give them are arbitrarily chosen, pick names that make sense to you.

This will write the output of the analysis into the subfolder cuffdiff.brain_vs_kidney (or whatever folder name you chose) and do a pairwise comparison of the samples (Note: the order in which you list the labels needs to be the same as the order of SAM/BAM files!). NB! To reduce the computing time we do not use the flag -b that correct for genome sequence (see manual for details), but we resolve issues arising from reads mapping to multiple loci in the genome using the flag -u.

The main file of interest to us is gene_exp.diff. It includes the analysis of differential expression. A quick way to find the cases of interest is to filter the file for genes that show evidence of differential expression (identified by the tag ‘yes’ in the ‘significant’ column).

<verbatim>$ head -n1 gene_exp.diff > results.txt
$ grep yes gene_exp.diff >> results.txt</verbatim>

(This copies the header of the output file as well as all rows tagged as significant into a new text file - open this file in a text editor or spread sheet program).

Using their _EnsEMBL_ accession numbers, you can go to <a href="http://www.ensembl.org/" target="_top">http://www.ensembl.org</a> to retrieve information on the function of these genes and see whether you can draw any conclusions as to why these genes would be differentially expressed between samples.
---++ Where to go next

So now you have analyzed the expression of genes between two samples. However, usually the work does not end here. For example, you may want to perform a thorough analysis of your output, visualze distributions and obtain statictics. This can be done either through clever scripting in R, or by use of a recently developed software suite called <a href="http://compbio.mit.edu/cummeRbund/" target="_top">CummeRbund</a>. It reads the native output from Cuffdiff, parses it into a database and provide ample options for in-depth analysis of the data. This package offer a lot of efficient parsing of the output files created by cuffdiff, however a recent update to Rsqlite package has broken the procedure whereby this package reads the data into R.
---++ Closing remarks

This tutorial has introduced you to a very straight-forward, but somewhat simplified pipeline for the analysis of RNA-seq data by use of a reference genome to study transcription. Both Cufflinks and Tophat come with additional parameters that we have not touched upon to avoid unnecessary confusion. Likewise, the read data we have used was strand-unspecific. This has certain drawbacks, specifically with respect to accuracy in the isoform analysis. Or perhaps you are not interested in comparing expression between pairs of samples but in a time series. For this reason as well as others, you may need to adjust one or several parameters to get the best results - depending on the nature of your data. We therefore highly recommend you to carefully read both manuals (and possible the original publications) so as to familiarize yourself with these additional options.
</div>