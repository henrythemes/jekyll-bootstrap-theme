=====================
Quality control
=====================

In this tutorial we will go through some of the key steps in performing a quality control on your samples. We will start with the read based quality control, using FastQC, and continue with mapping based QC using RseqQC.  

All the data you need for this lab is available in the folder:
``/proj/b2013006/webexport/downloads/courses/RNAseqWorkshop/QC/data/``

Or via web-browser at:
https://export.uppmax.uu.se/b2013006/downloads/courses/RNAseqWorkshop/QC/data

This folder contains:

* Two FASTQ files, with mate pair libraries for sample 1 that is used in several of the other exercises (see `data set summary page <intro.html>`_).
* BAM file and BAM-file index (bam.bai) for that sample mapped to the human genome using STAR
* count_table.txt - a table with number of reads per gene, using Ensembl annotations, created with HTseq-count


Before mapping - FastQC statistics on the reads.
================================================

FastQC aims to provide a simple way to do some quality control checks on raw sequence data coming from high throughput sequencing pipelines. It provides a modular set of analyses which you can use to get a quick impression of whether your data has any problems of which you should be aware before doing any further analysis.

The main functions of FastQC are:

* Import of data from BAM, SAM or FASTQ files (any variant)
* Providing a quick overview to tell you in which areas there may be problems
* Summary graphs and tables to quickly assess your data
* Export of results to an HTML-based permanent report
* Offline operation to allow automated generation of reports without running the interactive application

You can read more about the program and have a look at example reports at http://www.bioinformatics.babraham.ac.uk/projects/fastqc/.

Note: This program can be used for any type of NGS data, not only RNA-seq.

To run FastQC on uppmax you first need to load the module: ::

   module load bioinfo-tools
   module load FastQC/0.11.1

   # To see help information on the FastQC package:
   fastqc --help

   # Run for one FASTQ file:
   fastqc -o outdir fastqfile

   # Run on multiple FASTQ files:
   fastqc -o outdir fastqfile1 fastqfile2 etc.

   # You can use wildcards to run on all FASTQ files in a directory:
   fastqc -o outdir /proj/b2013006/INBOX/FASTQ/*fastq

In this case, only run FastQC on one file and take a look at the output. We have already prepared the outputs for all of the other samples. These can be viewed via a web-browser at:
https://export.uppmax.uu.se/b2013006/downloads/courses/RNAseqWorkshop/QC/fastQC/

Take a look at a few other files and see if they look similar in quality.

Mapping of reads
================
In this exercise we will actually not do any mapping (that will be done in a later exercise). But the commands used to run the mapping and postprocessing are: ::

   # NOTE: do not run, just an example of how we ran the program

   module load bioinfo-tools
   module load star/2.3.1o
   module load samtools

   cd /proj/b2013006/INBOX
   cd FASTQ_hg19_Gencode14.overhang75
   cd 7_111116_AD0341ACXX_137_10_index10__hg19_Gencode14.overhang75

   STAR --genomeDir /proj/b2013006/INBOX/hg19_Gencode14.overhang75 \
    --readFilesIn \
    /proj/b2013006/INBOX/FASTQ/7_111116_AD0341ACXX_137_10_index10_1.fastq \
    /proj/b2013006/INBOX/FASTQ/7_111116_AD0341ACXX_137_10_index10_2.fastq \
    --runThreadN 16 --outSAMstrandField intronMotif

   samtools view -bSh -o Aligned.out.bam Aligned.out.sam
   samtools sort Aligned.out.bam Aligned.out.sorted
   samtools index Aligned.out.sorted.bam


Copy Aligned.out.sorted.bam and Aligned.out.sorted.bam.bai (the BAM index file) from the data directory into your own folder or use the full path to the original BAM file in all your commands. Most programs require a BAM index file, but it is always the name of the BAM file that is provided on the command line, and the programs automatically look for a file with the same name and the .bai extension.

After mapping: map logs
=======================
The first step after you have finished your mapping is to get a general feel of how the mapping went. Most mapping programs produce some sort of summary output, either to a file or to standard out. For example, if using the mapper Bowtie you need to pipe that output to a file to see the summary statistics. In this case the samples were mapped with STAR, that by default creates a file called Log.final.out in the mapping directory. Here is one example of Log.final.out content: :: 

                                 Started job on |       Oct 16 20:21:39
                             Started mapping on |       Oct 16 20:27:04
                                    Finished on |       Oct 16 20:29:14
       Mapping speed, Million of reads per hour |       366.35

                          Number of input reads |       13229276
                      Average input read length |       202
                                    UNIQUE READS:
                   Uniquely mapped reads number |       11913568
                        Uniquely mapped reads % |       90.05%
                          Average mapped length |       198.41
                       Number of splices: Total |       9523918
            Number of splices: Annotated (sjdb) |       9443434
                       Number of splices: GT/AG |       9432792
                       Number of splices: GC/AG |       71488
                       Number of splices: AT/AC |       10675
               Number of splices: Non-canonical |       8963
                      Mismatch rate per base, % |       0.33%
                         Deletion rate per base |       0.01%
                        Deletion average length |       1.75
                        Insertion rate per base |       0.01%
                       Insertion average length |       1.39
                             MULTI-MAPPING READS:
        Number of reads mapped to multiple loci |       356839
             % of reads mapped to multiple loci |       2.70%
        Number of reads mapped to too many loci |       2102
             % of reads mapped to too many loci |       0.02%
                                  UNMAPPED READS:
       % of reads unmapped: too many mismatches |       0.00%
                 % of reads unmapped: too short |       7.21%
                     % of reads unmapped: other |       0.02%


The most important parts to look at are the proportion of uniquely mapping, multi-mapping and unmapped reads. We ideally want the uniquely mapping reads to be as high as possible. Multi-mapping or unmapped reads could indicate poor quality of the reads, adapter contamination or other reasons for low quality scores.

Another key point is the mismatch and indel rates. If they are very high, this could indicate that there has been some problems during the sequencing or during the library preparation.


After mapping: RseQC
====================

The RseQC package is one of many tools to get basic mapping statistics from your BAM files. This package provides a number of useful modules that can comprehensively evaluate high throughput sequence data, especially RNA-seq data. Some basic modules quickly inspect sequence quality, nucleotide composition bias, PCR bias and GC bias, while RNA-seq specific modules evaluate sequencing saturation, mapped reads distribution, coverage uniformity, strand specificity, etc. You can read more about the package at: http://rseqc.sourceforge.net/

The RseQC package contains many steps that are equivalent to FastQC analysis, e.g. read quality, sequence composition (NVC), GC-bias etc, but the results may be different since many of the low quality reads may not map to the genome and therefore will not be included in the BAM file.

Running all the QC steps takes a long time, so to save time, we only run the QC on a random selection of 10% of the reads. Random selection of reads can be performed with many different programs. Here we will use samtools: ::

    samtools view -b -s 0.1 Aligned.out.sorted.bam > Aligned.out.0.1.bam
    # then index the bamfile
    # (it is already sorted since you extracted reads from a sorted BAM file)
    samtools index Aligned.out.0.1.bam
 
The RseQC package is allready installed at Uppmax. Load the package: ::

    module add bioinfo-tools
    module add rseqc/2.4

Some steps of the RseQC package require a file with gene annotations in BED format. These can be downloaded from various sources. Some of the more common ones are UCSC, RefSeq and Ensembl. In this case, the RseQC team have already created annotation files in some common formats that can be downloaded from their website, but if you have data for a less studied organism you may need to create a BED-file on your own. 

Two annotation files have already been downloaded into ``/proj/b2013006/webexport/downloads/courses/RNAseqWorkshop/QC/annotation`` for you to use. These are: hg19.HouseKeepingGenes.bed  and hg19_RefSeq.bed. The folder also contains a reduced annotation file hg19_RefSeq_top1000.bed to speed things up. 

In this tutorial we will not run all the different parts of the RseQC package, only the most relevant ones for this experiment. The different scripts in the RseQC package are well described at their website (http://rseqc.sourceforge.net/), so read the instructions there and specify the input/output files to fit your file names and folder structure. 

The steps that we are going to run are:

1. geneBody_coverage.py
2. inner_distance.py
3. junction_saturation.py
4. read_distribution.py

Note: The geneBody_coverage.py script takes a very long time to run, so we have created a subsection of annotations to run it on. Use the file hg19_RefSeq_top1000.bed. This file was created with the command: ::

      # head -n 1000 hg19_RefSeq.bed > hg19_RefSeq_top1000.bed

Also note: When running read_distribution.py, an outfile cannot be specified. Instead you need to pipe (">") the output to a file, or look at the output in the terminal.


Run RseQC for one sample and have a look at your output. 

* Do most of your reads map to genes? 
* Do you have even coverage along the genes? 
* Do the reads cover most splice junctions? 
* Based on the inner distance plots, what do you think the average fragment size of the libraries was?

We have run the QC for all the samples and compiled summary files `here <https://export.uppmax.uu.se/b2013006/downloads/courses/RNAseqWorkshop/QC/output/>`_.
This folder contains one table that summarizes all the Log.final.out files from all the samples (summary_starlog.txt), and one pdf file with a few different plots to summarize those statistics (summary_starqc.pdf). There are also plots from the four RseqQC modules listed above including all the samples.

What is your conclusion, do your samples look good? Is there anything that looks strange in any sample, or do you feel comfortable using all the samples in your analysis?


Outlier detection and general overview of data
==============================================

One of the first steps once you have your libraries mapped to the genome and have filtered out low quality samples is to get a general overview of the samples. Logical first steps are to look for pairwise correlations between the samples, do some simple clustering and run principal component analysis (PCA). With these steps you can easily find out what the variation within your sample groups looks like and detect possible outliers or mixed up samples. We will run this analysis with a few simple R commands, but there are of course other options for how to run this analysis. 

For this exercise we have pre-calculated read counts per gene (according to Ensembl annotations) with commands like: ::

    # NOTE: Only given for reference
    #       Not supposed to be executed during the lab
    samtools view accepted_hits_137_1.bam | \
     sort > accepted_hits_prehtseq_137_1.sam
    htseq-count -s no -q accepted_hits_prehtseq_137_1.sam \
     Homo_sapiens.GRCh37.71.gtf > 137_1.counts

This was run for each of the samples and the counts were combined into a single table. You can get the count table from the data directory. You can run R on UPPMAX, or download the file to your local computer and do the analysis locally if you prefer.

The code to run in R: ::

  # read in the data
  counts <- read.delim("count_table.txt")
  head(counts)

As you can see, the samples are ordered with the 3 replicates from each group next to each other. So when we are to define colors for the samples we only have to repeat each color 3 times (this may not always be the case!) ::

  # define colors:
  col.def<-c("red","blue","green","magenta")
  sample.def<-c("ctrl", "t2h", "t6h", "t24h")
  colors <- rep(col.def, each=3)


Start with a PCA to se the general distribution. PCA of RNA-seq data is usually performed in log-scale. We also add a pseudo-count of +1 to avoid logging zero (gives infinity). You need to transpose - t() - the data matrix, otherwise you will run PCA on the genes instead of samples. ::

  myPca <- prcomp(t(log2(counts+1)))

This creates a list that contains:

* the samples mapping to each PC in myPca$x
* PC contribution to variance in myPca$sdev
* PC loadings for each gene in myPca$rotation

Now some plotting. In R you can either plot into a default window or direct all your output to a "device", that can be pdf, png, tiff etc. To open a new pdf device: ::

  pdf('pca_plot.pdf')
  # once you have plotted all you want to put into that file,
  # close it with dev.off()

Let's first make a simple plot of the first two principal components (PC1 vs PC2): ::

  plot(myPca$x[,1],myPca$x[,2],col=colors,pch=1)
  legend("topright",sample.def,pch=1,col=col.def)
  dev.off()

Sometimes the first two PCs may not be the ones that will best separate the sample groups, so it is a good idea to look at more PCs.
Here is one example that shows how to plot the top 5 PCs: ::

  pdf('pca_plot_5pc.pdf')
  tmpPcaData <- as.data.frame(myPca$x[,1:5])
  plot(tmpPcaData, col=colors,pch=1)
  dev.off()


Another thing to look at is the pairwise correlation between all the samples and see how they group based on correlation. Let's create one matrix with all pairwise Pearson correlations (again in log-space). ::

  nSamples<-ncol(counts)
  C<-mat.or.vec(nSamples,nSamples)
  for (i in 1:nSamples) {
     for (j in i:nSamples){
        if (i==j){ C[i,j]<-NA }
        else {
             c<-cor(log2(counts[,i]+1),log2(counts[,j]+1),method="pearson")
             C[i,j] = c
             C[j,i] = c
        }
     }
  }
  colnames(C)<-colnames(counts)
  rownames(C)<-colnames(counts)

This can also be calculated as one command with the R apply function, but to clarify what is being calculated we included a more descriptive code. Another way to do the same thing would be: ::

  C<-apply(log2(counts+1),2,cor,log2(counts+1),method="pearson")
  diag(C)<-NA


Now you will plot a heatmap with the correlations: ::

  pdf('correlation_heatmap.pdf')
  heatmap(C,symm=TRUE)
  dev.off()

Do the clusterings agree with what you expect? 
Which different sample groups are more similar? Are some sample groups more dissimilar compared to the others?


