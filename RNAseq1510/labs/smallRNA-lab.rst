============================
Small RNA computer excercise
============================

In this exercise we will analyze a few small RNA libraries, from Drosophila melanogaster (fruit fly) embryos and two cell lines (KC167 cells derived from whole embryos, and ML-DmD32 cells derived from adult wing discs). This is a subset of a much larger data set used to study microRNAs and other small RNAs in Drosophila:

- "emb_0-1". Small RNA from 0-1 hour embryos. (SRA id: SRR069840)
- "emb_6-10". Small RNA from 6-10 hour embryos. (SRA id:SRR069839) 
- "kc167". Small RNA from KC167 cells. (SRA id:SRR028728)
- "ml-DmD32_r1". Small RNA from ML-DmD32 cells replicate1. (SRA id: SRR069508)
- "ml-DmD32_r2". Small RNA from ML-DmD32 cells replicate2. (SRA id: SRR488696)

These data sets are described more in this paper: http://genome.cshlp.org/content/24/7/1236.full

The aim of the exercise is to show a simple way to process small RNA data and to quantify the expression of microRNAs, and to make some plots of global expression patterns in R. The exercise consists of two parts: First you will use a genome browser to get a feeling for what the small RNA reads mapped to the genome look like. Next, you will preprocess the small RNA reads, map them to the known microRNA loci, and quantify the expression of the microRNAs.



Preliminaries
=============

Running this on UPPMAX, start by loading the modules you will need: ::

	module load bioinfo-tools
	module load cutadapt
	module load bowtie

And create a directory to work in: ::

	# create a subfolder called "smallRNA"
	mkdir ~/glob/RNAseqCourse/smallRNA

All data and scripts required for this exercise can be found in 
``/proj/b2013006/webexport/downloads/courses/RNAseqWorkshop/smallRNA`` on UPPMAX and through this `URL <https://export.uppmax.uu.se/b2013006/downloads/courses/RNAseqWorkshop/smallRNA/>`_ .


This includes:

- 6 fastq files with the raw reads from the small RNA sequencing (in the subdirectory fastq).
 
- A fasta file with the sequence of all microRNA loci and a gff file with the coordinates of all microRNA loci in the Drosophila genome (in in the subdirectory mirbase).
 
- A script, "sam2expTable.pl" (in the subdirectory scripts), to count the reads mapping to each microRNA.
 
- Precomputed bam files with the sequencing data mapped to the entire Drosophila genome, which can be used for browsing in IGV (in the subdirectory mapped_to_genome).

Copy these files in the directory you will use for this exercise. On UPPMAX you can use the following command :: 

	cp -r /proj/b2013006/webexport/downloads/courses/RNAseqWorkshop/smallRNA dest

where dest is the destination directory. (This might take a while since the fastq files are quite large)

Browse small RNA reads 
======================

We will start by browsing how the small RNA reads look when mapped to the Drosophila genome. For this we will use pre-computed files, which can be viewed with IGV or some other genome browser. (Here only instructions for IGV are given.) Start IGV and load the files emb_0_1.sorted.bam and ml-DmD32_r2.sorted.bam. Also load the file with all microRNA annotations, dme_mirbase.gff3.

To load a file you first select the correct genome ("D. melanogaster r5.22") in the top left menu. Then go to the File menu, and select "Load from file", and select the files described above.

Type the name of a microRNA, e.g "mir-124", to go to that locus. You can see that the read mapping patterns are very distinct: (Almost) only the processed microRNAs end in the sequencing libraries, and you can see reads from both arms of the hairpin structure. 

While many microRNAs occur alone in the genome, other are arranged in clusters. Type "let-7" to browse such a cluster. How many microRNAs do you see in this region?

(To see something weird, go to "3R:18,118,436-18,118,767". Do you have any idea what this could be?)

Adaptor trimming
================

When sequencing small RNAs we are working with very short RNA fragments, typically shorter than the reads. This means that most reads will contain parts of adaptor sequences that were inserted during library preparation. These are found at the (3') end of the reads. Before we can do anything else with the data we have to remove these sequences. 

Look at any fastq file, e.g. using less: ::

	less in.fastq

Just by looking the nucleotide sequence, can you guess what the adaptor is? (This is actually a useful exercise. Many datasets are poorly annotated, and no information is given about the adaptor sequence.  When analyzing such data the only option is to infer the adaptor from the sequence data.)

There are many programs available for trimming adaptors. We will use a program called cutadapt (https://code.google.com/p/cutadapt/). You can run it with the following command: ::

	cutadapt -a adaptor --trimmed-only in.fastq --minimum-length=17 > trimmed.fastq

This trims the sequence given in adaptor from the reads in in.fastq and prints the results to a new file trimmed.fastq. It also applies the following filters: only reads where the adaptor was trimmed are printed to the output, and only reads that are at least 17 nucleotides after trimming are printed. For the data in this exercise, use the adaptor sequence CTGTAGGCACCATC.

Run this program on each of the 6 fastq files. This takes a few minutes per file. Make sure you give the resulting files good names (e.g. kc167_TRIM.fastq) so you can keep track of all files.

How many reads were removed because they didn't have the adaptor sequence or because they were to short?

Mapping
=======

The next step is to align (map) the reads to the genome sequence around the microRNA loci. We will use the program bowtie (http://bowtie-bio.sourceforge.net/index.shtml) to do this. We will map the reads against the microRNA loci in mirBase (http://www.mirbase.org), which is "the official" data base of microRNAs in many different species. To be able to map millions of reads very fast, bowtie creates an index of the sequence we map against. You can create the index using: ::

	bowtie-build seq.fastq index.name

Here, seq.fastq is the file with sequences we want to map against (in our case dme_mirs.fa) and index.name is the path and name of the bowtie index we create (e.g. "mydirectory/dme_mirs").

Now we can map all reads. We do this using the following command: ::

	bowtie -q -v 0 -k 10 -S -t index.name small_rna.fastq out.sam

Here, index.name is the bowtie index created above, small_rna.fastq is the file with the small RNA data (after trimming!) and out.sam is the resulting file. This maps the reads with the following settings: input is a fastq file (-q), no mismatches are allowed (-v 0), max 10 hits are reported for each read (-k 10), output is a sam file (-S) and the time the mapping took is printed to the screen (-t).  Run this command once for every file with trimmed reads.

(If you feel like it, try mapping one of the fastq files where the adaptor was not trimmed, and see what happens.)


Quantification of microRNAs
===========================

We can now summarize the mapped reads to see which microRNAs are expressed in the different samples, and to do some global comparisons. For this, we use the sam files created by bowtie. If you have not seen a sam file before,  have a look at one of the files, for example by running: ::

	less out.sam

Press space to scroll down into the file and q to exit the viewer. 

In the folder with all files for this exercise you will find a script sam2expTable.pl. This script reads all sam files in a folder, and for each file counts the reads mapping to each sequence (i.e. to each microRNA). It returns a table with one row per microRNA locus and one column for each sam file, where each element is the number of reads mapping to a specific microRNA from a specific sam file. Copy this script to somewhere in your folder, and do: ::

	chmod a+x sam2expTable.pl

to make the script executable. Then run it with: ::

	./sam2expTable.pl sam.dir > out.table

Here, sam.dir is the directory with all sam files and out.table the file to which the output is printed. This might take a few minutes.

Once the reads mapping to each microRNA have been counted, we can analyze the microRNA expression levels using R. Start R by typing: ::

	R

You will see a different prompt, since you are now typing commands to R. You can always exit R with quit(). Load the expression table you just created into R: ::

	exp.data <- read.table("out.table", header=TRUE, row.names=1, sep="\t")

Here out.table is the full path to the file with the expression table. You can look at the first 20 rows of the table by typing: ::

	exp.data[1:20,]

If some microRNAs are very similar, the same reads might map to them. See for examples dme-mir2b-1 and dme-mir2b-2. In this exercise we don't handle such cases in any special way. Can this be a problem? If so, how would you deal with it?

Since the log transformation we will do later cannot handle cases with zero reads, we add a dummy value of 1 read to each microRNA: ::

	exp.data <- exp.data + 1

To compare expression lkevels from different libraries, the read counts have to be normalized to compensate for different sequencing depths etc. For this we will use the TMM normalization. This normalization method uses a trimmed mean of M- values (TMM) between each pair of samples to find a set of scaling factors for the library sizes that minimize the log-fold changes between the samples for most genes (if you are interested in the details, see http://genomebiology.com/2010/11/3/r25). To use this method we need to load the edgeR module. edgeR is an R module with many useful functions for normalizing RNA-seq data and finding differentially expressed genes. Here we will only use one of the normalization functions. ::

	library(edgeR)

If you get an error message that the edgeR module is not installed on the computer you are using, you can download and install it with: ::

	source("http://bioconductor.org/biocLite.R")
	biocLite("edgeR")

In the normalization, we start by computing the factors by which the read counts from each library are rescaled: ::

	lib.size <- apply(exp.data,2,sum)
	scale.factors <- calcNormFactors(exp.data, method="TMM") 

Next, we apply the rescaling to the read counts for each library: ::

	norm.data <- t(t(exp.data)/(scale.factors*lib.size))

Finally, we log transform all values. This makes the analysis less sensitive to microRNAs with a huge number of reads: ::

	norm.data <- log(norm.data)

We can use principal component analysis (PCA) to get a global look of how similar the microRNA expression profiles are in the different libraries: ::

	mir.pca <- prcomp(t(norm.data))     ## compute principal components
	plot(mir.pca$x[,1], mir.pca$x[,2])  ## plot  PC1 and PC2
	text(mir.pca$x[,1], mir.pca$x[,2], rownames(mir.pca$x), cex=0.7, pos=4, col="red")

You should now see a plot on the screen. What can we learn from looking at the PCA plot?

We can also look at the loadings, i.e. how much each microRNA contributes to each principal component. To see which microRNAs are highly expressed in samples with high PC1, type: ::

	head(sort(mir.pca$rotation[,1], decreasing=TRUE))

To see which microRNAs are highly expressed in samples with low PC1, type: ::

	head(sort(mir.pca$rotation[,1], decreasing=FALSE))

(Some background about some specific microRNAs: bantam is known to prevent apoptosis by repressing pro-apoptosis genes, so it makes sense that it is  highly expressed in cell lines. The function of mir-184 is not known but it is interesting that it is also higher in cell lines than in normal tissue. mir-124 is a nervous system specific microRNA. It is  not surprising that it is higher expressed in embryos than in (non-neural) cell lines. mir-iab-4 is a developmental regulator. It is located in the Hox cluster and regulates Hox genes.) ::

Another way to get a global overview of the data is to use clustering and plot heatmaps. You can do this with the following command: ::

	heatmap(norm.data, scale="none", cexCol=0.2)

In the resulting plot each library is a column and each microRNA is a row. The color indicates the expression levels, with red being no reads and more yellow indicating higher expression. The dendrogram at the top shows how the libraries cluster together. What can you learn from looking at this plot? 

(There are some problems displaying plots etc. on UPPMAX when running in interactive mode. If you have trouble viewing the PCA plots and heatmaps, try:

- Log out of UPPMAX
- Log into UPPMAX again
- Do not go into interactive mode, just start R
- Type in all R commands again. )
