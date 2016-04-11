---
layout: default
title:  'Differential Expression'
---

### Differential expression analysis of RNA-seq data

## Overview


There are many software packages for differential expression analysis of RNA-seq data. We will first look at **Cuffdiff/CummeRbund** and then try the alternative method **DESeq**.

Several tools, such as DESeq and edgeR, start from read counts per gene and use the discrete nature of the data to do statistical tests appropriate for this kind of data. It can be argued that that such counts will never give quite the correct results because the presence of alernative isoforms confound the read counting. Cuffdiff therefore combines isoform-level quantification and differential expression testing into one framework and claim that they achieve better results because they are able to take into account the uncertainty of isoform quantification.

## Data set

In this exercise we are going to look at RNA-seq data from the A431 cell line. A431 is an epidermoid carcinoma cell line which is often used to study cancer and the cell cycle, and as a sort of positive control of epidermal growth factor receptor (EGFR) expression. A431 cells express very high levels of EGFR, in contrast to normal human
fibroblasts. In the experiment we are looking at today, A431 cells were treated with gefinitib, which is an EGFR inhibitor, and is used (under the trade name Iressa) as a drug to treat cancers with mutated and overactive EGFR. In the experiment, RNA was extracted at four time points: before the gefinitib treatment (t=0), and two, six and twenty-four hours after treatment (t=2, t=6, t=24, respectively), and sequenced using an Illumina HiSeq instrument in triplicates (thus there are 3x4=12 samples). Since gefinitib suppresses the expression of EGFR, the effects on the cell line are quite severe, as you will see when comparing the 24h time point to the control.

Some of the files that we will use have been pre-constructed by us ahead of the lab because it would take too much time to do it during the lab exercise. In those cases, we have indicated the commands that we used to generate the files, so that you should be able to reproduce the results on your own.


## CuffDiff and cummeRbund

** ATTENTION!! Some of the steps do not work on uppmax! We recommend that you run CummRBund this on your local computer **



We will start by looking at Cuffdiff and CummeRbund. These tools are part of the "Tuxedo" suite, which also includes Bowtie, TopHat and Cufflinks.

Cuffdiff, which you have already tried in an earlier exercise, is a command-line program that does the actual differential expression testing, and cummeRbund is an R package that reads Cuffdiff output and visualizes it in various nice ways.

Running Cuffdiff on the A431 RNA-seq data set would take too long for the purposes of this exercise and we have therefore executed it in advance so that you can access the results. For reference, these were the commands used:

	# Just for reference; do not run this
 	ALI=/path/to/TopHat/alignments
 	cuffdiff -o cuffdiff_all --labels 0h,2h,6h,24h -p 8 \
 	Homo_sapiens.GRCh37.71.gtf \
 	$ALI/sample_1/accepted_hits.bam,$ALI/sample_2/accepted_hits.bam,$ALI
 	sample_3/accepted_hits.bam \
 	$ALI/sample_4/accepted_hits.bam,$ALI/sample_5/accepted_hits.bam,$ALI
 	sample_6/accepted_hits.bam \
 	$ALI/sample_7/accepted_hits.bam,$ALI/sample_8/accepted_hits.bam,$ALI
 	sample_9/accepted_hits.bam \
 	$ALI/sample_10/accepted_hits.bam,$ALI/sample_11/accepted_hits.bam,$ALI
 	sample_12/accepted_hits.bam

It may be advantageous to download the output to your own computer, install the cummeRbund package (after installing R if you don't have it) and do the visualizations there. Uppmax has R installed. If you haven't installed R ahead of the course, please find the appropriate download for your [here](http://ftp.sunet.se/pub/lang/CRAN/). You might as well take the latest version (3.2.2).

If you want to work on Uppmax (note this is not necessary if you are running R on your own computer), login as we did yesterday. 

You can find a (gzipped archive of a) directory with CuffDiff output for a pairwise  differential expression analysis between all the time points on Uppmax at  `/proj/b2013006/webexport/downloads/courses/RNAseqWorkshop/diffExp/cuffdiff_all.tar.gz`. 

Or you can download it [here](https://export.uppmax.uu.se/b2013006/downloads/courses/RNAseqWorkshop/diffExp/cuffdiff_all.tar.gz)

Remember to untar and unzip it:

    tar -zxf cuffdiff_all.tar.gz
 
Then start R (either on your own system, or on Uppmax by simply typing
``R``) and install cummeRbund with these commands:

    source("http://bioconductor.org/biocLite.R")
    biocLite("cummeRbund")

Load the package:

    library(cummeRbund)

If you were running R on your local computer, you could bring up a tutorial with the command ``vignette("cummeRbund-manual")``. However, this does not work over the remote connection we are using for this course. Go to the [CummeRbund web page](http://compbio.mit.edu/cummeRbund/) to access the manual.

As you can see in the manual, CummeRbund provides lots of ways to visualize Cuffdiff output. We will just look at some of the possible visualizations, but it should be easy for you to try other ones from the manual after you have gone through the examples here.

Start by loading the Cuffdiff results. The command below should work if your R working directory contains the folder with Cuffdiff output (``cuffdiff_all``); otherwise you will have to specify the path to this folder.

    cuff <- readCufflinks("cuffdiff_all") 

This will take a while the first time you run it. For subsequent runs on the same data, it will be faster because cummeRbund will have stored an SQL database containing the information it needs to disk. ()If you have problems getting CummeRbund to work, updating Bioconductor, RSQLite and cummeRbund might help:

	source("http://bioconductor.org/biocLite.R")
	biocLite("BiocUpgrade")
	biocLite("RSQLite")
	biocLite("cummeRbund")

As mentioned above, the ``cuffdiff_all`` directory contains information about differential expression for genes and isoforms between all pairs of time points. We could start by looking at an overview of how many genes that are differentially expressed at a false discovery rate (FDR) of 1%:

    mySigMat <- sigMatrix(cuff, level="genes", alpha=0.01) # constructs the plot
    mySigMat # displays the plot

It would work to just give the ``sigMatrix(cuff, level="genes", alpha=0.01)`` command but it might be nice to keep the plot as an R object as done above.

The resulting plot shows how many genes that are differentially expressed between each pair of time points. It is easy to show numbers corresponding to isoforms instead; this is generally the case for all cummeRbund functions:

    sigMatrix(cuff, level="isoforms", alpha=0.01)

One way to see which genes that are differentially expressed in different comparisons is the following:

    mySigTable <- getSigTable(cuff,alpha=0.01,level="genes")

So for the 24h vs control comparison, you could select the ``X0hvsX24h`` (R doesn't want variable names to start with numbers so it adds an X in those cases) column from ``mySigTable`` and filter to those genes that have a 1 in that column:

    diff_0_24 <- mySigTable[,"X0hvsX24h"]
    cuffgenes.24h.ctrl <- names(diff_0_24[which(diff_0_24==1)])

Now you have a list of differentially expressed genes between 24h and control. This list could be saved to file using e g:

    write.table(cuffgenes.24h.ctrl, file="cuffgenes.24h.ctrl.txt", quote=F, row.names=F, col.names=F)

If you would just like to plot the expression of one specific gene across the time points, you could do something like the following. Let's use the *EGFR* gene since we know that it should be affected by the treatment in the experiment:

    myGene <- getGene(cuff, "EGFR")
    expressionPlot(myGene) # Will collapse replicates and only show gene level FPKM
    expressionPlot(myGene, replicates=T) # Will show replicate FPKMs
    expressionBarplot(myGene, replicates=T) # Show as bar plot instead

We can also plot FPKMs for all isoforms:

    expressionPlot(isoforms(myGene), replicates=T)
    expressionBarplot(isoforms(myGene), replicates=T) # Might be quite cluttered

Perhaps we are interested in genes that behave in the same way as our gene of interest. In that case we could use:

    egfr.similar <- findSimilar(cuff, "EGFR", n=20) # Will find the 20 genes that are most similar to EGFR 
    egfr.similar.expr <- expressionPlot(egfr.similar,logMode=T,showErrorbars=F)
    egfr.similar.expr

(This can take a while to run.) If we are not interested in a particular pattern but would like to find out which patterns seem to be present in the data, we can use cummeRbund's built-in clustering method. Let's say we assume that there are 10 different temporal patterns (the choice is up to you of course) and restrict the clustering to genes that are considered to be differentially expressed between at least one pair of time points. This will probably take a while to run:

    sig.gene.ids <- getSig(cuff,alpha=0.01,level='genes')
    sig.genes <- getGenes(cuff, sig.gene.ids)
    cl <- csCluster(sig.genes,k=10)
    clplot <- csClusterPlot(cl)
    clplot

Finally, we can look at some plots that visualize how similar the samples are to each other. In a real project we would probably have started by doing this at the beginning to check if the data look OK.

    MDSplot(genes(cuff),replicates=T)    

The MDS (multidimensional scaling) plot attempts to visualize the high-dimensional data (tens of thousands of genes) in a two-dimensional plot so that the distances between each sample are preserved as faithfully as possible.

There is also a heatmap function for samples. Of course there are also various heatmaps for gene lists; you can read more about that in the manual.

    csDistHeat(genes(cuff)) 

## DESeq

As mentioned above, the DESeq approach identifies differentially expressed genes based on counts of the number of reads mapped to each gene. DESeq is not limited to RNA-seq, but can be used for comparions of other count-based data, such gene expression profiles from tag sequencing or data from ChIP-seq experiments.

The DESeq method is implemented in the R packages DESeq and DESeq2. The latter is more recent, and recommended. The DESeq2 package is also available in several versions, tied to different versions of R (this applies to all Bioconductor packages). To use the most recent version of DESeq2, make sure you have the most recent R version installed. Also note that DESeq2 strictly requires R version
3.0 or above.

For this exercise we have pre-calculated read counts per gene (according to Ensembl annotations) with commands like:

     # Only given for reference, not supposed to be executed during the lab
     samtools view accepted_hits_137_1.bam | sort > accepted_hits_prehtseq_137_1.sam
     htseq-count -s no -q accepted_hits_prehtseq_137_1.sam Homo_sapiens.GRCh37.71.gtf > 137_1.counts

and combined the counts into a single table. You will import this table into R and use DESeq2 to get a list of differentially expressed genes. You can get the count table on UPPMAX here: ``/proj/b2013006/webexport/downloads/courses/RNAseqWorkshop/diffExp/count_table.txt``

Or you can download it [here](https://export.uppmax.uu.se/b2013006/downloads/courses/RNAseqWorkshop/diffExp/count_table.txt)

Start R and load the DESeq2 package:

    library(DESeq2)

If this does not work, you may need to go through the usual drill to install
the package:

     source("http://bioconductor.org/biocLite.R")
     biocLite("DESeq2")

The actual analysis is rather simple, after you have set up the data you are going to feed to DESeq2. Start by reading the file ``count_table.txt``. Of course you need to be in the same directory as the file for the following command to run cleanly:

    counts <- read.delim("count_table.txt")

You may want to look at the table with commands like ``head(counts)``. Next, you need to create a table with information about the samples:

    samples <- data.frame(timepoint = rep(c("ctrl", "t2h", "t6h", "t24h"), each=3))

Look at the content of the data frame that you created. (Type the name of the object, ``samples``, and press enter.) Note that the table only has one column, which indicates the time point for each of the 12 samples. For this simple experimental design, this is all we need: the timepoints define the groups that we wish to to compare. It doesn't really matter what you call the groups, as long as the names are distinct. For example, we could have used "t0h" instead of "ctrl" for the first time point.

If instead of coming from a cell line, these samples were (say) tumor samples from different patients, such that for example samples 0h_1, 2h_1, 6h_1 and 24h_1 were all from the same person at different time points, the sample description could be extended by one column:

    # Not to be used in the lab - just an example!
    samples <- data.frame(timepoint = rep(c("ctrl", "t2h", "t6h", "t24h"), each=3), individual=rep(1:3, 4))

This would facilitate a so-called factorial design and specifying it for DESeq2 would potentially give more statistical power than just comparing groups "blindly". However, we are not going to do this here. Now it's time to construct the data set object that DESeq2 needs to perform the analysis:

    ds <- DESeqDataSetFromMatrix(countData=counts, colData=samples, design=~timepoint)

This function call constructs a DESeq2 data set object using the arguments we provide: (1) count table; (2) sample description, and (3) experimental design.  All three arguments are mandatory.  The design is specified as a *formula* (another type of R object). In this case the formula is easy: the timepoints are really the only thing we can compare to each other. If we had had an additional factor as described above, we could have chosen to test for differences between timepoints while correcting for variability arising from individual differences, or to test for differences between individuals, while correcting for variation
arising from the timepoints. For example:

    # Just an example
    ds <- DESeqDataSetFromMatrix(countData=counts, colData=expr.desc, design=~timepoint + individual) # to test for differences between individuals    
    ds <- DESeqDataSetFromMatrix(countData=counts, colData=expr.desc, design=~individual + timepoint) #	to test	for differences	between	timepoints 

It can be useful to include the sample names in the data set object:

    colnames(ds) <- colnames(counts)

Now that we are set, we can proceed with the differential expression testing:

    ds <- DESeq(ds)

This very simple function call does all the hard work. Briefly, this
function performs three things:

- Compute a scaling factor for each sample to account for differences
  in read depth and complexity between samples
- Estimate the variance among samples
- Test for differences in expression among groups (time points, in our case)

For more details, see the manual page for the function:

  ? DESeq

You can also have a look at the manual for the DESeq2 package, which can be found on the [DESeq2 web page](http://www.bioconductor.org/packages/release/bioc/html/DESeq2.html). If you were running R locally, you would also be able to bring up the manual with the command ``vignette("DESeq2")``.

Now we just need to extract the results. Recall that we have expression data for four different time points. We can use the function ``results()`` to see the results of comparing two time points. Choose two time points that you would like to compare and give a command like:

    res <- results(ds, c("timepoint","t24h","ctrl"))

The object *res* is now a table with test results for each gene in the original count table, i.e. all annotated genes, both protein-coding and non-coding.  Use the function ``head()`` to inspect the results table.

Do you understand what the columns mean? You can see information such as the "base mean" (an average of the normalized mean counts per group), the log2 fold change between the groups, and the P-values (both "raw" and adjusted for multiple comparisons). If you are unsure how to interpret these data, discuss with one of the instructors.

We are not interested in results for all genes. For example, genes with zero or very low counts across all samples cannot be tested for differences in expression. Those will have a P-value of NA (not applicable). There are about 32,000 such genes, which we remove:

    nrow(res)
    sum( is.na(res$pvalue) )
    res <- res[ ! is.na(res$pvalue), ]
    nrow(res)

You probably want to focus on genes that are significant according to some criterion, such as false discovery rate (FDR) or log fold change. Filtering on the adjusted P-value (column *padj*) is equivalent to choosing a desired false discovery rate. For example, we can filter the results such that 1% are expected to be false positives (genes with no actual difference in expression):

    sig <- res[ which(res$padj < 0.01), ]

(Note: in some versions of DESeq, you may need to use ``res$FDR``
instead of ``res$padj``.)

How many significantly differentially expressed genes do you get? (Hint: try the ``dim()`` and ``nrow()`` functions).

To see the top "hits", we can sort the filtered result list by
statistical significance:

    sig <- sig[ order(sig$padj), ]

We convert the table to a data frame, so that we can manipulate and
view it more easily:

    sig <- as.data.frame(sig)
    head(sig, n=20)

If the table wraps over several lines, you can try to change some R
options before viewing the table:

    options(width=120)  ## Display width (number of characters)
    options(digits=5)   ## Number of digits to show for numbers
    head(sig, n=20)

You might want to compare the results from CuffDiff and DESeq2. The identifiers of the significant genes from DESeq2 can be easily obtained by:

    deseqgenes.24h.ctrl <- rownames(sig)

If you still have the list of significant genes between 0h and 24h from the CuffDiff/cummeRbund analysis in your R session, or if you have saved it to file, you can check how many of them that were picked up by both programs::

    common.24h.ctrl <- intersect(deseqgenes.24h.ctrl, cuffgenes.24h.ctrl)

The gene identifiers we work with above are Ensembl gene IDs. These are useful as unique identifiers, but does not tell us anything about what the genes do. One way to find out more about individual genes is to look up the identifiers at the [Ensembl web site](http://www.ensembl.org). We can also use the R package ``biomaRt`` to download a table of corresponding gene symbols (i.e. short gene names) from Ensembl. First we load the package:

    library(biomaRt)

If this does not work, you may need to install the package:

    source("http://bioconductor.org/biocLite.R")
    biocLite("biomaRt")
    library(biomaRt)

When you have successfully loaded the package, run the following:


    ensembl  <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL",  dataset="hsapiens_gene_ensembl")
    genemap <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                     filters = "ensembl_gene_id",
                     values = rownames(sig),
                     mart = ensembl )

(This step can also take a few minutes to run.) The data frame *genemap* now contains a mapping of Ensembl gene IDs to
gene symbols:

    head(genemap)

Let's add these gene symbols to the result table:

    symbols <- tapply(genemap$hgnc_symbol, genemap$ensembl_gene_id, paste, collapse="; ")
    sig$symbol <- symbols[ rownames(sig) ]
    head(sig)

(The ``tapply()`` function call above is needed to deal with cases where
there are multiple symbols for the same gene. This call maps each
Ensembl gene ID to a string of one more more gene symbols separated by
semi-colon.)

The DESeq2 package contains a function plotMA() that can be used to
visualize the differences in gene expression:

    plotMA(ds)

Do you understand what this plot shows? Look at the manual page for
the function, and run it again with the argument alpha set to
different values. Discuss with an instructor if you are unsure how to
interpret the plot.

If time allows, have a look at the [RNA-seq analysis workflow example](http://www.bioconductor.org/help/workflows/rnaseqGene/) on the BioConductor web site.
There is a section called [Visually exploring the dataset](http://www.bioconductor.org/help/workflows/rnaseqGene/#eda) on exploratory analysis of count data after regularized log transformation. This section shows how to make several plots that are useful for exploring RNA-seq data sets. If you don't have time to go through it now, try these commands and admire the resulting plots:

	## Apply regularized-log transform to counts
	rld <- rlogTransformation(ds)
	
	## Principal component analysis
	plotPCA(rld, intgroup="timepoint")
	
	## Heatmap of sample distances
	library("gplots")   # If this fails, run: install.packages("gplots")
	library("RColorBrewer")
	sampleDists <- dist(t(assay(rld)))
	sampleDistMatrix <- as.matrix( sampleDists )
	colours <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
	heatmap.2(sampleDistMatrix, trace="none", col=colours)
	
	## Heatmap of 35 most variable genes
	library("genefilter")
	topVarGenes <- head(order(rowVars(assay(rld)), decreasing=TRUE), 35)
	heatmap.2(assay(rld)[topVarGenes, ], scale="row",
	trace="none", dendrogram="column", margins=c(5, 10),
	col=colorRampPalette(rev(brewer.pal(9, "RdBu")))(255))

You may also want to try some of the examples from the cummeRbund manual. For example, the function s csDensity and fpkmSCVPlot produce nice looking and informative plots.

## Further reading

The algorithms used by Cuffdiff and DESeq are described in the papers by
[Trapnell et al. (2013)](http://www.ncbi.nlm.nih.gov/pubmed/23222703),
[Anders and Huber (2010)](http://genomebiology.com/content/11/10/R106) and
[Love et al. (2014)](http://biorxiv.org/content/early/2014/05/27/002832).

Both groups have also published descriptions of how to use their tools in *Nature Protocols*: [Trapnell et al. (2012)](http://www.ncbi.nlm.nih.gov/pubmed/22383036), [Anders et al. (2013)](http://www.ncbi.nlm.nih.gov/pubmed/23975260)

For a recent review and evaluation of a range of methods for
normalization and differential expression analysis of RNA-seq data,
see [Rapaport et al. (2013)](http://genomebiology.com/content/14/9/R95)
