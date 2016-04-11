---
layout: default
title:  'Differential Expression'
---




# CuffDiff and cummeRbund


Some of the files that we will use have been pre-constructed by us ahead of the lab because it would take too much time to do it during the lab exercise. In those cases, we have indicated the commands that we used to generate the files, so that you should be able to reproduce the results on your own.

** ATTENTION!! Some of the steps do not work on uppmax! We recommend that you run CummRBund this on your local computer **



We will work with Cuffdiff and CummeRbund. These tools are part of the "Tuxedo" suite, which also includes Bowtie, TopHat and Cufflinks.

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

It may be advantageous to download the output to your own computer, install the cummeRbund package (after installing R if you don't have it) and do the visualizations there. Uppmax has R installed. If you haven't installed R ahead of the course, please find the appropriate download for your [here](http://ftp.sunet.se/pub/lang/CRAN/). 


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

You may also want to try some of the examples from the cummeRbund manual. For example, the functions csDensity and fpkmSCVPlot produce nice looking and informative plots.

## Further reading

The algorithms used by Cuffdiff are described in the papers by
[Trapnell et al. (2013)](http://www.ncbi.nlm.nih.gov/pubmed/23222703),

They have also published descriptions of how to use their tools in *Nature Protocols*: [Trapnell et al. (2012)](http://www.ncbi.nlm.nih.gov/pubmed/22383036)

For a recent review and evaluation of a range of methods for
normalization and differential expression analysis of RNA-seq data,
see [Rapaport et al. (2013)](http://genomebiology.com/content/14/9/R95)
