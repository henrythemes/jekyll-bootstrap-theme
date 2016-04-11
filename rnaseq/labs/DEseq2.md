---
layout: default
title:  'Differential Expression'
---

# Differential expression analysis of RNA-seq data using DEseq2



## Data set

In this exercise we are going to look at RNA-seq data from the A431 cell line. A431 is an epidermoid carcinoma cell line which is often used to study cancer and the cell cycle, and as a sort of positive control of epidermal growth factor receptor (EGFR) expression. A431 cells express very high levels of EGFR, in contrast to normal human
fibroblasts. In the experiment we are looking at today, A431 cells were treated with gefinitib, which is an EGFR inhibitor, and is used (under the trade name Iressa) as a drug to treat cancers with mutated and overactive EGFR. In the experiment, RNA was extracted at four time points: before the gefinitib treatment (t=0), and two, six and twenty-four hours after treatment (t=2, t=6, t=24, respectively), and sequenced using an Illumina HiSeq instrument in triplicates (thus there are 3x4=12 samples). Since gefinitib suppresses the expression of EGFR, the effects on the cell line are quite severe, as you will see when comparing the 24h time point to the control.

Some of the files that we will use have been pre-constructed by us ahead of the lab because it would take too much time to do it during the lab exercise. In those cases, we have indicated the commands that we used to generate the files, so that you should be able to reproduce the results on your own.

## DESeq

The DESeq approach identifies differentially expressed genes based on counts of the number of reads mapped to each gene. DESeq is not limited to RNA-seq, but can be used for comparisons of other count-based data, such gene expression profiles from tag sequencing or data from ChIP-seq experiments.

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


## Further reading

The algorithms used by DESeq2 are described in the papers by
[Anders and Huber (2010)](http://genomebiology.com/content/11/10/R106) and
[Love et al. (2014)](http://biorxiv.org/content/early/2014/05/27/002832).

They have also published descriptions of how to use their tools in *Nature Protocols*: [Anders et al. (2013)](http://www.ncbi.nlm.nih.gov/pubmed/23975260)

