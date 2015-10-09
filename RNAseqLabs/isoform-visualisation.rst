====================================
Identifying isoforms in RNA-seq data
====================================

One of the advantages of RNA-seq data compared to microarrays is that you get 
isoform-level information 'for free' in addition to gene-level information. 
In this exercise, we are going to look at some ways to identify and visualize isoforms.

As you have hopefully read on the introduction page, we will use RNA-seq and quantitative 
mass-spectrometry (MS) data sets from the A431 cell line. These data were measured by a 
group at SciLifeLab Stockholm in a 'proteogenomics' study where the aim was to discover 
new genes or gene variants by deep proteomic profiling using an MS method, and mapping 
the obtained peptides back to the genome. 
The RNA-seq data was obtained to see if there was RNA-level support for the predicted novel 
peptides and transcript variants. We will look at one loci that were flagged by the research 
group as being interesting, and see what the RNA-seq data look like for that gene.


Strategies for using the RNA-seq data
=====================================

There are different ways to find out how the RNA-seq data shows the RAB11FIP5 gene to 
be expressed. Roughly speaking, we can talk about three different strategies:

- Mapping the reads to a reference genome and quantifying the gene and isoform FPKMs

- Mapping the reads to a reference genome and assembling transcripts based on the mapping (reference guided assembly)

- Assembling the reads *de novo*, ignoring the reference genome for now

In order to make these steps computationally feasible during the lab, we have extracted 
only those sequences that mapped to the RAB11FIP5 gene in each sample. These "sub-FASTQ" 
files can be found in ``/proj/b2013006/webexport/downloads/courses/RNAseqWorkshop/isoform/RAB11FIP5_fastqFiles``.
 

To do the reference guided assembly yourself go to `Reference guided assembly using Cufflinks 
<https://export.uppmax.uu.se/b2013006/courses/RNAseq201410/build/html/courseSource/isoform-lab.html>`_. 
This link also contains information on how to quantify already annotated genes and isoforms.

To do the *de novo* assembly yourself go to `Isoform detection using RNA seq *de novo* Assembly 
<https://export.uppmax.uu.se/b2013006/courses/RNAseq201410/build/html/courseSource/isoform-denovo.html>`_.


============================
Visualise isoform data
============================

For the gene **RAB11FIP5** you will now hopefully have generated your own data that you can look at. 
If everything worked you will now have :

 * One BAM file with short reads mapped to the genome 

 * One GTF file  containig differnt isoform of **RAB11FIP5** based on the short reads that mapped.
 
 * One BAM file with **Trinity** assembled transcripts mapped to the genome

Since it takes time to generate all data we have already created other files that you can also download and view in your browser.

Importing reference based isoform info to the **RAB11FIP5** gene
================================================================
This contains results files from the subset of reads that map to the **RAB11FIP5** gene that has then been used for 
reference based assembly of isoforms. 

You can download all BAM files and GTF files for all samples using a webinterface from `here. 
<https://export.uppmax.uu.se/b2013006/downloads/courses/RNAseqWorkshop/isoform/otherData/refBasedAssembly/RAB11FIP5>`_

Importing de novo assembled transcripts mapped to the **RAB11FIP5** gene
========================================================================
This contains results files from the subset of reads that map to the **RAB11FIP5** gene that has then been used for 
de novo assembled transcriptome and then mapped back on the genome. 

You can download all BAM files and GTF files for all samples using a webinterface from `here
<https://export.uppmax.uu.se/b2013006/downloads/courses/RNAseqWorkshop/isoform/otherData/deNovo/BAMfiles>`_


Importing reference based isoform info to the genome
====================================================
This contains results files from all the reads that map to the genome, which then been used for 
reference based assembly of isoforms in the genome. There is also a GTF file with the final merged isoform  
information from all 12 samples.  

You can download all BAM files and GTF files for all samples using a webinterface from `here
<https://export.uppmax.uu.se/b2013006/downloads/courses/RNAseqWorkshop/isoform/otherData/refBasedAssembly/Genome>`_

Importing the peptide track to the **RAB11FIP5** gene and the genome                                                           
===================================================================
As mentioned above, we will compare some identified peptides from a mass-spectrometry 
experiment with RNA-seq data from the same cell line. Let's start by importing the track 
with identified peptides from the MS experiment. 

You can download the BED file containing all mapped peptides to the genomeusing a webinterface from `here
<https://export.uppmax.uu.se/b2013006/downloads/courses/RNAseqWorkshop/isoform/otherData/>`_


From the BED file extension, IGV will automatically know to color the track according to peptide status
(green for annotated peptides, red for novel peptides).


Importing the Pac bio reads mapped to the genome                                                         
================================================
Unfortunately there are no pacbio reads that mapped to the **RAB11FIP5** gene but if you want to compare pacbio  reads from 
full length transcripts and compare them with reference based transcripts 

You can download the BAM file containing all mapped pac bio reads to the genome using a webinterface from `here 
<https://export.uppmax.uu.se/b2013006/downloads/courses/RNAseqWorkshop/isoform/otherData/>`_
























