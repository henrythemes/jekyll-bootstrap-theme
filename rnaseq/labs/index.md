---
layout: default
title:  'Index'
---

# Welcome to NBIS RNA-seq tutorial packages

This page contains links to different tutorials that are used in the RNA-seq course. Some of the tutorials 
are well documented and should be easy to follow. We also supply more beta versions of labs that requires more 
from the user and may contain errors. 


### Introduction

In the links below there are information about tools and data that we will use during the other labs. Please make sure you know the data and how to use **R** and **IGV** before you proceed with the other labs. 

*	[Introduction to the RNA seq data provided](intro)  
*	[Short introduction  to R](R_intro)  
*	[Short introduction to IGV](IGV) 

### Mapping reads 

This contains information on how to map reads to a reference sequence. In the tutorial you will learn how to use both **STAR** and **HiSAT2**
 
*	[Tutorial for mapping reads to a reference and converting them to the BAM format](mapping_reads) 

### Transcript assembly

This contains information regarding how to assemble short reads into transcripts

*	[Tutorial for reference guided assembly](isoform-lab)  
*	[Tutorial for *de novo* assembly](isoform-denovo)

### Visualise mapped reads and assembled transcripts on reference

When reads have been mapped to a reference and/or assembled to transcripts it is always a good idea to check on a reference what the results look like.
 
*	[Tutorial for isoform-visualisation](isoform-visualisation)  

### Quality control laboratory
Before doing any other analysis on mapped RNA-seq reads it is always important to do quality control of your mapped reads and that you do not have any obvious errors in your RNA-seq data 

*	[Tutorial for RNA seq Quality Control](QC_lab)   

### Small RNA analysis
When working with small RNA RNA-seq reads, this case miRNA, there are some analysis that are different  This will be covered in this labs.  

*	[Tutorial for small RNA analysis](smallRNA-lab)


### Differential expression analysis


There are many software packages for differential expression analysis of RNA-seq data.

Several tools, such as DESeq and edgeR, start from read counts per gene and use the discrete nature of the data to do statistical tests appropriate for this kind of data. It can be argued that that such counts will never give quite the correct results because the presence of alernative isoforms confound the read counting. Cuffdiff therefore combines isoform-level quantification and differential expression testing into one framework and claim that they achieve better results because they are able to take into account the uncertainty of isoform quantification. 

*	[Tutorial for differential expression analysis using DEseq2](DEseq2)
*	[Tutorial for differential expression analysis using Sleuth](kallisto)
*	[Tutorial for differential expression analysis using CuffDiff and CummRBund](CuffiDff)
*	[Tutorial for differential expression analysis using multi variate analysis in SIMCA](Simca_tutorial)
 
## Beta labs 
There are some labs that are more close to the cutting edge of analysis and therefore are not as well tested as the ones above. These are tools that have high potential and will most likely, if they hold, will be moved to the mature labs.
 
*	[Differential expression analysis using kallisto](kallisto)
*	[Single cell RNA PCA and clustering](Single_cell_RNA_PCA_and_Clustering)
 
 
 
## Caveat

We will try to keep these tutorials up to date. If you find any errors or things that you think should be updated please contact Johan (johan.reimegard@scilifelab.se) 
  		
