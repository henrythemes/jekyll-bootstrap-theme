---
layout: default
title:  'Exercise - Functional annotation'
---

# Functional annotation

Functional annotation is the process during which we try to put names to faces - what do the genes do that we have annotated and curated? Basically all existing approaches accomplish this by means of similarity. If a translation product has strong similarity to a protein that has previously been assigned a function, the rationale is that the function in this newly annotated transcript is probably the same. Of course, this thinking is a bit problematic (where do other functional annotations come from...?) and the method will break down the more distant a newly annotated genome is to existing reference data. A complementary strategy is to scan for more limited similarity - specifically to look for the motifs of functionally characterized protein domains. It doesn't directy tell you what the protein is doing exactly, but it can provide some first indication.

In this exercise we will use an approach that combines the search for full-sequence simliarity by means of 'Blast' against large public databases with more targeted characterization of functional elements through the InterproScan pipeline. Interproscan is a meta-search engine that can compare protein queries against numerous databases, including the protein family databases PFam and ProDom. The output from Blast and Interproscan can then be used to add some information to our annotation.

##Prepare the input data
Since we do not wish to spend too much time on this, we will again limit our analysis to chromosome 4. It is also robably best to choose the analysis with ab-initio predictions enabled (unless you found the other build to be more convincing). Maker produces a protein fasta file (called "annotations.proteins.fa") together with the annotation and this file should be located in your maker directory.

create a new folder for the functional annotation:  
mkdir functional\_annotation  
cd functional\_annotation

Now link the annotations.proteins.fa file you want to use into your folder.

## BLAST approach

A 'full' Blast analysis can run for several days and consume several GB of Ram. Consequently, for a huge amount of data it is recommended to parallelize this step doing analysis of chunks of tens or hundreds proteins.

### Perform Blast searches from the command line on Uppmax:

To run Blast on your data, use the Ncbi Blast+ package against a Drosophila-specific database (included in the folder we have provided for you, under blastdb/blastplus) - of course, any other NCBI database would also work:

*module load blast/2.2.29+*
*blastp -db /path/to/blastdb -query annotations.proteins.fa -outfmt 6 -out blast.out -num_threads 8*

Agains the Drosophila-specific database, the blast search takes about 5secs per protein request - depending on how many sequences you have submitted, you can make a fairly educted guess regarding the running time.

### Process the blast outout with Annie
The Blast outputs must be processed by [annie](http://genomeannotation.github.io/Annie/).
First download annie:  
*git clone https://github.com/genomeannotation/Annie.git*

then you should load python:
*module load python/2.7.6*

Annie is a convenient way to attach information to newly annotated transcripts. You can use it to visualize pathways, check for protein domains and other known sequence features. Next, you could write scripts of your own to merge annie output into your annotation. 


## Interproscan approach
Blast searches provide an indication about potential homology to known proteins. Interproscan, on the otherhand, combines a number of searches for conserved motifs and curated data sets of protein clusters etc. This step may also take fairly long time. It is recommended to paralellize it for huge amount of data by doing analysis of chunks of tens or hundreds proteins.

### Perform InterproScan analysis
InterproScan can be run through a website or from the command line on a linux server. Here we are interested in the command line approach:

Annotation -&gt; Interproscan - &gt; Run Interproscan

This will open a list of databases against which to search - the common choices here are HmmPfam, SuperFamil y and perhaps SignalPHmm. Since we can limit the number of sequences to run, it is fine to use all these databases (but perhaps only on 20-30 sequences).

Interproscan uses an internal database that related entries in public databases to established GO terms. By running the built-in option, we can add this information to our data set:

### load the retrieved functional information in your annotation file:
Next, you could write scripts of your own to merge interproscan output into your annotation. Incidentially, Maker comes with a small utility script that can take InterProscan output and add it to a Maker annotation file.

If you now copy this file back to Uppmax and into the same folder where the corresponding Maker gene annotation lives:

*ipr\_update\_gff maker.gff interproscan.tsv &gt; maker.with\_interpro.gff*

Where a match is found, the new file will now include a feature called Dbxref in the transcript feature field.

## What's next?

Because of Makers' compatibility with GMOD standards, an annotation augmented in one or both of this way can be loaded into e.g. WebApollo and will save annotators a lot of work when e.g. adding meta data to transcript models.

