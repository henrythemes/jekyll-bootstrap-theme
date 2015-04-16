---
layout: default
title:  'Exercise - Functional annotation'
---

# Functional annotation

Functional annotation is the process during which we try to put names to faces - what do the genes do that we have annotated and curated? Basically all existing approaches accomplish this by means of similarity. If a translation product has strong similarity to a protein that has previously been assigned a function, the rationale is that the function in this newly annotated transcript is probably the same. Of course, this thinking is a bit problematic (where do other functional annotations come from...?) and the method will break down the more distant a newly annotated genome is to existing reference data. A complementary strategy is to scan for more limited similarity - specifically to look for the motifs of functionally characterized protein domains. It doesn't directy tell you what the protein is doing exactly, but it can provide some first indication.

In this exercise we will use an approach that combines the search for full-sequence simliarity by means of 'Blast' against large public databases with more targeted characterization of functional elements through the InterproScan pipeline. Interproscan is a meta-search engine that can compare protein queries against numerous databases. The output from Blast and Interproscan can then be used to add some information to our annotation.

##Prepare the input data
Since we do not wish to spend too much time on this, we will again limit our analysis to chromosome 4. It is also probably best to choose the analysis with ab-initio predictions enabled (unless you found the other build to be more convincing). Maker produces a protein fasta file (called "annotations.proteins.fa") together with the annotation and this file should be located in your maker directory.

create a new folder for the functional annotation:  
*cd ~/*  
*mkdir practical4*  
*cd practical4*  

Now link the annotations.proteins.fa file you want to use into your folder. The command will looks like:
*ln -s ../practical2/maker_with_abinitio/annotations/*  

## Interproscan approach
 Interproscan combines a number of searches for conserved motifs and curated data sets of protein clusters etc. This step may take fairly long time. It is recommended to paralellize it for huge amount of data by doing analysis of chunks of tens or hundreds proteins.

### Perform [InterproScan](https://code.google.com/p/interproscan/wiki/DevDocIntroduction) analysis
InterproScan can be run through a website or from the command line on a linux server. Here we are interested in the command line approach.
<u>Interproscan allows to look up pathways, families, domains, sites, repeats, structural domains and other sequence features.</u>  

Launch Interproscan without any option if you want have a look about all the parameters.

- The '-app' option allows defining the database used. Here we will use the PfamA,ProDom,SuperFamily and PIRSF databases.  
- Interproscan uses an internal database that related entries in public databases to established GO terms. By running the '-goterms' option, we can add this information to our data set.
- If you enable the InterPro lookup ('-iprlookup'), you can also get the InterPro identifier corresponding to each motif retrieved: for example, the same motif is known as PF01623 in Pfam and as IPR002568 in InterPro. 
- The option '-pa' provides mappings from matches to pathway information (MetaCyc,UniPathway,KEGG,Reactome).

*module load InterProScan/5.10-50.0*  
*interproscan.sh --input annotations.proteins.fa --seqtype p -dp -pa -appl PfamA-27.0,ProDom-2006.1,SuperFamily-1.75,PIRSF-3.01 -goterms -iprlookup*

The analysis shoud take 2-3 secs per protein request - depending on how many sequences you have submitted, you can make a fairly educted guess regarding the running time.  
You will obtain 3 result files with the following extension '.gff3', '.tsv' and '.xml'. Explanation of these output are availabke [>>here<<](https://code.google.com/p/interproscan/wiki/OutputFormats).


### load the retrieved functional information in your annotation file:
Next, you could write scripts of your own to merge interproscan output into your annotation. Incidentially, Maker comes with utility scripts that can take InterProscan output and add it to a Maker annotation file.  

- ipr\_update\_gff: adds searchable tags to the gene and mRNA features in the GFF3 files.  
- iprscan2gff3: adds physical viewable features for daomains that can be displayed in JBrowse, Gbrowse, and Web Apollo.

If you now copy the .tsv file into the same folder where the corresponding Maker gene annotation lives:

*$SCRIPT\_PATH/ipr\_update\_gff maker.gff interproscan.tsv &gt; maker.with\_interpro.gff*

Where a match is found, the new file will now include features called Dbxref and/or Ontology_term in the gene and transcript feature field (9th column).


## BLAST approach
Blast searches provide an indication about potential homology to known proteins.
A 'full' Blast analysis can run for several days and consume several GB of Ram. Consequently, for a huge amount of data it is recommended to parallelize this step doing analysis of chunks of tens or hundreds proteins. This approach can be used to give a name to the genes and a function to the transcripts.

### Perform Blast searches from the command line on Uppmax:

To run Blast on your data, use the Ncbi Blast+ package against a Drosophila-specific database (included in the folder we have provided for you, under **blastdb/uniprot\_dmel/uniprot\_dmel.fa**) - of course, any other NCBI database would also work:

*module load blast/2.2.29+*  
*blastp -db /path/to/blastdb -query annotations.proteins.fa -outfmt 6 -out blast.out -num_threads 8*

Agains the Drosophila-specific database, the blast search takes about 2 secs per protein request - depending on how many sequences you have submitted, you can make a fairly educted guess regarding the running time.

### Process the blast outout with Annie
The Blast outputs must be processed to retrieve the information of the closest protein (best e-value) found by Blast. This work will be done using [annie](http://genomeannotation.github.io/Annie/).  

First download annie:  
*git clone https://github.com/genomeannotation/Annie.git*  

Then you should load python:  
*module load python/2.7.6*  

Now launch annie:  
*Annie/annie.py sprot blast.out maker.gff /path/to/blastdb/ maker\_annotation.annie*  

Annie writes in a 3-column table format file, providing gene name and mRNA product information. The purpose of annie is relatively simple. It recovers the information in the sequence header of the uniprot fasta file, from the best sequence found by Blast (the lowest e-value).

### load the retrieved information in your annotation file:  

Before to use the script allowing to load the information from Annie output to your annotation file you have to load some PATH to your profile. To do that just launch the following script:  
*$SCRIPT\_PATH/install\_perllib\_missing.sh*  
*source ~/.bash_profile*  

Now you should be able to use the following script:  
*$SCRIPT\_PATH/maker\_gff3manager\_JD\_V6.pl -f maker.with\_interpro.gff -b maker_annotation.annie -o finalOutputDir*  

That will add the name attribute to the "gene" feature and the description attribute (corresponding to the product information) to the "mRNA" feature into you annotation file. This script may be used for other purpose like to modify the ID value by something more conveniant (i.e FLYG00000001 instead of maker-4-exonerate_protein2genome-gene-8.41).  
The improved annotation is a file named "AllFeatures.gff" inside the finalOutputDir.



## What's next?

Because of Makers' compatibility with GMOD standards, an annotation augmented in one or both of this way can be loaded into e.g. WebApollo and will save annotators a lot of work when e.g. adding meta data to transcript models.

