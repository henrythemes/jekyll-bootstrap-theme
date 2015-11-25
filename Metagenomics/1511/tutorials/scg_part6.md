---
layout: default
title:  'Part 6: Analysis of a novel single cell'
---

# Part 6: Analysis of a novel single cell genome
---

Now that you have worked with several single cell genome datasets, you will be able to play around with the analysis you did before for a mystery single cell (dataset3). 
With the knowledge that you have acquired during this course, you can try to:

- make a decent assembly
- assess the completeness of your single cell
- try to identify your cell using rnammer or MEGAN
- try to asses if the dataset suffers from contamination
- find out any interesting features that might be encoded by the genome

In this exercise, you will work with MiSeq data produced from a different Single-cell Amplified Genome (SAG) than G5 in the first part of the tutorial. You can do the assemblies in the same way as the examples shown in Part 3. There is no need to optimize the steps, instead think about what is your first choice of assembler, settings and preprocessing with the experience you have now.

An overview of the steps are listed below:

0. File/folder structure and naming convention
1. Preprocessing
2. Assembly optimization (tool, kmers, flags)
3. Assembly stats from *'Quast'* and *'micomplete'* 
4. Prodigal predicttion of ORFs
5. RNAmmer predicttion of rRNAs
6. Blastn and Blastp of the genes
7. Blastn of the rRNA genes against Silva database
8. Overview of the results


We suggest that you start with question 1 and think of the steps necessary to obtain the answer. 

If you have more time you can play around with the other optional questions and decide where you want to focus. There are also a couple of optional exercises you can try for the different steps you went through. Choose the part you are most interested in.

### Questions
---

**Q7.1** What is the taxonomic affiliation of the SAG?  

Optional exercises

Biologically oriented questions:

**Q7.2** Can you say anything about the metabolism based on the assembled data?  

Bioinformatically oriented questions:

**Q7.3** Did you notice any major differences between the assemblies using different assembler, setting or input data (trimmed/untrimmed reads)? To make this question more specific you can look into the optional exercises below where some assembly optimization options are suggested.  
**Q7.4** Do you think the number of reads are enough to obtain a good assembly or should more sequences be obtained?  

## Optional exercises

Here is a list of optional exercises collected in one place, from the various parts of the tutorial that you went through.

Preprocessing - visualization

To visualize the quality of the reads you can use [FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) which provide plots for several checks together with some guidelines on which results might be suspicious. You can see for example a plot of qualities along the read length, look at the duplication level, and so on. You have learned this yesterday and you can try it for a single cell data if you like.

Preprocessing - merging reads

Another step you can do if your library setup is such that sequencing reads should be overlapping, is merging them. An example of how to do that is described [here](scg_part3_merging). Considering that all of the assemblers we use can take in paired reads, and some of them (Spades) actually do not recommend using the qualities that the merging result in, we skip this for the main assembly comparison. It can still be a useful step for other purposes.

Ray assembly optimization:

If you have time you can investigate the influence of various kmer lengths on the assembly results. Try for example using kmers increasing in steps of 10 from 30 to 64, which is the hard-coded limit. Check Ray log output file to make sure about the kmer actually used. If number larger than the threshold is given Ray changes it to the maximum allowed, makes a not of it in the log and proceeds. 

Spades assembly optimization:

If you have time you can investigate the influence of the flags on the assembly time and results. '--careful' flag uses *'bowtie'* tool to map the reads back to the contigs and check for errors due to bad quality sequences and correct these errors. This results in longer assembly times but should improve the results, especially for reads that were not pre-processed. SPAdes can handle single-cell genomic data that is known to be highly biased in terms of sequence coverage along the length of the genome by using the '--sc' flag.

Now you can also explore the influence of k-mers choice. SPAdes in default mode runs with **k-mers of 21, 33, and 55**. 
Try setting **k-mers to 55, 77, and 99**. To set these k-mers, you need to provide this parameter when running SPAdes:

```
-k 55,77,99
```


<div>
 <span style="float:left"><a class="btn btn-primary" href="scg_part5_2"> Previous page</a></span>
 <span style="float:right"><a class="btn btn-primary" href="#"> Next page</a></span>
</div>

