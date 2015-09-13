---
layout: default
title:  'Resequencing Extra'
---


# Extra

Except for Extra 2 and Extra 3 that need to be run in that order, the different parts are independent and can be run in any order you want. 

##Extra1: View data in the UCSC browser, extract interesting annotation tracks

The UCSC browser is a collection of tools for viewing and manipulating genomic data. Here we will only use a very limited set of all the functionalities, but there is extensive documentation that allows you to find all the possibilities. So while you have the data loaded into the UCSC browser, take the opportunity to look around and explore. 

Browser: https://genome.ucsc.edu

Examples of documentation :

http://www.genome.ucsc.edu/training/index.html
http://www.genome.ucsc.edu/training/ucscGeneFishing.pdf

Here we will do two simple things, upload and view your data and any annotation you find interesting. Download a bed file with annotations for a region of the genome. 

Step 1: Go the webpage for the browser

Step 2: Click on genome browser on the upper left

Step 3: Choose genome and assembly and the position you are interested in

Step 4: Scroll down to see the annotation tracks that are available for this region. Add some, remove some from what is shown according to your liking. 

Step 5: Just under the viewer window is on option “manage custom tracks”. Click there and upload your vcf-file with variants. 

Step 6: Use the zoom in/out buttons to get the optimal view of your data. 

Step 7: Click the home-button. On the left choose Table browser. 

Step 8: Choose genome and assembly and the position you are interested in. Also choose which track you are interested in ( one option Group: Genes and gene prediction, track RefSeq Genes). Chose output format BED, supply a file name and press get output. Save the file and use it in Extra 1. Do it again if you are interested in other annotation tracks. 

## Extra 2: Select subsets with bedTools

Originally BEDTools was written to handle the tab delimited .bed format where each line describes a feature. The three first columns are required and are chromosome - start coord - end coord. The are 9 additional optional columns for more information about the feature. The format is more closely described in http://genome.ucsc.edu/FAQ/FAQformat. Currently BEDTools supports input data in more formats including .bam, .bed, .gff or vcf. The BEDTools suite includes a large number of different commands, but a good point to start is intersectBed that can be used to find the common regions between two files. To learn more about this and other commands use either the flag -h or read more on the online documentation http://bedtools.readthedocs.org/en/latest/.

Task: As input use the .vcf file with your variants and the annotation track from the UCSC browser you downloaded in Extra 2. Or use the bed- file with RefSeq gene annotation in /proj/g2015006/labs/gatk/other/RefSeq_genes.bed. Play around the different options and try to answer some questions. How many variants are located within genes? How many variants are not located within genes? How many genes have variants within them?  

Useful hints:

Notice the -v option that gives you regions that do NOT overlap

The command wc -l filename will give you the number of lines in a file

## Extra 3: Annotate variants with annovar

"ANNOVAR is an efficient software tool to utilize update-to-date information to functionally annotate genetic variants detected from diverse genomes"

Annovar have a extensive online documentation (http://annovar.openbioinformatics.org/en/latest/). In this exercis you will use the program to annotate the variants in you .vcf file with gene annotation - is the variant located within a gene? What is the potential effect of the variant on that gene? The program also supplies additional annotation tracks that can be used according to your preferences.

Step 1: Annovar is available as a module on uppmax with bioinfo-tools. Load this module ( module add bioinfo-tools annovar)

Step 2: Annovar uses downloaded annotation tracks to annotate the input file. The program it self has built in utilities download tracks. But in an initial step we will use the annotation tracks already downloaded in the uppmax installation of annovar. The command table_annovar.pl can be used to perform the actual annotation. It takes as input the file that you want to annotate, the folder with annotations, the build version (default hg18, we use hg19), which protocol(s) you want to annotate the file with and what kind of operations that are included. Type table_annovar.pl -h for a more detailer description about the parameters. 

```bash
table_annovar.pl MERGED.illumina.low_coverage.17q.filtered.vcf /sw/apps/bioinfo/annovar/2014.11.12/milou/humandb/ -buildver hg19 -outfile myanno -vcfinput --protocol refGene -operation g
```

As output you get a number of files with the base name as given for the -outfile option. Annovar have transformed the .vcf file to a annovar input file, .avinput. In the variant_function file, the first and second column annotate variant effects on gene structure and the genes that are affected, yet the other columns are reproduced from avinput file. In the exonic_variant_function file, the first, second and third column annotate variant line number in avinput file, the variant effects on coding sequences and the gene/transcript being affected, yet the other columns are reproduced from avinput file. In the vcf file the annotations are included in the info-field of the .vcf file. You also get similar information in the tab-delimited file. From these files can you extract the information about how many variants in the .vcf files that are located within exons? How many are located within genes? How many are nonsynonymous?  

Step 3: If you would like to annotate the data with more types of annotation read more on the annovar website and create a humandb/ folder of your own where you download what you are interested in and then run the annotations. If you download large tracks, be sure to remove them again when you do not need them any more to avoid filling up the space on uppmax. 

## Extra 4: Run pipeline as a script

What you have just been doing could be described as a mapping and variant calling pipeline. How would you do if you would do this again but more efficiently? If you would do if for a large number of input files?

Instead of running interactively on just on sample at a time, the common workflow is to use the slurm queuing system and submit separate batch jobs to the cluster. All your commands would then be predefined and run on a core / node depending on what you ask for. The results will be written to files according to your commands the same way as if you would run them interactively. 

http://www.uppmax.uu.se/slurm-user-guide

If you have just one samples, you could write the commands by hand in the run script. If you would like to automate it and be able to run the same thing for a large number of samples you can use different scripting languages to do this. A simple option is to write it in bash, perl och python. Try for your self and think about different ways to make it efficient.
