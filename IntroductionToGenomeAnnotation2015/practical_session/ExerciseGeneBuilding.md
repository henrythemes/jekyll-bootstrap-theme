---
layout: default
title:  'Exercise Gene Building'
---

# Running the Maker gene build pipeline
## Overview

**Maker2** is a computational pipeline to automatically generate annotations from a range of input data - including proteins, ESTs, RNA-seq transcripts and ab-initio gene predictions. During this exercise, you will learn how to use Maker with different forms of input data, and how to judge the quality of the resulting annotations.
## Getting the input data

The Maker pipeline can work with any combination of the following data sets:

* Proteins from the same species or related species  

* Proteins from more distantly related organisms (e.g. Uniprot/Swissprot)  

* EST sequences from the same species or very closely related species  

* RNA-seq data from the same or very closely related species - in the form of splice sites or assembled transcripts  

* Ab-initio predictions from one or more tools (directly supported are: Augustus, Snap, GeneMark, Fgenesh)  

At minimum, most annotation projects will run with a protein data set, possibly complemented by some RNA-seq data. Popular examples of this are most of the traditional model systems, including human. However, a potential shortcoming of such approaches is that the comprehensiveness of the annotation depends directly on the input data. This can become a problem if our genome of interest is taxonomically distant to well-sequenced taxonomic groups so that only few protein matches can be found. Likewise, not all genes will be expressed at all times, making the generation of a comprehensive RNA-seq data set for annotation challenging.

We will therefore first run our annotation project in the traditional way, with proteins and ESTs, and then repeat the process with a well-trained ab-initio gene predictor. You can then compare the output to get an idea of how crucial the use of a gene predictor is. However, before we get our hands dirty, we need to understand Maker a little better...

The data for this you have already linked to your folder, but if not - run this command in a folder you have created in your home folder for this course:

*cd ~/*  
<i>cd \<some\_folder\></i>  
*ln -s /proj/g2014065/course\_data data*  

This data folder is write-proteced, it is only a resource for you to obtain data from, but not where you are writing your own outputs to!
## Organizing the data

The data we are providing for the course is organized in the following way (assuming you have sym-linked it to a folder name 'data':

data/human

data/dmel

- chromosome_4/

  - bam/

  - chromosome/

  - evidence/

  - maker

  - raw_computes

- full_genome

  - evidence/

  - genome/

  - raw/

blastdb/

scripts/

The folder scripts contains two perl scripts that we will use to format some data (referred to as $SCRIPT_PATH). The Blastdb folder will be used for the functional annotation exercise tomorrow.
## Loading Maker

Maker strings together a range of different tools into a complex pipeline (e.g. blast, exonerate, repeatmasker, augustus...). On Uppmax, loading all these tools and Maker into your global PATH is done simply by typing:

<i>module load bioinfo-tools</i>  
<i>module load maker/2.31</i>

If you are trying to run Maker on your own computer or cluster, make sure that in fact all its various dependencies are loaded.
## Understanding Makers control files

Makers behaviour and information on input data are specified in one of three control files. These are:

- maker_opts.ctl  
- maker_bopts.ctl  
- maker_exe.ctl

What are these files for?

'maker_exe.ctl' holds information on the location of the various binaries required by Maker (including Blast, Repeatmasker etc). Normally, all information in this file will be extracted from $PATH, so if everything is set up correctly, you will never have to look into this file.

Next, 'maker_bopts.ctl' provides access to a number of settings that control the behaviour of evidence aligners (blast, exonerate). The default settings will usually be fine, but if you want to try to annotate species with greater taxonomic distance to well-sequenced species, it may become necessary to decrease stringency of the e.g. blast alignments.

Finally, 'maker_opts.ctl' holds information on the location of input files and some of the parameters controlling the decision making during the gene building.
## Running Maker - Human contig

The first exercise will be a very short one in which you will create a gene build for a small piece of the human genome. This is to familiarize yourself with all the settings available in Maker. In the next exercise, the heavy-lifting of generating the protein alignments etc will already been done by us.

[Using the included test data](ExerciseHumanTestdata)
## Running Maker - Drosophila genome
### Getting the pre-computes

Running Maker on a full genome, even of an invertebrate, can take a considerable amount of time - especially if only few processing cores are available. We have therefore generate the genome-wide raw computes prior to the course. You can find the in the folder you have symlinked earlier under

_ln -s data/dmel/raw_

Ab-initio guided or evidence-based

Next, we will annotate the genome of the fruit fly _Drosophila melanogaster_. First without ab-initio predictions to guide the gene build, and afterwards with. Make sure you work with the sample data stored in the sub-folder '<b>chromosome_4</b>'. Even a small genome like Drosophila would take too long to run within the time we have for this course.

[Running Maker without ab-initio predictions](ExcerciseMakerNoAbinit)

[Running Maker with ab-initio predictions](ExcerciseMakerAbinit)
### Inspecting the output

The running of an annotation pipeline like Maker is not actually very hard. But the complicated work is only beginning. How to we best inspect the gene builds? Count features? Visualize it? Most importantly, what steps do we need to take to create a 'finished' annotation that we can use for scientific analyses?

[Comparing and evaluating annotations](ExcerciseMakerCompareAnnot)
## Closing remarks

This concludes the gene building part. We have learned how to use the Maker annotation pipeline and have created gene builds with and without ab-initio predictions. Moreover, we have employed some measures to describe and judge these annotations. An essential part that we decided to leave out is the training of ab-initio gene finders. The reason for this omission was that there isn't really any one best way to do this and your mileage may vary a lot based on your organism and input data. Perhaps the most direct approach available at the moment is a combination of evidence-based annotation with Maker and to use the resulting, crude gene models to train SNAP. Since Maker can improve ab-initio predictions 'on the fly', it can tolerate a bit of noise from a less-than-perfect ab-initio profile. If you are setting out on an annotation project, the BILS annotation platform would be happy to discuss the best approach for your data with you.

With that being said, generating a gene build is only one part of an annotation project. Next, we will inspect the annotation in genome browser and make an attempt at functional inference for the predicted gene models.