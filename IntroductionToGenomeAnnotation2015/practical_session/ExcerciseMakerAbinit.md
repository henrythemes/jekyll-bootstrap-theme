---
layout: default
title:  'Exercise - Running Maker with ab-initio predictions'
---

# Running Maker with ab-initio predictions

The recommended way of running Maker is in combination with one or more ab-initio profile models. Maker natively supports input from several tools, including augustus, snap and genemark. The choice of tool depends a bit on the organism that you are annotating - for example, GeneMark -ES is mostly recommended for fungi, whereas augustus and snap have a more general use.

The biggest problem with ab-initio models is the process of training them. It is usually recommended to have somewhere around 500-1000 curated gene models for this purpose. Naturally, this is a bit of a contradiction for a not-yet annotated genome.

However, if one or more good ab-initio profiles are available, they can potentially greatly enhance the quality of an annotation by filling in the blanks left by missing evidence. Interestingly, Maker even works with ab-initio profiles from somewhat distantly related species since it can create so-called hints from the evidence alignments, which the gene predictor can take into account to fine-tune the predictions.
## Create a new Maker project

In order to compare the performance of Maker with and without ab-initio predictions in a real-world scenario, we have first run a gene build without ab-initio predictions. Now, we run a similar analysis but enable ab-initio predictions through augustus.

Create a new folder for this Maker run:

<i>mkdir maker\_with\_abinitio</i>  
<i>cd maker\_with\_abinitio</i>

Now link the raw computes you want to use into your folder:

 - repeatmasker.chr4.gff (coordinates of repeatmasked regions)  
 - cufflinks2maker.chr4.gff (EST hints created from assembled transcripts)

In addition, you will also need the genome sequence and a protein set. Sym-link it from the the data directory created earlier:

*ln -s /path/to/chromosome\_4/chromosome/4.fa*

This time, we do specify a reference species to be used by augustus, which will enable ab-initio gene finding:

*augustus\_species=fly* #Augustus gene prediction species model  
...  
<i>protein2genome=0</i>  
<i>est2genome=0</i>

With these settings, Maker will run augustus to predict gene loci, but inform these predictions with information from the protein and est alignments. Note that if you are tempted to run both ab-initio predictions and gene model calling from evidence, you are out of luck. While Maker won't complain when trying to do this, it will simply not create hints for the ab-initio predictions. And in most cases, this will result in problems.
## Run Maker with ab-initio predictions

With everything configured, run Maker as you did for the previous analysis:

*maker -c 8*

We probably expect this to take a little bit longer than before, since we have added another step to our analysis.

## Compile the output

When Maker has finished, compile the output:

<i>$SCRIPT\_PATH/maker\_merge\_outputs.pl --outdir maker\_a1\_p0\_c0</i>  

And again, it is probably best to link the resulting output (maker.gff) to a result folder, under a descriptive name.