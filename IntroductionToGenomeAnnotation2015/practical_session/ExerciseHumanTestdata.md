---
layout: default
title:  'Exercise Using Human test data'
---

# Using the included test data

Maker comes with some test data to try out the pipeline. This test set is sufficiently small that we can run all the alignment steps without having to worry about the runtime too much.

The human test data set is located in the folder you symlinked earlier, named 'data/human'.
## Create project folder

First, we create a new folder in which to store all the configuration and input files. To do so, type:

<i>mkdir maker\_human</i>  
<i>cd maker\_human</i>

If you haven't done so already, load the maker module:

<i>module load bioinfo-tools</i>  
<i>module load maker/2.31</i>

Next, we create symbolic links to sequence files we wish to use in this exercise (located in data/human):

hsap\_contig.fasta - a piece of the human genome  
hsap\_protein.fasta - Proteins that map to the genomic region  
hsap\_est.fasta - EST data that maps to the genomic region

##Configure your maker project

Next, we create the 3 maker control files:

_maker -CTL_

Off these, only maker_opts.ctl is of concern to us. Have a look at the following sections and fill in the information as shown:

#-----Genome (these are always required)  
genome=hsap_contig.fasta #genome sequence (fasta file or fasta embeded in GFF3 file)  
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

...

#-----EST Evidence (for best results provide a file for at least one)  
est=hsap_est.fasta #set of ESTs or assembled mRNA-seq in fasta format  
altest= #EST/cDNA sequence file in fasta format from an alternate organism  
est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file  
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

...

#-----Protein Homology Evidence (for best results provide a file for at least one)  
protein=hsap_protein.fasta #protein sequence file in fasta format (i.e. from mutiple oransisms)  
protein_gff= #aligned protein homology evidence from an external GFF3 file

...

#-----Repeat Masking (leave values blank to skip repeat masking)  
model_org=human #select a model organism for RepBase masking in RepeatMasker  
rmlib= #provide an organism specific repeat library in fasta format for RepeatMasker   
repeat_protein=/pica/sw/apps/bioinfo/maker/2.31/milou/data/te_proteins.fasta #provide a fasta file of transposable element proteins for RepeatRunner  
rm_gff= #pre-identified repeat elements from an external GFF3 file  
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no  
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

...

#-----Gene Prediction  
snaphmm= #SNAP HMM file  
gmhmm= #GeneMark HMM file  
augustus_species=human #Augustus gene prediction species model  
fgenesh_par_file= #FGENESH parameter file  
pred_gff= #ab-initio predictions from an external GFF3 file  
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)  
est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no  
protein2genome=0 #infer predictions from protein homology, 1 = yes, 0 = no  
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no  
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs  
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no
---++ Run Maker

Above, ew have specified to annotate a short piece of the human genome using the augustus gene finder to create ab-initio predictions, combined with evidence alignments from both protein and EST data. To run this configured analysis, just type:

_maker -c 8_

This will start Maker on 8 cores, if everything is configured correctly.

*In a real-world scenario, executing maker in this way will not be the preferred option. Instead of using Makers' powerful parallelization via MPI, we are simply telling it to run on one node with 8 cores. This works fine for a small test data set, but will be insufficient for a large genome. Unfortunately, the MPI mode of Maker doesn't currently work on Milou.

## Inspect the output

### Finding your way around

By default, Maker will write the output of its different analyses into a folder named:

&lt;name_of_genome_fasta&gt;.maker.output

In our case:

hsap_contig.maker.output

Within the main output directory, Maker keeps a copy of the config files, a database (here: hsap_contig.db), directories for the blast databases created from your evidence data and a file called hsa_contig_master_datastore_index.log.

Out of these files, only the master_datastore_index is really interesting to us. It includes a log of all the contigs included in the genome fasta file - together with their processing status (ideally: FINISHED) and the location of the output files. Since Maker can technically run in parallel on a large number of contigs, it creates separate folders for each of these input data. For larger genomes, this can generate a very deep and confusing folder tree. The master_datastore_index helps you make sense of it:

NT_010783.15 hsap_contig_datastore/80/99/NT_010783.15/ STARTED  
NT_010783.15 hsap_contig_datastore/80/99/NT_010783.15/ FINISHED

This meens the contig NT_010783.15 was started - and finished, with all data (annotation, protein predictions etc) written to the subfolder hsap_contig_datastore/80/99/NT_010783.15/.

If you look into that folder, you will find the finished Maker annotation for this contig.

rw-rw-r- 1 marc b2011210 472193 Mar 24 10:16 'NT_010783%2E15.gff  
rw-rw-r- 1 marc b2011210 3599 Mar 24 10:16 NT_010783%2E15.maker.augustus_masked.proteins.fasta  
rw-rw-r- 1 marc b2011210 10388 Mar 24 10:16 NT_010783%2E15.maker.augustus_masked.transcripts.fasta  
rw-rw-r- 1 marc b2011210 176 Mar 24 10:16 NT_010783%2E15.maker.non_overlapping_ab_initio.proteins.fasta  
rw-rw-r- 1 marc b2011210 328 Mar 24 10:16 NT_010783%2E15.maker.non_overlapping_ab_initio.transcripts.fasta  
rw-rw-r- 1 marc b2011210 3931 Mar 24 10:16 NT_010783%2E15.maker.proteins.fasta  
rw-rw-r- 1 marc b2011210 20865 Mar 24 10:16 NT_010783%2E15.maker.transcripts.fasta  
rw-rw-r- 1 marc b2011210 4248 Mar 24 10:15 run.log  
drwxrwsr-x 3 marc b2011210 4096 Mar 24 10:16 theVoid.NT_010783.15

The main annotation file is 'NT_010783%2E15.gff' - including both the finished gene models and all the raw compute data. The other files include fasta files for the different sequence features that have been annotated - based on ab-initio predictions through augustus as well as on the finished gene models. The folder 'theVoid' inclused all the raw computations that Maker has preformed to synthesize the evidence into gene models.

## Understanding a Maker annotation

You have two options now for gathering the output in some usable form - copy select files by hand to wherever you want them. Or you can use a script that does the job for you (we have included an example in the script folder).

From your Maker folder, run the script called 'maker_merge_outputs.pl' to create an output file for all annotations and protein files:

$SCRIPT_PATH/maker_merge_outputs.pl

This will create two files:

annotations.gff  
annotations.proteins.fa

If you use 'less' to read the annotation file ([GFF3 format](http://www.sequenceontology.org/gff3.shtml)), you will see a range of different features:

##gff-version 3  
NT_010783.15 . contig 1 201444 . . . ID=NT_010783.15;Name=NT_010783.15  
NT_010783.15 maker gene 141432 160962 . + . ID=maker-NT_010783.15-augustus-gene-0.11;Name=maker-NT_010783.15-augustus-gene-0.11  
NT_010783.15 maker mRNA 141432 160962 . + . ID=maker-NT_010783.15-augustus-gene-0.11-mRNA-1;Parent=maker-NT_010783.15-augustus-gene-0.11;Name=maker-NT_010783.15-augustus-gene-0.11-mRNA-1;_AED=0.00;_eAED=0.00;_QI=334|1|1|1|0|0|5|700|112

...

For example, the above lines read:

A new contig is being shown, with the id 'NT_010783.15' and a length of 201444 nucleotides  
On this contig, a gene feature is located from position 141432 to 160962, on the forward strand and with the id 'maker-NT_010783.15-augustus-gene-0.11'. In this example, the id was given by maker - initially predicted by augustus and then refined by maker.  
On this contig, belonging to the gene, is located a transcript from position 141432 to 160962, on the forward strand and with the id 'maker-NT_010783.15-augustus-gene-0.11-mRNA-1'. It's quality, or AED score, is 0.00 - which means that the evidence alignments are in perfect agreement with the transcript model.

And so on.

Next, we want to separate the annotation file into it's different types of information. This is useful for a number of applications, like visualizing it as separate tracks in a genome browser. Or to compute some intersting numbers from the gene models.

First, lets split the annotation file:

<i>mkdir -p annotations</i>  
<i style="background-color: transparent;">$SCRIPT_PATH/split_gff_by_source.pl --input annotations.gff -d annotations</i>

This should create a bunch of files, including 'maker.gff' - which contains the actual gene models.

Next, we load the GenomeTools package:

_module load genometools/1.5.3_

And run its counting method on the gene models:

_gt stat maker.gff_

We could now also visualise all this information using a genome browser, such as [IGV](http://www.broadinstitute.org/igv/home). IGV requires a genome fasta file and any number of annotation files in GTF or GFF3 format (note that GFF3 formatted file tend to look a bit weird in IGV sometimes).

-- %USERSIG{IntroGenAnno1403Teacher - 2014-03-27}%

---++