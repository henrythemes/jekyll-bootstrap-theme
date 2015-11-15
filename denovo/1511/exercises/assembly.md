---
layout: default
title:  'Exercise: De Novo Assembly'
---

##Exercise: De Novo Assembly

 In this exercise you will assemble genomes de novo using commonly used assembly software. You will work with Illumina data of Rhodobacter sphaerioides, data that was used in the GAGE-B comparison of assemblers. All settings used for the different programs are the ones used by the GAGE-B project. Due to time-restrictions we are forced to stick to assemblers that run quickly, but if you have an assembly project of your own you are encouraged to try other assemblers too. Remember that no assembler is best for all projects; if possible you need to try several.

**Questions**

If you have questions about the lab after the course, you are welcome to contact me: henrik.lantz@bils.se

### Setup

Before you start working with the assemblers you need to setup a folder structure in your home folder to have easy access to the data.

Login to UPPMAX using the earlier supplied instructions and be sure to log in to your reserved node.

To make certain you are in your home folder, type: `cd`

Make a directory for all of today’s exercises using `mkdir RhodoAssembly`

And enter it using `cd RhodoAssembly`

Now make a copy of the folder which includes links to all data using `cp -r /proj/g2015027/private/nobackup/assembly_workshop/data .`

You are now ready to proceed to working with the first assembly program, Velvet. You will start by using the HiSeq data for all exercises, and later if time allows, redo with the MiSeq data and compare results. 

### Part1, HiSeq data

#### Spades

Start by making a folder called spades in RhodoAssembly and enter it with `cd spades`.

Now load the spades module:

```
module load bioinfo-tools
module load spades
```

Spades is very easy to run in the basic configuration. Usually you want to run with the "--careful" flag, but this will take too long for this exercise. It will take around 16 mins anyway, grab a coffee.

```
spades -t 8 --pe1-1 ../data/Rhodo_Hiseq_trimmed_read1.fastq --pe1-2 ../data/Rhodo_Hiseq_trimmed_read2.fastq -o SpadesOut
```

In SpadesOut you will now have a number of files, including contigs.fasta and scaffolds.fasta.

Take a look at the files using `less`. Can you see any regions where contigs have been scaffolded together?

We then calculate some statistocs and generate plots using Quast:
```
module load module load quast/2.3
quast.py -o spades -l Spades_scaffolds,Spades_contigs -t 1 scaffolds.fasta contigs.fasta
```

Download the whole quast result-folder (spades) to your own computer using scp and click on the reports.html file. Any big differences between the scaffolds and contigs files?

Now go on the next assembly program:

#### Abyss

First go to your folder RhodoAssembly and make a new folder using `mkdir Abyss`.

Now load the necessary modules:

```
module load abyss/1.3.7
module load bowtie
```

If you have paired end data you can start Abyss using the abyss-pe script:

```
abyss-pe k=31 l=1 n=5 s=100 np=8 name=asm lib='reads' reads=' ../data/Rhodo_Hiseq_trimmed_read1.fastq ../data/Rhodo_Hiseq_trimmed_read2.fastq' aligner=bowtie
```

Once done you will have two files called asm-contigs.fa and asm.scaffolds.fa. Now load these files together with the earlier Spades contigs. Can you based on these numbers say which assembler does the best job? Note that this is a trick question!

The next assembler we'll try is

#### SOAP denovo

SoapDeNovo is one of very few program that can assemble large genomes using high coverage Illumina data.

First, load the module:

```
module load soapdenovo/2.04-r240
```

Now make and enter a folder called 'soap'

SoapDeNovo requires a config file that you need to create yourself using a text editor like nano:

```
nano soap.config
```

Enter this information:

```
[LIB]

avg_ins=220

reverse_seq=0

asm_flags=3

rank=1

q1=../data/Rhodo_Hiseq_trimmed_read1.fastq

q2=../data/Rhodo_Hiseq_trimmed_read2.fastq
```

Exit and save the file by ctrl-x (if using nano) and answer yes when asked to save.

Start SoapDeNovo by:

```
SOAPdenovo-63mer all –K 55 –F –R –E –w –u –s soap.config –o asm –p 8>> SOAPdenovo.log
```

Check the result-files asm.contig and asm.scafSeq for N50 size and number of contigs/scaffolds and compare with earlier results using Quast.

SoapDeNovo also comes with a GapCloser utlity that tries to improve the assemblies by closing gaps in the scaffolds. Try it out using:

```
GapCloser –b soap.config –a asm.scafSeq –o asm.new.scafSeq –t 8 >> SOAPdenovo.log
```

Any improvements? 

### Part 2, MiSeq data

There is also MiSeq data for the same organism. Links are in the data folder. You should now try running the three programs you already tried using MiSeq data to see if you get any improvements. If you feel adventurous, you can also try to run the programs with both HiSeq and MiSeq at the same time.

*IMPORTANT! Remember to work in different directories and/or change the name of output directories so that you do not overwrite your old data!*

The MiSeq data have longer reads, and you therefore need to change the following parameters (and remember to use the MiSeq files this time):

- Abyss-pe - change k to 49
- SoapDeNovo - change in the config file avg_ins to 540, use on the command line a K value of 79.

Compare with your HiSeq results. Differences? 

### Part 3, Not had enough?

This is an optional exercise. Here you should try to use Mira on the HiSeq data. Load the module using: `module load mira/4.0rc4`

Here you will receive no help. Try to figure out how to use the program by googling, in particular try to find the manual. Once you get it to run, kill it. No, I am serious. Press ctrl-c to kill the process. It simply takes to long to run, but now you have probably learned a lot by trying to get it to work. smile Then check the result-files in `/proj/g2015027/private/assembly_workshop/mira/`

Was the longer running-time worth it?

By **Henrik Lantz**, BILS/SciLife/Uppsala University 
