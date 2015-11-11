---
layout: default
title:  'Part 3: Single cell genome assembly using SPAdes'
---

# Part 3: Single cell genome assembly using SPAdes

In this part of the course, you will start doing assemblies of 'real' (but reduced) single cell genome datasets using the assembler SPAdes (which was introduced by Kasia). 
The idea is that you will be exploring how different settings of SPAdes and different pre-treatments ('trimming') of your datasets will affect your assembly quality ('assembly metrics'). 
Below you will find a table in which different settings are indicated for running SPAdes (different 'flags' and datasets trimming).
Given that assembly is relatively time-consuming (even with the reduced datasets used here during the tutorial), we suggest that you distribute the different assemblies amongst different groups to save time. 
There are a total of 12 SPAdes assemblies to run in this exercise and it would be good to form a group of 4 people and each can take care of 3 assemblies. 
If you can complete these 12 assemblies and have time remaining, you can do Part 7 as a bonus exercise.  

Note that Part 6 is for running assemblies using MiSeq data.  
Note that you will have to fill in the results from the exercises in Tables 1 to 4:  
<p>
<style type="text/css">
.tg  {border-collapse:collapse;border-spacing:0;border-color:#999;}
.tg td{font-family:Arial, sans-serif;font-size:14px;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#999;color:#444;background-color:#F7FDFA;}
.tg th{font-family:Arial, sans-serif;font-size:14px;font-weight:normal;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#999;color:#fff;background-color:#26ADE4;}
.tg .tg-7un6{background-color:#ffffff;color:#000000;text-align:center;vertical-align:top}
.tg .tg-3xqa{background-color:#26ade4;font-weight:bold;color:#ffffff;text-align:center;vertical-align:top}
.tg .tg-baqh{text-align:center;vertical-align:top}
.tg .tg-amwm{font-weight:bold;text-align:center;vertical-align:top}
.tg .tg-j0tj{background-color:#D2E4FC;text-align:center;vertical-align:top}
</style>
<table class="tg">
  <tr>
    <th class="tg-amwm">Table 1</th>
    <th class="tg-7un6" colspan="3">HiSeq data - Original data</th>
  </tr>
  <tr>
    <td class="tg-3xqa">(Dataset 1)</td>
    <td class="tg-3xqa">With --sc and --careful flags</td>
    <td class="tg-3xqa">without --sc flag</td>
    <td class="tg-3xqa">without --careful flag</td>
  </tr>
  <tr>
    <td class="tg-baqh">Number of reads</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
  </tr>
  <tr>
    <td class="tg-j0tj">Assembly time</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
  </tr>
  <tr>
    <td class="tg-baqh">Number of contigs</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
  </tr>
  <tr>
    <td class="tg-j0tj">Total assembly size</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
  </tr>
  <tr>
    <td class="tg-baqh">Largest contig</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
  </tr>
  <tr>
    <td class="tg-j0tj">N50</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
  </tr>
  <tr>
    <td class="tg-baqh">G+C%</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
  </tr>
  <tr>
    <td class="tg-j0tj">Number of ORFs</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
  </tr>
  <tr>
    <td class="tg-baqh">Completeness (%)</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
  </tr>
</table>
</p>


<p>
<style type="text/css">
.tg  {border-collapse:collapse;border-spacing:0;border-color:#999;}
.tg td{font-family:Arial, sans-serif;font-size:14px;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#999;color:#444;background-color:#F7FDFA;}
.tg th{font-family:Arial, sans-serif;font-size:14px;font-weight:normal;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#999;color:#fff;background-color:#26ADE4;}
.tg .tg-7un6{background-color:#ffffff;color:#000000;text-align:center;vertical-align:top}
.tg .tg-3xqa{background-color:#26ade4;font-weight:bold;color:#ffffff;text-align:center;vertical-align:top}
.tg .tg-baqh{text-align:center;vertical-align:top}
.tg .tg-amwm{font-weight:bold;text-align:center;vertical-align:top}
.tg .tg-j0tj{background-color:#D2E4FC;text-align:center;vertical-align:top}
</style>
<table class="tg">
  <tr>
    <th class="tg-amwm">Table 2</th>
    <th class="tg-7un6" colspan="3">MiSeq data - Original data</th>
  </tr>
  <tr>
    <td class="tg-3xqa">(Dataset 2)</td>
    <td class="tg-3xqa">With --sc and --careful flags</td>
    <td class="tg-3xqa">without --sc flag</td>
    <td class="tg-3xqa">without --careful flag</td>
  </tr>
  <tr>
    <td class="tg-baqh">Number of reads</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
  </tr>
  <tr>
    <td class="tg-j0tj">Assembly time</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
  </tr>
  <tr>
    <td class="tg-baqh">Number of contigs</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
  </tr>
  <tr>
    <td class="tg-j0tj">Total assembly size</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
  </tr>
  <tr>
    <td class="tg-baqh">Largest contig</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
  </tr>
  <tr>
    <td class="tg-j0tj">N50</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
  </tr>
  <tr>
    <td class="tg-baqh">G+C%</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
  </tr>
  <tr>
    <td class="tg-j0tj">Number of ORFs</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
  </tr>
  <tr>
    <td class="tg-baqh">Completeness (%)</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
  </tr>
</table>
</p>


<p>
<style type="text/css">
.tg  {border-collapse:collapse;border-spacing:0;border-color:#999;}
.tg td{font-family:Arial, sans-serif;font-size:14px;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#999;color:#444;background-color:#F7FDFA;}
.tg th{font-family:Arial, sans-serif;font-size:14px;font-weight:normal;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#999;color:#fff;background-color:#26ADE4;}
.tg .tg-7un6{background-color:#ffffff;color:#000000;text-align:center;vertical-align:top}
.tg .tg-3xqa{background-color:#26ade4;font-weight:bold;color:#ffffff;text-align:center;vertical-align:top}
.tg .tg-baqh{text-align:center;vertical-align:top}
.tg .tg-amwm{font-weight:bold;text-align:center;vertical-align:top}
.tg .tg-j0tj{background-color:#D2E4FC;text-align:center;vertical-align:top}
</style>
<table class="tg">
  <tr>
    <th class="tg-amwm">Table 3</th>
    <th class="tg-7un6" colspan="3">HiSeq data - Processed with SeqPrep</th>
  </tr>
  <tr>
    <td class="tg-3xqa">(Dataset 1)</td>
    <td class="tg-3xqa">With --sc and --careful flags</td>
    <td class="tg-3xqa">without --sc flag</td>
    <td class="tg-3xqa">without --careful flag</td>
  </tr>
  <tr>
    <td class="tg-baqh">Number of reads</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
  </tr>
  <tr>
    <td class="tg-j0tj">Assembly time</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
  </tr>
  <tr>
    <td class="tg-baqh">Number of contigs</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
  </tr>
  <tr>
    <td class="tg-j0tj">Total assembly size</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
  </tr>
  <tr>
    <td class="tg-baqh">Largest contig</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
  </tr>
  <tr>
    <td class="tg-j0tj">N50</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
  </tr>
  <tr>
    <td class="tg-baqh">G+C%</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
  </tr>
  <tr>
    <td class="tg-j0tj">Number of ORFs</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
  </tr>
  <tr>
    <td class="tg-baqh">Completeness (%)</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
  </tr>
</table>
</p>


<p>
<style type="text/css">
.tg  {border-collapse:collapse;border-spacing:0;border-color:#999;}
.tg td{font-family:Arial, sans-serif;font-size:14px;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#999;color:#444;background-color:#F7FDFA;}
.tg th{font-family:Arial, sans-serif;font-size:14px;font-weight:normal;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#999;color:#fff;background-color:#26ADE4;}
.tg .tg-7un6{background-color:#ffffff;color:#000000;text-align:center;vertical-align:top}
.tg .tg-3xqa{background-color:#26ade4;font-weight:bold;color:#ffffff;text-align:center;vertical-align:top}
.tg .tg-baqh{text-align:center;vertical-align:top}
.tg .tg-amwm{font-weight:bold;text-align:center;vertical-align:top}
.tg .tg-j0tj{background-color:#D2E4FC;text-align:center;vertical-align:top}
</style>
<table class="tg">
  <tr>
    <th class="tg-amwm">Table 4</th>
    <th class="tg-7un6" colspan="3">MiSeq data - Processed with SeqPrep</th>
  </tr>
  <tr>
    <td class="tg-3xqa">(Dataset 2)</td>
    <td class="tg-3xqa">With --sc and --careful flags</td>
    <td class="tg-3xqa">without --sc flag</td>
    <td class="tg-3xqa">without --careful flag</td>
  </tr>
  <tr>
    <td class="tg-baqh">Number of reads</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
  </tr>
  <tr>
    <td class="tg-j0tj">Assembly time</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
  </tr>
  <tr>
    <td class="tg-baqh">Number of contigs</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
  </tr>
  <tr>
    <td class="tg-j0tj">Total assembly size</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
  </tr>
  <tr>
    <td class="tg-baqh">Largest contig</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
  </tr>
  <tr>
    <td class="tg-j0tj">N50</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
  </tr>
  <tr>
    <td class="tg-baqh">G+C%</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
  </tr>
  <tr>
    <td class="tg-j0tj">Number of ORFs</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
    <td class="tg-j0tj">-</td>
  </tr>
  <tr>
    <td class="tg-baqh">Completeness (%)</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
    <td class="tg-baqh">-</td>
  </tr>
</table>
</p>


Actual tables to be filled in are provided in Google Docs and the links can be found below. 
Note that each group will consist of 4 people except Group 8 which will consist of 3 people. 
You should talk to each other to form the groups and make sure that you work in groups to discuss who will work on which assembly.

[Group 1:](https://docs.google.com/spreadsheet/ccc?key=0AuNHyFPCsxthdGRKMXJwdF9jVDMzX2lGMkdJSDdOcnc&usp=sharing)  
[Group 2:](https://docs.google.com/spreadsheet/ccc?key=0AuNHyFPCsxthdC0tdzFySDFyaDNIdEh4M01xMXFQb3c&usp=sharing)  
[Group 3:](https://docs.google.com/spreadsheet/ccc?key=0AuNHyFPCsxthdHhwdUdUWUxBbnd0eC15WkJhS29iV3c&usp=sharing)  
[Group 4:](https://docs.google.com/spreadsheet/ccc?key=0AuNHyFPCsxthdFZRcXBjN0lrMGV5NmNuRnUzV2RkT0E&usp=sharing)  
[Group 5:](https://docs.google.com/spreadsheet/ccc?key=0AuNHyFPCsxthdEhkV3hIejJaMDZrWDFqd29XNTZFbmc&usp=sharing)  
[Group 6:](https://docs.google.com/spreadsheet/ccc?key=0AuNHyFPCsxthdFZDQzhLamJocHI0M0ZBQ0dMUDRrSFE&usp=sharing)  
[Group 7:](https://docs.google.com/spreadsheet/ccc?key=0AuNHyFPCsxthdDU2OUk4ank4c1A1VVhhbjZPaldtN2c&usp=sharing)  
[Group 8:](https://docs.google.com/spreadsheets/d/1Q3QBvPYzQ1kFHjWu0O5Jf1wgSH-9sxMrcndoLlNW92Y/edit?usp=sharing)  


## 3.1a. Running SPAdes on untrimmed sequences
---

The following set of commands are to be typed in your compute node (for example mXX - look up using *jobinfo -u username* command). 
Make sure you are typing them in the compute node and not log in node. Go back to Part 1 to check how to log in to your compute node.
Before starting the exercises, you should make a folder in your home directory where the exercises will be run.

```sh
mkdir ~/metagenomics_exercises
cd ~/metagenomics_exercises/
```

Make a folder for SPAdes assemblies: 

```sh
mkdir spades_assemblies
cd spades_assemblies
```

Next, make symbolic links of sequences in that folder:

```sh
ln -s /proj/g2014180/nobackup/single_cell_exercises/sequences/dataset1/G5_Hiseq_[12].fastq G5_Hiseq_1.fastq
#ln -s /proj/g2014180/nobackup/single_cell_exercises/sequences/dataset1/G5_Hiseq_2.fastq G5_Hiseq_2.fastq
ln -s /proj/g2014180/nobackup/single_cell_exercises/sequences/dataset2/G5_Miseq_[12].fastq G5_Miseq_1.fastq
#ln -s /proj/g2014180/nobackup/single_cell_exercises/sequences/dataset2/G5_Miseq_2.fastq G5_Miseq_2.fastq
```

Now you are almost ready to run assemblies! But before you can start assemblies, you need to load SPAdes module first.  
To load the SPAdes assembler, type:

```sh
module load bioinfo-tools
module load spades/3.1.1
```

Now, you can run SPAdes. To run SPAdes with --sc and --careful flags (our default assembly setting), type:

```sh
time spades.py --sc --careful -t 8 -m 24 -1 G5_Hiseq_1.fastq -2 G5_Hiseq_2.fastq -o G5_Hiseq_sc_careful
```

This command will launch SPAdes assembly but also checks how long the assembly takes. After the assembly has completed, check the time it took for the program to run. You should look at the 'real' time. 
Record the time in the table for HiSeq data. Enter the time it takes for the assembly to finish in the Table 1 (printed for you). 
In general, typing the command 'time' before other commands will help you check how long the computation took.

## 3.1b. Using the '--sc' flag
---

SPAdes can handle single-cell genomic data that is known to be highly biased in terms of sequence coverage along the length of the genome. 
In order for SPAdes to be able to handle biased sequence coverage, you need to supply the *'--sc'* flag when running the assembly. 
In this exercise, you will omit the *'--sc'* flag to see how this affects the assembly. You should be able to compare results between assemblies with or without *'--sc'* flag. 

Type the following command:

```sh
time spades.py --careful -t 8 -m 24 -1 G5_Hiseq_1.fastq -2 G5_Hiseq_2.fastq -o G5_Hiseq_careful
```

## 3.1c. Using the '--careful' flag
---

Here, you will repeat the assembly in very much the same way as the previous exercise (3.1b) but this time you will omit the *'--careful'* flag instead of the *'--sc'* flag. 
The *'--careful'* flag uses *'bowtie'* tool to map the reads back to contigs and check for errors due to bad quality sequences and correct these errors. 
This results in longer assembly times than not using the *'--careful'* flag. This time, name the output directory as *'G5_Hiseq_sc'*. 
Follow previous examples to get the assembly metrics and tabulate the results. Notice how long it takes for the assembly to finish and record the time in Table 1.


## 3.1d. Merging reads with SeqPrep
---
To merge read pairs that have significant overlaps, we will use the tool called *'SeqPrep'*. First, load the module.

```sh
module load SeqPrep/2013-11-14
SeqPrep -f G5_Hiseq_1.fastq -r G5_Hiseq_2.fastq -1 G5_Hiseq_merge_1.fastq.gz -2 G5_Hiseq_merge_2.fastq.gz -s G5_Hiseq_merged.fastq.gz -q 30
```

Note that SeqPrep merges read pairs if overlaps between the read pairs are identified. Quality threshold of 30, for example, can be specified by '-q 30' flag 
for overlaps with some mismatches to be counted. It can also remove adapter sequences optionally.  
To run SPAdes assembly using merged reads, type the following command:

```sh
time spades.py --sc --careful -1 G5_Hiseq_merge_1.fastq.gz -2 G5_Hiseq_merge_2.fastq.gz -s G5_Hiseq_merged.fastq.gz -t 8 -m 24 -o G5_Hiseq_SeqPrep_sc_careful
```

You are now providing 3 input files to the assembler; 2 unmerged read pairs and 1 merged reads. After this assembly is done, you should repeat the assemblies again but omitting *'--sc'* or *'--careful'* flags as in the previous exercises. 
Name the next two assemblies as G5_Hiseq_SeqPrep_sc and G5_Hiseq_SeqPrep_careful.

## 3.2. Assessing assembly quality using Quast
---

After the assembly runs are done, you will use this tool called 'Quast' to calculate the basic metrics such as the length of largest contig, N50, etc.  
To run 'Quast', first go to the assembly directory and type the following commands:  

```sh
module load quast/2.3
quast.py -o assembly_metrics contigs.fasta
```

Results from 'Quast' will be found in the 'assembly_metrics' folder. Take a look at the files produced by the program. 
The summary file containing the stats is 'report.txt'. You should see the assembly metrics such as N50, G+C%, largest contig, total length (i.e., total length of all contigs added). 
Repeat Quast for the rest of the assemblies. Tabulate the results.

## 3.3. Gene prediction using Prodigal
---

Prodigal is a tool that can identify open reading frames (ORFs) in microbial genomes (bacteria or archaea). 
In this exercise, you will learn how to use Prodigal to predict ORFs and to prepare them for running Blastp later.  
Type the following commands to run Prodigal:  

```sh
module load prodigal/2.60
```

Go into a specific assembly folder (for example G5_Hiseq_careful_sc) and run:  

```sh
prodigal -i contigs.fasta -a contigs.prodigal.faa -o contigs.prodigal.gff
```

Take a look at the files produced by Prodigal. For example, type less *'contigs.prodigal.faa'* and see what the contents look like. 
Can you count how many genes are predicted by Prodigal? Use the example given in the Part 1 to count the number of genes (look for a repeating pattern in the file and search for that pattern). 
Tabulate the results.

## 3.4. Identifying your cell
---

First, you will run this tool called 'rnammer' to predict the regions within assembled contigs that contain ribosomal RNA sequences (rRNA), such as 16S, 23S, and 5S rRNA sequences.  
Go into the folder containing the assembled contigs and type this command:  

```sh
module load rnammer/1.2
rnammer -S bac -m lsu,ssu,tsu -gff contigs.rnammer.gff -f contigs.rnammer.fasta < contigs.fasta
```

Take a look at the file (*'contigs.rnammer.gff'*) produced by 'rnammer'.  
Can you identify the positive matches predicted by 'rnammer'?  
Can you interpret the result output?  
Next, you will run 'blastn' against the 'Silva' database. 
Depending on whether or not you have identified 16S or 23S sequences, you will need to run 'blastn' on a different database.  
First, load the Blast module:

```sh
module load blast/2.2.29+
```

If you have identified a 16S rRNA sequence (check the gff file), type:

```sh
blastn -query contigs.rnammer.fasta -db /proj/g2014180/nobackup/single_cell_exercises/databases/SSURef_NR99_115_tax_silva_trunc.dna.fasta -evalue 1e-6 -num_threads 8 -out contigs.rnammer.16S_silva.blastn
```

If you have identified 23S rRNA sequence, then type:

```sh
blastn -query contigs.rnammer.fasta -db /proj/g2014180/nobackup/single_cell_exercises/databases/LSURef_115_tax_silva_trunc.dna.fasta -evalue 1e-6 -num_threads 8 -out contigs.rnammer.23S_silva.blastn
```

After running Blastn, can you identify what organism G5 belongs to?

## 3.5. Running completeness estimates
---

Of course, you will be interested how 'complete' your assembly is. In order to determine completeness, we make use of an in-house script that checks for the presence of a number of universally conserved genes. 
Based on how many of these genes will be identified, an estimation of genome completeness will be made.

First, load the correct version of HMMER needed for this analysis:
You will need to unload the older version of HMMER that was loaded with rnammer program as the script you will use in this part requires the latest version of HMMER.

```sh
module unload hmmer/2.3.2-gcc
module load hmmer/3.1b1-gcc
```

Then go into the assembly directory and run the following perl script:

```sh
perl /proj/g2014180/nobackup/single_cell_exercises/scripts/micomplete.pl -h /proj/g2014180/nobackup/single_cell_exercises/scripts/Bact139.hmm -w /proj/g2014180/nobackup/single_cell_exercises/scripts/Bact139.weights -p contigs.prodigal.faa -t 8 -e -c 1e-10
```

If the script completed without any errors, you should see that the script printed something like this:

Completeness: ??  
N markers found: ?? out of 139

The script looks for unique marker genes that are present in single copies in most bacteria and estimates how complete the genome is based on the marker genes found. Completeness is shown as fractions and you should multiply this with 100 to get the percentage completeness.


## Questions:
---

**Q3.1:** Did you notice how many read pairs from HiSeq or MiSeq data were merged by SeqPrep?  
Is there a reason why a certain data set has higher merge rates than the other? How many reads were discarded in the process?  
**Q3.2:** Did you notice any differences between HiSeq and MiSeq data assembled using both --sc and --careful flags, i.e., default?  
**Q3.3:** Did you notice any differences in the quality of assembly when --sc and --careful were omitted in each data set?  
**Q3.4:** What do you think is the best way to assess the 'quality' of an assembly? (e.g. total size, N50, number of predicted ORFs, completeness)  
**Q3.5:** What do you think is the best way to assemble this particular SC dataset? Why?  
**Q3.6:** What is the identity of the organism based on the analyses you have performed? What phylum does it belong to and is there any closely related organisms in the databases?  
**Q3.7:** Try to find out in what type of environment you might find similar organisms in .

Note: If you are ahead of time and have completed the exercises (and have plenty of time) you can start doing Part 6, which is similar to this exercise but using MiSeq data.

<!--- 
Notes for next year's workshop
- Uppmax support at the beginning of the workshop (troubleshooting, node reservation, technical support - cables etc) => not enough ethernet socket
- Having MEGAN and Artemis installed locally instead of on Milou (15-30min for sorting out computer problems/ software installation, unless we have dedicated computers/terminals). 
Students are to pre-install the tools on their computers if they will use their own laptops.
- table of contents
- useful Linux commands such as:
pwd to check current working directory
tab completion
arrow up for previous command
- define learning objectives
- include pre-course materials about fasta/q formats, basic unix commands (count number of reads in fasta or fastq files as checks).
- Longer introductory lecture on basics such as NGS data
- flowchart of steps involved in generating single cell genome data to help visualize how we got the data (overview of what part they are involved in)
- check points (completeness estimates - distribution of marker genes, chimera generation, MEGAN summary, Artemis)
- extra day for binning (together with metagenomic workshop)
- guest lectures (national/international) prior to hands-on workshops
-->
