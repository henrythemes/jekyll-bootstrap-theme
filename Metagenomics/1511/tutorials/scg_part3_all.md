---
layout: default
title:  'Part 3: Single cell genome assembly using SPAdes'
---

# Part 3: Single cell genome assembly using SPAdes

<p>
In this part of the course, you will start doing assemblies of 'real' (but reduced) single cell genome datasets using the assembler SPAdes (which was introduced by Kasia). 
The idea is that you will be exploring how different settings of SPAdes and different pre-treatments ('trimming') of your datasets will affect your assembly quality ('assembly metrics'). 
Below you will find a table in which different settings are indicated for running SPAdes (different 'flags' and datasets).
Given that assembly is relatively time-consuming (even with the reduced datasets used here during the tutorial), we suggest that you distribute the different assemblies amongst different groups to save time. 
There are a total of 12 SPAdes assemblies to run in this exercise and it would be good to form a group of 4 people and each can take care of 3 assemblies.
</p>

<!--- 
If you can complete these 12 assemblies and have time remaining, you can do Part 7 as a bonus exercise. 
Note that Part 6 is for running assemblies using MiSeq data. 
Note that you will have to fill in the results from the exercises in Tables 1 to 4 
--->

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


## 3.1a. Preparing your data

The following set of commands are to be typed in your compute node (for example mXX - look up using ```jobinfo -u username``` command). 
Make sure you are typing them in the compute node and not log in node. Go back to Part 1 to check how to log in to your compute node.  
Before starting the exercises, you should make a folder named *single_cell_exercises* in your home directory where the exercises will be run.
Then create 2 folders *dataset1* and *dataset2* where the raw data will be linked

```sh
mkdir ~/single_cell_exercises
cd ~/single_cell_exercises/
mkdir dataset1 dataset2
```

Next, make symbolic links of sequences in those folders:

```sh
ln -s /proj/g2015028/nobackup/single_cell_exercises/sequences/dataset1/ dataset1/
ln -s /proj/g2015028/nobackup/single_cell_exercises/sequences/dataset2/ dataset2/
```
**Please, do not modify those files**

Check that the data is present in those 2 datasets folders

```sh
ls dataset1
ls dataset2
```

You should now see 2 files per dataset, a forward fastq file ```_R1_001.fastq``` and its reverse ```_R2_001.fastq```  

Later in some commands we use the variables *sample* and *merge*, the following commands set those variables. 
We will also load the softwares we need to work.  
**In case you loose your connection, you will need to redo this step again.**  

#### For *Hiseq* data without merging:
```sh
sample=Hiseq
merge=''
cd ~/single_cell_exercises/dataset1
source /proj/g2015028/nobackup/single_cell_exercises/modules_load
```

#### For *Hiseq* data with merging:
```sh
sample=Hiseq
merge=SeqPrep_
cd ~/single_cell_exercises/dataset1
source /proj/g2015028/nobackup/single_cell_exercises/modules_load
```

#### For *Miseq* data without merging:
```sh
sample=Miseq
merge=''
cd ~/single_cell_exercises/dataset2
source /proj/g2015028/nobackup/single_cell_exercises/modules_load
```

#### For *Miseq* data with merging:
```sh
sample=Miseq
merge=SeqPrep_
cd ~/single_cell_exercises/dataset2
source /proj/g2015028/nobackup/single_cell_exercises/modules_load
```

## 3.1b. Merging reads with SeqPrep

**This part is only for the ones that will work on Merged reads. Skip this if you do the assembly on the raw reads**  
To merge read pairs that have significant overlaps, we will use the tool called *'SeqPrep'*. 

First, create a folder where to store the output then merge the reads:

```sh
mkdir merged_reads
SeqPrep -q 30 -f G5_${sample}_1.fastq -r G5_${sample}_2.fastq \
-1 merged_reads/G5_${sample}_merge_1.fastq.gz \
-2 merged_reads/G5_${sample}_merge_2.fastq.gz \
-s merged_reads/G5_${sample}_merged.fastq.gz
```

Note that SeqPrep merges read pairs if overlaps between the read pairs are identified. Quality threshold of 30, for example, can be specified by ```-q 30``` flag 
for overlaps with some mismatches to be counted. It can also remove adapter sequences optionally.  



## 3.1c. Assemble your data Using SPAdes

Make a folder for SPAdes assemblies: 

```sh
mkdir spades_assemblies
#cd spades_assemblies
```

You have to run SPAdes with 3 different flags and report the assembly time in the spreadsheet:
* with ```--sc``` only
* with ```--careful``` only
* with both ```--sc --careful```

To run SPAdes with --sc and --careful flags (our default assembly setting), type:  

* if you work with the **raw** reads  
>```sh
>time spades.py --sc --careful -t 8 -m 24 \
>-1 G5_${sample}_R1_001.fastq \
>-2 G5_${sample}_R2_001.fastq \
>-o spades_assemblies/G5_${sample}_sc_careful
>```

* if you work with the **merged** reads  
>```sh
>time spades.py --sc --careful -t 8 -m 24  \
>-1 merged_reads/G5_${sample}_merge_1.fastq.gz \
>-2 merged_reads/G5_${sample}_merge_2.fastq.gz \
>-s G5_${sample}_merged.fastq.gz \
>-o spades_assemblies/G5_${sample}_SeqPrep_sc_careful
>```

>You are now providing 3 input files to the assembler; 2 unmerged read pairs and 1 merged reads. After this assembly is done, 

Repeat the assemblies again but omitting *'--sc'* or *'--careful'* flags as in the previous exercises. 
Name the next two assemblies as *G5_Hiseq_SeqPrep_careful* and *G5_Hiseq_SeqPrep_sc*.

*Notice: Those previous commands will launch SPAdes assembly but also check how long the assembly takes. After the assembly has completed, check the time it took for the program to run. You should look at the 'real' time. 
Record the time in the spreadsheet table.
In general, typing the command ```time``` before other commands will help you check how long the computation took.*

**Using '--sc' flag**  
SPAdes can handle single-cell genomic data that is known to be highly biased in terms of sequence coverage along the length of the genome. 
In order for SPAdes to be able to handle biased sequence coverage, you need to supply the *'--sc'* flag when running the assembly.


**Using '--careful' flag**  
This flag uses *'bowtie'* tool to map the reads back to the contigs and check for errors due to bad quality sequences and correct these errors. 
This results in longer assembly times than not using the *'--careful'* flag.


## 3.2. Assessing assembly quality using Quast


After assembling the reads into contigs, you will use this tool called 'Quast' to calculate the basic metrics such as the length of largest contig, N50, etc.  
To run 'Quast', first we will create a link to the different contig files in the assembly folder. 
Use the following commands to help yourself:  

```sh
cd spades_assemblies
ln -s `pwd`/G5_${sample}_${merge}sc_careful/contigs.fasta G5_${sample}_${merge}sc_careful.fasta
```

then run quast on all assemblies at once:

```sh
#module load quast/2.3
quast.py -o assembly_metrics *.fasta
cd ..
```

Results from 'Quast' will be found in the 'assembly_metrics' folder. Take a look at the files produced by the program. 
The summary file containing the stats is 'report.txt'. You should see the assembly metrics such as N50, G+C%, largest contig, total length (i.e., total length of all contigs added).  

Report the results in the spreadsheet.

## 3.3. Gene prediction using Prodigal


Prodigal is a tool that can identify open reading frames (ORFs) in microbial genomes (bacteria or archaea). 
In this exercise, you will learn how to use Prodigal to predict ORFs and to prepare them for running Blastp later.  
Type the following commands to run Prodigal:  

For each assembly file run the following command, adapt this command according to your file name:  

```sh
prodigal \
-i G5_${sample}_${merge}sc_careful.fasta \
-a G5_${sample}_${merge}sc_careful.prodigal.faa \
-o G5_${sample}_${merge}sc_careful.prodigal.gff
```

Take a look at the files produced by Prodigal. For example, type  
```less G5_${sample}_${merge}sc_careful.prodigal.faa``` and see what the contents look like.  
Can you count how many genes are predicted by Prodigal?  
*Hint: Use the example given in the Part 1 to count the number of genes (look for a repeating pattern in the file and search for that pattern).* 
Tabulate the results in the spreadsheet.

## 3.4. Running completeness estimates


Of course, you will be interested how 'complete' your assembly is. In order to determine completeness, we make use of an in-house script that checks for the presence of a number of universally conserved genes. 
Based on how many of these genes will be identified, an estimation of genome completeness will be made.

From the assembly directory and run the following perl script (repeat the command for all prodigal files):

```sh
perl /proj/g2015028/nobackup/single_cell_exercises/scripts/micomplete.pl \
-h /proj/g2015028/nobackup/single_cell_exercises/scripts/Bact139.hmm \
-w /proj/g2015028/nobackup/single_cell_exercises/scripts/Bact139.weights \
-p G5_${sample}_${merge}sc_careful.prodigal.faa \
-t 8 -e -c 1e-10
```

If the script completed without any errors, you should see that the script printed something like this:

*Completeness: ??  
N markers found: ?? out of 139*  
For all asseblies, report the completeness in the spreadsheet.

The script looks for unique marker genes that are present in single copies in most bacteria and estimates how complete the genome is based on the marker genes found. Completeness is shown as fractions and you should multiply this with 100 to get the percentage completeness.



## 3.5. Identifying your cell

First, you will run this tool called 'rnammer' to predict the regions within assembled contigs that contain ribosomal RNA sequences (rRNA), such as 16S, 23S, and 5S rRNA sequences.  
Go into the folder containing the assembled contigs and type this command:  

```sh
cd ..
mkdir 16s_blast
cp spades_assemblies/*.fasta ./16s_blast
cd 16s_blast
rnammer -S bac -m lsu,ssu,tsu \
-gff G5_${sample}_${merge}sc_careful.rnammer.gff \
-f G5_${sample}_${merge}sc_careful.rnammer.fasta \
< G5_${sample}_${merge}sc_careful.fasta
```

Take a look at the file (*'contigs.rnammer.gff'*) produced by 'rnammer'.  
Can you identify the positive matches predicted by 'rnammer'?  
Can you interpret the result output?  
Next, you will run 'blastn' against the 'Silva' database. 
Depending on whether or not you have identified 16S or 23S sequences, you will need to run 'blastn' on a different database.  

If you have identified a 16S rRNA sequence (check the gff file), type:

```sh
DB=/proj/g2015028/nobackup/single_cell_exercises/databases/SSURef_NR99_115_tax_silva_trunc.dna.fasta
blastn -db $DB -evalue 1e-6 -num_threads 8 -query G5_${sample}_${merge}sc_careful.rnammer.fasta \
-out G5_${sample}_${merge}sc_careful.rnammer.16S_silva.blastn
```

If you have identified 23S rRNA sequence, then type:

```sh
DB=/proj/g2015028/nobackup/single_cell_exercises/databases/LSURef_115_tax_silva_trunc.dna.fasta
blastn -db $DB -evalue 1e-6 -num_threads 8 -query G5_${sample}_${merge}sc_careful.rnammer.fasta \
-out G5_${sample}_${merge}sc_careful.rnammer.23S_silva.blastn
```

After running Blastn, can you identify what organism G5 belongs to?

For the next part you need to come back in the 


## Questions:


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
Next, make symbolic links of sequences in that folder:

```sh
ln -s /proj/g2015028/nobackup/single_cell_exercises/sequences/dataset1/G5_Hiseq_R1_001.fastq .
ln -s /proj/g2015028/nobackup/single_cell_exercises/sequences/dataset1/G5_Hiseq_R2_001.fastq .
ln -s /proj/g2015028/nobackup/single_cell_exercises/sequences/dataset2/G5_Miseq_R1_001.fastq .
ln -s /proj/g2015028/nobackup/single_cell_exercises/sequences/dataset2/G5_Miseq_R2_001.fastq .
```

Now you are almost ready to run assemblies! But before you can start assemblies, you need to load SPAdes module first.  
To load the SPAdes assembler, type:

```sh
module load bioinfo-tools
module load spades/3.1.1
```  
--->  


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
--->
