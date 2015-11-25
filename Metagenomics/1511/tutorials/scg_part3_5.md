---
layout: default
title:  'Part 3: Single cell genome assembly using SPAdes'
---

# Part 3: Single cell genome assembly

## 3.5. Gene prediction using Prodigal


Prodigal is a tool that can identify open reading frames (ORFs) in microbial genomes (bacteria or archaea). 
In this exercise, you will learn how to use Prodigal to predict ORFs and to prepare them for running Blastp later.  
Type the following commands to run Prodigal:  

For each assembly file run the following command. This command is based on the output file for *IDBA* assembler, do the same for Ray and SPAdes output.  
*Hint: use the command ```ls *.fasta``` to see all assembly files.*  
Please, adapt this command according to your file name:  

```sh
prodigal \
-i G5_${sample}${trim}_IDBA_contig.fasta \
-a G5_${sample}${trim}_IDBA_contig.prodigal.faa \
-o G5_${sample}${trim}_IDBA_contig.prodigal.gff
```

Take a look at the files produced by Prodigal. For example, type  
```less G5_${sample}${trim}_IDBA_contig.prodigal.faa``` and see what the contents look like.  
Can you count how many genes are predicted by Prodigal?  
*Hint: Use the example given in the Part 1 to count the number of genes (look for a repeating pattern in the file and search for that pattern using grep).* 
Tabulate the results in the spreadsheet.

## 3.6. Running completeness estimates


Of course, you will be interested to know how 'complete' is your assembly. In order to determine completeness, we make use of an in-house script that checks for the presence of a number of universally conserved genes. 
Based on how many of these genes will be identified, an estimation of genome completeness will be made.

From the assembly directory and run the following perl script (repeat the command for all prodigal files):

```sh
# Fix dep
module load hmmer/3.1b1-gcc
# End fix
perl /proj/g2015028/nobackup/single_cell_exercises/scripts/micomplete.pl \
-h /proj/g2015028/nobackup/single_cell_exercises/scripts/Bact139.hmm \
-w /proj/g2015028/nobackup/single_cell_exercises/scripts/Bact139.weights \
-p G5_${sample}${trim}_IDBA_contig.prodigal.faa \
-t 8 -e -c 1e-10
```

If the script completed without any errors, you should see that the script printed something like this:

*Completeness: ??  
N markers found: ?? out of 139*  
For all assemblies, report the completeness in the spreadsheet.

The script looks for unique marker genes that are present in single copies in most bacteria and estimates how complete the genome is based on the marker genes found. Completeness is shown as fractions and you should multiply this with 100 to get the percentage completeness.

<div>
 <span style="float:left"><a class="btn btn-primary" href="scg_part3_4"> Previous page</a></span>
 <span style="float:right"><a class="btn btn-primary" href="scg_part3_7"> Next page</a></span>
</div>
