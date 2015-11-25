---
layout: default
title:  'Part 4: Coverage plotting, chimera detection and inspection'
---

# Part 4: Coverage plotting, chimera detection and inspection

## 4.1. Reads mapping

Mapping can be done with various read-mappers such as *bwa, bowtie, gsnap* etc... You will be using bwa since it is fast and it can also map partial reads back to our reference.  

Make a folder called 'coverage_mapping' in the *~/single_cell_exercises/datasetX* folder and move into it.  
*Notice: Replace X in datasetX by the dataset number you were working on in the previous part.*  
Copy the contigs file from an example run:

```sh
cp /proj/g2015028/nobackup/single_cell_exercises/assembly_examples/contigs.fasta .
```

First you need to create an index for your contigs:

```sh
bwa index contigs.fasta #[time to run: < 1 sec]
```

This should produce 5 more files in your directory called contigs.fasta.[amb|ann|bwt|pac|sa]. 
Make sure you know where the fastq files are, in this example the location is in a directory above current location, modify '../' if needed.  
Then map your reads back to your contigs:

```sh
bwa mem -t 8 contigs.fasta ../G5_Hiseq_R1_001.fastq ../G5_Hiseq_R2_001.fastq > G5_vs_contigs.sam #[time to run: 2.5 sec]
```

This will produce a SAM file (Sequence Alignment/Map) which is a flat text file. 
These files are often quite large, so converting it to a smaller, binary, format is often a good idea. 
To do this, and to do some additional manipulations to the mapped reads, you will be using picard-tools.

```sh
java -jar /sw/apps/bioinfo/picard/1.92/milou/SamFormatConverter.jar INPUT=G5_vs_contigs.sam OUTPUT=G5_vs_contigs.bam #[time to run: 13 sec]
```

We can now safely remove our original SAM file in order to save space

```sh
rm G5_vs_contigs.sam
```

The BAM file (Binary Alignment/Map) is, as the name implies, a binary file, and can therefor not be viewed with normal tools such as less or more. 
If you still want to see the file as plain text, you can use picard-tools ViewSam. 
This tool prints the whole file to the screen. For easier viewing, pipe the output to 'less'.

```sh
java -jar /sw/apps/bioinfo/picard/1.92/milou/ViewSam.jar INPUT=G5_vs_contigs.bam | less
```

If you want to know more about SAM files, have a look at http://samtools.sourceforge.net/SAMv1.pdf. 
Especially section 1.4 is helpful for understanding the different columns. 
In order to more easily make sense of the information in the alignment file, you should sort it. 
If you viewed your bam file, you might have noticed that the aligned reads comes in the same order as in the input reads. 
This is not tremendously helpful. A more easily interpretative way would be to have the aligned reads show up in the order they appear along the various contigs.  
Or to put it another way, with increasing coordinates. This can be achieved by using picard-tools SortSam:

```sh
java -jar /sw/apps/bioinfo/picard/1.92/milou/SortSam.jar INPUT=G5_vs_contigs.bam OUTPUT=G5_vs_contigs_sorted.bam SORT_ORDER=coordinate #[time to run: 25 sec]
```

If you now view the sorted bam file with *ViewSam.jar*, you will see that the contigs appear in order after each other.


<div>
 <span style="float:left"><a class="btn btn-primary" href="scg_part4"> Previous page</a></span>
 <span style="float:right"><a class="btn btn-primary" href="scg_part4_2"> Next page</a></span>
</div>
