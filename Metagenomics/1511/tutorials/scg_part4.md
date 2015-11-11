---
layout: default
title:  'Part 4: Coverage plotting, chimera detection and inspection'
---

# Part 4: Coverage plotting, chimera detection and inspection

As explained during the lectures, the genome amplification process (MDA) results in a coverage bias across the genome and induces chimera formation. 
Here you will have a look at the 'best' assembly that you managed to obtain, and try to assess the level of bias and chimera formation.
In order to visualize how the coverage a normal single-cell assembly, you need to map your reads used in the assembly, back to your assembled contigs. 
This can be done with various read-mappers such as bwa, bowtie, gsnap etc. You will be using bwa since it is fast and it can also map partial reads back to our reference.
Make a folder called 'coverage_mapping' in the 'metagenomics_exercises' folder.  
Copy the contigs file from an example run:

```sh
cp /proj/g2014180/nobackup/single_cell_exercises/assembly_examples/contigs.fasta .
```

First you need to create an index for your contigs:

```sh
module load bioinfo-tools
module load bwa/0.7.5a
bwa index contigs.fasta #[time to run: < 1 sec]
```

This should produce 5 more files in your directory called contigs.fasta.[amb|ann|bwt|pac|sa]. 
Make sure you know where the fastq files are, in this example the location is in a directory above current location, modify '../' if needed.  
Then map your reads back to your contigs:

```sh
bwa mem -t 8 contigs.fasta ../spades_assemblies/G5_Hiseq_1.fastq ../spades_assemblies/G5_Hiseq_2.fastq > G5_vs_contigs.sam [time to run: 2.5 sec]
```

This will produce a SAM file (Sequence Alignment/Map) which is a flat text file. 
These files are often quite large, so converting it to a smaller, binary, format is often a good idea. 
To do this, and to do some additional manipulations to the mapped reads, you will be using picard-tools.

```sh
module load picard/1.92
java -jar /sw/apps/bioinfo/picard/1.92/milou/SamFormatConverter.jar INPUT=G5_vs_contigs.sam OUTPUT=G5_vs_contigs.bam [time to run: 13 sec]
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

## 4.1 Assessing coverage bias in Artemis
---

To view the reads mapped back to our contigs you can use artemis. In order to load your BAM file, you must first index it.  
You can do this with picard-tools BuildBamIndex:

```sh
java -jar /sw/apps/bioinfo/picard/1.92/milou/BuildBamIndex.jar INPUT=G5_vs_contigs_sorted.bam O=G5_vs_contigs_sorted.bam.bai
```

Open artemis with the following command:  

```sh
module load artemis/15.0.0
art &
```

Open your assembled contigs with *File -> Open -> contigs.fasta*
This will open all your contigs in a single window, concatenated one after the other. 
In this newly opened window, you can load your BAM file by choosing File -> Read BAM / VCF - Select. 
Make sure you select the sorted BAM file and click 'OK'.
You can now see all the reads mapped back to our assembled contigs. 
Inspect the various contigs and see how the coverage differs quite a lot on the various contigs.
**It is important to keep in mind that this is not a whole genome, but rather concatenated contigs.**
The order of the contigs here have no biological relevance, they are completely arbitrary.
By right-clicking in the top most are you can select *Graph -> Coverage* to show a coverage plot of the viewed region. 
One thing to note with this plot is that it is a sliding window mean, and the size of the sliding window is affected by how large region you are viewing at the time. 
In the top-right corner the maximum coverage is shown. By zooming in, using the middle slider, you can set the zoom level.
Have a look at the contigs called *'NODE_7'* and *'NODE_4'*. 
They show some interesting coverage patterns. You can easily find them by going *Goto -> Navigator (or Ctrl-G)*. 
Then fill in *'NODE_7_'* or *'NODE_4_'* under the 'Goto Feature with this qualifier value' -> Goto. 
If for any reason you can't type in the contig names in the search box, you should browse the list of contig names at the bottom of Artemis browser and double click on the contig name.


## 4.2 Detection and inspection of chimeric reads
---

The MDA reaction will produce chimeric DNA in several ways. 
If you want to get a more detailed picture of how this is happening, a very good starting point is the Lasken et al. paper from 2007. 
See http://www.biomedcentral.com/1472-6750/7/19
One issue with single-cell data is that chimeric DNA will be formed and sequenced just like any other DNA. 
This will be assembled into chimeric contigs which can be quite hard to distinguish from correct contigs. 
One way of getting around this is by using multiple SAGs from the same organism. 
If they have been amplified independently, chimeras will form, but there will be different chimeras in each SAG. 
By co-assembling n-1 SAGs together, and mapping the reads from the excluded SAG against this co-assembly, 
chimeric reads can be identified since they will be split and mapped to different contigs, or parts of contigs.
However, in most cases, like the one we are examining today, there is only 1 SAG. 
Identifying chimeras in a single SAG is quite difficult since the chimeric mapped reads will map perfect against the assembled chimeric contigs. 
We can look for traces of this though. If the chimeric reads were amplified late, resulting in very few reads actually being chimeric, 
one might still catch these. To find chimeric reads, we return to our SAM file. The sixth field in the SAM file holds the mapped reads CIGAR value. 
This value shows you how a read has been mapped in terms of matches and mismatches. 
The most common value you would see if you open your SAM file would be “100M”. This means that 100 bases of our 100 bp read was mapped against the contig.
We can also find a few that have a partial mapping. This usually looks something like “40S60M”, meaning that the first 40 bases were masked off, 
then the following 60 bases mapped against our contig. 
If the masked of part of our read is mapped against another coordinate we can tell by the twelfth field in the SAM file. 
It is prefixed with “SA:Z”. To find these kinds of reads we can use ‘grep’.

```sh
java -jar /sw/apps/bioinfo/picard/1.92/milou/ViewSam.jar I=G5_vs_contigs_sorted.bam | grep "SA:Z:"
```

If we also want to count the number of such reads we can use greps '-c' flag.

```sh
java -jar /sw/apps/bioinfo/picard/1.92/milou/ViewSam.jar I=G5_vs_contigs_sorted.bam | grep -c "SA:Z:"
# Or we could also use 'wc' (word count)
java -jar /sw/apps/bioinfo/picard/1.92/milou/ViewSam.jar I=G5_vs_contigs_sorted.bam | grep "SA:Z:" | wc -l
```

One thing to keep in mind is that some reads will be at the very ends of contigs, and that could be the reason why they are not mapped in full-length. 
To correct for this we modify our command with some awk-magic.

```sh
java -jar /sw/apps/bioinfo/picard/1.92/milou/ViewSam.jar I=G5_vs_contigs_sorted.bam | grep "SA:Z:" | awk '{split($3,array,"_"); if($4 > 200 && $4 < (array[4]-300)) print $0}' | wc -l
```

For those of you unfamiliar with awk, the command gets the length of each contig by looking in the contig name which looks something like
NODE_24_length_2520_cov_2.03245_ID_101896
Awk then splits this into parts separated by "_" and stores them in a data-structure called an array. 
We then tell awk to only print the line ($0) if the first position of our reads (fourth field in the SAM file, $4) 
is more than 200 bases from the start, or more than than 300 bases from the end. 
Keep in mind that all reads with an alternative alignment isn’t necessarily chimeric in origin. 
Some reads might be split due to assembly-errors or sequencing errors. If we want to do a better job at finding chimeras we would need more clonal SAGs.

## 4.3 Insert size
---

When you’ve done your read-mapping, you can also inspect the insert-sizes, size of the actual sequenced fragments, quite easily. 
Picard-tools have a ready made tool for this called CollectInsertSizeMetrics.

```sh
java -jar /sw/apps/bioinfo/picard/1.92/milou/CollectInsertSizeMetrics.jar HISTOGRAM_FILE=G5_vs_contigs_inssizePlot.pdf INPUT=G5_vs_contigs_sorted.bam OUTPUT=G5_vs_contigs_inssize.out #[time to run = 2.5 sec]
```

You can inspect the insert size distribution by opening the pdf file in firefox.

```sh
firefox G5_vs_contigs_inssizePlot.pdf
#or Okular.
okular G5_vs_contigs_inssizePlot.pdf
```

### Questions
---

**Q4.1** What was the highest coverage for any of your contigs?  
**Q4.2** Do you think that the completeness of your SC genome will improve with more sequence data? Why? Will your genome ever be ‘complete’?  
**Q4.3** How would you explain the coverage patterns you can see on the contigs called NODE_7 and NODE_4.  
**Q4.4** How many chimeric reads did you find in your original dataset? Do you think this value is a true representation of the total amount of chimeras?  
**Q4.5** What was the average insert size in this dataset?  
