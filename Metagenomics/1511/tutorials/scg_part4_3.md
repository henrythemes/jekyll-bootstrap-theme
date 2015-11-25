---
layout: default
title:  'Part 4: Coverage plotting, chimera detection and inspection'
---

# Part 4: Coverage plotting, chimera detection and inspection

## 4.3 Detection and inspection of chimeric reads
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


<div>
 <span style="float:left"><a class="btn btn-primary" href="scg_part4_2"> Previous page</a></span>
 <span style="float:right"><a class="btn btn-primary" href="scg_part4_4"> Next page</a></span>
</div>
