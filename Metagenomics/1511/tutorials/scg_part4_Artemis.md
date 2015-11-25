---
layout: default
title:  'Part 4: Coverage plotting, chimera detection and inspection'
---

## 4.x Assessing coverage bias in Artemis
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


