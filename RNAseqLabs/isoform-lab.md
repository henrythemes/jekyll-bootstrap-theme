---
layout: default
title:  'Isoform Reference'
---

#Reference guided assembly using Cufflinks

Assembly of genes and isoforms using an assembly usig cufflinks is a two step process. 
First you map all the reads from your experiment to the reference sequence. For information on how to map reads to a reference go to the [mapping reads tutorial](mapping_reads).   Then you run another step where you use the mapped reads to idenitify potenital genes including introns and exons.  




##Data available for exercise

There are also other files that have been pre-mapped and pre-analyzed using STAR and cufflinks and those can be found in
``/proj/b2013006/webexport/downloads/courses/RNAseqWorkshop/isoform/otherData`` on UPPMAX and through this `URL <https://export.uppmax.uu.se/b2013006/downloads/courses/RNAseqWorkshop/isoform/otherData>`_ .
 



##Reference guided assembly using Cufflinks

Cufflinks can do reference based assembly, which means 
that it tries to discover transcripts, disregarding gene annotation (actually there
is an option to use it as well but we will ignore that for now), just based on the 
mapped reads to the genome. This functionality works even on our small files.

Try to do a reference guided assembly. This is done simply by running Cufflinks 
without feeding it a GTF file with the -G flag::

     module load cufflinks/2.2.1
     cufflinks -o /path/to/outputDirectory sampleNN_RAB11FIP5.sorted.bam

Substitute the appropriate names for the BAM file and the output directory. When 
Cufflinks has finished (which should hardly take any time at all), the output 
directory will contain a file called ``transcripts.gtf``. Rename the file to a 
name that reflects what created the GTF file.  You can import that in 
the usual way into IGV as a track.

Was Cufflinks able to assemble your alignments into something that makes sense?
 
Other alternatives for reference-based assembly include 
`Scripture <http://www.broadinstitute.org/software/scripture>`_, 





##Quantification with Cufflinks

Cufflinks is a well-known software package for estimating gene and isoform 
abundances in a BAM or SAM files based on an annotation file in GTF format 
. However, we run 
into problems here because of the small file size. Cufflinks needs a certain amount 
of data to be able to do its estimations so it will not work very well on our small 
alignment files. Therefore we have run it on the full alignment files on each sample 
and provided the merged results at ``/proj/b2013006/webexport/downloads/courses/RNAseqWorkshop/isoform/otherData/isoform_fpkm_table.txt``
(isoform FPKM) and ``/proj/b2013006/webexport/downloads/courses/RNAseqWorkshop/isoform/otherData/RNAseqfpkm_table.txt``.
If you are interested, you can check in the isoform FPKM file how much each isoform 
is deemed to have been expressed by Cufflinks, but we will not go into that today. 
(These results will be revisited tomorrow, in the differential expression lab exercise.)

For reference, the commands we used were of the form::

     # Only for reference, does not need to be executed during the exercise
     cufflinks -p 8 -G /proj/b2013006/webexport/downloads/courses/RNAseqWorkshop/isoform/otherData/Homo_sapiens.GRCh38.77.fixed.gtf -o cufflinks_out_137_1 accepted_hits_137_1.bam

The ``-G`` option points to an annotation file in GTF format for which to calculate
FPKM values. The input here is a BAM file which is just a binary version of a SAM file.  

Other options for doing abundance estimation are `RSEM <http://deweylab.biostat.wisc.edu/rsem/>`_ 
or the flexible `RPKMforgenes.py script <http://sandberg.cmb.ki.se/media/data/rnaseq/instructions-rpkmforgenes.html>`_.








