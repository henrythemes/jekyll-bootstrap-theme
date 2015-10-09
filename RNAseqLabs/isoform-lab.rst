=========================================
Reference guided assembly using Cufflinks
=========================================
Assembly of genes and isoforms using an assembly usig cufflinks is a two step process. 
First you map all the reads from your experiment to the reference sequence then you run another 
step where you use the mapped reads to idenitify potenital genes including introns and exons.  


Data available for exercise
---------------------------

All fastqfiles that you will need for this exercise can be found in 
``/proj/b2013006/webexport/downloads/courses/RNAseqWorkshop/isoform/referenceBased/data`` on UPPMAX and through this `URL <https://export.uppmax.uu.se/b2013006/downloads/courses/RNAseqWorkshop/isoform/referenceBased/data>`_ .

If you want to map more files for practise you can continue with the files found in 
``/proj/b2013006/webexport/downloads/courses/RNAseqWorkshop/isoform/RAB11FIP5_fastqFiles`` on UPPMAX and through this `URL <https://export.uppmax.uu.se/b2013006/downloads/courses/RNAseqWorkshop/isoform/RAB11FIP5_fastqFiles>`_ .

There are also other files that have been pre mapped and pre analyzed using STAR and cufflinks and those can be found in
``/proj/b2013006/webexport/downloads/courses/RNAseqWorkshop/isoform/otherData`` on UPPMAX and through this `URL <https://export.uppmax.uu.se/b2013006/downloads/courses/RNAseqWorkshop/isoform/otherData>`_ .
 

Mapping short reads to a reference
----------------------------------

Here we will map the reads to the hg19 reference genome using a popular RNA-seq 
aligner, **STAR**. There are many many features that can be tweaked using STAR. 
Read below for the flags we use for this exercise. Remember to change names accordingly 
so that you can run your program properly and know which files you have used.

To load the STAR package on Uppmax, execute::

     module load bioinfo-tools
     module load star
     module load samtools

Now you can map the reads from one of the samples (or several; it's up to you 
which ones(s)) to map using a command such as the one below. ::
  
  mkdir outDir
    
  STAR  --genomeDir /proj/b2013006/webexport/downloads/courses/RNAseqWorkshop/reference/hg19_Gencode14.overhang75  --readFilesIn sample1_RAB11FIP5_1.fastq sample1_RAB11FIP5_2.fastq --runThreadN 2 --outSAMstrandField intronMotif --outFileNamePrefix outDir/
	
flags used are 

* ``--genomeDir /path/to/STARgenome`` is the path to a pre-rendered reference library that STAR uses to map reads to the genome. 

*  ``--readFilesIn /path/to/read1 [/path/to/read2 ]`` is where you should add your fastq files that you will map to the reference.

*  ``--runThreadN N`` is the number of threads that will be used by the program.

*  ``--outSAMstrandField intronMotif`` adds SAM information so that unstranded reads function with cufflinks 

*  ``--outFileNamePrefix outDir`` tells you were the output data ends up. 


  
This should run fairly quickly and put a file called ``Aligned.out.sam`` in 
the directory that you specified with ``--outFileNamePrefix``. 

Next step is to convert it to a BAM file and rename it to a more representable name. 
A good naming praxis is to name the file to correspond what you mapped. As an example if you mapped sample 12
you should rename the mapped SAM file to a file with the name ``sample12_RAB11FIP5.bam`` 
The renaming and BAMfile conversion can be done in one step. Then to view them you also have to sort the hits and index them: ::

  
  samtools view -bSh -o sampleNN_RAB11FIP5.bam /path/to/Aligned.out.sam
  samtools sort sampleNN_RAB11FIP5.bam  sampleNN_RAB11FIP5.sorted
  samtools index sampleNN_RAB11FIP5.sorted.bam


The sorted, indexed bam file can then be viewed in IGV. 


Reference guided assembly using Cufflinks
-----------------------------------------

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
`iReckon <http://compbio.cs.toronto.edu/ireckon/>`_ and 
`SLIDE <https://sites.google.com/site/jingyijli/>`_. These may require some 
annotation as input but they can discover (and quantify) new isoforms. 




Quantification with Cufflinks
=============================

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








