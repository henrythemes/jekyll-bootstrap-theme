====================================================
Viewing data in the Integrated Genomics Viewer (IGV)
====================================================


In these exercises we will use the  Integrated Genomics Viewer (IGV). 
IGV is installed on UPPMAX, but since you will then access it over a network connection, the graphics might be sluggish, 
so we recommend that you download and run the browsers locally and download the files that you will look at 
locally. 

If you want to try to run IGV on Uppmax, please refer to the 
`Uppmax page for IGV instructions <http://www.uppmax.uu.se/starting-integrative-genomics-viewer-igv-on-milou>`_ 
for advice.  

Visualising a BED file or GTF file                                                          
================================

In IGV, select File > Load from File ... and navigate to the BED/GTF file (on 
Uppmax according to above or to the location where you downloaded it locally). From 
the BED/GTF file extension, IGV will automatically treat the information in the file accordingly. 

Visualising a .sam or .bam file
===============================
If the file output is in the SAM format (file ends with .sam), which is a uncompressed test format, 
you need to convert it to the BAM format (file ends with .bam). A BAM file contains the same information 
as the SAM file but now it is in a binary compressed format unreadable for a human. 

To convert a SAM file to BAM format type: ::

     samtools view -bS fileName.sam >fileName.bam


Before the visualization, you need to sort it and then build a so-called 
BAM file index::

     samtools sort fileName.bam fileName.sorted
     samtools index fileName.sorted.bam

Then, in IGV, just choose File > Load from File ... and select the BAM file. 
If you are running IGV locally and did the mapping on Uppmax, you will have to 
download the BAM file and the corresponding index file (.bai) from your work folder 
to your own computer first.


Going to a locus
================

In IGV, click the textbox which has the word 
Go on its right hand side. This textbox will typically contain genomic coordinates for 
the locus you are presently looking at, but it can also be used to find gene locations. 
Type locus name in it and press Enter.




