---
layout: default
title:  'Part 6: Comparing MiSeq vs. HiSeq data'
---

# Part 7: Analysis of a novel single cell genome (bonus exercise)
---

Now that you have worked with several single cell genome datasets, you will be able to play around with a bonus dataset of a mystery single cell (dataset3). With the knowledge that you have acquired during this course, try to:

- make a decent assembly
- assess the completeness of your single cell
- try to identify your cell using rnammer or MEGAN
- try to asses if the dataset suffers from contamination
- find out any interesting features that might be encoded by the genome

In this bonus exercise, you will work with MiSeq data produced from a different Single-cell Amplified Genome (SAG) than G5. 
The assemblies here will be carried out in pretty much the same way as those examples shown in Part 3.1.
An overview of the steps are listed below:

1. Run an assembly using both *'--sc'* and *'--careful'* flags
2. Run an assembly using different k-mers set
3. Trim/merge the read pairs using SeqPrep and repeat 1 and 2
4. Run *'Quast'* and *'micomplete'* tools
5. Run Prodigal and RNAmmer to predict ORFs and rRNA regions
6. Run Blastn and Blastp to check contaminants
7. Blast the rRNA genes identified against Silva database to identify the organism
8. Tabulate the results

In this bonus exercise, a different parameter will be introduced, i.e., k-mers. SPAdes in default mode runs with **k-mers of 21, 33, and 55**. 
In this exercise, you will set the **k-mers to 55, 77, and 99**. To set these k-mers, you need to provide this parameter when running SPAdes:

```
-k 55,77,99
```

Record the results after the assemblies and assembly statistics are calculated.

<p>
<style type="text/css">
.tg  {border-collapse:collapse;border-spacing:0;border-color:#999;}
.tg td{font-family:Arial, sans-serif;font-size:14px;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#999;color:#444;background-color:#F7FDFA;}
.tg th{font-family:Arial, sans-serif;font-size:14px;font-weight:normal;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#999;color:#fff;background-color:#26ADE4;}
.tg .tg-yw4l{vertical-align:top}
.tg .tg-pxng{background-color:#ffffff;color:#000000;vertical-align:top}
.tg .tg-25to{background-color:#26ade4;font-weight:bold;color:#ffffff;vertical-align:top}
.tg .tg-6k2t{background-color:#D2E4FC;vertical-align:top}
</style>
<table class="tg">
  <tr>
    <th class="tg-yw4l">Table 5</th>
    <th class="tg-pxng" colspan="2">MiSeq data - Original data</th>
  </tr>
  <tr>
    <td class="tg-25to">(Dataset 3 - N21)</td>
    <td class="tg-25to">With --sc and --careful flags</td>
    <td class="tg-25to">using different k-mers (eg: 55,77,99)</td>
  </tr>
  <tr>
    <td class="tg-yw4l">Number of reads</td>
    <td class="tg-yw4l">-</td>
    <td class="tg-yw4l">-</td>
  </tr>
  <tr>
    <td class="tg-6k2t">Assembly time</td>
    <td class="tg-6k2t">-</td>
    <td class="tg-6k2t">-</td>
  </tr>
  <tr>
    <td class="tg-yw4l">Number of contigs</td>
    <td class="tg-yw4l">-</td>
    <td class="tg-yw4l">-</td>
  </tr>
  <tr>
    <td class="tg-6k2t">Total assembly size</td>
    <td class="tg-6k2t">-</td>
    <td class="tg-6k2t">-</td>
  </tr>
  <tr>
    <td class="tg-yw4l">Largest contig</td>
    <td class="tg-yw4l">-</td>
    <td class="tg-yw4l">-</td>
  </tr>
  <tr>
    <td class="tg-6k2t">N50</td>
    <td class="tg-6k2t">-</td>
    <td class="tg-6k2t">-</td>
  </tr>
  <tr>
    <td class="tg-yw4l">G+C%</td>
    <td class="tg-yw4l">-</td>
    <td class="tg-yw4l">-</td>
  </tr>
  <tr>
    <td class="tg-6k2t">Number of ORFs</td>
    <td class="tg-6k2t">-</td>
    <td class="tg-6k2t">-</td>
  </tr>
  <tr>
    <td class="tg-yw4l">Completeness (%)</td>
    <td class="tg-yw4l">-</td>
    <td class="tg-yw4l">-</td>
  </tr>
</table>
</p>

<p>
<style type="text/css">
.tg  {border-collapse:collapse;border-spacing:0;border-color:#999;}
.tg td{font-family:Arial, sans-serif;font-size:14px;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#999;color:#444;background-color:#F7FDFA;}
.tg th{font-family:Arial, sans-serif;font-size:14px;font-weight:normal;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#999;color:#fff;background-color:#26ADE4;}
.tg .tg-yw4l{vertical-align:top}
.tg .tg-pxng{background-color:#ffffff;color:#000000;vertical-align:top}
.tg .tg-25to{background-color:#26ade4;font-weight:bold;color:#ffffff;vertical-align:top}
.tg .tg-6k2t{background-color:#D2E4FC;vertical-align:top}
</style>
<table class="tg">
  <tr>
    <th class="tg-yw4l">Table 6</th>
    <th class="tg-pxng" colspan="2">MiSeq data - Processed with SeqPrep<br></th>
  </tr>
  <tr>
    <td class="tg-25to">(Dataset 3 - N21)</td>
    <td class="tg-25to">With --sc and --careful flags</td>
    <td class="tg-25to">using different k-mers (eg: 55,77,99)</td>
  </tr>
  <tr>
    <td class="tg-yw4l">Number of reads</td>
    <td class="tg-yw4l">-</td>
    <td class="tg-yw4l">-</td>
  </tr>
  <tr>
    <td class="tg-6k2t">Assembly time</td>
    <td class="tg-6k2t">-</td>
    <td class="tg-6k2t">-</td>
  </tr>
  <tr>
    <td class="tg-yw4l">Number of contigs</td>
    <td class="tg-yw4l">-</td>
    <td class="tg-yw4l">-</td>
  </tr>
  <tr>
    <td class="tg-6k2t">Total assembly size</td>
    <td class="tg-6k2t">-</td>
    <td class="tg-6k2t">-</td>
  </tr>
  <tr>
    <td class="tg-yw4l">Largest contig</td>
    <td class="tg-yw4l">-</td>
    <td class="tg-yw4l">-</td>
  </tr>
  <tr>
    <td class="tg-6k2t">N50</td>
    <td class="tg-6k2t">-</td>
    <td class="tg-6k2t">-</td>
  </tr>
  <tr>
    <td class="tg-yw4l">G+C%</td>
    <td class="tg-yw4l">-</td>
    <td class="tg-yw4l">-</td>
  </tr>
  <tr>
    <td class="tg-6k2t">Number of ORFs</td>
    <td class="tg-6k2t">-</td>
    <td class="tg-6k2t">-</td>
  </tr>
  <tr>
    <td class="tg-yw4l">Completeness (%)</td>
    <td class="tg-yw4l">-</td>
    <td class="tg-yw4l">-</td>
  </tr>
</table>
</p>

### Questions
---

**Q7.1** Did you notice the major differences between the two assemblies (untrimmed reads)?  
What about the differences between two assemblies for trimmed reads?  
**Q7.2** What do you think happened by setting a longer k-mer sizes than default settings?  
**Q7.3** Do you think the number of reads are enough to obtain a good assembly or if more sequences should be obtained?  
**Q7.4** What is the taxonomic affiliation of the SAG?  
**Q7.5** Can you say anything about the metabolism based on the assembled data?  
