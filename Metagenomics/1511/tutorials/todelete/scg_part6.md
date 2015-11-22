---
layout: default
title:  'Part 6: Comparing MiSeq vs. HiSeq data'
---

# Part 6: Comparing MiSeq vs. HiSeq data
---

Thus far, you have been working with paired end HiSeq data. Depending on your needs and project design, 
you might be working with NGS data that has been generated on a MiSeq sequencer. 
Here you will be working with a paired end (2x250 bp) MiSeq dataset (dataset2) from the same sequencing library (hence the same ‘cell’), 
allowing you to compare assemblies with different datasets generated from the same libraries. 
As before, you will be assembling with SPAdes, exploring assemblies generated with different flags and filtering settings, see Table 2. 
Again, given the time consuming nature of the assemblies, try to distribute the jobs between different groups.
Here, you will need to figure out the commands needed to successfully run assemblies using MiSeq data. 
The commands are very similar and you should be able to figure out how to run the assembly correctly using HiSeq commands as examples. 
Check the examples given in Part 3.1. Note that you will run assemblies on untrimmed/unmerged MiSeq data (original data) and 
also on merged (data processed with SeqPrep).

Run ‘Quast’ and ‘micomplete’ tools to calculate metrics and completeness estimates and record these in the table.

