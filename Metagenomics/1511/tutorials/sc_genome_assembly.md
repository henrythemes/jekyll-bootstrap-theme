---
layout: default
title:  'Single Cell Genomics Tutorial'
---

# Introduction: Single Cell Genomics Tutorial

<p>Single cell genomics is an emerging technology that allows one to explore the genome sequence of individual cells. 
During this tutorial you will work with real single cell genome data from a single cell that was isolated from a hot spring in Yellowstone National Park (USA). 
The data you will work with is part of a larger project ('PUZZLE_CELL') that aims to identify and genomically probe novel prokaryotic lineages, and to gain insight in the origin and evolution of life.  </p>
The data that you will work with is paired-end Illumina reads. 
We have chosen to have you work with both HiSeq ( 2x100 bp - dataset 1 ) and MiSeq ( 2x250 bp - dataset 2 ) datasets. 
Both datasets were generated from the same single cell, hence allowing you to develop a feeling for what you might want to used in any potential future SCG project. 
In addition, there is a third MiSeq dataset (2x250 bp - dataset 3) that contains a completely new organism, which might be a bit more challenging to work with. Depending on how much time you have at the end of the tutorial you can choose what you want to do. It would be good to at least take some time to think about what kind of analysis would give you answer to the question 'What is this cell?' There are several more suggestions for optional exercises that you can choose if you happen to have some extra time at the end of the tutorial.

## Overview of steps in this exercise

We will have a lunch break. We will split the discussion of results into two parts, one when you have the assembly results and then another one to summarize the whole day. 
Here is a schematic [Workflow](../slides/scg_workflow.pdf) of what you will be doing.  

1. [Connecting to UPPMAX](connectToUppmax)  
2. [Familiarizing with data](scg_part2)  
3. [Single-cell genome assembly](scg_part3)  
3.1. [Organize working folder](scg_part3_1)  
3.2. [Pre-processing](scg_part3_2)  
3.3. [Assembly](scg_part3_3)  
3.4. [Assessing assembly quality using Quast](scg_part3_4)  
3.5. [Gene prediction using Prodigal](scg_part3_5)  
3.6. [Running completeness estimates](scg_part3_5)  
3.7. [Identifying ribosomal RNAs](scg_part3_7)  
4. [Assessing read coverage and chimera checking](scg_part4)  
4.1 [Reads mapping](scg_part4_1)  
4.2 [Assessing coverage bias](scg_part4_2)  
4.3 [Detection and inspection of chimeric reads](scg_part4_3)  
4.4 [Insert size] (scg_part4_4)  
5. [Exploring your single cell genome](scg_part5)  
5.1 [Contamination analysis in MEGAN](scg_part5_1)  
5.2 [Functional analysis in MEGAN](scg_part5_2)  
6. [Analysis of a novel single-cell genome](scg_part6) 

<div>
 <span style="float:left"><a class="btn btn-primary" href="../schedule"> Previous page</a></span>
 <span style="float:right"><a class="btn btn-primary" href="connectToUppmax"> Next page</a></span>
</div> 
