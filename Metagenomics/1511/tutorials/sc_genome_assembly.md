---
layout: default
title:  'Single Cell Genomics Tutorial'
---

# Introduction: Single Cell Genomics Tutorial

Single cell genomics is an emerging technology that allows one to explore the genome sequence of individual cells. 
During this tutorial you will work with real single cell genome data from a single cell that was isolated from a hot spring in Yellowstone National Park (USA). 
The data you will work with is part of a larger project ('PUZZLE_CELL') that aims to identify and genomically probe novel prokaryotic lineages, and to gain insight in the origin and evolution of life. 
The data that you will work with is paired-end Illumina reads. 
We have chosen to have you work with both HiSeq (2x100 bp - dataset 1) and MiSeq (2x250 bp - dataset 2) datasets. 
Both datasets were generated from the same single cell, hence allowing you to develop a feeling for what you might want to used in any potential future SCG project. 
In addition, there is a third MiSeq dataset (2x250 bp - dataset 3) that contains a completely new organism, which might be a bit more challenging to work with. 
You can chose to work with the latter dataset if you happen to have some extra time at the end of the tutorial.

## Overview of steps in this exercise

### Exercises are split into 2 sessions:

**Session 1**

1. [Connecting to UPPMAX](connectToUppmax)
2. [Familiarizing with data](scg_part2)
3. [Single-cell genome assembly using SPAdes (HiSeq data)](scg_part3)

**Session 2**

4. [Assessing read coverage and chimera checking](scg_part4)
5. [Checking for contaminants](scg_part5)
6. [Single-cell genome assembly using SPAdes (MiSeq data) and comparison between HiSeq vs. MiSeq data](scg_part6)
7. [Analysis of a novel single-cell genome (bonus exercise)](scg_part7)