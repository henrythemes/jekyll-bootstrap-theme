---
layout: default
title:  'Part 5: Exploring your single cell genome'
---

# Part 5: Exploring your single cell genome

## 5.2 Functional analysis of predicted protein content in MEGAN

In addition to being a great tool to detect contaminants, you can also use MEGAN to get a rough estimate of which proteins are encoded in your SAG. 
For this, MEGAN utilizes the KEGG (Kyoto Encyclopedia of Genes and Genomes) Pathways database in combination with a BLASTp file in which 
you searched NCBI's 'nr' database for protein sequences similar to the predicted protein sequences in your SAG (see 3c). 
MEGAN inspects the BLASTp outfile and 'maps' all BLASTp-hits to all KEGG Pathways. 
This allows you to get a clear overview of which biological pathways and processes are present in your SAG. 

Note however that the KEGG Pathways database does not contain all proteins and 
several proteins you might be interested in could be missed by this analysis.  
Execute the BLASTp search:  

```sh
blastp -query contigs.prodigal.faa -db /proj/g2015028/nobackup/single_cell_exercises/databases/db -evalue 1e-5 -num_threads 8 -out contigs.prodigal_vs_db.blastp
```

The database that you are searching now is a modified version of NCBI's 'nr'. 
We have severely reduced the size to reduce the search time to 10-15 minutes, while still finding all relevant proteins. 
Open MEGAN and import the BLASTp file as before but this time you need to select the Content tab and check Analyse KEGG content, before hitting Apply. 
Name your MEGAN file *'contigs.prodigal_vs_db.rma'*.
Explore the taxonomic classification using functions described in *5.1*.

### Questions

**Q5.6** Can you find the taxa (or taxon) that you identified as contamination in 5.1?  
Why (not)?  

**Q5.7** Apparently the fraction of ‘No hits’ is smaller in the BLASTp search compared to the BLASTn search. Can you think of a reason that explains this?  
Look for the KEGG button in the toolbar and hit it. A new window will pop up. Like with the taxonomy window, the KEGG Pathways are organized in a tree-like structure, making it easy to browse through all the different pathways.
For each pathway there is a corresponding graphical overview. Double click on the name of a pathway to see the graphical overview. Basically, each overview consists of the illustrated pathway and a list of proteins that make up this pathway. In these overviews it is possible to see which particular proteins of a given pathway are present in your SAG and which ones are absent. If present, the box that contains the name of protein is marked with green. The amount of proteins present in your SAG that ‘map’ to a given protein determine the intensity of the green color.  
Explore the KEGG Pathways.  

**Q5.8** Which RNA Polymerase subunits can you find?  

**Q5.9** Does this SAG encode genes for flagellar assembly? If so, which ones?   

**Q5.10** How many genes encoding for ribosomal proteins are present in the SAG?  

**Q5.11** Look at metabolism pathways. Would you say this bacterium is able to grow aerobically?  

**Q5.12** What other potentially interesting proteins can you find?  
Explain for each one why do you think they are interesting.  

**Q5.13** Based on what you find in the KEGG window in MEGAN, what can you say about the lifestyle of the bacterium from which this SAG was retrieved?  
For example, would you say it is motile, does it grow aerobically or anaerobically, is it free-living or symbiotic?  
Try to say as much as possible.  

**Q5.14** Given the estimated completeness of the genome (see 3f), how confident are you about statements regarding the nature of this bacterium?  
Close MEGAN and return to the terminal.  


<div>
 <span style="float:left"><a class="btn btn-primary" href="scg_part5_1"> Previous page</a></span>
 <span style="float:right"><a class="btn btn-primary" href="scg_part6"> Next page</a></span>
</div>

