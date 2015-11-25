---
layout: default
title:  'Part 3: Single cell genome assembly using SPAdes'
---

# Part 3: Single cell genome assembly

<p>
In this part of the course, you will start doing assemblies of 'real' (but reduced) single cell genome datasets. We will compare two single cell specific assemblers, namely Spades and IDBA-UD, and one 'general-purpose' assembler called Ray (which were introduced by Kasia). The idea is that you will be able to compare the results of these different assemblers on two kinds of datasets (HiSeq and MiSeq), as well as different pre-treatments ('trimming'). You will also have a chance to explore how to decide which assembly is the best ('assembly metrics'), as there is no simple answer to this question.
Out of the total 12 assemblies we like you to compare, we suggest one person does only 3 on one of the datasets and pre-treatments. This way you can focus and skip handling too many folders and files. Assembly is also relatively time-consuming (although we have prepared reduced datasets for the tutorial to keep the times reasonable). So if you work as a group of 4 you can collect the results in one summary table we will look at. You will find a list of these tables below.
</p>

Actual tables to be filled in are provided in Google Docs and the links can be found below. 
You should talk to each other to form the groups and split the work. Do not worry if you miss something, we will collect result from all groups and discuss it together.  

[Group 1:](https://docs.google.com/spreadsheet/ccc?key=0AuNHyFPCsxthdGRKMXJwdF9jVDMzX2lGMkdJSDdOcnc&usp=sharing)  
[Group 2:](https://docs.google.com/spreadsheet/ccc?key=0AuNHyFPCsxthdC0tdzFySDFyaDNIdEh4M01xMXFQb3c&usp=sharing)  
[Group 3:](https://docs.google.com/spreadsheet/ccc?key=0AuNHyFPCsxthdHhwdUdUWUxBbnd0eC15WkJhS29iV3c&usp=sharing)  
[Group 4:](https://docs.google.com/spreadsheet/ccc?key=0AuNHyFPCsxthdFZRcXBjN0lrMGV5NmNuRnUzV2RkT0E&usp=sharing)  
[Group 5:](https://docs.google.com/spreadsheet/ccc?key=0AuNHyFPCsxthdEhkV3hIejJaMDZrWDFqd29XNTZFbmc&usp=sharing)  
[Group 6:](https://docs.google.com/spreadsheet/ccc?key=0AuNHyFPCsxthdFZDQzhLamJocHI0M0ZBQ0dMUDRrSFE&usp=sharing)  
[Group 7:](https://docs.google.com/spreadsheet/ccc?key=0AuNHyFPCsxthdDU2OUk4ank4c1A1VVhhbjZPaldtN2c&usp=sharing)  
[Group 8:](https://docs.google.com/spreadsheets/d/1Q3QBvPYzQ1kFHjWu0O5Jf1wgSH-9sxMrcndoLlNW92Y/edit?usp=sharing)  


3.1. [Organize working folder](scg_part3_1)  
3.2. [Pre-processing](scg_part3_2)  
3.3. [Assembly](scg_part3_3)  
3.4. [Assessing assembly quality using Quast](scg_part3_4)  
3.5. [Gene prediction using Prodigal](scg_part3_5)  
3.6. [Running completeness estimates](scg_part3_5)  
3.7. [Identifying ribosomal RNAs](scg_part3_7)  


These steps will help you obtain results to think about the following questions:


## Questions:


**Q3.1:** Did you notice how many reads were discarded in the pre-processing? Do the numbers differ between the Miseq and Hiseq datasets? You can use the group google doc to see the results for all the datasets.  

**Q3.2:** Did you notice any differences between HiSeq and MiSeq data and the different assemblers?  

**Q3.3:** Did you notice the impact of trimming, is it the same for all assemblers?  

**Q3.4:** What do you think is the best way to assess the 'quality' of an assembly? (e.g. total size, N50, number of predicted ORFs, completeness)  

**Q3.5:** What do you think is the best way to assemble this particular SC dataset? Why?  

**Q3.6:** What is the identity of the organism based on the analyses you have performed? What phylum does it belong to and is there any closely related organisms in the databases?  

**Q3.7:** Try to find out in what type of environment you might find similar organisms in.


<div>
 <span style="float:left"><a class="btn btn-primary" href="scg_part2"> Previous page</a></span>
 <span style="float:right"><a class="btn btn-primary" href="scg_part3_1"> Next page</a></span>
</div>


<!--- 
Notes for next year's workshop
- Uppmax support at the beginning of the workshop (troubleshooting, node reservation, technical support - cables etc) => not enough ethernet socket
- Having MEGAN and Artemis installed locally instead of on Milou (15-30min for sorting out computer problems/ software installation, unless we have dedicated computers/terminals). 
Students are to pre-install the tools on their computers if they will use their own laptops.
- table of contents
- useful Linux commands such as:
pwd to check current working directory
tab completion
arrow up for previous command
- define learning objectives
- include pre-course materials about fasta/q formats, basic unix commands (count number of reads in fasta or fastq files as checks).
- Longer introductory lecture on basics such as NGS data
- flowchart of steps involved in generating single cell genome data to help visualize how we got the data (overview of what part they are involved in)
- check points (completeness estimates - distribution of marker genes, chimera generation, MEGAN summary, Artemis)
- extra day for binning (together with metagenomic workshop)
- guest lectures (national/international) prior to hands-on workshops
--->
