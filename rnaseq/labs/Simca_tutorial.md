---
layout: default
title:  'Multivariate analysis'
---



# Multivariate analysis cookbook


## Software used in this exercise is Simca from Umetrics. 


## About the data set for MVA excercise
Thoracic aortic aneurysm (TAA) is a pathological widening of the aorta resulting from degeneration of the extracellular matrix and loss of smooth muscle cells in the tunica media. TAA is an asymptomatic disease before the actual rupture of the aorta, a condition that is lethal if not treated on time. There are several different etiologies of TAA that predispose individuals to TAA, the most common one being aneurysm associated with bicuspid aortic valve (BAV), and idopathic causes of TAA. Patients with both bicuspid and tricuspid aortic valves are prone to get dilatation, however. The dataset presented in this exercise is probe set expression levels in TGFβ signaling pathway (Affymetrix Exon Array platform) in patients that have dilated and non dilated aorta tissues with BAV or TAV. The question we wanted to find out is if there are different mechanisms of dilatation in patients with TAV compared to those with BAV. The results are published in the following reference:

**Kurtovic**, Paloschi, Folkersen, Gottfries, Franco-Cereceda, and Eriksson (2011) **"Diverging alternative splicing fingerprints in TGFβ signaling pathway identified in thoracic aortic aneurysms"**; Molecular Medicine, 17; 665-675] that you can find [here](http://www.ncbi.nlm.nih.gov/pubmed/21448509)

This is a cookbook for handling multivariate analysis using SIMCA. 


By following the steps in this cookbook with the excel sheet downloaded [here](https://export.uppmax.uu.se/b2013006/downloads/courses/RNAseqWorkshop/multivariateAnalysis/Multivariate_excercise_data_set_151020.xlsb)

## Data import
•	Open Simca  
•	File -> New -> Regular project -> Multivariate_excercise_data_set_151020 -> click on the second column (“CvsD”) and choose as Y variable -> Finish -> Save  


## Principal component analysis

### TAV

•	Edit M1 -> Variables -> CvsD -> Exclude -> Observations -> Find (write CB) -> Exclude -> Find (write DB) -> Exclude -> OK -> Autofit -> Scores -> Right click Score plot -> Properties -> Color -> Identifiers -> Length: 2  

### BAV

•	New -> As M1 -> Variables -> CvsD -> Exclude -> Observations -> ctrl A -> Include-> Find (write CT) -> Exclude -> Find (write DT) -> Exclude -> OK -> Autofit -> Scores -> Right click Score plot -> Properties -> Color -> Identifiers -> Length: 2

## Orthogonal partial least squares to latent structures discriminant analysis (OPLS)

### TAV

•	New -> As M1 -> Variables -> CvsD -> Y -> Model type: OPLS -> OK -> Autofit -> Scores -> Loadings -> Pred X-Y -> Right click loading plot -> Sort ascending -> Right click loading plot -> Create list 


### BAV

•	New -> As M3 -> Variables -> CvsD -> Y -> Model type: OPLS -> OK -> Autofit -> Scores -> Loadings -> Pred X-Y -> Right click loading plot -> Sort ascending -> Right click loading plot -> Create list   

