---
layout: default
title:  'PCA and clustering'
---

# PCA and clustering on a single cell RNA-seq dataset

Here are some examples on how to run PCA/Clustering on a single cell RNA-seq dataset. These methods can also be applied to any other type of dataset, such as RNA-seq or other high throuput data. 
The dataset used is single-cell RNA-seq data from mouse embryonic development from Deng. et al. Science 2014, Vol. 343 no. 6167 pp. 193-196, "Single-Cell RNA-Seq Reveals Dynamic, Random Monoallelic Gene Expression in Mammalian Cells"

All data you need is available in the folder: 
/proj/b2013006/webexport/downloads/courses/RNAseqWorkshop/pca_clust

Copy all that data to your folder of choice (on uppmax or your own computer) and start R in that folder. 

You will need some R packages, gplots and plotrix, these can be installed with the command:
    	install.packages(c("gplots","plotrix")) 

All the commands that are run in this example can also be found in the file: run_PCA_clust.R


## Data processing

First read in the data and define colors/symbols for plotting

	DATA<-read.table("rpkms_Deng2014_preinplantation.txt")
	nS<-ncol(DATA) #number of samples
	nG<-nrow(DATA) #number of genes
	
	
We want to create a vector of colors for each embryonic stage, but also symbols for each individual embryo. To get the stage and embry definitions for each sample we split the names either by "_" or "."
	
	sample.names<-colnames(DATA)
	stage<-unlist(lapply(strsplit(sample.names,"_"),function(x) x[1]))
	embryo<-unlist(lapply(strsplit(sample.names,"\\."),function(x) x[1]))
	
	stages<-unique(stage)
	# have 11 different stages, define 11 colors for those.
	coldef.stage<-c("black","red","green","blue","cyan","magenta","yellow","pink","gray","brown","orange")
	col.stage<-mat.or.vec(1,nS)
	for (i in 1:length(stages)){
	    idx<-which(stage==stages[i])
	    col.stage[idx]<-coldef.stage[i]
	}
	
	embryos<-unique(embryo)
	# have 42 different embryos, use the first default 18 symbols rotated and we will never have 
	# a combination of the same color/symbol
	pchdef.embryo<-c(1:18,1:18,1:6)
	pch.embryo<-mat.or.vec(1,nS)
	for (i in 1:length(embryos)){
	    idx<-which(embryo==embryos[i])
	    pch.embryo[idx]<-pchdef.embryo[i]
	}

## PCA

There are some custom functions in PCA_RNAseq_functions.R that are called run.pca, pca.plot, pca.contribution.plot and pca.loadings.plot. Have a look at the file for documentation of the scripts.

	# load the custom PCA functions
	source("PCA_RNAseq_functions.R")
	
	# first run pca, should take about 1 min.
	PC<-run.pca(DATA)
	# should give you a prcomp object (PC)
	
	# for plotting, open a PDF
	pdf("pca_all_genes.pdf")
	
	# now we can plot the first 2 PCs
	pca.plot(PC,col=col.stage,pch=pch.embryo,main="first PCA")
	legend("topleft",stages,col=coldef.stage,pch=16,cex=0.5,bty='n')
	
	# and plot with the first 5 PCs
	pca.plot(PC,col=col.stage,pch=pch.embryo,main="first PCA",selpc=1:5)
	
	# plot PC contribution
	pca.contribution.plot(PC)
	
	# and top gene loadings for the first 5 components
	pca.loadings.plot(PC)
	
	# close the PDF
	dev.off()
	
Have a look at the file you created, what type of variance does the different PCs capture? 
How many PCs are informative? Do the genes that contribute to each PC make sense? 

## Coloring in PCA

We can also add in coloring by any color we want, lets use the expression of the top genes for PC1 & PC2. Here I have used the function color.scale from the "plotrix" package to define a color range with green-yellow-red scale.

	# load the library
	library(plotrix)
	
	pdf("pca_top_loading_genes_genes.pdf")
	# define plotting of 4 plots (2 rows, 2 columns)
	par(mfrow=c(2,2))
	n1<-names(sort(PC$rotation[,1]))
	n2<-names(sort(PC$rotation[,2]))
	# get first/last gene from PC1 & 2
	plotgenes<-c(tail(n1,1),n1[1],tail(n2,1),n2[1])
	for (gene in plotgenes) {
	    idx<-match(gene, rownames(DATA))
	    expr<-log2(as.numeric(DATA[idx,]+1))
	    # color scale with red for high values, yellow intermeidate,  green for low.
	    col<-color.scale(expr,c(0,1,1),c(1,1,0),0)
	    pca.plot(PC,col=col,pch=16,main=gene)
	}
	dev.off()
	
## Different settings in PCA

The PCA that was run automatically transforms the rpkm-values to log-space and does centering of the data, test doing it without logging and witout centering and compare the results.
	
	PC.nolog<-run.pca(DATA,log.vals=FALSE)
	PC.nocenter<-run.pca(DATA,center=FALSE)
	PC.nocenterlog<-run.pca(DATA,center=FALSE,log.vals=FALSE)
	
	pdf("pca_comparisons.pdf")
	par(mfrow=c(2,2))
	pca.plot(PC,col=col.stage,pch=pch.embryo,main="first PCA")
	legend("topleft",stages,col=coldef.stage,pch=16,cex=0.3,bty='n')
	pca.plot(PC.nolog,col=col.stage,pch=pch.embryo,main="no logging")
	pca.plot(PC.nocenter,col=col.stage,pch=pch.embryo,main="no centering")
	pca.plot(PC.nocenterlog,col=col.stage,pch=pch.embryo,main="no centering, no logging")
	dev.off()
	
What are the main differences, why is that do you think? Why should the RPKM-values be logged? 

## PCA based on blastocyst stages only

The different embryonic stages separated out quite well in the first PCA, but at the blastocyst stage the cells do not separate by timepoint. Try running a PCA with only those cells and see if you can get them to separate.

	# get the index for all the blastocyst cells
	blasto<-grep("blast",colnames(DATA))
	PC.blast<-run.pca(DATA[,blasto])
	
	# run and plot PCA
	pdf("pca_blastocyst.pdf")
	pca.plot(PC.blast,col=col.stage[blasto],pch=pch.embryo[blasto],main="Blastocyst PCA")
	legend("topleft",stages,col=coldef.stage,pch=16,cex=0.5,bty='n')
	pca.plot(PC.blast,col=col.stage[blasto],pch=pch.embryo[blasto],main="Blastocyst PCA",selpc=1:5)
	pca.contribution.plot(PC.blast)
	pca.loadings.plot(PC)
	dev.off()
	
Have a look at the PCA with blastocyst cells, do you see clear separation of the timepoints at any of the PCs?


## Clustering

Now, lets try some different clustering methods. Quite often, clustering is based on pairwise correlations. So let's start with calculating pairwise correlations for all samples. 

OBS! it takes a long time to run if your computer is slow, so there is a file prepared for you to read in if it is taking too long

     	   # OBS! This shows how the correlations were calculated, you may instead load the file
	   C<-mat.or.vec(nS,nS)
	   for (i in 1:nS) {
	       for (j in 1:nS){
	           if (i==j){ C[i,j]<-NA }
	           else {  C[i,j] = cor(log2(DATA[,i]+1),log2(DATA[,j]+1),method="pearson") }
	       }
	   }
	   colnames(C)<-colnames(DATA)
	   rownames(C)<-colnames(DATA)
	   write.table(C,file="pairwise_pearson_correlations.txt")

Or load the file directly:

	   C<-read.table("pairwise_pearson_correlations.txt")
	   C<-as.matrix(C)
	

Run clustering based on the correlations

	dist.corr<-as.dist(1-C) 
	hcl.corr<-hclust(dist.corr,method="ward.D2")
	
For comparison a test with a different clustering method, average linkage:

	hcl.corr2<-hclust(dist.corr,method="average")
	
Another option is to do clustering based on euklidean distances:

	dist.euk<-dist(log2(t(DATA)+1))
	hcl.euk<-hclust(dist.euk,method="ward.D2")
	
Lets plot a heatmap with the correlations and the results from the different clustering methods:
	
	library(gplots)
	pdf("heatmap_pairwise_correlations.pdf")
	par(xpd=T)
	heatmap.2(C,ColSideColors=col.stage,RowSideColors=col.stage,Colv=as.dendrogram(hcl.corr),Rowv=as.dendrogram(hcl.corr),scale="none",trace="none",main="correlation, Ward")
	legend("topright",stages,fill=coldef.stage,cex=0.5,bty='n',inset=c(0,-0.15,0,0))
	heatmap.2(C,ColSideColors=col.stage,RowSideColors=col.stage,Colv=as.dendrogram(hcl.corr2),Rowv=as.dendrogram(hcl.corr2),scale="none",trace="none",main="correlation, average")
	legend("topright",stages,fill=coldef.stage,cex=0.5,bty='n',inset=c(0,-0.15,0,0))
	heatmap.2(C,ColSideColors=col.stage,RowSideColors=col.stage,Colv=as.dendrogram(hcl.euk),Rowv=as.dendrogram(hcl.euk),scale="none",trace="none",main="euklidean distance, Ward")
	legend("topright",stages,fill=coldef.stage,cex=0.5,bty='n',inset=c(0,-0.15,0,0))
	dev.off()
	
Another common clustering method is K-means clustering, lets try that as well, with a few different settings for k:

	km7<-kmeans(log2(t(DATA)+1),7)
	km10<-kmeans(log2(t(DATA)+1),10)
	km15<-kmeans(log2(t(DATA)+1),15)
	
## Plot the clusters from hierarchical clustering in PCA-space

To get clusters from a hierarchical clustering we have to cut the branches of the dendrogram, this is done with the function "cutree", either with desired number of final clusters, or the height for cutting.

	
Split the different hclust objects into 7 clusters:

	clusters.corr<-cutree(hcl.corr,7)
	clusters.corr2<-cutree(hcl.corr2,7)
	clusters.euk<-cutree(hcl.euk,7)
	
Now, lets plot them onto PCA-space, with PC1+PC2:

	pdf("pca_clustering_methods.pdf")
	par(mfrow=c(3,2),mar=c(1,1,4,1))
	pca.plot(PC,col=clusters.corr,main="clusters from correlation,Ward",selpc=1:2)
	pca.plot(PC,col=clusters.corr2,main="clusters from correlation,average",selpc=1:2)
	pca.plot(PC,col=clusters.euk,main="clusters from euklidean dist,Ward",selpc=1:2)
	pca.plot(PC,col=km7$cluster,main="clusters from k-means, 7",selpc=1:2)
	pca.plot(PC,col=km10$cluster,main="clusters from k-means, 10",selpc=1:2)
	pca.plot(PC,col=km15$cluster,main="clusters from k-means, 15",selpc=1:2)
	dev.off()
	
	
Or we can plot the first 5 PCs, to check if the splitting of clusters is captured by any of the lower PCs

	pdf("pca_clustering_methods_5PC.pdf")
	pca.plot(PC,col=clusters.corr,main="clusters from correlation,Ward",selpc=1:5)
	pca.plot(PC,col=clusters.corr2,main="clusters from correlation,average",selpc=1:5)
	pca.plot(PC,col=clusters.euk,main="clusters from euklidean dist,Ward",selpc=1:5)
	pca.plot(PC,col=km7$cluster,main="clusters from k-means, 7",selpc=1:5)
	pca.plot(PC,col=km10$cluster,main="clusters from k-means, 10",selpc=1:5)
	pca.plot(PC,col=km15$cluster,main="clusters from k-means, 15",selpc=1:5)
	dev.off()
	
	
Using the hcl.corr object (based on pairwise correlations and Ward distance) how many clusters do you think are optimal? How many clusters do we need to cut the dendrogram into to separate out mid/late 2-cell? Does this clustering make sense? 
	
Test some different cutoffs with:

	pdf("pca_hcl_corr_n_test.pdf")
	par(mfrow=c(3,3),mar=c(1,1,4,1))
	pca.plot(PC,col=col.stage,pch=pch.embryo,main="first PCA")
	legend("topleft",stages,col=coldef.stage,pch=16,cex=0.3,bty='n')
	for (ncl in c(5,7,10,15,20,25,30,40)){
	    clusters.corr<-cutree(hcl.corr,ncl)
	    pca.plot(PC,col=clusters.corr,main=sprintf("%d clusters",ncl),selpc=1:2)
	}
	dev.off()
	



## PCA or MDS

Classical MDS (based on euklidean distances) should be identical to PCA, so let's run both and compare

	#classical MDS using the distance matrix we created before:
	fit.euk <- cmdscale(dist.euk, eig = TRUE, k = 2)
	# or based on the correlation distance matrix:
	fit.corr <- cmdscale(dist.corr, eig = TRUE, k = 2)
	
We can also run MDS using the limma function for MDS that has a different distance measure and a gene selection step first.
Description:
 This function is a variation on the usual multdimensional scaling (or principle coordinate) plot, in that a distance measure particularly appropriate for the microarray context is used. The distance between each pair of samples (columns) is the root-mean-square deviation (Euclidean distance) for the top ‘top’ genes. Distances on the plot can be interpreted as _leading log2-fold-change_, meaning the typical (root-mean-square) log2-fold-change between the samples for the genes that distinguish those samples.


	library(limma)
	mds.limma <- plotMDS(log2(DATA+1))
	

Now lets plot all of them and compare:
	
	pdf("mds_pca_comparison.pdf")
	par(mfrow=c(2,2),mar=c(1,1,4,1))
	pca.plot(PC,col=col.stage,pch=pch.embryo,main="first PCA")
	plot(fit.euk$points,col=col.stage,pch=pch.embryo,main="MDS euklidean")
	plot(fit.corr$points,col=col.stage,pch=pch.embryo,main="MDS correlation")
	plot(mds$cmdscale.out,col=col.stage,pch=pch.embryo,main="MDS limma")
	dev.off()
	
