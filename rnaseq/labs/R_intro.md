---
layout: default
title:  'R-intro'
---

# Very quick introduction to R

In several exercises in this course, we use the statistical
programming environment R. Here follows a very quick introduction,
primarily for those who have never used R before. You can go through
this now, or return to it when you first encounter R in the exercises.

First start R on your computer. How to do this depends on your
operating system. If you are on a Linux or Mac OS X system, you would
typically execute the command ``R`` at the shell prompt.

In this course, we are running R on a Linux server. If you would like
to install R on your own computer, you can find the appropriate
download for your system [here](http://ftp.sunet.se/pub/lang/CRAN/).

When R starts, you should see a text message stating the version of
R you are running and some further information, followed by a
command-line “prompt” ( ``>``, a greater-than sign). The prompt means that
R is waiting for you to type a command. To get a feel for how this
works, try out some arithmetic::

	 > 2 + 3
	 [1] 5
	 > 5 * 4 + 10
	 [1] 30
And some function calls

	> abs(-5)
	[1] 5
	> sum(1,5,10)
	[1] 16

To find out how to use a function, type its name preceded by a question mark::

	> ? sin

This will bring up some help documentation for the function. You can
use the arrow keys to scroll the help text up and down. Press q to get
back to the R prompt.

Now try this:
	
	> a <- 10

This command created an object called a. Objects are an important
concept in R (as in many other programming languages), and we will be
creating more of them in the RNA-seq exercises. We can inspect an
object by just typing its name:
	
	> a
	[1] 10

We can also change the value of an object that we’ve created:

	> a <- 2 * a	
	> a
	[1] 20

The object *a* created above is a vector with a single element. To
create a vector with several elements, you can use the function *c*::


	> b <- c(1, 2, 10)
	> b
	[1]  1  2 10

Or the colon operator::

	> 1:10
	[1]  1  2  3  4  5  6  7  8  9 10

A matrix can be created with the function *cbind*::

	> b <- cbind(1:10, 101:110)
	> b
       [,1] [,2]
	[1,]    1  101
	[2,]    2  102
	[3,]    3  103
	[4,]    4  104
	[5,]    5  105
	[6,]    6  106
	[7,]    7  107
	[8,]    8  108
	[9,]    9  109
	[10,]   10  110

We can then use indices to access selected elements of the matrix:

	> b[1,]
	[1]   1 101
	> b[, 2]
	[1] 101 102 103 104 105 106 107 108 109 110
	> b[c(5,8), 2]
	[1] 105 108

You can find manuals for R and more information on the [R web site](http://www.r-project.org/)

# Installing packages in R and loading them

## Installing packages from BioConductor
To install packages from BioConductor you first have to load a function biocLite.R from bioconductor and then use that function to install the package. That is done in two lines of codes

	## try http if https is not available
	source("https://bioconductor.org/biocLite.R")
	biocLite("Name of your package")
	
	##Example installing DEseq2
	
	## try http if https is not available
	source("https://bioconductor.org/biocLite.R")
	biocLite("DESeq2")
	 
## Installing packages from CRAN

To install a CRAN package in R, use the install.packages() function. This simple command downloads the package from a specified repository (by default, CRAN) and installs it on your machine:

	install.packages("Name of your CRAN package")
	
	##An example
	
	install.packages("ggplot2")

	
## Installing packages from github

Sometimes you want to install packages that are still under development. Those are mostly found on github. This example show how you can switch between a stable version and the development version. If you only want to install from github ignore the dev_mode() part

	install.packages("devtools")

	library(devtools)

	dev_mode(on=T)

	install_github("hadley/ggplot2")

	# use dev ggplot2 now

	# when finished do:

	dev_mode(on=F)  #and you are back to having stable ggplot2

## Load package/library

Use library() to load the package after installing it

	#Example of loading DEseq and ggplot2
	library(DEseq)
	library(ggplot2)

# Quit r
 
	q() or quit()	 

	 