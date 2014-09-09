# Version of January 19, 2013.
# John K. Kruschke  
# johnkruschke@gmail.com
# http://www.indiana.edu/~kruschke/BEST/
#
# This program is believed to be free of errors, but it comes with no guarantee! 
# The user bears all responsibility for interpreting the results. 
# Please check the webpage above for updates or corrections.
#
################################################################################
### To run this program, please prepare your computer as follows.
### 1. Install the general-purpose programming language R from  
###      http://www.r-project.org/
###    Install the version of R appropriate for your computer's operating
###    system (Windows, MacOS, or Linux).   
### 2. Install the Bayesian MCMC sampling program JAGS from
###      http://mcmc-jags.sourceforge.net/
###    Again, intall the version appropriate for your operating system.
### 3. Install the R editor, RStudio, from
###      http://rstudio.org/
###    This editor is not necessary, but highly recommended.
### 4. Invoke RStudio (or R by itself). At the command prompt, type
###      install.packages("rjags")
###    You will be asked to select an internet server; choose one. The
###    package rjags communicates with JAGS from R.
### 5. Make sure that the following programs are all
###    in the same folder as this file:
###      BESTexample.R (this file)
###      BEST.R
###      plotPost.R
###      HDIofMCMC.R
###      HDIofICDF.R
###      BESTexamplePower.R 
###      openGraphSaveGraph.R
### 6. Make sure that R's working directory is the folder in which those 
###    files reside. In RStudio, use menu tabs Tools -> Set Working Directory.
###    If working in R, use menu tabs File -> Change Dir.
### 7. After the above actions are accomplished, this program should
###    run as-is in R. You may "source" it to run the whole thing at once,
###    or run lines one at a time, in order.
### Problems? Please use the search term "installing" at  
### http://doingbayesiandataanalysis.blogspot.com/
################################################################################

# OPTIONAL: Clear R's memory and graphics:
rm(list=ls())  # Careful! This clears all of R's memory!
graphics.off() # This closes all of R's graphics windows.

# Get the functions loaded into R's working memory:
source("BEST.R")

# Specify data as vectors: 
# (Replace with your own data as needed. R can read many formats of data files,
# see the commands "scan" or "read.table" and its variants.)
y1 = c(101,100,102,104,102,97,105,105,98,101,100,123,105,103,100,95,102,106,
       109,102,82,102,100,102,102,101,102,102,103,103,97,97,103,101,97,104,
       96,103,124,101,101,100,101,101,104,100,101)
y2 = c(99,101,100,101,102,100,97,101,104,101,102,102,100,105,88,101,100,
       104,100,100,100,101,102,103,97,101,101,100,101,99,101,100,100,
       101,100,99,101,100,102,99,100,99)

# Run the Bayesian analysis:
mcmcChain = BESTmcmc( y1 , y2 ) 

# Plot the results of the Bayesian analysis:
postInfo = BESTplot( y1 , y2 , mcmcChain , pairsPlot=TRUE )
# Show detailed summary info on console:
show( postInfo ) 
# You can save the plot(s) using the pull-down menu in the R graphics window,
# or by using the following:
# saveGraph( file="BESTexample" , type="eps" )
# saveGraph( file="BESTexample" , type="jpeg" )

# Save the data and results for future use:
save( y1, y2, mcmcChain, postInfo, file="BESTexampleMCMC.Rdata" )
# To re-load the saved data and MCMC chain, type: load( "BESTexampleMCMC.Rdata" ) 

#-------------------------------------------------------------------------------
