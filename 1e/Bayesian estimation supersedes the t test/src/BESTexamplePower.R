# Version of May 26, 2012.
# John K. Kruschke  
# johnkruschke@gmail.com
# http://www.indiana.edu/~kruschke/BEST/
#
# This program is believed to be free of errors, but it comes with no guarantee! 
# The user bears all responsibility for interpreting the results. 
# Please check the webpage above for updates or corrections.
#
### ***************************************************************
### ******** SEE FILE BESTexample.R FOR INSTRUCTIONS **************
### ***************************************************************

# OPTIONAL: Clear R's memory and graphics:
rm(list=ls())  # Careful! This clears all of R's memory!
graphics.off() # This closes all of R's graphics windows.

# Get the functions loaded into R's working memory:
source("BEST.R")

#-------------------------------------------------------------------------------
# RETROSPECTIVE POWER ANALYSIS.
# !! This section assumes you have already run BESTexample.R !!
# Re-load the saved data and MCMC chain from the previously conducted 
# Bayesian analysis. This re-loads the variables y1, y2, mcmcChain, etc.
load( "BESTexampleMCMC.Rdata" )
power = BESTpower( mcmcChain , N1=length(y1) , N2=length(y2) , 
                  ROPEm=c(-0.1,0.1) , ROPEsd=c(-0.1,0.1) , ROPEeff=c(-0.1,0.1) , 
                  maxHDIWm=2.0 , maxHDIWsd=2.0 , maxHDIWeff=0.2 , nRep=1000 ,
                  mcmcLength=10000 , saveName = "BESTexampleRetroPower.Rdata" ) 

#-------------------------------------------------------------------------------
# PROSPECTIVE POWER ANALYSIS, using fictitious strong data.
# Generate large fictitious data set that expresses hypothesis:
prospectData = makeData( mu1=108, sd1=17, mu2=100, sd2=15, nPerGrp=1000, 
                         pcntOut=10, sdOutMult=2.0, rnd.seed=NULL )
y1pro = prospectData$y1 # Merely renames simulated data for convenience below.
y2pro = prospectData$y2 # Merely renames simulated data for convenience below.
# Generate Bayesian posterior distribution from fictitious data:
# (uses fewer than usual MCMC steps because it only needs nRep credible 
# parameter combinations, not a high-resolution representation)
mcmcChainPro = BESTmcmc( y1pro , y2pro , numSavedSteps=2000 )  
postInfoPro = BESTplot( y1pro , y2pro , mcmcChainPro , pairsPlot=TRUE )  
save( y1pro, y2pro, mcmcChainPro, postInfoPro, 
      file="BESTexampleProPowerMCMC.Rdata" )
# Now compute the prospective power for planned sample sizes:
N1plan = N2plan = 50 # specify planned sample size
powerPro = BESTpower( mcmcChainPro , N1=N1plan , N2=N2plan , showFirstNrep=5 ,
                   ROPEm=c(-1.5,1.5) , ROPEsd=c(-0.0,0.0) , ROPEeff=c(-0.0,0.0) , 
                   maxHDIWm=15.0 , maxHDIWsd=10.0 , maxHDIWeff=1.0 , nRep=1000 ,
                   mcmcLength=10000 , saveName = "BESTexampleProPower.Rdata" ) 

#-------------------------------------------------------------------------------
