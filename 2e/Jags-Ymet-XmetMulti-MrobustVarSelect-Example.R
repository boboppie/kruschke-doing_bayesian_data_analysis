# Example for Jags-Ymet-XmetMulti-MrobustVarSelect.R 
#------------------------------------------------------------------------------- 
# Optional generic preliminaries:
graphics.off() # This closes all of R's graphics windows.
rm(list=ls())  # Careful! This clears all of R's memory!
#------------------------------------------------------------------------------- 
# Load data file and specity column names of x (predictors) and y (predicted):
myData = read.csv( file="Guber1999data.csv" )

# UNCOMMENT ONE OF THE FOLLOWING SECTIONS (In RStudio, select and ctrl-shift-C)
#.............................................................................
# Two predictors:
# myData = read.csv( file="Guber1999data.csv" )
# yName = "SATT" ; xName = c("Spend","PrcntTake")
# fileNameRoot = "Guber1999data-Jags-VarSelect-" 
# numSavedSteps=15000 ; thinSteps=25
#.............................................................................
# # Two predictors with redundant predictor:
# myData = read.csv( file="Guber1999data.csv" )
# PropNotTake = (100-myData[,"PrcntTake"])/100
# myData = cbind( myData , PropNotTake )
# yName = "SATT" ; xName = c("Spend","PrcntTake","PropNotTake")
# fileNameRoot = "Guber1999data-Jags-Redund-VarSelect-" 
# numSavedSteps=15000 ; thinSteps=30
#.............................................................................
# # Two predictors with two redundant predictors:
# myData = read.csv( file="Guber1999data.csv" )
# PropNotTake = (100-myData[,"PrcntTake"])/100
# Partic = myData[,"PrcntTake"]/10 
# myData = cbind( myData , PropNotTake , Partic )
# yName = "SATT" ; xName = c("Spend","PrcntTake","PropNotTake","Partic")
# fileNameRoot = "Guber1999data-Jags-Redund2-VarSelect-" 
# numSavedSteps=15000 ; thinSteps=15
#.............................................................................
# # Four predictors:
myData = read.csv( file="Guber1999data.csv" )
yName = "SATT" ; xName = c("Spend","PrcntTake","StuTeaRat","Salary")
# Specify name of folder (i.e., directory) in which to save output files:
dirName="Guber1999data-Jags-4X-VarSelect"
# Create the folder:
if(!dir.exists(dirName)){dir.create(dirName)}
# Specify prefix (i.e., filename root) for names of saved output files:
fileNameRoot = paste0(dirName,"/Guber1999data-Jags-4X-VarSelect-")
numSavedSteps=15000 ; thinSteps=20
#.............................................................................
# # Append columns of random predictors:
# myData = read.csv( file="Guber1999data.csv" )
# standardize = function(v){(v-mean(v))/sd(v)}
# Ny=nrow(myData)
# NxRand = 12
# set.seed(47405)
# for ( xIdx in 1:NxRand ) {
#   xRand = standardize(rnorm(Ny))
#   myData = cbind( myData , xRand )
#   colnames(myData)[ncol(myData)] = paste0("xRand",xIdx)
# }
# yName = "SATT" ; xName = c("Spend","PrcntTake", paste0("xRand",1:NxRand) )
# fileNameRoot = "Guber1999data-Jags-RandX-VarSelect-" 
# numSavedSteps=15000 ; thinSteps=5
#.............................................................................
# # Two predictors with interaction.
# myData = read.csv( file="Guber1999data.csv" )
# # Append named column with interaction product:
# myData = cbind( myData , SpendXPrcnt=myData[,"Spend"]*myData[,"PrcntTake"] )
# yName = "SATT" ; xName = c("Spend","PrcntTake","SpendXPrcnt")
# fileNameRoot = "Guber1999data-Jags-Inter-VarSelect-" 
# numSavedSteps=15000 ; thinSteps=25
#.............................................................................
graphFileType = "png" 
#------------------------------------------------------------------------------- 
# Load the relevant model into R's working memory:
source("Jags-Ymet-XmetMulti-MrobustVarSelect.R")
#------------------------------------------------------------------------------- 
# Generate the MCMC chain:
#startTime = proc.time()
mcmcCoda = genMCMC( data=myData , xName=xName , yName=yName , 
                    numSavedSteps=numSavedSteps , thinSteps=thinSteps , 
                    saveName=fileNameRoot )
#stopTime = proc.time()
#duration = stopTime - startTime
#show(duration)
#------------------------------------------------------------------------------- 
# Display diagnostics of chain, for specified parameters:
parameterNames = varnames(mcmcCoda) # get all parameter names
for ( parName in parameterNames ) {
  diagMCMC( codaObject=mcmcCoda , parName=parName , 
            saveName=fileNameRoot , saveType=graphFileType )
}

#------------------------------------------------------------------------------- 
# Get summary statistics of chain:
summaryInfo = smryMCMC( mcmcCoda , 
                        saveName=fileNameRoot )
show(summaryInfo)
# Display posterior information:
plotMCMC( mcmcCoda , data=myData , xName=xName , yName=yName , 
          pairsPlot=TRUE , showCurve=FALSE ,
          saveName=fileNameRoot , saveType=graphFileType )
#------------------------------------------------------------------------------- 
