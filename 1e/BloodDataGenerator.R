# Fictitious blood data. The correlations, means, and SDs used here are
# fabricated for pedagogical purposes and may have no resemblance to real data.
# Kruschke, J. K. (2010). Doing Bayesian data analysis:
# A Tutorial with R and BUGS. Academic Press / Elsevier Science.

# Specify the names of the predictors:
xNames = c("Systolic","Diastolic","Weight","Cholesterol","Height","Age")
nX = length(xNames) # number of predictors
# SPECIFY THE CORRELATIONS BETWEEN PREDICTORS:
if ( T ) { # zero correlations everywhere
rMat = matrix( c(   1 ,  0 ,  0 ,  0 ,  0 ,  0 ,
                    0 ,  1 ,  0 ,  0 ,  0 ,  0 ,
                    0 ,  0 ,  1 ,  0 ,  0 ,  0 ,
                    0 ,  0 ,  0 ,  1 ,  0 ,  0 ,
                    0 ,  0 ,  0 ,  0 ,  1 ,  0 ,
                    0 ,  0 ,  0 ,  0 ,  0 ,  1  ) , ncol=nX ) }
if ( F ) { # first two predictors strongly correlated
rMat = matrix( c(   1 , .95,  0 ,  0 ,  0 ,  0 ,
                   .95,  1 ,  0 ,  0 ,  0 ,  0 ,
                    0 ,  0 ,  1 ,  0 ,  0 ,  0 ,
                    0 ,  0 ,  0 ,  1 ,  0 ,  0 ,
                    0 ,  0 ,  0 ,  0 ,  1 ,  0 ,
                    0 ,  0 ,  0 ,  0 ,  0 ,  1  ) , ncol=nX ) }
if ( F ) { # first two uncorrelated, but other predictors correlated
rMat = matrix( c(  1  ,  0  , .6 , .2 , .1 , .1 ,
                   0  ,  1  , .6 , .2 , .1 , .1 ,
                  .6  , .6  ,  1 , .4 , .2 , .2 ,
                  .2  , .2  , .4 ,  1 ,  0 , .3 ,
                  .1  , .1  , .2 ,  0 ,  1 ,  0 ,
                  .1  , .1  , .2 , .3 ,  0 ,  1 ) , ncol=nX ) }
if ( F ) { # first two strongly correlated with others also correlated
rMat = matrix( c(  1  , .95 , .6 , .2 , .1 , .1 ,
                  .95 ,  1  , .6 , .2 , .1 , .1 ,
                  .6  , .6  ,  1 , .4 , .2 , .2 ,
                  .2  , .2  , .4 ,  1 ,  0 , .3 ,
                  .1  , .1  , .2 ,  0 ,  1 ,  0 ,
                  .1  , .1  , .2 , .3 ,  0 ,  1 ) , ncol=nX ) }
# SPECIFY THE NUMBER OF DATA POINTS:
nSubj = 200
mVec = rep(0,nX) # means of predictors
require(MASS) # package needed for mvrnorm() function in next line
set.seed(47405)
xMat = mvrnorm( n=nSubj , mu=mVec , Sigma=rMat ) #
# SPECIFY THE REGRESSION COEFFICIENTS:
betaVec = c( 0 , 2 , 2 , 1 , 0 , 0.5 )
# SPECIFY THE PROPORTION OF PREDICTED VALUES THAT ARE 1's, WHICH IN TURN
# DETERMINES THE THRESHOLD (i.e., negative intercept). THIS IS ACCURATE ONLY
# FOR LARGE REGRESSION COEFFICIENTS; OTHERWISE SPECIFY MORE EXTREME PROPORTION:
proportionOnes = 0.5  # e.g., about .05 for actual .10
heartAttackLinear = xMat %*% betaVec
threshold = quantile( heartAttackLinear , 1-proportionOnes )
heartAttackProb = 1 / ( 1 + exp( -1 * ( heartAttackLinear - threshold ) ) )
y = 0*heartAttackProb
for ( sIdx in 1:nSubj ) {
    y[sIdx] = sample( x=c(0,1) , size=1 , 
                   prob=c( 1-heartAttackProb[sIdx] , heartAttackProb[sIdx] ) )
}
cat("Generated proportion of 1's in data: ",mean(y),"\n")
# Convert to "real world" scale values (multiply by SD, add mean):
xMat[,1] = xMat[,1] * 17 + 125 # systolic
xMat[,2] = xMat[,2] * 11 + 80  # diastolic
xMat[,3] = xMat[,3] * 30 + 150 # weight
xMat[,4] = xMat[,4] * 30 + 130 # cholest
xMat[,5] = xMat[,5] * 3 + 65 # height
xMat[,6] = xMat[,6] * 15 + 50 # age
xMat = round( xMat )
# Assemble the values into a matrix:
dataMat = cbind( y , xMat )
colnames(dataMat) = c( "HeartAttack" , xNames )
# Write the matrix to a table to be loaded by other programs:
write.table( dataMat , file="BloodDataGeneratorOutput.txt" , row.names=F , col.names=T )
