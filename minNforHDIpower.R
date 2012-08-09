minNforHDIpower = function( genPriorMean , genPriorN ,
                     HDImaxwid=NULL , nullVal=NULL , ROPE=c(nullVal,nullVal) ,
                     desiredPower=0.8 , audPriorMean=0.5 , audPriorN=2 ,
                     HDImass=0.95 , initSampSize=20 , verbose=T ) {
   if ( !xor( is.null(HDImaxwid) , is.null(nullVal) ) ) {
      stop("One and only one of HDImaxwid and nullVal must be specified.")
   }
   source("HDIofICDF.R")
   # Convert prior mean and N to a,b parameter values of beta distribution.
   genPriorA = genPriorMean * genPriorN
   genPriorB = ( 1.0 - genPriorMean ) * genPriorN
   audPriorA = audPriorMean * audPriorN
   audPriorB = ( 1.0 - audPriorMean ) * audPriorN
   # Initialize loop for incrementing sampleSize
   sampleSize = initSampSize
   notPowerfulEnough = TRUE
   # Increment sampleSize until desired power is achieved.
   while( notPowerfulEnough ) {
      zvec = 0:sampleSize # All possible z values for N flips.
      # Compute probability of each z value for data-generating prior.
      pzvec = exp( lchoose( sampleSize , zvec )
                   + lbeta( zvec + genPriorA , sampleSize-zvec + genPriorB )
                   - lbeta( genPriorA , genPriorB ) )
      # For each z value, compute HDI. hdiMat is min, max of HDI for each z.
      hdiMat = matrix( 0 , nrow=length(zvec) , ncol=2 )
      for ( zIdx in 1:length(zvec) ) {
         z = zvec[zIdx]
         hdiMat[zIdx,] = HDIofICDF( qbeta , 
                                    shape1 = z + audPriorA ,
                                    shape2 = sampleSize - z + audPriorB ,
                                    credMass = HDImass )
      }
      hdiWid = hdiMat[,2] - hdiMat[,1]
      if ( !is.null( HDImaxwid ) ) {
         powerHDI = sum( pzvec[ hdiWid < HDImaxwid ] )
      }
      if ( !is.null( nullVal ) ) {
         powerHDI = sum( pzvec[ hdiMat[,1] > ROPE[2] | hdiMat[,2] < ROPE[1] ] )
      }
      if ( verbose ) {
         cat( " For sample size = ", sampleSize , ", power = " , powerHDI ,
              "\n" , sep="" ) ; flush.console()
      }
      if ( powerHDI > desiredPower ) {
         notPowerfulEnough = FALSE 
      } else {
         sampleSize = sampleSize + 1
      }
   } # End while( notPowerfulEnough )    
   # Return the minimal sample size that achieved the desired power.
   return( sampleSize )
} # end of function