HDIofGrid = function( probMassVec , credMass=0.95 ) {
   # Arguments:
   #   probMassVec is a vector of probability masses at each grid point.
   #   credMass is the desired mass of the HDI region.
   # Return value:
   #   A list with components:
   #   indices is a vector of indices that are in the HDI
   #   mass is the total mass of the included indices
   #   height is the smallest component probability mass in the HDI
   # Example of use: For determining HDI of a beta(30,12) distribution
   #   approximated on a grid:
   #   > probDensityVec = dbeta( seq(0,1,length=201) , 30 , 12 )
   #   > probMassVec = probDensityVec / sum( probDensityVec )
   #   > HDIinfo = HDIofGrid( probMassVec )
   #   > show( HDIinfo )
   sortedProbMass = sort( probMassVec , decreasing=T )
   HDIheightIdx = min( which( cumsum( sortedProbMass ) >= credMass ) )
   HDIheight = sortedProbMass[ HDIheightIdx ]
   HDImass = sum( probMassVec[ probMassVec >= HDIheight ] )
   return( list( indices = which( probMassVec >= HDIheight ) ,
                 mass = HDImass , height = HDIheight ) )
}