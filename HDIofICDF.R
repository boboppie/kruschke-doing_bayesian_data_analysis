HDIofICDF = function( ICDFname , credMass=0.95 , tol=1e-8 , ... ) {
   # Arguments:
   #   ICDFname is R's name for the inverse cumulative density function
   #     of the distribution.
   #   credMass is the desired mass of the HDI region.
   #   tol is passed to R's optimize function.
   # Return value:
   #   Highest density iterval (HDI) limits in a vector.
   # Example of use: For determining HDI of a beta(30,12) distribution, type
   #   HDIofICDF( qbeta , shape1 = 30 , shape2 = 12 )
   #   Notice that the parameters of the ICDFname must be explicitly named;
   #   e.g., HDIofICDF( qbeta , 30 , 12 ) does not work.
   # Adapted and corrected from Greg Snow's TeachingDemos package.
   incredMass =  1.0 - credMass
   intervalWidth = function( lowTailPr , ICDFname , credMass , ... ) {
      ICDFname( credMass + lowTailPr , ... ) - ICDFname( lowTailPr , ... )
   }
   optInfo = optimize( intervalWidth , c( 0 , incredMass ) , ICDFname=ICDFname ,
                        credMass=credMass , tol=tol , ... )
   HDIlowTailPr = optInfo$minimum
   return( c( ICDFname( HDIlowTailPr , ... ) ,
	            ICDFname( credMass + HDIlowTailPr , ... ) ) )
}