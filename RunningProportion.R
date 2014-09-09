# Goal: Toss a coin N times and compute the running proportion of heads.
N = 500	# Specify the total number of flips, denoted N.
# Generate a random sample of N flips for a fair coin (heads=1, tails=0);
# the function "sample" is part of R:
#set.seed(47405) # Uncomment to set the "seed" for the random number generator.
flipsequence = sample( x=c(0,1) , prob=c(.5,.5) , size=N , replace=TRUE )
# Compute the running proportion of heads:
r = cumsum( flipsequence ) # The function "cumsum" is built in to R.
n = 1:N                    # n is a vector.
runprop = r / n            # component by component division.
# Graph the running proportion:
# To learn about the parameters of the plot function,
# type help('par') at the R command prompt.
# Note that "c" is a function in R.
plot( n , runprop , type="o" , log="x" ,
	  xlim=c(1,N) , ylim=c(0.0,1.0) , cex.axis=1.5 ,
	  xlab="Flip Number" , ylab="Proportion Heads" , cex.lab=1.5 ,
	  main="Running Proportion of Heads" , cex.main=1.5 )
# Plot a dotted horizontal line at y=.5, just as a reference line:
lines( c(1,N) , c(.5,.5) , lty=3 )
# Display the beginning of the flip sequence. These string and character
# manipulations may seem mysterious, but you can de-mystify by unpacking
# the commands starting with the innermost parentheses or brackets and
# moving to the outermost.
flipletters = paste( c("T","H")[ flipsequence[ 1:10 ] + 1 ] , collapse="" )
displaystring = paste( "Flip Sequence = " , flipletters , "..." , sep="" )
text( 5 , .9 , displaystring , adj=c(0,1) , cex=1.3 )
# Display the relative frequency at the end of the sequence.
text( N , .3 , paste("End Proportion =",runprop[N]) , adj=c(1,0) , cex=1.3 )
# To save graphs, please see update at
# http://doingbayesiandataanalysis.blogspot.com/2013/01/uniform-r-code-for-opening-saving.html