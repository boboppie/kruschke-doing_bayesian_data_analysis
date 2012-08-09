x = seq( from = -2 , to = 2 , by = 0.1 )   # Specify vector of x values.
y = x^2                                    # Specify corresponding y values.
plot( x , y , type = "l" )                 # Make a graph of the x,y points.
dev.copy2eps( file = "SimpleGraph.eps" )   # Save the plot to an EPS file.