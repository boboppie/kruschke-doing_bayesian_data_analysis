# For NHST exercise regarding two-tier testing.

N1 = 30       # Number of flips for first test. Try 17.
N2 = 15       # Number of _additional_ flips for second test. Try 27 or 50.

theta = .5    # Hypothesized bias of coin.
FAmax = .05   # False Alarm maximum for a single test.
NT = N1 + N2  # Total number of flips.

# Determine critical values for N1:
# EXPLAIN what each function does and why, including
# dbinom, cumsum, which, max, and (0:N)[...]
loCritN1 = (0:N1)[ max( which( cumsum( dbinom(0:N1,N1,theta) ) <= FAmax/2 ) ) ]
hiCritN1 = (N1:0)[ max( which( cumsum( dbinom(N1:0,N1,theta) ) <= FAmax/2 ) ) ]
# Compute actual false alarm rate for those critical values.
# EXPLAIN what this does and why.
FA1 = sum( ( 0:N1 <= loCritN1 | 0:N1 >= hiCritN1 ) * dbinom(0:N1,N1,theta) )
cat( "N1:",N1 , ", lo:",loCritN1 , ", hi:",hiCritN1 , ", FA:",FA1 , "\n" )

# Determine critical values for NT:
# EXPLAIN what each function does and why, including
# dbinom, cumsum, which, max, and (0:N)[...]
loCritNT = (0:NT)[ max( which( cumsum( dbinom(0:NT,NT,theta) ) <= FAmax/2 ) ) ]
hiCritNT = (NT:0)[ max( which( cumsum( dbinom(NT:0,NT,theta) ) <= FAmax/2 ) ) ]
# Compute actual false alarm rate for those critical values.
# EXPLAIN what this does and why.
FAT = sum( ( 0:NT <= loCritNT | 0:NT >= hiCritNT ) * dbinom(0:NT,NT,theta) )
cat( "NT:",NT , ", lo:",loCritNT , ", hi:",hiCritNT , ", FA:",FAT , "\n" )

# Determine actual false alarm rate for the two-tier test:
# EXPLAIN each of the matrices below --- what is in each one?
Z1mat = matrix( 0:N1 , nrow=N2+1 , ncol=N1+1 , byrow=TRUE )
ZTmat = outer( 0:N2 , 0:N1 , "+" )
pZTmat = outer( dbinom( 0:N2 , N2 , theta ) , dbinom( 0:N1 , N1 , theta ) )
# EXPLAIN the matrices in computation below.
FA1or2 = sum( ( ( ZTmat <= loCritNT | ZTmat >= hiCritNT ) # double dagger matrix
              | ( Z1mat <= loCritN1 | Z1mat >= hiCritN1 ) # single dagger matrix
              ) * pZTmat )
cat( "Two tier FA:" , FA1or2 , "\n" )
