z_obs = 30 ; N_obs = 46
nulltheta = .5
tail_prob = 0  # Zero initial value for accumulation over possible N.
for ( N in 1 : (3*N_obs) ) {  # Start at 1 to avoid /0. 3*N_obs is arbitrary.
  # Create vector of z values such that z/N >= z_obs/N_obs
  zvec = (0:N)[ (0:N)/N >= z_obs/N_obs ]
  tail_prob = tail_prob + (
                dpois( N , N_obs ) * sum( dbinom( zvec , N , nulltheta ) ) )
}
show( tail_prob )
