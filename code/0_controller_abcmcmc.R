# ABC_MCMC function for MDRVX model MCMC model runs on single thread
# Requirements: Input: - Vector of parameter values Output: - Vector of summary
# statistics NB Summary statistics are used in the rho/distance calculation for
# the ABC rejection step.

## NOTES HEADER
# 1. Constraints - RH/GK files
# 2. Fit targets
# -- Including MDR fit code and targets
# 3. Summary statistics return
# 4. Plot results
# 5. Hopkins SGE


# Controller file for ABC ABC

# Read initial parameter values for ABC MCMC function and generate parameter vector
sampling_params <- read.csv("some_param_init.csv", header = T)
sampling_params <- unlist(sampling_params, use.names = T)

# NB: The parameter input vector must be NAMED
mdrvx_abcmcmc <- function(parameters) {

  # Convert parameter vector to list, then assign list contents to internal
  # operating environment of function.  Variables are now accessible by name within
  # the function.
  list2env(as.list(parameters), envir = environment())

  # Set up parameter modifications and constraints


  # Source mainloop function
  << Full text of mainloop>>

  # Execute mainloop

  # Assign summary statistics vector


  # Return summary stats vector as required by EasyABC
  return(summary_statistics)

}
