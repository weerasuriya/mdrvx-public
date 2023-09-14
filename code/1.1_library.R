# Timekeeper function - returns the current year and substep, given the current,
# start year and dt
timekeeper <- function(i, dt, year1) {
  # Calculate number of steps per annum
  if (dt > 1)
    stop("Timestep too big for timekeeper function")
  steps_pa = as.integer((1/dt))

  # Special case - dt=1
  if (dt == 1) {
    yr = year1 + i - 1
    sbstep = 1
    assign("substep", 1, envir = environment())
    assign("current_year", yr, envir = environment())
    return(c(sbstep, yr))
  }

  # If stepnumber / steps_pa divides without remainder, then is last step in the
  # year
  if (i%%steps_pa == 0) {
    sbstep <- steps_pa
    yr <- year1 - 1 + (i%/%steps_pa)
  } else {
    # If stepnumber / steps_pa has a remainder, this remainder is the step of the
    # year
    sbstep <- i%%steps_pa
    yr <- year1 + (i%/%steps_pa)
  }
  yrindex <- yr - 1899
  assign("substep", sbstep, envir = environment())
  assign("current_year", yr, envir = environment())
  assign("cyindex", yrindex, envir = environment())
  return(c(sbstep, yr))
}

yrfinder <- function(i, steps_pa = 4, year1 = 1900) {
  j <- ifelse(test = (i%%steps_pa == 0),
              yes = year1 - 1 + (i%/%steps_pa),
              no = year1 + (i%/%steps_pa))
  return(j)
}

# Convert annual risk to sub-annual risk (per dt timestep)
risk2risk <- function(arisk, dt) {
  # if (any(arisk >= 1))
  #   stop("Annual risk too high. Must be <1")
  stepspa <- as.integer(1/dt)
  asurv   <- 1 - arisk
  strisk  <- 1 - (asurv^(1/stepspa))
  rm(arisk, asurv)
  return(strisk)
}

# Convert matrix from time steps to years, **take mean** over the timestep values
# and return yearly value
annual_mean <- function(table, dt = inputs$dt, year1 = inputs$year1, yearend = inputs$yearend) {
  if (1%%dt != 0)
    stop("Number of steps in a year not clean")
  if (length(dim(table)) > 2)
    stop("Table has too many dimensions")
  ann   <- c(rep(year1:yearend, each = (1/dt)))
  table <- cbind(ann, table)
  table <- aggregate(table, by = list(table[, 1]), mean)
  table <- table[, -1]
  return(table)
}

# Convert matrix from time steps to years, **sum** over the timestep values and
# return yearly value
annual_sum   <- function(table, dt = inputs$dt, year1 = inputs$year1, yearend = inputs$yearend) {
  if (1%%dt != 0)
    stop("Number of steps in a year not clean")
  if (length(dim(table)) > 2)
    stop("Table has too many dimensions")
  ann        <- c(rep(year1:yearend, each = (1/dt)))
  table      <- cbind(ann, table)
  table      <- aggregate(table, by = list(table[, 1]), sum)
  table      <- table[, -1]
  table[, 1] <- table[, 1]/(1/dt)
  return(table)
}

# Function to generate per-100K population rates over arbitrary age ranges
# Input arguments are a 2-dimensional table (rows = ages, columns = arbitrary)
psadj <- function(table, age_range, denom = 1e+05) {
  # Adjust for year column + zero-age
  age_range <- age_range + 2
  age_totals <- rowSums(ir$ann_agewise[, ..age_range])
  var_totals <- rowSums(table[, ..age_range])
  adj_result <- (var_totals/age_totals) * denom
  return(adj_result)
}

# GLF-CDR generator
glf <- function(A = 0, K = 1, C = 1, B = 1, Q = 1, v = 1, M = 2003, tstart = 1970, tend = 2050, et = 2099, ef = 1900) {
  climb = A + ((K - A)/((C + (Q * exp(-B * ((tstart:tend) - M)))^(1/v))))
  rbind(cbind((ef:(tstart - 1)), min(climb)), cbind((tstart:tend), climb), cbind(((tend + 1):et), max(climb)))
}


# CDR scaling function
CDR_scaler <- function(scl, cdr) {
  if (scl < 0) {
    pmin(((1 + scl) * cdr), 1)
  } else {
    pmin(((scl * (1 - cdr)) + cdr), 1)
  }
}


# Testing vapply subfunction
vsf <- function(agegrp, cmp = NULL, stp = current_step) {
  sum(la[cmp, stp, agegrp])
}

# Stepfinder
sf <- function(yr, dt = 0.25, yrstart = year1) {
   ((yr - yrstart) * (1 / dt)) + (1:(1 / dt))
}
