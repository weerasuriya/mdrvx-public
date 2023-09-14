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

# # Historical population mortality scaling factor applicator function
# rm_adj <- function(rmort, mort) {
#   if (rmort < 0) {
#     return(as.vector((rmort * mort) + mort))
#   } else {
#     return(as.vector((rmort * (1 - mort)) + mort))
#   }
# }

# Convert annual risk to sub-annual risk (per dt timestep)
risk2risk <- function(arisk, dt) {
  if (any(arisk >= 1))
    stop("Annual risk too high. Must be <1")
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
annual_sum <- function(table, dt = inputs$dt, year1 = inputs$year1, yearend = inputs$yearend) {
  if (1%%dt != 0)
    stop("Number of steps in a year not clean")
  if (length(dim(table)) > 2)
    stop("Table has too many dimensions")
  ann <- c(rep(year1:yearend, each = (1/dt)))
  table <- cbind(ann, table)
  table <- aggregate(table, by = list(table[, 1]), sum)
  table <- table[, -1]
  table[, 1] <- table[, 1]/(1/dt)
  return(table)
}

# Function to generate per-100K population rates over arbitrary age ranges
# Input arguments are a 2-dimensional table (rows = ages, columns = arbitrary)
psadj <- function(table, age_range, extra = "") {
  if (extra == "m") {
    table <- annual_mean(table, dt, year1, yearend)
  } else if (extra == "s") {
    table <- annual_sum(table, dt, year1, yearend)
  }
  age_totals <- rowSums(calc$ann_agewise[, (1 + age_range)])
  var_totals <- rowSums(table[, 1 + age_range])
  adj_result <- (var_totals/age_totals) * 1e+05
  return(adj_result)
}


## INITIALISING MATRICES ###

# Initialise Timing #
steps <- (yearend - year1 + 1) %/% dt

# Prevalent compartment vector
p_compartments <- c(
  "S",
  "NTDS_L",
  "NTDS_II",
  "NTDS_IN",
  "NTDS_T",
  "NTDS_R",
  "NTDR_L",
  "NTDR_II",
  "NTDR_IN",
  "NTDR_T_IIphi",
  "NTDR_T_IIpsi",
  "NTDR_T_INphi",
  "NTDR_T_INpsi",
  "NTDR_R",
  "PTDS_II",
  "PTDS_IN",
  "PTDS_T",
  "PTDS_R",
  "PTDR_L",
  "PTDR_II",
  "PTDR_IN",
  "PTDR_T_IIphi",
  "PTDR_T_IIpsi",
  "PTDR_T_INphi",
  "PTDR_T_INpsi",
  "PTDR_R",
  "NTDR_mdt",
  "NTDR_nid",
  "PTDR_mdt",
  "PTDR_nid"
)

# Initialise the "Large Array" (la) - the axes of this 3D array are: compartment name, timestep, age
la <- array(0, dim = c(length(p_compartments), steps, max_age), dimnames = list(p_compartments, c(as.character(1:steps)), c(as.character(1:max_age))))

# Background mortality matrix
bg_mort <- array(0, dim = c(length(p_compartments), steps, max_age), dimnames = list(p_compartments, c(as.character(1:steps)), c(as.character(1:max_age))))

# List of Transits (i.e. flows)
transit <- c(
  "DSTB_deaths",
  "DRTB_deaths",
  "DSTB_inc",
  "DRTB_inc",
  # "DSTB_new",
  # "DRTB_new",
  # "DSTB_relapse",
  # "DRTB_relapse",
  # "DSTB_T_failure",
  # "DRTB_T_failure",
  "DSTB_onRx",
  "DRTB_onRx",
  "bgmort",
  # "MDR_invivo",
  # "MDR_exvivo",
  # "DSTB_rxexit",
  # "DSTB_rxexit_prop",
  # "DSTB.II.exit",
  # "DSTB.II.exit.prop",
  # "DSTB.IN.exit",
  # "DSTB.IN.exit.prop",
  # "TB_new",
  # "TB_relapse",
  "TB_inc",
  "DSTB_initRx",
  "DRTB_initRx",
  "DS_nt_intx",
  "DS_nt_intx_f",
  "DS_pt_intx",
  "DS_pt_intx_f",
  "DR_nt_intx",
  "DR_nt_intx_f",
  "DR_pt_intx",
  "DR_pt_intx_f",
  "All_nt_intx_f",
  "All_pt_intx_f",
  "All_nt_intx",
  "All_pt_intx",
  "DRTB_initRxLab",
  "DS_nt_inc",
  "DS_nt_inc_f",
  "DS_pt_inc",
  "DS_pt_inc_f",
  "DR_nt_inc",
  "DR_nt_inc_f",
  "DR_pt_inc",
  "DR_pt_inc_f",
  "All_nt_inc_f",
  "All_pt_inc_f",
  "All_nt_inc",
  "All_pt_inc",
  "total_cost"
)

# Initialise the 'Transit Array' (ta) - the axes of this 3D array are:
# compartment name (i.e. flow), timestep, age
ta <- array(0, dim = c(length(transit), steps, max_age), dimnames = list(transit, c(as.character(1:steps)), c(as.character(1:max_age))))

cost_types <- c(
  "ds_dx",
  "dr_dx",
  "dst",
  "ds_tx",
  "dr_tx",
  "tbrx"
)

# Initialise the 'Transit Array' (ta) - the axes of this 3D array are:
# compartment name (i.e. flow), timestep, age
ca <- array(0, dim = c(length(cost_types), steps, max_age), dimnames = list(cost_types, c(as.character(1:steps)), c(as.character(1:max_age))))

# Special arrays
# DSTB incidence array
dsia           <- matrix(0, steps, max_age)
rownames(dsia) <- paste("Step", as.character(1:steps))
colnames(dsia) <- paste("Age", as.character(1:max_age))

# DRTB incidence array
dria           <- matrix(0, steps, max_age)
rownames(dria) <- paste("Step", as.character(1:steps))
colnames(dria) <- paste("Age", as.character(1:max_age))

# AllTB incidence array
tbia           <- matrix(0, steps, max_age)
rownames(tbia) <- paste("Step", as.character(1:steps))
colnames(tbia) <- paste("Age", as.character(1:max_age))

# DSTB TB mortality array
dsma           <- matrix(0, steps, max_age)
rownames(dsma) <- paste("Step", as.character(1:steps))
colnames(dsma) <- paste("Age", as.character(1:max_age))

# Initialise Step1 Compartment Populations ####

for (age in 1:length(total_population[1,])) {
  for (compartment in dimnames(la)[1]) {
    la[compartment, 1, age] <- total_population[1, age] * init_tb_prop[compartment, prop_col]
  }
}

# Contact Matrix components ####
# psize and I for contact matirx
psizematrix <- matrix(0, steps, 4)
DS_Imatrix <- matrix(0, steps, 4)
DR_Imatrix <- matrix(0, steps, 4)

# Denominator matrix
psizematrix[1, 1] <- sum(la[, 1, 1:6])
psizematrix[1, 2] <- sum(la[, 1, 7:20])
psizematrix[1, 3] <- sum(la[, 1, 21:65])
psizematrix[1, 4] <- sum(la[, 1, 66:max_age])

## Total Infectious Cases by contact matrix age classes - DSTB
DS_Imatrix[1, 1] <- sum(la["PTDS_II", 1, 1:6], la["NTDS_II", 1, 1:6])
DS_Imatrix[1, 2] <- sum(la["PTDS_II", 1, 7:20], la["NTDS_II", 1, 7:20])
DS_Imatrix[1, 3] <- sum(la["PTDS_II", 1, 21:65], la["NTDS_II", 1, 21:65])
DS_Imatrix[1, 4] <- sum(la["PTDS_II", 1, 66:max_age], la["NTDS_II", 1, 66:max_age])

## Total Infectious Cases by contact matrix age classes - DRTB
DR_Imatrix[1, 1] <- sum(la["PTDR_II", 1, 1:6], la["NTDR_II", 1, 1:6], la["PTDR_T_IIphi", 1, 1:6], la["NTDR_T_IIphi", 1, 1:6])
DR_Imatrix[1, 2] <- sum(la["PTDR_II", 1, 7:20], la["NTDR_II", 1, 7:20], la["PTDR_T_IIphi", 1, 7:20], la["NTDR_T_IIphi", 1, 7:20])
DR_Imatrix[1, 3] <- sum(la["PTDR_II", 1, 21:65], la["NTDR_II", 1, 21:65], la["PTDR_T_IIphi", 1, 21:65], la["NTDR_T_IIphi", 1, 21:65])
DR_Imatrix[1, 4] <- sum(la["PTDR_II", 1, 66:max_age], la["NTDR_II", 1, 66:max_age], la["PTDR_T_IIphi", 1, 66:max_age], la["NTDR_T_IIphi", 1, 66:max_age])

# Force Of Infection ####
## Drug Sensitive
NTDS_lambda <- matrix(0, steps, max_age)
PTDS_lambda <- NTDS_lambda

## Drug Resistant
NTDR_lambda <- matrix(0, steps, max_age)
PTDR_lambda <- NTDR_lambda


# Apply the DR_ts scaling factor
DR_neta <- DS_neta * DR_ts

## Initialise the first step force of infection
# Calculate transmission terms - Lambda
NTDS_lambda_raw <- colSums(-(myneta[1:4, 1:max_age]) * z * ((DS_Imatrix[1, 1:4])/(psizematrix[1, 1:4])))
NTDS_lambda[1, 1:max_age] <- t(DS_neta * (1 - exp(NTDS_lambda_raw)))
PTDS_lambda <- NTDS_lambda

NTDR_lambda_raw <- colSums(-(myneta[1:4, 1:max_age]) * z * ((DR_Imatrix[1, 1:4])/(psizematrix[1, 1:4])))
NTDR_lambda[1, 1:max_age] <- t(DR_neta * (1 - exp(NTDR_lambda_raw)))
PTDR_lambda <- NTDR_lambda

# Scale time-dependent annual rates to relevant time step
ui          <- risk2risk(ui, dt)
uni         <- risk2risk(uni, dt)
ut          <- risk2risk(ut, dt)
rchild      <- risk2risk(rchild, dt)
radult      <- risk2risk(radult, dt)
relderly    <- risk2risk(relderly, dt)
vchild      <- risk2risk(vchild, dt)
vadult      <- risk2risk(vadult, dt)
velderly    <- risk2risk(velderly, dt)
n           <- risk2risk(n, dt)
nelderly    <- risk2risk(nelderly, dt)
xi_init     <- risk2risk(xi_init, dt)
death_rate  <- cbind(death_rate[, 1], risk2risk(death_rate[, -1], dt))
# ptdr_r <- risk2risk(ptdr_r, dt)
# ptdr_relderly <- risk2risk(ptdr_relderly, dt)

# #############################################################
# ######## Testing Zone: Parameter value assignment ###########
#
# nt_dst_prob <- 0.1
# nt_nid_prob <- 0.1
# pt_nid_prob <- pt_dst_prob * 0.5
#
# #############################################################
# #############################################################

# Initialise DST and NID probability vectors
nt_dst_p     <- vector(mode = "numeric", length = (yearend - year1 + 1))
pt_dst_p     <- vector(mode = "numeric", length = (yearend - year1 + 1))
ntii_emp_p   <- vector(mode = "numeric", length = (yearend - year1 + 1))
ptii_emp_p   <- vector(mode = "numeric", length = (yearend - year1 + 1))
ntin_emp_p   <- vector(mode = "numeric", length = (yearend - year1 + 1))
ptin_emp_p   <- vector(mode = "numeric", length = (yearend - year1 + 1))
ntii_tx_corr <- vector(mode = "numeric", length = (yearend - year1 + 1))
ptii_tx_corr <- vector(mode = "numeric", length = (yearend - year1 + 1))
ntin_tx_corr <- vector(mode = "numeric", length = (yearend - year1 + 1))
ptin_tx_corr <- vector(mode = "numeric", length = (yearend - year1 + 1))

# pt_dst_p[1:(dst_start_year - year1)] <- 0
# nt_dst_p[1:(dst_start_year - year1)] <- 0
# nt_dst_p[(dst_start_year - year1 + 1):(yearend - year1 + 1)] <- nt_dst_prob
# pt_dst_p[(dst_start_year - year1 + 1):(yearend - year1 + 1)] <- pt_dst_prob
#
# # Initialise DRTB Empirical treatment probability
# ntii_emp_p[1:(dst_start_year - year1)] <- 0
# ntin_emp_p[1:(dst_start_year - year1)] <- 0
# ptii_emp_p[1:(dst_start_year - year1)] <- 0
# ptin_emp_p[1:(dst_start_year - year1)] <- 0
# ntii_emp_p[(dst_start_year - year1 + 1):(yearend - year1 + 1)] <- (1 - nt_dst_p[(dst_start_year - year1 + 1):(yearend - year1 + 1)]) * ntii_emp_tx_p
# ntin_emp_p[(dst_start_year - year1 + 1):(yearend - year1 + 1)] <- ntin_emp_tx_p
# ptii_emp_p[(dst_start_year - year1 + 1):(yearend - year1 + 1)] <- (1 - pt_dst_p[(dst_start_year - year1 + 1):(yearend - year1 + 1)]) * ptii_emp_tx_p
# ptin_emp_p[(dst_start_year - year1 + 1):(yearend - year1 + 1)] <- ptin_emp_tx_p
#
# # Intialise DRTB Correct Treatment Probability
# ntii_tx_corr[1:(dst_start_year - year1)] <- 0
# ntin_tx_corr[1:(dst_start_year - year1)] <- 0
# ptii_tx_corr[1:(dst_start_year - year1)] <- 0
# ptin_tx_corr[1:(dst_start_year - year1)] <- 0
# ntii_tx_corr[(dst_start_year - year1 + 1):(yearend - year1 + 1)] <- nt_dst_p[(dst_start_year - year1 + 1):(yearend - year1 + 1)] + ntii_emp_p[(dst_start_year - year1 + 1):(yearend - year1 + 1)]
# ntin_tx_corr[(dst_start_year - year1 + 1):(yearend - year1 + 1)] <- ntin_emp_p[(dst_start_year - year1 + 1):(yearend - year1 + 1)]
# ptii_tx_corr[(dst_start_year - year1 + 1):(yearend - year1 + 1)] <- pt_dst_p[(dst_start_year - year1 + 1):(yearend - year1 + 1)] + ptii_emp_p[(dst_start_year - year1 + 1):(yearend - year1 + 1)]
# ptin_tx_corr[(dst_start_year - year1 + 1):(yearend - year1 + 1)] <- ptin_emp_p[(dst_start_year - year1 + 1):(yearend - year1 + 1)]

# Case Detection and Treatment initiation - holding yearwise arrays

hcdr        <- matrix(NA, yearend - year1 + 1, max_age)
hntdr_cdr   <- matrix(NA, yearend - year1 + 1, max_age)
hptdr_cdr   <- matrix(NA, yearend - year1 + 1, max_age)
hkappa      <- matrix(NA, yearend - year1 + 1, max_age)
hntdr_kappa <- matrix(NA, yearend - year1 + 1, max_age)
hptdr_kappa <- matrix(NA, yearend - year1 + 1, max_age)
hcdrscaling <- matrix(NA, yearend - year1 + 1, 3)
hn          <- matrix(NA, yearend - year1 + 1, max_age)
hui         <- matrix(NA, yearend - year1 + 1, max_age)

# Array to hold calculated SDR
sdr <- vector(mode = "numeric", length = steps)
sdr_start_yr_steps <- ((sdr_start_year - year1) * (1/dt)) + seq(1:(1/dt))
sdr[sdr_start_yr_steps] <- sdr_base

# Garbage collection
rm(age, compartment, sdr_start_yr_steps)
# The parameters fed into the function

NTDS_II_n        <- n
NTDS_II_nelderly <- nelderly

NTDS_IN_n        <- n
NTDS_IN_nelderly <- nelderly

PTDS_II_n        <- NTDS_II_n
PTDS_II_nelderly <- NTDS_II_nelderly

PTDS_IN_n        <- NTDS_IN_n
PTDS_IN_nelderly <- NTDS_IN_nelderly

NTDR_II_n        <- n
NTDR_II_nelderly <- nelderly

NTDR_IN_n        <- n
NTDR_IN_nelderly <- nelderly

PTDR_II_n        <- NTDR_II_n
PTDR_II_nelderly <- NTDR_II_nelderly

PTDR_IN_n        <- NTDR_IN_n
PTDR_IN_nelderly <- NTDR_IN_nelderly

NTDS_fchild      <- fchild
PTDS_fchild      <- fchild
NTDR_fchild      <- fchild
PTDR_fchild      <- fchild

NTDS_fadult      <- fadult
PTDS_fadult      <- fadult
NTDR_fadult      <- fadult
PTDR_fadult      <- fadult

NTDS_felderly    <- felderly
PTDS_felderly    <- felderly
NTDR_felderly    <- felderly
PTDR_felderly    <- felderly

NTDS_rchild      <- rchild
PTDS_rchild      <- rchild
NTDR_rchild      <- rchild
PTDR_rchild      <- rchild

NTDS_radult      <- radult
PTDS_radult      <- radult
NTDR_radult      <- radult
PTDR_radult      <- radult

NTDS_relderly    <- relderly
PTDS_relderly    <- relderly
NTDR_relderly    <- relderly
PTDR_relderly    <- relderly

NTDS_pchild      <- pchild
PTDS_pchild      <- pchild
NTDR_pchild      <- pchild
PTDR_pchild      <- pchild

NTDS_padult      <- padult
PTDS_padult      <- padult
NTDR_padult      <- padult
PTDR_padult      <- padult

NTDS_pelderly    <- pelderly
PTDS_pelderly    <- pelderly
NTDR_pelderly    <- pelderly
PTDR_pelderly    <- pelderly

NTDS_vchild      <- vchild
PTDS_vchild      <- vchild
NTDR_vchild      <- vchild
PTDR_vchild      <- vchild

NTDS_vadult      <- vadult
PTDS_vadult      <- vadult
NTDR_vadult      <- vadult
PTDR_vadult      <- vadult

NTDS_velderly    <- velderly
PTDS_velderly    <- velderly
NTDR_velderly    <- velderly
PTDR_velderly    <- velderly

NTDS_omega       <- omega
PTDS_omega       <- omega
NTDR_omega       <- omega
PTDR_omega       <- omega

DS_CDRscale      <- CDRscale
DR_CDRscale      <- CDRscale

DS_CDRscaleE     <- CDRscaleE
DR_CDRscaleE     <- CDRscaleE

NTDS_x           <- x
PTDS_x           <- x
NTDR_x           <- x
PTDR_x           <- x

NTDR_T_IIpsi_tau <- tau_psi
NTDR_T_INpsi_tau <- tau_psi
PTDR_T_IIpsi_tau <- tau_psi
PTDR_T_INpsi_tau <- tau_psi

NTDR_T_IIphi_tau <- tau_phi
NTDR_T_INphi_tau <- tau_phi
PTDR_T_IIphi_tau <- tau_phi
PTDR_T_INphi_tau <- tau_phi

# Create age-wise parameter vectors for key natural history (nh) outcome parameters
# Params: p,f,r,v,n

# Natural history 'Older adult' (age 55-65) parameter calculations
# Older adult parameters for natural cure rate, n
NTDS_II_noldadu <- (NTDS_II_n + NTDS_II_nelderly)/2
NTDS_IN_noldadu <- (NTDS_IN_n + NTDS_IN_nelderly)/2
PTDS_II_noldadu <- (PTDS_II_n + PTDS_II_nelderly)/2
PTDS_IN_noldadu <- (PTDS_IN_n + PTDS_IN_nelderly)/2

NTDR_II_noldadu <- (NTDR_II_n + NTDR_II_nelderly)/2
NTDR_IN_noldadu <- (NTDR_IN_n + NTDR_IN_nelderly)/2
PTDR_II_noldadu <- (PTDR_II_n + PTDR_II_nelderly)/2
PTDR_IN_noldadu <- (PTDR_IN_n + PTDR_IN_nelderly)/2

# Older adult parameters for reactivation Rate, v
NTDS_voldadu    <- (NTDS_vadult + NTDS_velderly)/2
NTDR_voldadu    <- (NTDR_vadult + NTDR_velderly)/2
PTDS_voldadu    <- (PTDS_vadult + PTDS_velderly)/2
PTDR_voldadu    <- (PTDR_vadult + PTDR_velderly)/2

# Older adult parameters for reactivation Rate, r
NTDS_roldadu    <- (NTDS_radult + NTDS_relderly)/2
NTDR_roldadu    <- (NTDR_radult + NTDR_relderly)/2
PTDS_roldadu    <- (PTDS_radult + PTDS_relderly)/2
PTDR_roldadu    <- (PTDR_radult + PTDR_relderly)/2

# Age-wise parameter vectors - Natural History Parameters (p,f,r,v,n)
NTDS_p = c((rep(NTDS_pchild, l = chiyrs)), (rep(NTDS_padult, l = aduyrs)), (rep(NTDS_pelderly, l = eldyrs)))
NTDS_f = c((rep(NTDS_fchild, l = chiyrs)), (rep(NTDS_fadult, l = aduyrs)), (rep(NTDS_felderly, l = eldyrs)))
NTDS_v = c((rep(NTDS_vadult, l = (yaduyrs + chiyrs))), (rep(NTDS_voldadu, l = (aduyrs - yaduyrs))), (rep(NTDS_velderly, l = eldyrs)))
NTDS_r = c((rep(NTDS_rchild, l = chiyrs)), (rep(NTDS_radult, l = yaduyrs)), (rep(NTDS_roldadu, l = (aduyrs - yaduyrs))), (rep(NTDS_relderly, l = eldyrs)))

NTDS_II_n = c((rep(NTDS_II_n, l = chiyrs)), (rep(NTDS_II_n, l = yaduyrs)), (rep(NTDS_II_noldadu, l = (aduyrs - yaduyrs))), (rep(NTDS_II_nelderly, l = eldyrs)))
NTDS_IN_n = c((rep(NTDS_IN_n, l = chiyrs)), (rep(NTDS_II_n, l = yaduyrs)), (rep(NTDS_IN_noldadu, l = (aduyrs - yaduyrs))), (rep(NTDS_IN_nelderly, l = eldyrs)))

NTDR_p = c((rep(NTDR_pchild, l = chiyrs)), (rep(NTDR_padult, l = aduyrs)), (rep(NTDR_pelderly, l = eldyrs)))
NTDR_f = c((rep(NTDR_fchild, l = chiyrs)), (rep(NTDR_fadult, l = aduyrs)), (rep(NTDR_felderly, l = eldyrs)))
NTDR_v = c((rep(NTDR_vadult, l = (yaduyrs + chiyrs))), (rep(NTDR_voldadu, l = (aduyrs - yaduyrs))), (rep(NTDR_velderly, l = eldyrs)))
NTDR_r = c((rep(NTDR_rchild, l = chiyrs)), (rep(NTDR_radult, l = yaduyrs)), (rep(NTDR_roldadu, l = (aduyrs - yaduyrs))), (rep(NTDR_relderly, l = eldyrs)))

NTDR_II_n = c((rep(NTDR_II_n, l = chiyrs)), (rep(NTDR_II_n, l = yaduyrs)), (rep(NTDR_II_noldadu, l = (aduyrs - yaduyrs))), (rep(NTDR_II_nelderly, l = eldyrs)))
NTDR_IN_n = c((rep(NTDR_IN_n, l = chiyrs)), (rep(NTDR_II_n, l = yaduyrs)), (rep(NTDR_IN_noldadu, l = (aduyrs - yaduyrs))), (rep(NTDR_IN_nelderly, l = eldyrs)))

PTDS_p = c((rep(PTDS_pchild, l = chiyrs)), (rep(PTDS_padult, l = aduyrs)), (rep(PTDS_pelderly, l = eldyrs)))
PTDS_f = c((rep(PTDS_fchild, l = chiyrs)), (rep(PTDS_fadult, l = aduyrs)), (rep(PTDS_felderly, l = eldyrs)))
PTDS_v = c((rep(PTDS_vadult, l = (yaduyrs + chiyrs))), (rep(PTDS_voldadu, l = (aduyrs - yaduyrs))), (rep(PTDS_velderly, l = eldyrs)))
PTDS_r = c((rep(PTDS_rchild, l = chiyrs)), (rep(PTDS_radult, l = yaduyrs)), (rep(PTDS_roldadu, l = (aduyrs - yaduyrs))), (rep(PTDS_relderly, l = eldyrs)))

PTDS_II_n = c((rep(PTDS_II_n, l = chiyrs)), (rep(PTDS_II_n, l = yaduyrs)), (rep(PTDS_II_noldadu, l = (aduyrs - yaduyrs))), (rep(PTDS_II_nelderly, l = eldyrs)))
PTDS_IN_n = c((rep(PTDS_IN_n, l = chiyrs)), (rep(PTDS_II_n, l = yaduyrs)), (rep(PTDS_IN_noldadu, l = (aduyrs - yaduyrs))), (rep(PTDS_IN_nelderly, l = eldyrs)))

# PTDR_p = c((rep(PTDR_pchild, l = chiyrs)), (rep(PTDR_padult, l = aduyrs)), (rep(PTDR_pelderly, l = eldyrs)))
# Disable PTDR_L compartment
PTDR_p = c(rep(1, 100))
PTDR_f = c((rep(PTDR_fchild, l = chiyrs)), (rep(PTDR_fadult, l = aduyrs)), (rep(PTDR_felderly, l = eldyrs)))
PTDR_v = c((rep(PTDR_vadult, l = (yaduyrs + chiyrs))), (rep(PTDR_voldadu, l = (aduyrs - yaduyrs))), (rep(PTDR_velderly, l = eldyrs)))
PTDR_r = c((rep(PTDR_rchild, l = chiyrs)), (rep(PTDR_radult, l = yaduyrs)), (rep(PTDR_roldadu, l = (aduyrs - yaduyrs))), (rep(PTDR_relderly, l = eldyrs)))

PTDR_II_n = c((rep(PTDR_II_n, l = chiyrs)), (rep(PTDR_II_n, l = yaduyrs)), (rep(PTDR_II_noldadu, l = (aduyrs - yaduyrs))), (rep(PTDR_II_nelderly, l = eldyrs)))
PTDR_IN_n = c((rep(PTDR_IN_n, l = chiyrs)), (rep(PTDR_II_n, l = yaduyrs)), (rep(PTDR_IN_noldadu, l = (aduyrs - yaduyrs))), (rep(PTDR_IN_nelderly, l = eldyrs)))

# PTDR_r <- c(rep(ptdr_r, max_age))
# Apply uiscaleC
if (uiscaleC < 0 ) {
  uichild  <- ui + (uiscaleC * ui)
  unichild <- uni + (uiscaleC * uni)
  utchild  <- ut + (uiscaleC * ut)
} else if (uiscaleC >= 0) {
  uichild  <- ui + (uiscaleC * (1 - ui))
  unichild <- uni + (uiscaleC * (1 - uni))
  utchild  <- ut + (uiscaleC * (1 - ut))
}

# Apply uiscaleA
if (uiscaleA < 0 ) {
  uiadult  <- ui + (uiscaleA * ui)
  uniadult <- uni + (uiscaleA * uni)
  utadult  <- ut + (uiscaleA * ut)
} else if (uiscaleA >= 0) {
  uiadult  <- ui + (uiscaleA * (1 - ui))
  uniadult <- uni + (uiscaleA * (1 - uni))
  utadult  <- ut + (uiscaleA * (1 - ut))
}

# Apply uiscaleE
if (uiscaleE < 0 ) {
  uielderly  <- ui + (uiscaleE * ui)
  unielderly <- uni + (uiscaleE * uni)
  utelderly  <- ut + (uiscaleE * ut)
} else if (uiscaleE >= 0) {
  uielderly  <- ui + (uiscaleE * (1 - ui))
  unielderly <- uni + (uiscaleE * (1 - uni))
  utelderly  <- ut + (uiscaleE * (1 - ut))
}

# Construct age-wise mortality vectors with scaled mortality rates
NTDS_II_u  <- c((rep(uichild, l = chiyrs)), (rep(uiadult, l = aduyrs)), (rep(uielderly, l = eldyrs)))
NTDS_IN_u  <- c((rep(unichild, l = chiyrs)), (rep(uniadult, l = aduyrs)), (rep(unielderly, l = eldyrs)))
NTDS_T_u   <- c((rep(utchild, l = chiyrs)), (rep(utadult, l = aduyrs)), (rep(utelderly, l = eldyrs)))

PTDS_II_u  <- c((rep(uichild, l = chiyrs)), (rep(uiadult, l = aduyrs)), (rep(uielderly, l = eldyrs)))
PTDS_IN_u  <- c((rep(unichild, l = chiyrs)), (rep(uniadult, l = aduyrs)), (rep(unielderly, l = eldyrs)))
PTDS_T_u   <- c((rep(utchild, l = chiyrs)), (rep(utadult, l = aduyrs)), (rep(utelderly, l = eldyrs)))

NTDR_II_u  <- c((rep(uichild, l = chiyrs)), (rep(uiadult, l = aduyrs)), (rep(uielderly, l = eldyrs)))
NTDR_IN_u  <- c((rep(unichild, l = chiyrs)), (rep(uniadult, l = aduyrs)), (rep(unielderly, l = eldyrs)))
NTDR_T_u   <- c((rep(utchild, l = chiyrs)), (rep(utadult, l = aduyrs)), (rep(utelderly, l = eldyrs)))

PTDR_II_u  <- c((rep(uichild, l = chiyrs)), (rep(uiadult, l = aduyrs)), (rep(uielderly, l = eldyrs)))
PTDR_IN_u  <- c((rep(unichild, l = chiyrs)), (rep(uniadult, l = aduyrs)), (rep(unielderly, l = eldyrs)))
PTDR_T_u   <- c((rep(utchild, l = chiyrs)), (rep(utadult, l = aduyrs)), (rep(utelderly, l = eldyrs)))

NTDR_mdt_u <- c((rep(uichild, l = chiyrs)), (rep(uiadult, l = aduyrs)), (rep(uielderly, l = eldyrs)))
PTDR_mdt_u <- c((rep(uichild, l = chiyrs)), (rep(uiadult, l = aduyrs)), (rep(uielderly, l = eldyrs)))
NTDR_nid_u <- c((rep(unichild, l = chiyrs)), (rep(uniadult, l = aduyrs)), (rep(unielderly, l = eldyrs)))
PTDR_nid_u <- c((rep(unichild, l = chiyrs)), (rep(uniadult, l = aduyrs)), (rep(unielderly, l = eldyrs)))
# Case Detection Rate
# Create age and year-wise CDR matrices
DS_CDR <- matrix(0, (yearend - year1 + 1), max_age)

# Assign the unscaled CDRs to all ages.
DS_CDR[, 1:max_age] <- uDS_CDR

# Age-wise scaling of CDR DS_CDRScaled, DS_CDRScaledE, DR_CDRScaled,
# Older adult CDRScaling
DS_CDRscaleO <- (DS_CDRscale + DS_CDRscaleE)/2

# ## DSTB
# ## Apply scaling to children and younger (15-55) adults
if (DS_CDRscale < 0) {
  DS_CDR[, 1:(chiyrs + yaduyrs)] <- sapply(((DS_CDRscale * DS_CDR[, 1:(chiyrs + yaduyrs)]) + DS_CDR[, 1:(chiyrs + yaduyrs)]), min, 1)
} else {
  DS_CDR[, 1:(chiyrs + yaduyrs)] <- sapply(((DS_CDRscale * (1 - DS_CDR[, 1:(chiyrs + yaduyrs)])) + DS_CDR[, 1:(chiyrs + yaduyrs)]), min, 1)
}

# ## Apply scaling to older adults (55-65)
if (DS_CDRscaleO < 0) {
  DS_CDR[, (chiyrs + yaduyrs + 1):(chiyrs + aduyrs)] <- sapply(((DS_CDRscaleO * DS_CDR[, (chiyrs + yaduyrs + 1):(chiyrs + aduyrs)]) + DS_CDR[, (chiyrs + yaduyrs + 1):(chiyrs + aduyrs)]), min, 1)
} else {
  DS_CDR[, (chiyrs + yaduyrs +1 ):(chiyrs + aduyrs)] <- sapply(((DS_CDRscaleO * (1 - DS_CDR[, (chiyrs + yaduyrs + 1):(chiyrs + aduyrs)])) + DS_CDR[, (chiyrs + yaduyrs + 1):(chiyrs + aduyrs)]), min, 1)
}
#
# ## Apply scaling to elderly (65+)
if (DS_CDRscaleE < 0) {
  DS_CDR[, (chiyrs + aduyrs + 1):(max_age)] <- sapply(((DS_CDRscaleE * DS_CDR[, (chiyrs + aduyrs + 1):(max_age)]) + DS_CDR[, (chiyrs + aduyrs + 1):(max_age)]), min, 1)
} else {
  DS_CDR[, (chiyrs + aduyrs + 1):(max_age)] <- sapply((DS_CDRscaleE * (1 - DS_CDR[, (chiyrs + aduyrs + 1):(max_age)]) + DS_CDR[, (chiyrs + aduyrs + 1):(max_age)]), min, 1)
}

# Scale down to per time step
DS_CDR <- 1 - (1 - DS_CDR)^(dt)

# Add a column of years
DS_CDR <- cbind(year1:yearend, DS_CDR)

# Create matrix to hold new treatment initiation rates
NTDS_II_kap_main <- matrix(0, yearend - year1 + 1, max_age)
PTDS_II_kap_main <- matrix(0, yearend - year1 + 1, max_age)
NTDS_IN_kap_main <- matrix(0, yearend - year1 + 1, max_age)
PTDS_IN_kap_main <- matrix(0, yearend - year1 + 1, max_age)

NTDR_II_kap_main <- matrix(0, yearend - year1 + 1, max_age)
PTDR_II_kap_main <- matrix(0, yearend - year1 + 1, max_age)
NTDR_IN_kap_main <- matrix(0, yearend - year1 + 1, max_age)
PTDR_IN_kap_main <- matrix(0, yearend - year1 + 1, max_age)

NTDR_II_kap_mdt  <- matrix(0, yearend - year1 + 1, max_age)
PTDR_II_kap_mdt  <- matrix(0, yearend - year1 + 1, max_age)
NTDR_IN_kap_nid  <- matrix(0, yearend - year1 + 1, max_age)
PTDR_IN_kap_nid  <- matrix(0, yearend - year1 + 1, max_age)


######################### Treatment Probabilities for DRTB ##################################
# Create linear scale up of empirical treatment probability until DST start year1

pt_dst_p[1:(dst_start_year - year1)] <- 0
nt_dst_p[1:(dst_start_year - year1)] <- 0
nt_dst_p[(dst_start_year - year1 + 1):(yearend - year1 + 1)] <- nt_dst_prob
pt_dst_p[(dst_start_year - year1 + 1):(yearend - year1 + 1)] <- pt_dst_prob

# Initialise DRTB Empirical treatment probability
# Empirical treatment probability = 0 until 1 year before when MDR acquisition begins
ntii_emp_p[1:(mdrac_yr - year1)] <- 0
ntin_emp_p[1:(mdrac_yr - year1)] <- 0
ptii_emp_p[1:(mdrac_yr - year1)] <- 0
ptin_emp_p[1:(mdrac_yr - year1)] <- 0

# Internal representation of empirical ntii/ptii treatment probability
internal_ntii_emp_tx_p <- (1 - nt_dst_prob) * ntii_emp_tx_p
internal_ptii_emp_tx_p <- (1 - pt_dst_prob) * ptii_emp_tx_p

# Linear increase in empirical treatment probability until dst_start_year
for (year in (mdrac_yr - year1 + 1):(dst_start_year - year1)) {
  index <- year - (mdrac_yr - year1)
  ntin_emp_p[year] <- index * ntin_emp_tx_p / (dst_start_year - mdrac_yr)
  ptin_emp_p[year] <- index * ptin_emp_tx_p / (dst_start_year - mdrac_yr)
  ntii_emp_p[year] <- index * internal_ntii_emp_tx_p / (dst_start_year - mdrac_yr)
  ptii_emp_p[year] <- index * internal_ptii_emp_tx_p / (dst_start_year - mdrac_yr)
}

# Static empirical treatment probability after dst_start_year
ntii_emp_p[(dst_start_year - year1 + 1):(yearend - year1 + 1)] <- (1 - nt_dst_p[(dst_start_year - year1 + 1):(yearend - year1 + 1)]) * ntii_emp_tx_p
ntin_emp_p[(dst_start_year - year1 + 1):(yearend - year1 + 1)] <- ntin_emp_tx_p
ptii_emp_p[(dst_start_year - year1 + 1):(yearend - year1 + 1)] <- (1 - pt_dst_p[(dst_start_year - year1 + 1):(yearend - year1 + 1)]) * ptii_emp_tx_p
ptin_emp_p[(dst_start_year - year1 + 1):(yearend - year1 + 1)] <- ptin_emp_tx_p

# Intialise DRTB Correct Treatment Probability
# ntii_tx_corr[1:(dst_start_year - year1)] <- 0
# ntin_tx_corr[1:(dst_start_year - year1)] <- 0
# ptii_tx_corr[1:(dst_start_year - year1)] <- 0
# ptin_tx_corr[1:(dst_start_year - year1)] <- 0
# ntii_tx_corr[(dst_start_year - year1 + 1):(yearend - year1 + 1)] <- nt_dst_p[(dst_start_year - year1 + 1):(yearend - year1 + 1)] + ntii_emp_p[(dst_start_year - year1 + 1):(yearend - year1 + 1)]
# ntin_tx_corr[(dst_start_year - year1 + 1):(yearend - year1 + 1)] <- ntin_emp_p[(dst_start_year - year1 + 1):(yearend - year1 + 1)]
# ptii_tx_corr[(dst_start_year - year1 + 1):(yearend - year1 + 1)] <- pt_dst_p[(dst_start_year - year1 + 1):(yearend - year1 + 1)] + ptii_emp_p[(dst_start_year - year1 + 1):(yearend - year1 + 1)]
# ptin_tx_corr[(dst_start_year - year1 + 1):(yearend - year1 + 1)] <- ptin_emp_p[(dst_start_year - year1 + 1):(yearend - year1 + 1)]
ntii_tx_corr <- ntii_emp_p + nt_dst_p
ntin_tx_corr <- ntin_emp_p
ptii_tx_corr <- ptii_emp_p + pt_dst_p
ptin_tx_corr <- ptin_emp_p

############################################################################################

# Calculate treatment initiation rates

for (year in 1:(yearend - year1 + 1)) {
  #DSTB
  NTDS_II_kap_main[year, ] <- ((NTDS_II_u + NTDS_II_n) * DS_CDR[year, -1])/(1 - DS_CDR[year, -1])
  PTDS_II_kap_main[year, ] <- ((PTDS_II_u + PTDS_II_n) * DS_CDR[year, -1])/(1 - DS_CDR[year, -1])
  NTDS_IN_kap_main[year, ] <- ((NTDS_IN_u + NTDS_IN_n) * (DS_CDR[year, -1] * e))/(1 - (DS_CDR[year, -1] * e))
  PTDS_IN_kap_main[year, ] <- ((PTDS_IN_u + PTDS_IN_n) * (DS_CDR[year, -1] * e))/(1 - (DS_CDR[year, -1] * e))
  # #DRTB
  # NTDR_II_kap_main[year, ] <- ((NTDR_II_u + NTDR_II_n) * NTDR_CDR[year, -1])/(1 - NTDR_CDR[year, -1])
  # PTDR_II_kap_main[year, ] <- ((PTDR_II_u + PTDR_II_n) * PTDR_CDR[year, -1])/(1 - PTDR_CDR[year, -1])
  # NTDR_IN_kap_main[year, ] <- ((NTDR_IN_u + NTDR_IN_n) * (NTDR_CDR[year, -1] * e))/(1 - (NTDR_CDR[year, -1] * e))
  # PTDR_IN_kap_main[year, ] <- ((PTDR_IN_u + PTDR_IN_n) * (PTDR_CDR[year, -1] * e))/(1 - (PTDR_CDR[year, -1] * e))

  NTDR_II_kap_main[year, ] <- NTDS_II_kap_main[year, ] * ntii_tx_corr[year]
  NTDR_II_kap_mdt[year, ]  <- NTDS_II_kap_main[year, ] * (1 - ntii_tx_corr[year])
  PTDR_II_kap_main[year, ] <- PTDS_II_kap_main[year, ] * ptii_tx_corr[year]
  PTDR_II_kap_mdt[year, ]  <- PTDS_II_kap_main[year, ] * (1 - ptii_tx_corr[year])
  NTDR_IN_kap_main[year, ] <- NTDS_IN_kap_main[year, ] * ntin_tx_corr[year]
  NTDR_IN_kap_nid[year, ]  <- NTDS_IN_kap_main[year, ] * (1 - ntin_tx_corr[year])
  PTDR_IN_kap_main[year, ] <- PTDS_IN_kap_main[year, ] * ptin_tx_corr[year]
  PTDR_IN_kap_nid[year, ]  <- PTDS_IN_kap_main[year, ] * (1 - ptin_tx_corr[year])
}

# Rebind a column of years - DSTB
NTDS_II_kap_main <- cbind(DS_CDR[, 1], NTDS_II_kap_main)
PTDS_II_kap_main <- cbind(DS_CDR[, 1], PTDS_II_kap_main)
NTDS_IN_kap_main <- cbind(DS_CDR[, 1], NTDS_IN_kap_main)
PTDS_IN_kap_main <- cbind(DS_CDR[, 1], PTDS_IN_kap_main)

# Rebind a column of years - DRTB
NTDR_II_kap_main <- cbind(DS_CDR[, 1], NTDR_II_kap_main)
PTDR_II_kap_main <- cbind(DS_CDR[, 1], PTDR_II_kap_main)
NTDR_IN_kap_main <- cbind(DS_CDR[, 1], NTDR_IN_kap_main)
PTDR_IN_kap_main <- cbind(DS_CDR[, 1], PTDR_IN_kap_main)

# Rebind column of years - DRTB misdiagnosis
NTDR_II_kap_mdt <- cbind(DS_CDR[, 1], NTDR_II_kap_mdt)
PTDR_II_kap_mdt <- cbind(DS_CDR[, 1], PTDR_II_kap_mdt)
NTDR_IN_kap_nid <- cbind(DS_CDR[, 1], NTDR_IN_kap_nid)
PTDR_IN_kap_nid <- cbind(DS_CDR[, 1], PTDR_IN_kap_nid)

# Test if any sum of natural cure, mortality and treatment initiation rate exceeds 1
if (any(
  (NTDS_II_kap_main[, -1] + NTDS_II_n + NTDS_II_u > 1) |
  (PTDS_II_kap_main[, -1] + PTDS_II_n + PTDS_II_u > 1) |
  (NTDS_IN_kap_main[, -1] + NTDS_IN_n + NTDS_IN_u > 1) |
  (PTDS_IN_kap_main[, -1] + PTDS_IN_n + PTDS_IN_u > 1) |
  (NTDR_II_kap_main[, -1] + NTDR_II_n + NTDR_II_u + NTDR_II_kap_mdt[, -1] > 1) |
  (PTDR_II_kap_main[, -1] + PTDR_II_n + PTDR_II_u + PTDR_II_kap_mdt[, -1] > 1) |
  (NTDR_IN_kap_main[, -1] + NTDR_IN_n + NTDR_IN_u + NTDR_IN_kap_nid[, -1] > 1) |
  (PTDR_IN_kap_main[, -1] + PTDR_IN_n + PTDR_IN_u + PTDR_IN_kap_nid[, -1] > 1)
)) {
  return("Reject-ntu")
}
