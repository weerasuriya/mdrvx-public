mainloop <- function(data, mode = 0, para_static, para_variable, startcon, vaccine, transmission = 1, mdrac = 1, rx = 1) {
  time.start <- Sys.time()
  require(here)
  setwd(here("code"))

  list2env(data, environment(), hash = TRUE)
  list2env(para_static, environment(), hash = TRUE)
  list2env(para_variable, environment(), hash = TRUE)
  list2env(startcon, environment(), hash = TRUE)

  # Set up age classes - children, adults, younger adults, elderly
  chiyrs    <- 15
  aduyrs    <- 50
  yaduyrs   <- 40
  eldyrs    <- (max_age - chiyrs - aduyrs)

  # Global MDR switch. If mdrac==0, then disable MDR acquisition in this run.
  if (mdrac == 0) {
    if (exists("xi") == TRUE) {
      rm(xi)
    }
    assign("NT_xi", 0, envir = environment())
    assign("PT_xi", 0, envir = environment())
    mdrac_yr <- year1
  } else (
    mdrac_yr <- mdrac
  )

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

  # Case detection switch
  if (rx == 0) {
    DS_CDR[, 2:101] <- 0
    # DR_CDR[, 2:101] <- 0
  }

  ## Main iterator
  for (i in 2:steps) {

    current_step <- i
    # Set the substep and current_year variables for this iteration
    current_year <- timekeeper(i, dt, year1)[2]
    substep      <- timekeeper(i, dt, year1)[1]

    # # Deaths - select age-specific background death rates for current year
    if (bgu == 1) {
      if (substep == 1) {
        u <- death_rate[match(current_year - 1, death_rate[, 1]), 2:(max_age + 1)]
      } else if (substep != 1) {
        u <- death_rate[match(current_year, death_rate[, 1]), 2:(max_age + 1)]
      }
    } else if (bgu == 0) {
      u <- rep(0, 100)
    }

    # If MDR enabled, then start resistance acquisition in this year
    if (mdrac != 0) {
      if (current_year < mdrac) {
        assign("NT_xi", 0, envir = environment())
        assign("PT_xi", 0, envir = environment())
      } else if (current_year >= mdrac) {
        assign("NT_xi", xi_init, envir = environment())
        assign("PT_xi", xi_init, envir = environment())
      }
    }

    # Set treatment initiation rate - DSTB
    NTDS_II_kappa <- NTDS_II_kap_main[match(current_year, NTDS_II_kap_main[, 1]), 2:(max_age + 1)]
    PTDS_II_kappa <- PTDS_II_kap_main[match(current_year, PTDS_II_kap_main[, 1]), 2:(max_age + 1)]
    NTDS_IN_kappa <- NTDS_IN_kap_main[match(current_year, NTDS_IN_kap_main[, 1]), 2:(max_age + 1)]
    PTDS_IN_kappa <- PTDS_IN_kap_main[match(current_year, PTDS_IN_kap_main[, 1]), 2:(max_age + 1)]

    # Set treatment initiation rate - DRTB
    NTDR_II_kappa <- NTDR_II_kap_main[match(current_year, NTDR_II_kap_main[, 1]), 2:(max_age + 1)]
    PTDR_II_kappa <- PTDR_II_kap_main[match(current_year, PTDR_II_kap_main[, 1]), 2:(max_age + 1)]
    NTDR_IN_kappa <- NTDR_IN_kap_main[match(current_year, NTDR_IN_kap_main[, 1]), 2:(max_age + 1)]
    PTDR_IN_kappa <- PTDR_IN_kap_main[match(current_year, PTDR_IN_kap_main[, 1]), 2:(max_age + 1)]

    # Set misdiagnosed treatment initiation rate
    NTDR_II_mdt <- NTDR_II_kap_mdt[match(current_year, NTDR_II_kap_mdt[, 1]), 2:(max_age + 1)]
    PTDR_II_mdt <- PTDR_II_kap_mdt[match(current_year, PTDR_II_kap_mdt[, 1]), 2:(max_age + 1)]
    NTDR_IN_nid <- NTDR_IN_kap_nid[match(current_year, NTDR_IN_kap_nid[, 1]), 2:(max_age + 1)]
    PTDR_IN_nid <- PTDR_IN_kap_nid[match(current_year, PTDR_IN_kap_nid[, 1]), 2:(max_age + 1)]

    # Select treatment success/failure rates
    # DSTB
    NTDS_psi    <- psi[match(current_year, psi[, 1]), 2] * (1 - exp(-(1/0.5) * dt)) #Exit timing constant
    NTDS_phi    <- (1 - psi[match(current_year, psi[, 1]), 2]) * (1 - exp(-(1/0.5) * dt))

    PTDS_psi    <- psi[match(current_year, psi[, 1]), 3] * (1 - exp(-(1/0.5) * dt))
    PTDS_phi    <- (1 - psi[match(current_year, psi[, 1]), 3]) * (1 - exp(-(1/0.5) * dt))

    # DRTB
    NTDR_II_psi <- psi[match(current_year, psi[, 1]), 3]
    NTDR_II_phi <- 1 - NTDR_II_psi

    NTDR_IN_psi <- psi[match(current_year, psi[, 1]), 3]
    NTDR_IN_phi <- 1 - NTDR_IN_psi

    PTDR_II_psi <- psi[match(current_year, psi[, 1]), 4]
    PTDR_II_phi <- 1 - PTDR_II_psi

    PTDR_IN_psi <- psi[match(current_year, psi[, 1]), 4]
    PTDR_IN_phi <- 1 - PTDR_IN_psi

    NTDR_mdt_exit <- (1 - exp(-(1/0.5) * dt))
    PTDR_mdt_exit <- (1 - exp(-(1/0.5) * dt))
    NTDR_nid_exit <- (1 - exp(-(1/0.5) * dt))
    PTDR_nid_exit <- (1 - exp(-(1/0.5) * dt))

    if (current_year < dst_start_year) {
      nt_mdt_loss  <- 1
      pt_mdt_loss  <- 1
      nt_nid_loss  <- 1
      pt_nid_loss  <- 1
    } else if (current_year >= dst_start_year) {
      nt_mdt_loss  <- 1 - ntii_tx_corr[current_year - year1 + 1]
      pt_mdt_loss  <- 1 - ptii_tx_corr[current_year - year1 + 1]
      nt_nid_loss  <- 1 - ntin_tx_corr[current_year - year1 + 1]
      pt_nid_loss  <- 1 - ptin_tx_corr[current_year - year1 + 1]
    }

    # Equation age range initialisation
    a1    <- NULL
    a0    <- NULL
    surv  <- NULL
    surv4 <- NULL
    surv3 <- NULL

    if (substep == 1) {

      # Start Step-1-of-year calculations

      # Birth model - add absolute number of births for this year into population
      # vector
      if (fert == 1) {
        la["S", i, 1] <- births[match(current_year, births[, 1]), 2]
      } else if (fert == 0) {
        la["S", i, 1] <- 0
      }

      # Calculate transmission terms - Lambda
      NTDS_lambda_raw <- colSums(-(myneta[1:4, 1:max_age]) * z * ((DS_Imatrix[i - 1, 1:4])/(psizematrix[i - 1, 1:4])))
      NTDS_lambda[i - 1, 1:max_age] <- t(DS_neta * (1 - exp(NTDS_lambda_raw)))
      PTDS_lambda <- NTDS_lambda

      NTDR_lambda_raw <- colSums(-(myneta[1:4, 1:max_age]) * z * ((DR_Imatrix[i - 1, 1:4])/(psizematrix[i - 1, 1:4])))
      NTDR_lambda[i - 1, 1:max_age] <- t(DR_neta * (1 - exp(NTDR_lambda_raw)))
      PTDR_lambda <- NTDR_lambda

      # Source and run the epi Equations
      rm(a1, a0)
      a1 <- 2:(max_age - 1)
      a0 <- 1:(max_age - 2)
      a4 <- 100
      a3 <- 99


rm(surv)
surv <- 1 - (u[a0])

# Susceptibles
la["S", i, a1] <-   (surv * la["S", i - 1, a0]) +
                   -(NTDS_lambda[i - 1, a0] * la["S", i - 1, a0]) +
                   -(NTDR_lambda[i - 1, a0] * la["S", i - 1, a0])

rm(surv)
surv <- 1 - (u[a0])

## Never Treated Drug Sensitive TB
## NTDS-Latent
la["NTDS_L", i, a1] <- (surv * la["NTDS_L", i - 1, a0]) +
                       ((1 - NTDS_p[a0]) * NTDS_lambda[i - 1, a0] * la["S", i - 1, a0]) +
                       ((1 - NTDS_p[a0]) * NTDS_lambda[i - 1, a0] * NTDS_x * la["NTDR_L", i - 1, a0]) +
                      -(NTDS_v[a0]  * la["NTDS_L", i - 1, a0]) +
                      -(NTDS_x * NTDS_lambda[i - 1, a0] * NTDS_p[a0] * la["NTDS_L", i - 1, a0]) +
                      -(NTDR_x * NTDR_lambda[i - 1, a0] * la["NTDS_L", i - 1, a0])

## Infectious NTDS
la["NTDS_II", i, a1] <- (surv * la["NTDS_II", i - 1, a0]) +
                        ## Conversion of Non-infectious to Infectious
                        (NTDS_omega  * la["NTDS_IN", i - 1, a0]) +
                        ## Deconstructed
                        (NTDS_p[a0] * NTDS_f[a0] * NTDS_lambda[i - 1, a0] * la["S", i - 1, a0]) +
                        (NTDS_p[a0] * NTDS_f[a0] * NTDS_lambda[i - 1, a0] * NTDS_x * la["NTDR_L", i - 1, a0]) +
                        (NTDS_p[a0] * NTDS_f[a0] * NTDS_lambda[i - 1, a0] * NTDS_x * la["NTDS_L", i - 1, a0]) +
                        ## End deconstruct
                        ## Reactivation from Latent
                        (NTDS_v[a0]  * NTDS_f[a0] * la["NTDS_L", i - 1, a0]) +
                        ## Reinfection of Resolved - DS
                        (NTDS_lambda[i - 1, a0] * NTDS_p[a0] * NTDS_f[a0] * NTDS_x * la["NTDS_R", i - 1, a0]) +
                        ## Reactivation from Resolved
                        (NTDS_f[a0] * NTDS_r[a0]  * la["NTDS_R", i - 1, a0]) +
                        ## Reinfection of Resolved - DR
                        (NTDS_lambda[i - 1, a0] * NTDS_p[a0] * NTDS_f[a0] * NTDS_x * la["NTDR_R", i - 1, a0]) +
                        ## Natural Cure
                       -(NTDS_II_n[a0]  * la["NTDS_II", i - 1, a0]) +
                        ## Case detection
                       -(NTDS_II_kappa[a0] * la["NTDS_II", i - 1, a0]) +
                        ## TB Death
                       -(NTDS_II_u[a0] * la["NTDS_II", i - 1, a0])

## Non-Infectious NTDS
la["NTDS_IN", i,a1] <- (surv * la["NTDS_IN", i-1,a0]) +
                       ## Deconstructed
                       (NTDS_p[a0] * (1 - NTDS_f[a0]) * NTDS_lambda[i - 1, a0] * la["S", i - 1, a0]) +
                       (NTDS_p[a0] * (1 - NTDS_f[a0]) * NTDS_lambda[i - 1, a0] * NTDS_x * la["NTDR_L", i - 1, a0]) +
                       (NTDS_p[a0] * (1 - NTDS_f[a0]) * NTDS_lambda[i - 1, a0] * NTDS_x * la["NTDS_L", i - 1, a0]) +
                       ## End deconstruct
                       ## Reactivation from Latent
                       (NTDS_v[a0]  * (1 - NTDS_f[a0]) * la["NTDS_L", i - 1, a0]) +
                       ## Reinfection of Resolved - DS
                       (NTDS_lambda[i - 1, a0] * NTDS_p[a0] * (1 - NTDS_f[a0]) * NTDS_x * la["NTDS_R", i - 1,a0]) +
                       ## Reactivation from Resolved
                       ((1 - NTDS_f[a0]) * NTDS_r[a0]  * la["NTDS_R", i-1,a0]) +
                       ## Reinfection of Resolved - DR
                       (NTDS_lambda[i - 1, a0] * NTDS_p[a0] * (1 - NTDS_f[a0]) * NTDS_x * la["NTDR_R", i - 1, a0]) +
                       ## Natural Cure
                      -(NTDS_IN_n[a0]  * la["NTDS_IN", i - 1, a0]) +
                       ## Case Detection
                      -(NTDS_IN_kappa[a0] * la["NTDS_IN", i - 1, a0]) +
                       ## TB Death
                      -(NTDS_IN_u[a0] * la["NTDS_IN", i - 1, a0]) +
                       ## Conversion from Non-infectious to Infectious
                      -(NTDS_omega  * la["NTDS_IN", i - 1, a0])

## NTDS on treatment
la["NTDS_T", i,a1]  <- (surv * la["NTDS_T", i-1,a0]) +
                       ## Treatment initiation from Infectious
                       (NTDS_II_kappa[a0]*la["NTDS_II", i-1,a0]) +
                       ## Treatment initiation from Non-infectious
                       (NTDS_IN_kappa[a0]*la["NTDS_IN", i-1,a0]) +
                       ## Resistance acquisition
                      -(NT_xi * la["NTDS_T", i-1,a0]) +
                       ## Treatment failures
                      -(NTDS_phi * la["NTDS_T", i-1,a0]) +
                       ## Treatment successes
                      -(NTDS_psi * la["NTDS_T", i-1,a0]) +
                       ## Death on treatment
                      -(NTDS_T_u[a0] * la["NTDS_T", i-1,a0])

# NTDS Resolved
la["NTDS_R", i, a1] <- (surv * la["NTDS_R", i - 1, a0]) +
                       ## Natural cure Infectious
                       (NTDS_II_n[a0]  * la["NTDS_II", i - 1, a0]) +
                       ## Natural cure Non-infectious
                       (NTDS_IN_n[a0]  * la["NTDS_IN", i - 1, a0]) +
                       ## Reinfection of Resolved - DS
                      -(NTDS_lambda[i - 1, a0] * NTDS_p[a0] * NTDS_x * la["NTDS_R", i - 1, a0]) +
                       ## Reactivation of Resolved
                      -(NTDS_r[a0]  * la["NTDS_R", i - 1, a0]) +
                       ## Reinfection of Resolved - DR
                      -(NTDR_lambda[i - 1, a0] * NTDR_p[a0] * NTDR_x * la["NTDS_R", i - 1, a0])

rm(surv4, surv3)
surv4 <- 1 - (u[a4] )
surv3 <- 1 - (u[a3] )

# Susceptibles
la["S", i, a4] <-      (surv4 * la["S", i - 1, a4]) +
                      -(NTDS_lambda[i - 1, a4] * la["S", i - 1, a4]) +
                      -(NTDR_lambda[i - 1, a4] * la["S", i - 1, a4]) +
# Accumulate
                       (surv3 * la["S", i - 1, a3]) +
                      -(NTDS_lambda[i - 1, a3] * la["S", i - 1, a3]) +
                      -(NTDR_lambda[i - 1, a3] * la["S", i - 1, a3])

# Never Treated Drug Sensitive TB
# NTDS-latent
la["NTDS_L", i, a4] <- (surv4 * la["NTDS_L", i - 1, a4]) +
                       ((1 - NTDS_p[a4]) * NTDS_lambda[i - 1, a4] * la["S", i - 1, a4]) +
                       ((1 - NTDS_p[a4]) * NTDS_lambda[i - 1, a4] * NTDS_x * la["NTDR_L", i - 1, a4]) +
                      -(NTDS_v[a4] * la["NTDS_L", i - 1, a4]) +
                      -(NTDS_x * NTDS_lambda[i - 1, a4] * NTDS_p[a4] * la["NTDS_L", i - 1, a4]) +
                      -(NTDR_x * NTDR_lambda[i - 1, a4] * la["NTDS_L", i - 1, a4]) +
                      #Accumulate
                       (surv3 * la["NTDS_L", i - 1, a3]) +
                       ((1 - NTDS_p[a3]) * NTDS_lambda[i - 1, a3] * la["S", i - 1, a3]) +
                       ((1 - NTDS_p[a3]) * NTDS_lambda[i - 1, a3] * NTDS_x * la["NTDR_L", i - 1, a3]) +
                      -(NTDS_v[a3] * la["NTDS_L", i - 1, a3]) +
                      -(NTDS_x * NTDS_lambda[i - 1, a3] * NTDS_p[a3] * la["NTDS_L", i - 1, a3]) +
                      -(NTDR_x * NTDR_lambda[i - 1, a3] * la["NTDS_L", i - 1, a3])

### Infectious DS-TB
la["NTDS_II", i, a4] <- (surv4 * la["NTDS_II", i - 1, a4]) +
                        (NTDS_omega * la["NTDS_IN", i - 1, a4]) +
                        (NTDS_p[a4] * NTDS_f[a4] * NTDS_lambda[i - 1, a4] * la["S", i - 1, a4]) +
                        (NTDS_p[a4] * NTDS_f[a4] * NTDS_lambda[i - 1, a4] * NTDS_x * la["NTDR_L", i - 1, a4]) +
                        (NTDS_p[a4] * NTDS_f[a4] * NTDS_lambda[i - 1, a4] * NTDS_x * la["NTDS_L", i - 1, a4]) +
                        (NTDS_v[a4] * NTDS_f[a4] * la["NTDS_L", i - 1, a4]) +
                        (NTDS_lambda[i - 1, a4] * NTDS_p[a4] * NTDS_f[a4] * NTDS_x * la["NTDS_R", i - 1, a4]) +
                        (NTDS_f[a4] * NTDS_r[a4] * la["NTDS_R", i - 1, a4]) +
                        (NTDS_lambda[i - 1, a4] * NTDS_p[a4] * NTDS_f[a4] * NTDS_x * la["NTDR_R", i - 1, a4]) +
                       -(NTDS_II_n[a4] * la["NTDS_II", i - 1, a4]) +
                       -(NTDS_II_kappa[a4] * la["NTDS_II", i - 1, a4]) +
                       -(NTDS_II_u[a4] * la["NTDS_II", i - 1, a4]) +
                       # Accumulate
                        (surv3 * la["NTDS_II", i - 1, a3]) +
                       (NTDS_omega * la["NTDS_IN", i - 1, a3]) +
                        (NTDS_p[a3] * NTDS_f[a3] * NTDS_lambda[i - 1, a3] * la["S", i - 1, a3]) +
                        (NTDS_p[a3] * NTDS_f[a3] * NTDS_lambda[i - 1, a3] * NTDS_x * la["NTDR_L", i - 1, a3]) +
                        (NTDS_p[a3] * NTDS_f[a3] * NTDS_lambda[i - 1, a3] * NTDS_x * la["NTDS_L", i - 1, a3]) +
                        (NTDS_v[a3] * NTDS_f[a3] * la["NTDS_L", i - 1, a3]) +
                        (NTDS_lambda[i - 1, a3] * NTDS_p[a3] * NTDS_f[a3] * NTDS_x * la["NTDS_R", i - 1, a3]) +
                        (NTDS_f[a3] * NTDS_r[a3] * la["NTDS_R", i - 1, a3]) +
                        (NTDS_lambda[i - 1, a3] * NTDS_p[a3] * NTDS_f[a3] * NTDS_x * la["NTDR_R", i - 1, a3]) +
                       -(NTDS_II_n[a3] * la["NTDS_II", i - 1, a3]) +
                       -(NTDS_II_kappa[a3] * la["NTDS_II", i - 1, a3]) +
                       -(NTDS_II_u[a3] * la["NTDS_II", i - 1, a3])

## Non-infectious DS-TB
la["NTDS_IN", i,a4] <-  (surv4 * la["NTDS_IN", i-1,a4]) +
                        (NTDS_p[a4] * (1 - NTDS_f[a4]) * NTDS_lambda[i - 1, a4] * la["S", i - 1, a4]) +
                        (NTDS_p[a4] * (1 - NTDS_f[a4]) * NTDS_lambda[i - 1, a4] * NTDS_x * la["NTDR_L", i - 1, a4]) +
                        (NTDS_p[a4] * (1 - NTDS_f[a4]) * NTDS_lambda[i - 1, a4] * NTDS_x * la["NTDS_L", i - 1, a4]) +
                        (NTDS_v[a4] * (1 - NTDS_f[a4]) * la["NTDS_L", i - 1, a4]) +
                        (NTDS_lambda[i - 1, a4] * NTDS_p[a4] * (1 - NTDS_f[a4]) * NTDS_x * la["NTDS_R", i-1,a4]) +
                        ((1 - NTDS_f[a4]) * NTDS_r[a4] * la["NTDS_R", i-1,a4]) +
                        (NTDS_lambda[i - 1, a4] * NTDS_p[a4] * (1 - NTDS_f[a4]) * NTDS_x * la["NTDR_R", i - 1, a4]) +
                       -(NTDS_IN_n[a4] * la["NTDS_IN", i - 1, a4]) +
                       -(NTDS_IN_kappa[a4] * la["NTDS_IN", i - 1, a4]) +
                       -(NTDS_IN_u[a4] * la["NTDS_IN", i - 1, a4]) +
                       -(NTDS_omega * la["NTDS_IN", i - 1, a4]) +
                       # Accumulate
                        (surv3 * la["NTDS_IN", i-1,a3]) +
                        (NTDS_p[a3] * (1 - NTDS_f[a3]) * NTDS_lambda[i - 1, a3] * la["S", i - 1, a3]) +
                        (NTDS_p[a3] * (1 - NTDS_f[a3]) * NTDS_lambda[i - 1, a3] * NTDS_x * la["NTDR_L", i - 1, a3]) +
                        (NTDS_p[a3] * (1 - NTDS_f[a3]) * NTDS_lambda[i - 1, a3] * NTDS_x * la["NTDS_L", i - 1, a3]) +
                        (NTDS_v[a3] * (1 - NTDS_f[a3]) * la["NTDS_L", i - 1, a3]) +
                        (NTDS_lambda[i - 1, a3] * NTDS_p[a3] * (1 - NTDS_f[a3]) * NTDS_x * la["NTDS_R", i-1,a3]) +
                        ((1 - NTDS_f[a3]) * NTDS_r[a3] * la["NTDS_R", i-1,a3]) +
                        (NTDS_lambda[i - 1, a3] * NTDS_p[a3] * (1 - NTDS_f[a3]) * NTDS_x * la["NTDR_R", i - 1, a3]) +
                       -(NTDS_IN_n[a3] * la["NTDS_IN", i - 1, a3]) +
                       -(NTDS_IN_kappa[a3] * la["NTDS_IN", i - 1, a3]) +
                       -(NTDS_IN_u[a3] * la["NTDS_IN", i - 1, a3]) +
                       -(NTDS_omega * la["NTDS_IN", i - 1, a3])

## NTDS On Treatment
la["NTDS_T", i,a4]  <- (surv4 * la["NTDS_T", i-1,a4]) +
                      (NTDS_II_kappa[a4]*la["NTDS_II", i-1,a4]) +
                      (NTDS_IN_kappa[a4]*la["NTDS_IN", i-1,a4]) -
                      (NT_xi * la["NTDS_T", i-1,a4]) -
                      (NTDS_phi * la["NTDS_T", i-1,a4]) -
                      (NTDS_psi * la["NTDS_T", i-1,a4]) -
                      (NTDS_T_u[a4] * la["NTDS_T", i-1,a4]) +
                      #Accumulate
                      (surv3 * la["NTDS_T", i-1,a3]) +
                      (NTDS_II_kappa[a3]*la["NTDS_II", i-1,a3]) +
                      (NTDS_IN_kappa[a3]*la["NTDS_IN", i-1,a3]) -
                      (NT_xi * la["NTDS_T", i-1,a3]) -
                      (NTDS_phi * la["NTDS_T", i-1,a3]) -
                      (NTDS_psi * la["NTDS_T", i-1,a3]) -
                      (NTDS_T_u[a3] * la["NTDS_T", i-1,a3])

# NTDS Resolved
la["NTDS_R", i, a4] <-  (surv4 * la["NTDS_R", i - 1, a4]) +
                        (NTDS_II_n[a4] * la["NTDS_II", i - 1, a4]) +
                        (NTDS_IN_n[a4] * la["NTDS_IN", i - 1, a4]) +
                       -(NTDS_lambda[i - 1, a4] * NTDS_p[a4] * NTDS_x * la["NTDS_R", i - 1, a4]) +
                       -(NTDS_r[a4] * la["NTDS_R", i - 1, a4]) +
                       -(NTDR_lambda[i - 1, a4] * NTDR_p[a4] * NTDR_x * la["NTDS_R", i - 1, a4]) +
                        # Accumulate
                        (surv3 * la["NTDS_R", i - 1, a3]) +
                        (NTDS_II_n[a3] * la["NTDS_II", i - 1, a3]) +
                        (NTDS_IN_n[a3] * la["NTDS_IN", i - 1, a3]) +
                       -(NTDS_lambda[i - 1, a3] * NTDS_p[a3] * NTDS_x * la["NTDS_R", i - 1, a3]) +
                       -(NTDS_r[a3] * la["NTDS_R", i - 1, a3]) +
                       -(NTDR_lambda[i - 1, a3] * NTDR_p[a3] * NTDR_x * la["NTDS_R", i - 1, a3])

rm(surv)
surv <- 1 - (u[a0] )

## Never Treated Drug Resistant

## NTDR-Latent
la["NTDR_L", i, a1] <- (surv * la["NTDR_L", i - 1, a0]) +
                       ## New Infections from Susceptible
                       ((1 - NTDR_p[a0]) * NTDR_lambda[i - 1, a0] * la["S", i - 1, a0]) +
                       ## New Infections from NTDS-Latent
                       ((1 - NTDR_p[a0]) * NTDR_lambda[i - 1, a0] * NTDR_x * la["NTDS_L", i - 1, a0]) +
                       ## Reactivation from Latent
                      -(NTDR_v[a0] * la["NTDR_L", i - 1, a0]) +
                       ## DRTB Reinfection of NTDR-Latent
                      -(NTDR_x * NTDR_lambda[i - 1, a0] * NTDR_p[a0] * la["NTDR_L", i - 1, a0]) +
                       ## DSTB Reinfection of NTDR-Latent
                      -(NTDR_x * NTDS_lambda[i - 1, a0] * la["NTDR_L", i - 1, a0])

## Infectious NTDR
la["NTDR_II", i,a1] <- (surv * la["NTDR_II", i - 1, a0]) +
                       ## Conversion of non-infectious to infectious
                       (NTDR_omega * la["NTDR_IN", i - 1, a0]) +
                       ## Deconstructed compound term - new infections in susceptble, NTDR-Latent and NTDS-Latent
                       (NTDR_p[a0] * NTDR_f[a0] * NTDR_lambda[i - 1, a0] * la["S", i - 1, a0]) +
                       (NTDR_p[a0] * NTDR_f[a0] * NTDR_lambda[i - 1, a0] * NTDR_x * la["NTDR_L", i - 1, a0]) +
                       (NTDR_p[a0] * NTDR_f[a0] * NTDR_lambda[i - 1, a0] * NTDR_x * la["NTDS_L", i - 1, a0]) +
                       ## End deconstruct
                       ## Reactivation from Latent
                       (NTDR_v[a0] * NTDR_f[a0] * la["NTDR_L", i - 1, a0]) +
                       ## Reinfection of Resolved - DR
                       (NTDR_lambda[i - 1, a0] * NTDR_p[a0] * NTDR_f[a0] * NTDR_x * la["NTDR_R", i - 1, a0]) +
                       ## Reactivation from Resolved
                       (NTDR_f[a0] * NTDR_r[a0] * la["NTDR_R", i - 1, a0]) +
                       ## Reinfection of Resolved - DS
                       (NTDR_lambda[i - 1, a0] * NTDR_p[a0] * NTDR_f[a0] * NTDR_x * la["NTDS_R", i - 1, a0]) +
                       ## Natural cure
                      -(NTDR_II_n[a0] * la["NTDR_II", i - 1, a0]) +
                       ## Case detection
                      -(NTDR_II_kappa[a0] * la["NTDR_II", i - 1, a0]) +
                       ## Misdiagnosis and Treatment
                      -(NTDR_II_mdt[a0] * la["NTDR_II", i - 1, a0]) +
                       ## TB death
                      -(NTDR_II_u[a0] * la["NTDR_II", i - 1, a0])

## Non-Infectious NTDR
la["NTDR_IN", i,a1] <- (surv * la["NTDR_IN", i - 1, a0]) +
                       # Deconstructed
                       (NTDR_p[a0] * (1 - NTDR_f[a0]) * NTDR_lambda[i - 1, a0] * la["S", i - 1, a0]) +
                       (NTDR_p[a0] * (1 - NTDR_f[a0]) * NTDR_lambda[i - 1, a0] * NTDR_x * la["NTDR_L", i - 1, a0]) +
                       (NTDR_p[a0] * (1 - NTDR_f[a0]) * NTDR_lambda[i - 1, a0] * NTDR_x * la["NTDS_L", i - 1, a0]) +
                       # End deconstruct
                       ## Reactivation from Latent
                       (NTDR_v[a0] * (1 - NTDR_f[a0]) * la["NTDR_L", i - 1, a0]) +
                       ## Reinfection of Resolved - DR
                       (NTDR_lambda[i - 1, a0] * NTDR_p[a0] * (1 - NTDR_f[a0]) * NTDR_x * la["NTDR_R", i - 1, a0]) +
                       ## Reactivation from Resolved
                       ((1 - NTDR_f[a0]) * NTDR_r[a0] * la["NTDR_R", i - 1, a0]) +
                       ## Reinfection of Resolved - DS
                       (NTDR_lambda[i - 1, a0] * NTDR_p[a0] * (1 - NTDR_f[a0]) * NTDR_x * la["NTDS_R", i - 1, a0]) +
                       ## Natural cure
                      -(NTDR_IN_n[a0] * la["NTDR_IN", i - 1, a0]) +
                       ## Case detection
                      -(NTDR_IN_kappa[a0] * la["NTDR_IN", i - 1, a0]) +
                       ## Misdiagnosis and Treatment
                      -(NTDR_IN_nid[a0] * la["NTDR_IN", i - 1, a0]) +
                       ## TB death
                      -(NTDR_IN_u[a0] * la["NTDR_IN", i - 1, a0]) +
                       ## Conversion from non-infectious to infectious
                      -(NTDR_omega * la["NTDR_IN", i - 1, a0])

## Treatment - NTDR
## Treatment Succeeding - Originating from Infectious DR-TB
la["NTDR_T_IIpsi", i, a1] <- (surv * la["NTDR_T_IIpsi", i - 1, a0]) +
                             (NTDR_II_kappa[a0] * NTDR_II_psi * la["NTDR_II", i - 1, a0]) +
                            -(NTDR_T_IIpsi_tau * la["NTDR_T_IIpsi", i - 1, a0]) +
                            -(NTDR_T_u[a0] * la["NTDR_T_IIpsi", i - 1, a0])

## Treatment Failing - Originating from Infectious DR-TB
la["NTDR_T_IIphi", i, a1] <- (surv * la["NTDR_T_IIphi", i - 1, a0]) +
                             (NTDR_II_kappa[a0] * NTDR_II_phi * la["NTDR_II", i - 1, a0]) +
                            -(NTDR_T_IIphi_tau * la["NTDR_T_IIphi", i - 1, a0]) +
                            -(NTDR_II_u[a0] * la["NTDR_T_IIphi", i - 1, a0])

## Treatment Succeeding - Originating from Non-Infectious DR-TB
la["NTDR_T_INpsi", i, a1] <- (surv * la["NTDR_T_INpsi", i - 1, a0]) +
                             (NTDR_IN_kappa[a0] * NTDR_IN_psi * la["NTDR_IN", i - 1, a0]) +
                            -(NTDR_T_INpsi_tau * la["NTDR_T_INpsi", i - 1, a0]) +
                            -(NTDR_T_u[a0] * la["NTDR_T_INpsi", i - 1, a0])

## Treatment Failing - Originating from Non-Infectious DR-TB
la["NTDR_T_INphi", i, a1] <- (surv * la["NTDR_T_INphi", i - 1, a0]) +
                             (NTDR_IN_kappa[a0] * NTDR_IN_phi * la["NTDR_IN", i - 1, a0]) +
                            -(NTDR_T_INphi_tau * la["NTDR_T_INphi", i - 1, a0]) +
                            -(NTDR_IN_u[a0] * la["NTDR_T_INphi", i - 1, a0])

## Resolved - DR TB
la["NTDR_R", i, a1] <- (surv * la["NTDR_R", i - 1, a0]) +
                       ## Natural cure - Infectious
                       (NTDR_II_n[a0] * la["NTDR_II", i - 1, a0]) +
                       ## Natural cure - Non-infectious
                       (NTDR_IN_n[a0] * la["NTDR_IN", i - 1, a0]) +
                       ## Reinfection of Resolved - to DR
                      -(NTDR_lambda[i - 1, a0] * NTDR_p[a0] * NTDR_x * la["NTDR_R", i - 1, a0]) +
                       ## Reactivation of Resolved
                      -(NTDR_r[a0] * la["NTDR_R", i - 1, a0]) +
                       ## Reinfection of Resolved - to DS
                      -(NTDS_lambda[i - 1, a0] * NTDS_p[a0] * NTDS_x * la["NTDR_R", i - 1, a0])

## Misdiagnosed and Treated - Infectious
la["NTDR_mdt", i, a1] <-  (surv * la["NTDR_mdt", i - 1, a0]) +
                          # Misdiagnosed and treated
                          (NTDR_II_mdt[a0] * la["NTDR_II", i - 1, a0]) +
                          # Exit outwards
                          -(la["NTDR_mdt", i - 1, a0] * NTDR_mdt_exit) +
                          -(NTDR_mdt_u[a0] * la["NTDR_mdt", i - 1, a0])

## Misdiagnosed and Treated - Npn-Infectious
la["NTDR_nid", i, a1] <-   (surv * la["NTDR_nid", i - 1, a0]) +
                           # Misdiagnosed and treated
                          (NTDR_IN_nid[a0] * la["NTDR_IN", i - 1, a0]) +
                           # Exit outwards
                          -(la["NTDR_nid", i - 1, a0] * NTDR_nid_exit) +
                          -(NTDR_nid_u[a0] * la["NTDR_nid", i - 1, a0])

rm(surv4, surv3)
surv4 <- 1 - (u[a4])
surv3 <- 1 - (u[a3])

## Never Treated Drug Resistant

## NTDR-Latent
la["NTDR_L", i, a4] <-  (surv4 * la["NTDR_L", i - 1, a4]) +
                        (surv3 * la["NTDR_L", i - 1, a3]) +
                        ## New Infections from Susceptible
                        ((1 - NTDR_p[a4]) * NTDR_lambda[i - 1, a4] * la["S", i - 1, a4]) +
                        ((1 - NTDR_p[a3]) * NTDR_lambda[i - 1, a3] * la["S", i - 1, a3]) +
                        ## New Infections from NTDS-Latent
                        ((1 - NTDR_p[a4]) * NTDR_lambda[i - 1, a4] * NTDR_x * la["NTDS_L", i - 1, a4]) +
                        ((1 - NTDR_p[a3]) * NTDR_lambda[i - 1, a3] * NTDR_x * la["NTDS_L", i - 1, a3]) +
                        ## Reactivation from Latent
                       -(NTDR_v[a4] * la["NTDR_L", i - 1, a4]) +
                       -(NTDR_v[a3] * la["NTDR_L", i - 1, a3]) +
                        ## DRTB Reinfection of NTDR-Latent
                       -(NTDR_x * NTDR_lambda[i - 1, a4] * NTDR_p[a4] * la["NTDR_L", i - 1, a4]) +
                       -(NTDR_x * NTDR_lambda[i - 1, a3] * NTDR_p[a3] * la["NTDR_L", i - 1, a3]) +
                        ## DSTB Reinfection of NTDR-Latent
                       -(NTDR_x * NTDS_lambda[i - 1, a4] * la["NTDR_L", i - 1, a4]) +
                       -(NTDR_x * NTDS_lambda[i - 1, a3] * la["NTDR_L", i - 1, a3])

### Infectious NTDR
la["NTDR_II", i,a4] <-  (surv4 * la["NTDR_II", i - 1, a4]) +
                        (surv3 * la["NTDR_II", i - 1, a3]) +
                        ## Conversion of non-infectious to infectious
                        (NTDR_omega * la["NTDR_IN", i - 1, a4]) +
                        (NTDR_omega * la["NTDR_IN", i - 1, a3]) +
                        ## Deconstructed
                        (NTDR_p[a4] * NTDR_f[a4] * NTDR_lambda[i - 1, a4] * la["S", i - 1, a4]) +
                        (NTDR_p[a4] * NTDR_f[a4] * NTDR_lambda[i - 1, a4] * NTDR_x * la["NTDR_L", i - 1, a4]) +
                        (NTDR_p[a4] * NTDR_f[a4] * NTDR_lambda[i - 1, a4] * NTDR_x * la["NTDS_L", i - 1, a4]) +
                        (NTDR_p[a3] * NTDR_f[a3] * NTDR_lambda[i - 1, a3] * la["S", i - 1, a3]) +
                        (NTDR_p[a3] * NTDR_f[a3] * NTDR_lambda[i - 1, a3] * NTDR_x * la["NTDR_L", i - 1, a3]) +
                        (NTDR_p[a3] * NTDR_f[a3] * NTDR_lambda[i - 1, a3] * NTDR_x * la["NTDS_L", i - 1, a3]) +
                        ## End deconstruct
                        ## Reactivation from Latent
                        (NTDR_v[a4] * NTDR_f[a4] * la["NTDR_L", i - 1, a4]) +
                        (NTDR_v[a3] * NTDR_f[a3] * la["NTDR_L", i - 1, a3]) +
                        ## Reinfection of Resolved - DR
                        (NTDR_lambda[i - 1, a4] * NTDR_p[a4] * NTDR_f[a4] * NTDR_x * la["NTDR_R", i - 1, a4]) +
                        (NTDR_lambda[i - 1, a3] * NTDR_p[a3] * NTDR_f[a3] * NTDR_x * la["NTDR_R", i - 1, a3]) +
                        ## Reactivation from Resolved
                        (NTDR_f[a4] * NTDR_r[a4] * la["NTDR_R", i - 1, a4]) +
                        (NTDR_f[a3] * NTDR_r[a3] * la["NTDR_R", i - 1, a3]) +
                        ## Reinfection of Resolved - DS
                        (NTDR_lambda[i - 1, a4] * NTDR_p[a4] * NTDR_f[a4] * NTDR_x * la["NTDS_R", i - 1, a4]) +
                        (NTDR_lambda[i - 1, a3] * NTDR_p[a3] * NTDR_f[a3] * NTDR_x * la["NTDS_R", i - 1, a3]) +
                        ## Natural cure
                       -(NTDR_II_n[a4] * la["NTDR_II", i - 1, a4]) +
                       -(NTDR_II_n[a3] * la["NTDR_II", i - 1, a3]) +
                        ## Case detection
                       -(NTDR_II_kappa[a4] * la["NTDR_II", i - 1, a4]) +
                       -(NTDR_II_kappa[a3] * la["NTDR_II", i - 1, a3]) +
                        ## Misdiagnosis and Treatment
                       -(NTDR_II_mdt[a4] * la["NTDR_II", i - 1, a4]) +
                       -(NTDR_II_mdt[a3] * la["NTDR_II", i - 1, a3]) +
                        ## TB death
                       -(NTDR_II_u[a4] * la["NTDR_II", i - 1, a4]) +
                       -(NTDR_II_u[a3] * la["NTDR_II", i - 1, a3])

### Non-Infectious NTDR
la["NTDR_IN", i,a4] <-  (surv4 * la["NTDR_IN", i - 1, a4]) +
                        (surv3 * la["NTDR_IN", i - 1, a3]) +
                        # Deconstructed
                        (NTDR_p[a4] * (1 - NTDR_f[a4]) * NTDR_lambda[i - 1, a4] * la["S", i - 1, a4]) +
                        (NTDR_p[a4] * (1 - NTDR_f[a4]) * NTDR_lambda[i - 1, a4] * NTDR_x * la["NTDR_L", i - 1, a4]) +
                        (NTDR_p[a4] * (1 - NTDR_f[a4]) * NTDR_lambda[i - 1, a4] * NTDR_x * la["NTDS_L", i - 1, a4]) +
                        (NTDR_p[a3] * (1 - NTDR_f[a3]) * NTDR_lambda[i - 1, a3] * la["S", i - 1, a3]) +
                        (NTDR_p[a3] * (1 - NTDR_f[a3]) * NTDR_lambda[i - 1, a3] * NTDR_x * la["NTDR_L", i - 1, a3]) +
                        (NTDR_p[a3] * (1 - NTDR_f[a3]) * NTDR_lambda[i - 1, a3] * NTDR_x * la["NTDS_L", i - 1, a3]) +
                        # End deconstruct
                        ## Reactivation from Latent
                        (NTDR_v[a4] * (1 - NTDR_f[a4]) * la["NTDR_L", i - 1, a4]) +
                        (NTDR_v[a3] * (1 - NTDR_f[a3]) * la["NTDR_L", i - 1, a3]) +
                        ## Reinfection of Resolved - DR
                        (NTDR_lambda[i - 1, a4] * NTDR_p[a4] * (1 - NTDR_f[a4]) * NTDR_x * la["NTDR_R", i - 1, a4]) +
                        (NTDR_lambda[i - 1, a3] * NTDR_p[a3] * (1 - NTDR_f[a3]) * NTDR_x * la["NTDR_R", i - 1, a3]) +
                        ## Reactivation from Resolved
                        ((1 - NTDR_f[a4]) * NTDR_r[a4] * la["NTDR_R", i - 1, a4]) +
                        ((1 - NTDR_f[a3]) * NTDR_r[a3] * la["NTDR_R", i - 1, a3]) +
                        ## Reinfection of Resolved - DS
                        (NTDR_lambda[i - 1, a4] * NTDR_p[a4] * (1 - NTDR_f[a4]) * NTDR_x * la["NTDS_R", i - 1, a4]) +
                        (NTDR_lambda[i - 1, a3] * NTDR_p[a3] * (1 - NTDR_f[a3]) * NTDR_x * la["NTDS_R", i - 1, a3]) +
                        ## Natural cure
                       -(NTDR_IN_n[a4] * la["NTDR_IN", i - 1, a4]) +
                       -(NTDR_IN_n[a3] * la["NTDR_IN", i - 1, a3]) +
                        ## Case detection
                       -(NTDR_IN_kappa[a4] * la["NTDR_IN", i - 1, a4]) +
                       -(NTDR_IN_kappa[a3] * la["NTDR_IN", i - 1, a3]) +
                        ## Misdiagnosis and Treatment
                       -(NTDR_IN_nid[a4] * la["NTDR_IN", i - 1, a4]) +
                       -(NTDR_IN_nid[a3] * la["NTDR_IN", i - 1, a3]) +
                        ## TB death
                       -(NTDR_IN_u[a4] * la["NTDR_IN", i - 1, a4]) +
                       -(NTDR_IN_u[a3] * la["NTDR_IN", i - 1, a3]) +
                        ## Conversion from non-infectious to infectious
                       -(NTDR_omega * la["NTDR_IN", i - 1, a4]) +
                       -(NTDR_omega * la["NTDR_IN", i - 1, a3])

## Treatment - NTDR
## Treatment Succeeding - Originating from Infectious DR-TB
la["NTDR_T_IIpsi", i, a4] <- (surv4 * la["NTDR_T_IIpsi", i - 1, a4]) +
                             (surv3 * la["NTDR_T_IIpsi", i - 1, a3]) +
                             (NTDR_II_kappa[a4] * NTDR_II_psi * la["NTDR_II", i - 1, a4]) +
                             (NTDR_II_kappa[a3] * NTDR_II_psi * la["NTDR_II", i - 1, a3]) +
                            -(NTDR_T_IIpsi_tau * la["NTDR_T_IIpsi", i - 1, a4]) +
                            -(NTDR_T_IIpsi_tau * la["NTDR_T_IIpsi", i - 1, a3]) +
                            -(NTDR_T_u[a4] * la["NTDR_T_IIpsi", i - 1, a4]) +
                            -(NTDR_T_u[a3] * la["NTDR_T_IIpsi", i - 1, a3])

## Treatment Failing - Originating from Infectious DR-TB
la["NTDR_T_IIphi", i, a4] <- (surv4 * la["NTDR_T_IIphi", i - 1, a4]) +
                             (surv3 * la["NTDR_T_IIphi", i - 1, a3]) +
                             (NTDR_II_kappa[a4] * NTDR_II_phi * la["NTDR_II", i - 1, a4]) +
                             (NTDR_II_kappa[a3] * NTDR_II_phi * la["NTDR_II", i - 1, a3]) +
                            -(NTDR_T_IIphi_tau * la["NTDR_T_IIphi", i - 1, a4]) +
                            -(NTDR_T_IIphi_tau * la["NTDR_T_IIphi", i - 1, a3]) +
                            -(NTDR_II_u[a4] * la["NTDR_T_IIphi", i - 1, a4]) +
                            -(NTDR_II_u[a3] * la["NTDR_T_IIphi", i - 1, a3])

## Treatment Succeeding - Originating from Non-Infectious DR-TB
la["NTDR_T_INpsi", i, a4] <- (surv4 * la["NTDR_T_INpsi", i - 1, a4]) +
                             (surv3 * la["NTDR_T_INpsi", i - 1, a3]) +
                             (NTDR_IN_kappa[a4] * NTDR_IN_psi * la["NTDR_IN", i - 1, a4]) +
                             (NTDR_IN_kappa[a3] * NTDR_IN_psi * la["NTDR_IN", i - 1, a3]) +
                            -(NTDR_T_INpsi_tau * la["NTDR_T_INpsi", i - 1, a4]) +
                            -(NTDR_T_INpsi_tau * la["NTDR_T_INpsi", i - 1, a3]) +
                            -(NTDR_T_u[a4] * la["NTDR_T_INpsi", i - 1, a4]) +
                            -(NTDR_T_u[a3] * la["NTDR_T_INpsi", i - 1, a3])

## Treatment Failing - Originating from Non-Infectious DR-TB
la["NTDR_T_INphi", i, a4] <- (surv4 * la["NTDR_T_INphi", i - 1, a4]) +
                             (surv3 * la["NTDR_T_INphi", i - 1, a3]) +
                             (NTDR_IN_kappa[a4] * NTDR_IN_phi * la["NTDR_IN", i - 1, a4]) +
                             (NTDR_IN_kappa[a3] * NTDR_IN_phi * la["NTDR_IN", i - 1, a3]) +
                            -(NTDR_T_INphi_tau * la["NTDR_T_INphi", i - 1, a4]) +
                            -(NTDR_T_INphi_tau * la["NTDR_T_INphi", i - 1, a3]) +
                            -(NTDR_IN_u[a4] * la["NTDR_T_INphi", i - 1, a4]) +
                            -(NTDR_IN_u[a3] * la["NTDR_T_INphi", i - 1, a3])

## Resolved - DR TB
la["NTDR_R", i, a4] <- (surv4 * la["NTDR_R", i - 1, a4]) +
                       (surv3 * la["NTDR_R", i - 1, a3]) +
                       ## Natural cure - Infectious
                       (NTDR_II_n[a4] * la["NTDR_II", i - 1, a4]) +
                       (NTDR_II_n[a3] * la["NTDR_II", i - 1, a3]) +
                       ## Natural cure - Non-infectious
                       (NTDR_IN_n[a4] * la["NTDR_IN", i - 1, a4]) +
                       (NTDR_IN_n[a3] * la["NTDR_IN", i - 1, a3]) +
                       ## Reinfection of Resolved - to DR
                      -(NTDR_lambda[i - 1, a4] * NTDR_p[a4] * NTDR_x * la["NTDR_R", i - 1, a4]) +
                      -(NTDR_lambda[i - 1, a3] * NTDR_p[a3] * NTDR_x * la["NTDR_R", i - 1, a3]) +
                       ## Reactivation of Resolved
                      -(NTDR_r[a4] * la["NTDR_R", i - 1, a4]) +
                      -(NTDR_r[a3] * la["NTDR_R", i - 1, a3]) +
                       ## Reinfection of Resolved - to DS
                      -(NTDS_lambda[i - 1, a4] * NTDS_p[a4] * NTDS_x * la["NTDR_R", i - 1, a4]) +
                      -(NTDS_lambda[i - 1, a3] * NTDS_p[a3] * NTDS_x * la["NTDR_R", i - 1, a3])

## Misdiagnosed and Treated - Infectious
la["NTDR_mdt", i, a4] <-  (surv4 * la["NTDR_mdt", i - 1, a4]) +
                          # Misdiagnosed and treated
                          (NTDR_II_mdt[a4] * la["NTDR_II", i - 1, a4]) +
                          # Exit outwards
                          -(la["NTDR_mdt", i - 1, a4] * NTDR_mdt_exit) +
                          -(NTDR_mdt_u[a4] * la["NTDR_mdt", i - 1, a4]) +
                          (surv3 * la["NTDR_mdt", i - 1, a3]) +
                          # Misdiagnosed and treated
                          (NTDR_II_mdt[a3] * la["NTDR_II", i - 1, a3]) +
                          # Exit outwards
                          -(la["NTDR_mdt", i - 1, a3] * NTDR_mdt_exit) +
                          -(NTDR_mdt_u[a3] * la["NTDR_mdt", i - 1, a3])

## Misdiagnosed and Treated - Non-Infectious
la["NTDR_nid", i, a4] <-  (surv4 * la["NTDR_nid", i - 1, a4]) +
                          # Misdiagnosed and treated
                          (NTDR_IN_nid[a4] * la["NTDR_IN", i - 1, a4]) +
                          # Exit outwards
                          -(la["NTDR_nid", i - 1, a4] * NTDR_nid_exit) +
                          -(NTDR_nid_u[a4] * la["NTDR_nid", i - 1, a4]) +
                          (surv3 * la["NTDR_mdt", i - 1, a3]) +
                          # Misdiagnosed and treated
                          (NTDR_IN_nid[a3] * la["NTDR_IN", i - 1, a3]) +
                          # Exit outwards
                          -(la["NTDR_nid", i - 1, a3] * NTDR_nid_exit) +
                          -(NTDR_nid_u[a3] * la["NTDR_nid", i - 1, a3])

rm(surv)
surv <- 1 - (u[a0] )

## Previously Treated Drug Sensitive TB

## Infectious PTDS
la["PTDS_II", i, a1] <- (surv * la["PTDS_II", i - 1, a0]) +
                        ## Conversion of non-infectious to infectious
                        (PTDS_omega * la["PTDS_IN", i - 1, a0]) +
                        ## Reinfection of Resolved - DS
                        (PTDS_lambda[i - 1, a0] * PTDS_p[a0] * PTDS_f[a0] * PTDS_x * la["PTDS_R", i - 1, a0]) +
                        ## Reactivation from Resolved
                        (PTDS_f[a0] * PTDS_r[a0] * la["PTDS_R", i - 1, a0]) +
                        ## Reinfection of Resolved - DR
                        (PTDS_lambda[i - 1, a0] * PTDS_p[a0] * PTDS_f[a0] * PTDS_x * la["PTDR_R", i - 1, a0]) +
                        ## Natural cure
                       -(PTDS_II_n[a0] * la["PTDS_II", i - 1, a0]) +
                        ## Case detection
                       -(PTDS_II_kappa[a0] * la["PTDS_II", i - 1, a0]) +
                        ## TB death
                       -(PTDS_II_u[a0] * la["PTDS_II", i - 1, a0])

## Non-Infectious PTDS
la["PTDS_IN", i, a1] <- (surv * la["PTDS_IN", i - 1, a0]) +
                        ## NTDS Treatment Failures
                        (NTDS_phi * la["NTDS_T", i - 1, a0]) +
                        ## PTDS Treatment Failures
                        (PTDS_phi * la["PTDS_T", i - 1, a0]) +
                        ## Reinfection of Resolved - DS
                        (PTDS_lambda[i - 1, a0] * PTDS_p[a0] * (1 - PTDS_f[a0]) * PTDS_x * la["PTDS_R", i - 1, a0]) +
                        ## Reactivation from Resolved
                        ((1 - PTDS_f[a0]) * PTDS_r[a0] * la["PTDS_R", i - 1, a0]) +
                        ## Reinfection of Resolved - DR
                        (PTDS_lambda[i - 1, a0] * PTDS_p[a0] * (1 - PTDS_f[a0]) * PTDS_x * la["PTDR_R", i - 1, a0]) +
                        ## Natural cure
                       -(PTDS_IN_n[a0] * la["PTDS_IN", i - 1, a0]) +
                        ## Case detection
                       -(PTDS_IN_kappa[a0] * la["PTDS_IN", i - 1, a0]) +
                        ## TB death
                       -(PTDS_IN_u[a0] * la["PTDS_IN", i - 1, a0]) +
                        ## Conversion from non-infectious to infectious
                       -(PTDS_omega * la["PTDS_IN", i - 1, a0])

## PTDS on treatment
la["PTDS_T", i, a1] <-  (surv * la["PTDS_T", i - 1, a0]) +
                        ## Treatment initiation from Infectious
                        (PTDS_II_kappa[a0] * la["PTDS_II", i - 1, a0]) +
                        ## Treatment initiation from Non-infectious
                        (PTDS_IN_kappa[a0] * la["PTDS_IN", i - 1, a0]) +
                        ## Resistance acquisition
                       -(PT_xi * la["PTDS_T", i - 1, a0]) +
                        ## Treatment failures
                       -(PTDS_phi * la["PTDS_T", i - 1, a0]) +
                        ## Treatment successes
                       -(PTDS_psi * la["PTDS_T", i - 1, a0]) +
                        ## Death on treatment
                       -(PTDS_T_u[a0] * la["PTDS_T", i - 1, a0])

## PTDS Resolved
la["PTDS_R", i, a1] <- (surv * la["PTDS_R", i - 1, a0]) +
                       ## Natural cure - Infectious
                       (PTDS_II_n[a0] * la["PTDS_II", i - 1, a0]) +
                       ## Natural cure - Non-infectious
                       (PTDS_IN_n[a0] * la["PTDS_IN", i - 1, a0]) +
                       ## NTDS Treatment Success
                       (NTDS_psi * la["NTDS_T", i - 1, a0]) +
                       ## PTDS Treatment Success
                       (PTDS_psi * la["PTDS_T", i - 1, a0]) +
                       ## Reinfection of Resolved - to DS
                      -(PTDS_lambda[i - 1, a0] * PTDS_p[a0] * PTDS_x * la["PTDS_R", i - 1, a0]) +
                       ## Reactivation of Resolved
                      -(PTDS_r[a0] * la["PTDS_R", i - 1, a0]) +
                       ## Reinfection of Resolved - to DR
                      -(PTDR_lambda[i - 1, a0] * PTDR_p[a0] * PTDR_x * la["PTDS_R", i - 1, a0])

rm(surv4, surv3)
surv4 <- 1 - (u[a4] )
surv3 <- 1 - (u[a3] )

## Previously Treated for TB : Drug Sensitive
### Infectious DS-TB
la["PTDS_II", i, a4] <- (surv4 * la["PTDS_II", i - 1, a4]) +
                        (PTDS_omega * la["PTDS_IN", i - 1, a4]) +
                        (PTDS_lambda[i - 1, a4] * PTDS_p[a4] * PTDS_f[a4] * PTDS_x * la["PTDS_R", i - 1, a4]) +
                        (PTDS_f[a4] * PTDS_r[a4] * la["PTDS_R", i - 1, a4]) +
                        (PTDS_lambda[i - 1, a4] * PTDS_p[a4] * PTDS_f[a4] * PTDS_x * la["PTDR_R", i - 1, a4]) +
                       -(PTDS_II_n[a4] * la["PTDS_II", i - 1, a4]) +
                       -(PTDS_II_kappa[a4] * la["PTDS_II", i - 1, a4]) +
                       -(PTDS_II_u[a4] * la["PTDS_II", i - 1, a4]) +
                        # Accumulate
                        (surv3 * la["PTDS_II", i - 1, a3]) +
                        (PTDS_omega * la["PTDS_IN", i - 1, a3]) +
                        (PTDS_lambda[i - 1, a3] * PTDS_p[a3] * PTDS_f[a3] * PTDS_x * la["PTDS_R", i - 1, a3]) +
                        (PTDS_f[a3] * PTDS_r[a3] * la["PTDS_R", i - 1, a3]) +
                        (PTDS_lambda[i - 1, a3] * PTDS_p[a3] * PTDS_f[a3] * PTDS_x * la["PTDR_R", i - 1, a3]) +
                       -(PTDS_II_n[a3] * la["PTDS_II", i - 1, a3]) +
                       -(PTDS_II_kappa[a3] * la["PTDS_II", i - 1, a3]) +
                       -(PTDS_II_u[a3] * la["PTDS_II", i - 1, a3])

### Non-Infectious DS-TB
la["PTDS_IN", i, a4] <- (surv4 * la["PTDS_IN", i - 1, a4]) +
                        (NTDS_phi * la["NTDS_T", i - 1, a4]) +
                        (PTDS_phi * la["PTDS_T", i - 1, a4]) +
                        (PTDS_lambda[i - 1, a4] * PTDS_p[a4] * (1 - PTDS_f[a4]) * PTDS_x * la["PTDS_R", i - 1, a4]) +
                        ((1 - PTDS_f[a4]) * PTDS_r[a4] * la["PTDS_R", i - 1, a4]) +
                        (PTDS_lambda[i - 1, a4] * PTDS_p[a4] * (1 - PTDS_f[a4]) * PTDS_x * la["PTDR_R", i - 1, a4]) +
                       -(PTDS_IN_n[a4] * la["PTDS_IN", i - 1, a4]) +
                       -(PTDS_IN_kappa[a4] * la["PTDS_IN", i - 1, a4]) +
                       -(PTDS_IN_u[a4] * la["PTDS_IN", i - 1, a4]) +
                       -(PTDS_omega * la["PTDS_IN", i - 1, a4]) +
                       # Accumulate
                        (surv3 * la["PTDS_IN", i - 1, a3]) +
                        (NTDS_phi * la["NTDS_T", i - 1, a3]) +
                        (PTDS_phi * la["PTDS_T", i - 1, a3]) +
                        (PTDS_lambda[i - 1, a3] * PTDS_p[a3] * (1 - PTDS_f[a3]) * PTDS_x * la["PTDS_R", i - 1, a3]) +
                        ((1 - PTDS_f[a3]) * PTDS_r[a3] * la["PTDS_R", i - 1, a3]) +
                        (PTDS_lambda[i - 1, a3] * PTDS_p[a3] * (1 - PTDS_f[a3]) * PTDS_x * la["PTDR_R", i - 1, a3]) +
                       -(PTDS_IN_n[a3] * la["PTDS_IN", i - 1, a3]) +
                       -(PTDS_IN_kappa[a3] * la["PTDS_IN", i - 1, a3]) +
                       -(PTDS_IN_u[a3] * la["PTDS_IN", i - 1, a3]) +
                       -(PTDS_omega * la["PTDS_IN", i - 1, a3])

### Treatment DS-TB
la["PTDS_T", i, a4] <- (surv4 * la["PTDS_T", i - 1, a4]) +
                      (PTDS_II_kappa[a4] * la["PTDS_II", i - 1, a4]) +
                      (PTDS_IN_kappa[a4] * la["PTDS_IN", i - 1, a4]) -
                      (PT_xi * la["PTDS_T", i - 1, a4]) -
                      (PTDS_phi * la["PTDS_T", i - 1, a4]) -
                      (PTDS_psi * la["PTDS_T", i - 1, a4]) -
                      (PTDS_T_u[a4] * la["PTDS_T", i - 1, a4]) +
                      #Accumulate
                      (surv3 * la["PTDS_T", i - 1, a3]) +
                      (PTDS_II_kappa[a3] * la["PTDS_II", i - 1, a3]) +
                      (PTDS_IN_kappa[a3] * la["PTDS_IN", i - 1, a3]) -
                      (PT_xi * la["PTDS_T", i - 1, a3]) -
                      (PTDS_phi * la["PTDS_T", i - 1, a3]) -
                      (PTDS_psi * la["PTDS_T", i - 1, a3]) -
                      (PTDS_T_u[a3] * la["PTDS_T", i - 1, a3])

### PTDS Resolved
la["PTDS_R", i, a4] <-  (surv4 * la["PTDS_R", i - 1, a4]) +
                        (PTDS_II_n[a4] * la["PTDS_II", i - 1, a4]) +
                        (PTDS_IN_n[a4] * la["PTDS_IN", i - 1, a4]) +
                        ((NTDS_psi * la["NTDS_T", i - 1, a4]) +
                        (PTDS_psi * la["PTDS_T", i - 1, a4])) +
                       -(PTDS_lambda[i - 1, a4] * PTDS_p[a4] * PTDS_x * la["PTDS_R", i - 1, a4]) +
                       -(PTDS_r[a4] * la["PTDS_R", i - 1, a4]) +
                       -(PTDR_lambda[i - 1, a4] * PTDR_p[a4] * PTDR_x * la["PTDS_R", i - 1, a4]) +
                        # Accumulate
                        (surv3 * la["PTDS_R", i - 1, a3]) +
                        (PTDS_II_n[a3] * la["PTDS_II", i - 1, a3]) +
                        (PTDS_IN_n[a3] * la["PTDS_IN", i - 1, a3]) +
                        (NTDS_psi * la["NTDS_T", i - 1, a3]) +
                        (PTDS_psi * la["PTDS_T", i - 1, a3]) +
                       -(PTDS_lambda[i - 1, a3] * PTDS_p[a3] * PTDS_x * la["PTDS_R", i - 1, a3]) +
                       -(PTDS_r[a3] * la["PTDS_R", i - 1, a3]) +
                       -(PTDR_lambda[i - 1, a3] * PTDR_p[a3] * PTDR_x * la["PTDS_R", i - 1, a3])

rm(surv)
surv <- 1 - (u[a0] )

## Previously Treated Drug Resistant
## PTDR-Latent
la["PTDR_L", i, a1] <- (surv * la["PTDR_L", i - 1, a0]) +
                       ((PT_xi * (1 - PTDR_p[a0])) * (la["PTDS_T", i - 1, a0])) +
                       ((NT_xi * (1 - PTDR_p[a0])) * (la["NTDS_T", i - 1, a0])) +
                      -(PTDR_v[a0] * la["PTDR_L", i - 1, a0])

## Infectious PTDR
la["PTDR_II", i, a1] <- (surv * la["PTDR_II", i - 1, a0]) +
                        ## Conversion of non-infectious to infectious
                        (PTDR_omega * la["PTDR_IN", i - 1, a0]) +
                        ## Reactivation from Latent
                        (PTDR_v[a0] * PTDR_f[a0] * la["PTDR_L", i - 1, a0]) +
                        ## Reinfection of Resolved - DR
                        (PTDR_lambda[i - 1, a0] * PTDR_p[a0] * PTDR_f[a0] * PTDR_x * la["PTDR_R", i - 1, a0]) +
                        ## Reactivation from Resolved
                        (PTDR_f[a0] * (PTDR_r[a0] ) * la["PTDR_R", i - 1, a0]) +
                        ## Reinfection of Resolved - DS
                        (PTDR_lambda[i - 1, a0] * PTDR_p[a0] * PTDR_f[a0] * PTDR_x * la["PTDS_R", i - 1, a0]) +
                        ## MDR Acquisition
                        (NT_xi * PTDR_p[a0] * PTDR_f[a0] * la["NTDS_T", i - 1, a0] * nt_mdt_loss) +
                        (PT_xi * PTDR_p[a0] * PTDR_f[a0] * la["PTDS_T", i - 1, a0] * pt_mdt_loss) +
                        ## PTDR Treatment FAILURES - Infectious
                        (PTDR_T_IIphi_tau * la["PTDR_T_IIphi", i - 1, a0]) +
                        ## NTDR Treatment FAILURES - Infectious
                        (NTDR_T_IIphi_tau * la["NTDR_T_IIphi", i - 1, a0]) +
                        # Misdiagnosed - re-entry from Infectious MDT compartment
                        (la["NTDR_mdt", i - 1, a0] * NTDR_mdt_exit * nt_mdt_loss) +
                        (la["PTDR_mdt", i - 1, a0] * PTDR_mdt_exit * pt_mdt_loss) +
                        ## Natural cure
                       -(PTDR_II_n[a0] * la["PTDR_II", i - 1, a0]) +
                        ## Case detection
                       -(PTDR_II_kappa[a0] * la["PTDR_II", i - 1, a0]) +
                        ## Misdiagnosis and treatment
                       -(PTDR_II_mdt[a0] * la["PTDR_II", i - 1, a0]) +
                        ## TB death
                       -(PTDR_II_u[a0] * la["PTDR_II", i - 1, a0])

## Non-Infectious PTDR
la["PTDR_IN", i, a1] <- (surv * la["PTDR_IN", i - 1, a0]) +
                        ## Reactivation from Latent
                        (PTDR_v[a0] * (1 - PTDR_f[a0]) * la["PTDR_L", i - 1, a0]) +
                        ## Reinfection of Resolved - DS
                        (PTDR_lambda[i - 1, a0] * PTDR_p[a0] * (1 - PTDR_f[a0]) * PTDR_x * la["PTDS_R", i - 1, a0]) +
                        ## Reactivation from Resolved
                        ((1 - PTDR_f[a0]) * PTDR_r[a0] * la["PTDR_R", i - 1, a0]) +
                        ## Reinfection of Resolved - DR
                        (PTDR_lambda[i - 1, a0] * PTDR_p[a0] * (1 - PTDR_f[a0]) * PTDR_x * la["PTDR_R", i - 1, a0]) +
                        ## MDR Acquisition
                        (PT_xi * PTDR_p[a0] * (1 - PTDR_f[a0]) * la["PTDS_T", i - 1, a0] * pt_nid_loss) +
                        (NT_xi * PTDR_p[a0] * (1 - PTDR_f[a0]) * la["NTDS_T", i - 1, a0] * nt_nid_loss) +
                        ## PTDR Treatment FAILURES - Non-Infectious
                        (PTDR_T_INphi_tau * la["PTDR_T_INphi", i - 1, a0]) +
                        ## NTDR Treatment FAILURES - Non-Infectious
                        (NTDR_T_INphi_tau * la["NTDR_T_INphi", i - 1, a0]) +
                        # Misdiagnosed - re-entry from Non-Infectious NID compartment
                        (la["NTDR_nid", i - 1, a0] * NTDR_nid_exit * nt_nid_loss) +
                        (la["PTDR_nid", i - 1, a0] * PTDR_nid_exit * pt_nid_loss) +
                        ## Natural cure
                       -(PTDR_IN_n[a0] *la["PTDR_IN", i - 1, a0]) +
                        ## Case detection
                       -(PTDR_IN_kappa[a0] * la["PTDR_IN", i - 1, a0]) +
                        ## Misdiagnosis and treatment
                       -(PTDR_IN_nid[a0] * la["PTDR_IN", i - 1, a0]) +
                        ## TB death
                       -(PTDR_IN_u[a0] * la["PTDR_IN", i - 1, a0]) +
                        ## Conversion from non-infectious to infectious
                       -(PTDR_omega * la["PTDR_IN", i - 1, a0])

## PTDR on treatment
## Treatment Succeeding - Originating from Infectious DR-TB
la["PTDR_T_IIpsi", i, a1] <-  (surv * la["PTDR_T_IIpsi", i - 1, a0]) +
                              (PTDR_II_kappa[a0] * PTDR_II_psi * la["PTDR_II", i - 1, a0]) +
                              # Misdiagnosed - re-entry from MDR compartment
                              (la["NTDR_mdt", i - 1, a0] * NTDR_mdt_exit * (1 - nt_mdt_loss) * PTDR_II_psi) +
                              (la["PTDR_mdt", i - 1, a0] * PTDR_mdt_exit * (1 - pt_mdt_loss) * PTDR_II_psi) +
                              (NT_xi * PTDR_p[a0] * PTDR_f[a0] * la["NTDS_T", i - 1, a0] * (1 - nt_mdt_loss) * PTDR_II_psi) +
                              (PT_xi * PTDR_p[a0] * PTDR_f[a0] * la["PTDS_T", i - 1, a0] * (1 - pt_mdt_loss) * PTDR_II_psi) +
                             -(PTDR_T_IIpsi_tau * la["PTDR_T_IIpsi", i - 1, a0]) +
                             -(PTDR_T_u[a0] * la["PTDR_T_IIpsi", i - 1, a0])

## Treatment Failing - Originating from Infectious DR-TB
la["PTDR_T_IIphi", i, a1] <-  (surv * la["PTDR_T_IIphi", i - 1, a0]) +
                              (PTDR_II_kappa[a0] * PTDR_II_phi * la["PTDR_II", i - 1, a0]) +
                              # Misdiagnosed - re-entry from MDR compartment
                              (la["NTDR_mdt", i - 1, a0] * NTDR_mdt_exit * (1 - nt_mdt_loss) * PTDR_II_phi) +
                              (la["PTDR_mdt", i - 1, a0] * PTDR_mdt_exit * (1 - pt_mdt_loss) * PTDR_II_phi) +
                              (NT_xi * PTDR_p[a0] * PTDR_f[a0] * la["NTDS_T", i - 1, a0] * (1 - nt_mdt_loss) * PTDR_II_phi) +
                              (PT_xi * PTDR_p[a0] * PTDR_f[a0] * la["PTDS_T", i - 1, a0] * (1 - pt_mdt_loss) * PTDR_II_phi) +
                             -(PTDR_T_IIphi_tau * la["PTDR_T_IIphi", i - 1, a0]) +
                             -(PTDR_II_u[a0] * la["PTDR_T_IIphi", i - 1, a0])

## Treatment Succeeding - Originating from Non-Infectious DR-TB
la["PTDR_T_INpsi", i, a1] <-  (surv * la["PTDR_T_INpsi", i - 1, a0]) +
                              (PTDR_IN_kappa[a0] * PTDR_IN_psi * la["PTDR_IN", i - 1, a0]) +
                              # Misdiagnosed - re-entry from MDR compartment
                              (la["NTDR_nid", i - 1, a0] * NTDR_nid_exit * (1 - nt_nid_loss) * PTDR_IN_psi) +
                              (la["PTDR_nid", i - 1, a0] * PTDR_nid_exit * (1 - pt_nid_loss) * PTDR_IN_psi) +
                              (NT_xi * PTDR_p[a0] * (1 - PTDR_f[a0]) * la["NTDS_T", i - 1, a0] * (1 - nt_nid_loss) * PTDR_IN_psi) +
                              (PT_xi * PTDR_p[a0] * (1 - PTDR_f[a0]) * la["PTDS_T", i - 1, a0] * (1 - pt_nid_loss) * PTDR_IN_psi) +
                             -(PTDR_T_INpsi_tau * la["PTDR_T_INpsi", i - 1, a0]) +
                             -(PTDR_T_u[a0] * la["PTDR_T_INpsi", i - 1, a0])

## Treatment Failing - Originating from Non-Infectious DR-TB
la["PTDR_T_INphi", i, a1] <-  (surv * la["PTDR_T_INphi", i - 1, a0]) +
                              (PTDR_IN_kappa[a0] * PTDR_IN_phi * la["PTDR_IN", i - 1, a0]) +
                              # Misdiagnosed - re-entry from MDR compartment
                              (la["NTDR_nid", i - 1, a0] * NTDR_nid_exit * (1 - nt_nid_loss) * PTDR_IN_phi) +
                              (la["PTDR_nid", i - 1, a0] * PTDR_nid_exit * (1 - pt_nid_loss) * PTDR_IN_phi) +
                              (NT_xi * PTDR_p[a0] * (1 - PTDR_f[a0]) * la["NTDS_T", i - 1, a0] * (1 - nt_nid_loss) * PTDR_IN_phi) +
                              (PT_xi * PTDR_p[a0] * (1 - PTDR_f[a0]) * la["PTDS_T", i - 1, a0] * (1 - pt_nid_loss) * PTDR_IN_phi) +
                             -(PTDR_T_INphi_tau * la["PTDR_T_INphi", i - 1, a0]) +
                             -(PTDR_II_u[a0] * la["PTDR_T_INphi", i - 1, a0])

## PTDR Resolved
la["PTDR_R", i, a1] <- (surv * la["PTDR_R", i - 1, a0]) +
                       ## Natural cure - Infectious
                       (PTDR_II_n[a0] * la["PTDR_II", i - 1, a0]) +
                       ## Natural cure - Non-infectious
                       (PTDR_IN_n[a0] * la["PTDR_IN", i - 1, a0]) +
                       ## NTDR Treatment Success - Infectious
                       (NTDR_T_IIpsi_tau * la["NTDR_T_IIpsi", i - 1, a0]) +
                       ## NTDR Treatment Success - Non-infectious
                       (NTDR_T_INpsi_tau * la["NTDR_T_INpsi", i - 1, a0]) +
                       ## PTDR Treatment Success - Infectious
                       (PTDR_T_INpsi_tau * la["PTDR_T_INpsi", i - 1, a0]) +
                       ## PTDR Treatment Success - Non-infectious
                       (PTDR_T_IIpsi_tau * la["PTDR_T_IIpsi", i - 1, a0]) +
                       ## Reinfection of Resolved - to DS
                      -(PTDS_lambda[i - 1, a0] * PTDS_p[a0] * PTDS_x * la["PTDR_R", i - 1, a0]) +
                       ## Reactivation of Resolved
                      -(PTDR_r[a0] * la["PTDR_R", i - 1, a0]) +
                      ## Reinfection of Resolved - to DR
                      -(PTDR_lambda[i - 1, a0] * PTDR_p[a0] * PTDR_x * la["PTDR_R", i - 1, a0])

## Misdiagnosed and Treated - Infectious
la["PTDR_mdt", i, a1] <-  (surv * la["PTDR_mdt", i - 1, a0]) +
                          # Misdiagnosed and treated
                          (PTDR_II_mdt[a0] * la["PTDR_II", i - 1, a0]) +
                          # Exit outwards
                          -(la["PTDR_mdt", i - 1, a0] * PTDR_mdt_exit) +
                          -(PTDR_mdt_u[a0] * la["PTDR_mdt", i - 1, a0])

## Misdiagnosed and Treated - Npn-Infectious
la["PTDR_nid", i, a1] <-   (surv * la["PTDR_nid", i - 1, a0]) +
                           # Misdiagnosed and treated
                          (PTDR_IN_nid[a0] * la["PTDR_IN", i - 1, a0]) +
                           # Exit outwards
                          -(la["PTDR_nid", i - 1, a0] * PTDR_nid_exit) +
                          -(PTDR_nid_u[a0] * la["PTDR_nid", i - 1, a0])

rm(surv4, surv3)
surv4 <- 1 - (u[a4] )
surv3 <- 1 - (u[a3] )

## Previously Treated for TB: Drug Resistant
### DR-LTBI
la["PTDR_L", i, a4] <-  (surv4 * la["PTDR_L", i - 1, a4]) +
                        (PT_xi * (1 - PTDR_p[a4]) * la["PTDS_T", i - 1, a4]) +
                        (NT_xi * (1 - PTDR_p[a4]) * la["NTDS_T", i - 1, a4]) +
                       -(PTDR_v[a4] * la["PTDR_L", i - 1, a4]) +
                        # Accumulate
                        (surv3 * la["PTDR_L", i - 1, a3]) +
                        (PT_xi * (1 - PTDR_p[a3]) * la["PTDS_T", i - 1, a3]) +
                        (NT_xi * (1 - PTDR_p[a3]) * la["NTDS_T", i - 1, a3]) +
                       -(PTDR_v[a3] * la["PTDR_L", i - 1, a3])

### Infectious DR-TB
la["PTDR_II", i, a4] <- (surv4 * la["PTDR_II", i - 1, a4]) +
                        (PTDR_omega * la["PTDR_IN", i - 1, a4]) +
                        (PTDR_v[a4] * PTDR_f[a4] * la["PTDR_L", i - 1, a4]) +
                        (PTDR_lambda[i - 1, a4] * PTDR_p[a4] * PTDR_f[a4] * PTDR_x * la["PTDR_R", i - 1, a4]) +
                        (PTDR_f[a4] * (PTDR_r[a4] ) * la["PTDR_R", i - 1, a4]) +
                        (PTDR_lambda[i - 1, a4] * PTDR_p[a4] * PTDR_f[a4] * PTDR_x * la["PTDS_R", i - 1, a4]) +
                        (NT_xi * PTDR_p[a4] * PTDR_f[a4] * la["NTDS_T", i - 1, a4] * nt_mdt_loss) +
                        (PT_xi * PTDR_p[a4] * PTDR_f[a4] * la["PTDS_T", i - 1, a4] * pt_mdt_loss) +
                        (PTDR_T_IIphi_tau * la["PTDR_T_IIphi", i - 1, a4]) +
                        (NTDR_T_IIphi_tau * la["NTDR_T_IIphi", i - 1, a4]) +
                        (la["NTDR_mdt", i - 1, a4] * NTDR_mdt_exit * nt_mdt_loss) +
                        (la["PTDR_mdt", i - 1, a4] * PTDR_mdt_exit * pt_mdt_loss) +
                       -(PTDR_II_n[a4] * la["PTDR_II", i - 1, a4]) +
                       -(PTDR_II_kappa[a4] * la["PTDR_II", i - 1, a4]) +
                       -(PTDR_II_mdt[a4] * la["PTDR_II", i - 1, a4]) +
                       -(PTDR_II_u[a4] * la["PTDR_II", i - 1, a4]) +
                        # Accumulate
                        (surv3 * la["PTDR_II", i - 1, a3]) +
                        (PTDR_omega * la["PTDR_IN", i - 1, a3]) +
                        (PTDR_v[a3] * PTDR_f[a3] * la["PTDR_L", i - 1, a3]) +
                        (PTDR_lambda[i - 1, a3] * PTDR_p[a3] * PTDR_f[a3] * PTDR_x * la["PTDR_R", i - 1, a3]) +
                        (PTDR_f[a3] * (PTDR_r[a3] ) * la["PTDR_R", i - 1, a3]) +
                        (PTDR_lambda[i - 1, a3] * PTDR_p[a3] * PTDR_f[a3] * PTDR_x * la["PTDS_R", i - 1, a3]) +
                        (NT_xi * PTDR_p[a3] * PTDR_f[a3] * la["NTDS_T", i - 1, a3] * nt_mdt_loss) +
                        (PT_xi * PTDR_p[a3] * PTDR_f[a3] * la["PTDS_T", i - 1, a3] * pt_mdt_loss) +
                        (PTDR_T_IIphi_tau * la["PTDR_T_IIphi", i - 1, a3]) +
                        (NTDR_T_IIphi_tau * la["NTDR_T_IIphi", i - 1, a3]) +
                        (la["NTDR_mdt", i - 1, a3] * NTDR_mdt_exit * nt_mdt_loss) +
                        (la["PTDR_mdt", i - 1, a3] * PTDR_mdt_exit * pt_mdt_loss) +
                       -(PTDR_II_n[a3] * la["PTDR_II", i - 1, a3]) +
                       -(PTDR_II_kappa[a3] * la["PTDR_II", i - 1, a3]) +
                       -(PTDR_II_mdt[a3] * la["PTDR_II", i - 1, a3]) +
                       -(PTDR_II_u[a3] * la["PTDR_II", i - 1, a3])

### Non-Infectious DR-TB
la["PTDR_IN", i, a4] <- (surv4 * la["PTDR_IN", i - 1, a4]) +
                        (PTDR_v[a4] * (1 - PTDR_f[a4]) * la["PTDR_L", i - 1, a4]) +
                        (PTDR_lambda[i - 1, a4] * PTDR_p[a4] * (1 - PTDR_f[a4]) * PTDR_x * la["PTDR_R", i - 1, a4]) +
                        ((1 - PTDR_f[a4]) * PTDR_r[a4] * la["PTDR_R", i - 1, a4]) +
                        (PTDR_lambda[i - 1, a4] * PTDR_p[a4] * (1 - PTDR_f[a4]) * PTDR_x * la["PTDS_R", i - 1, a4]) +
                        (PT_xi * PTDR_p[a4] * (1 - PTDR_f[a4]) * la["PTDS_T", i - 1, a4] * pt_nid_loss) +
                        (NT_xi * PTDR_p[a4] * (1 - PTDR_f[a4]) * la["NTDS_T", i - 1, a4] * nt_nid_loss) +
                        (PTDR_T_INphi_tau * la["PTDR_T_INphi", i - 1, a4]) +
                        (NTDR_T_INphi_tau * la["NTDR_T_INphi", i - 1, a4]) +
                        (la["NTDR_nid", i - 1, a4] * NTDR_nid_exit * nt_nid_loss) +
                        (la["PTDR_nid", i - 1, a4] * PTDR_nid_exit * pt_nid_loss) +
                       -(PTDR_IN_n[a4] *la["PTDR_IN", i - 1, a4]) +
                       -(PTDR_IN_kappa[a4] * la["PTDR_IN", i - 1, a4]) +
                       -(PTDR_IN_nid[a4] * la["PTDR_IN", i - 1, a4]) +
                       -(PTDR_IN_u[a4] * la["PTDR_IN", i - 1, a4]) +
                       -(PTDR_omega * la["PTDR_IN", i - 1, a4]) +
                        # Accumulate
                        (surv3 * la["PTDR_IN", i - 1, a3]) +
                        (PTDR_v[a3] * (1 - PTDR_f[a3]) * la["PTDR_L", i - 1, a3]) +
                        (PTDR_lambda[i - 1, a3] * PTDR_p[a3] * (1 - PTDR_f[a3]) * PTDR_x * la["PTDR_R", i - 1, a3]) +
                        ((1 - PTDR_f[a3]) * PTDR_r[a3] * la["PTDR_R", i - 1, a3]) +
                        (PTDR_lambda[i - 1, a3] * PTDR_p[a3] * (1 - PTDR_f[a3]) * PTDR_x * la["PTDS_R", i - 1, a3]) +
                        (PT_xi * PTDR_p[a3] * (1 - PTDR_f[a3]) * la["PTDS_T", i - 1, a3] * pt_nid_loss) +
                        (NT_xi * PTDR_p[a3] * (1 - PTDR_f[a3]) * la["NTDS_T", i - 1, a3] * nt_nid_loss) +
                        (PTDR_T_INphi_tau * la["PTDR_T_INphi", i - 1, a3]) +
                        (NTDR_T_INphi_tau * la["NTDR_T_INphi", i - 1, a3]) +
                        (la["NTDR_nid", i - 1, a3] * NTDR_nid_exit * nt_nid_loss) +
                        (la["PTDR_nid", i - 1, a3] * PTDR_nid_exit * pt_nid_loss) +
                       -(PTDR_IN_n[a3] *la["PTDR_IN", i - 1, a3]) +
                       -(PTDR_IN_kappa[a3] * la["PTDR_IN", i - 1, a3]) +
                       -(PTDR_IN_nid[a3] * la["PTDR_IN", i - 1, a3]) +
                       -(PTDR_IN_u[a3] * la["PTDR_IN", i - 1, a3]) +
                       -(PTDR_omega * la["PTDR_IN", i - 1, a3])

### Treatment Succeeding - Originating from Infectious DR-TB
la["PTDR_T_IIpsi", i, a4] <- (surv4 * la["PTDR_T_IIpsi", i - 1, a4]) +
                      (PTDR_II_kappa[a4] * PTDR_II_psi * la["PTDR_II", i - 1, a4]) +
                      (la["NTDR_mdt", i - 1, a4] * NTDR_mdt_exit * (1 - nt_mdt_loss) * PTDR_II_psi) +
                      (la["PTDR_mdt", i - 1, a4] * PTDR_mdt_exit * (1 - pt_mdt_loss) * PTDR_II_psi) +
                      (NT_xi * PTDR_p[a4] * PTDR_f[a4] * la["NTDS_T", i - 1, a4] * (1 - nt_mdt_loss) * PTDR_II_psi) +
                      (PT_xi * PTDR_p[a4] * PTDR_f[a4] * la["PTDS_T", i - 1, a4] * (1 - pt_mdt_loss) * PTDR_II_psi) +
                     -(PTDR_T_IIpsi_tau * la["PTDR_T_IIpsi", i - 1, a4]) +
                     -(PTDR_T_u[a4] * la["PTDR_T_IIpsi", i - 1, a4]) +
                      #Accumulate
                      (surv3 * la["PTDR_T_IIpsi", i - 1, a3]) +
                      (PTDR_II_kappa[a3] * PTDR_II_psi * la["PTDR_II", i - 1, a3]) +
                      (la["NTDR_mdt", i - 1, a3] * NTDR_mdt_exit * (1 - nt_mdt_loss) * PTDR_II_psi) +
                      (la["PTDR_mdt", i - 1, a3] * PTDR_mdt_exit * (1 - pt_mdt_loss) * PTDR_II_psi) +
                      (NT_xi * PTDR_p[a3] * PTDR_f[a3] * la["NTDS_T", i - 1, a3] * (1 - nt_mdt_loss) * PTDR_II_psi) +
                      (PT_xi * PTDR_p[a3] * PTDR_f[a3] * la["PTDS_T", i - 1, a3] * (1 - pt_mdt_loss) * PTDR_II_psi) +
                     -(PTDR_T_IIpsi_tau * la["PTDR_T_IIpsi", i - 1, a3]) +
                     -(PTDR_T_u[a3] * la["PTDR_T_IIpsi", i - 1, a3])

### Treatment Failing - Originating from Infectious DR-TB
la["PTDR_T_IIphi", i, a4] <- (surv4 * la["PTDR_T_IIphi", i - 1, a4]) +
                      (PTDR_II_kappa[a4] * PTDR_II_phi * la["PTDR_II", i - 1, a4]) +
                      (la["NTDR_mdt", i - 1, a4] * NTDR_mdt_exit * (1 - nt_mdt_loss) * PTDR_II_phi) +
                      (la["PTDR_mdt", i - 1, a4] * PTDR_mdt_exit * (1 - pt_mdt_loss) * PTDR_II_phi) +
                      (NT_xi * PTDR_p[a4] * PTDR_f[a4] * la["NTDS_T", i - 1, a4] * (1 - nt_mdt_loss) * PTDR_II_phi) +
                      (PT_xi * PTDR_p[a4] * PTDR_f[a4] * la["PTDS_T", i - 1, a4] * (1 - pt_mdt_loss) * PTDR_II_phi) +
                     -(PTDR_T_IIphi_tau * la["PTDR_T_IIphi", i - 1, a4]) +
                     -(PTDR_II_u[a4] * la["PTDR_T_IIphi", i - 1, a4]) +
                      #Accumulate
                      (surv3 * la["PTDR_T_IIphi", i - 1, a3]) +
                      (PTDR_II_kappa[a3] * PTDR_II_phi * la["PTDR_II", i - 1, a3]) +
                      (la["NTDR_mdt", i - 1, a3] * NTDR_mdt_exit * (1 - nt_mdt_loss) * PTDR_II_phi) +
                      (la["PTDR_mdt", i - 1, a3] * PTDR_mdt_exit * (1 - pt_mdt_loss) * PTDR_II_phi) +
                      (NT_xi * PTDR_p[a3] * PTDR_f[a3] * la["NTDS_T", i - 1, a3] * (1 - nt_mdt_loss) * PTDR_II_phi) +
                      (PT_xi * PTDR_p[a3] * PTDR_f[a3] * la["PTDS_T", i - 1, a3] * (1 - pt_mdt_loss) * PTDR_II_phi) +
                     -(PTDR_T_IIphi_tau * la["PTDR_T_IIphi", i - 1, a3]) +
                     -(PTDR_II_u[a3] * la["PTDR_T_IIphi", i - 1, a3])

### Treatment Succeeding - Originating from Non-Infectious DR-TB
la["PTDR_T_INpsi", i, a4] <- (surv4 * la["PTDR_T_INpsi", i - 1, a4]) +
                      (PTDR_IN_kappa[a4] * PTDR_IN_psi * la["PTDR_IN", i - 1, a4]) +
                      (la["NTDR_nid", i - 1, a4] * NTDR_nid_exit * (1 - nt_nid_loss) * PTDR_IN_psi) +
                      (la["PTDR_nid", i - 1, a4] * PTDR_nid_exit * (1 - pt_nid_loss) * PTDR_IN_psi) +
                      (NT_xi * PTDR_p[a4] * (1 - PTDR_f[a4]) * la["NTDS_T", i - 1, a4] * (1 - nt_nid_loss) * PTDR_IN_psi) +
                      (PT_xi * PTDR_p[a4] * (1 - PTDR_f[a4]) * la["PTDS_T", i - 1, a4] * (1 - pt_nid_loss) * PTDR_IN_psi) +
                     -(PTDR_T_INpsi_tau * la["PTDR_T_INpsi", i - 1, a4]) +
                     -(PTDR_T_u[a4] * la["PTDR_T_INpsi", i - 1, a4]) +
                      #Accumulate
                      (surv3 * la["PTDR_T_INpsi", i - 1, a3]) +
                      (PTDR_IN_kappa[a3] * PTDR_IN_psi * la["PTDR_IN", i - 1, a3]) +
                      (la["NTDR_nid", i - 1, a3] * NTDR_nid_exit * (1 - nt_nid_loss) * PTDR_IN_psi) +
                      (la["PTDR_nid", i - 1, a3] * PTDR_nid_exit * (1 - pt_nid_loss) * PTDR_IN_psi) +
                      (NT_xi * PTDR_p[a3] * (1 - PTDR_f[a3]) * la["NTDS_T", i - 1, a3] * (1 - nt_nid_loss) * PTDR_IN_psi) +
                      (PT_xi * PTDR_p[a3] * (1 - PTDR_f[a3]) * la["PTDS_T", i - 1, a3] * (1 - pt_nid_loss) * PTDR_IN_psi) +
                     -(PTDR_T_INpsi_tau * la["PTDR_T_INpsi", i - 1, a3]) +
                     -(PTDR_T_u[a3] * la["PTDR_T_INpsi", i - 1, a3])

### Treatment Failing - Originating from Non-Infectious DR-TB
la["PTDR_T_INphi", i, a4] <- (surv4 * la["PTDR_T_INphi", i - 1, a4]) +
                      (PTDR_IN_kappa[a4] * PTDR_IN_phi * la["PTDR_IN", i - 1, a4]) +
                      (la["NTDR_nid", i - 1, a4] * NTDR_nid_exit * (1 - nt_nid_loss) * PTDR_IN_phi) +
                      (la["PTDR_nid", i - 1, a4] * PTDR_nid_exit * (1 - pt_nid_loss) * PTDR_IN_phi) +
                      (NT_xi * PTDR_p[a4] * (1 - PTDR_f[a4]) * la["NTDS_T", i - 1, a4] * (1 - nt_nid_loss) * PTDR_IN_phi) +
                      (PT_xi * PTDR_p[a4] * (1 - PTDR_f[a4]) * la["PTDS_T", i - 1, a4] * (1 - pt_nid_loss) * PTDR_IN_phi) +
                     -(PTDR_T_INphi_tau * la["PTDR_T_INphi", i - 1, a4]) +
                     -(PTDR_II_u[a4] * la["PTDR_T_INphi", i - 1, a4]) +
                      #Accumulate
                      (surv3 * la["PTDR_T_INphi", i - 1, a3]) +
                      (PTDR_IN_kappa[a3] * PTDR_IN_phi * la["PTDR_IN", i - 1, a3]) +
                      (la["NTDR_nid", i - 1, a3] * NTDR_nid_exit * (1 - nt_nid_loss) * PTDR_IN_phi) +
                      (la["PTDR_nid", i - 1, a3] * PTDR_nid_exit * (1 - pt_nid_loss) * PTDR_IN_phi) +
                      (NT_xi * PTDR_p[a3] * (1 - PTDR_f[a3]) * la["NTDS_T", i - 1, a3] * (1 - nt_nid_loss) * PTDR_IN_phi) +
                      (PT_xi * PTDR_p[a3] * (1 - PTDR_f[a3]) * la["PTDS_T", i - 1, a3] * (1 - pt_nid_loss) * PTDR_IN_phi) +
                     -(PTDR_T_INphi_tau * la["PTDR_T_INphi", i - 1, a3]) +
                     -(PTDR_II_u[a3] * la["PTDR_T_INphi", i - 1, a3])

### Resolved - DR TB
la["PTDR_R", i, a4] <- (surv4 * la["PTDR_R", i - 1, a4]) +
                       (PTDR_II_n[a4] * la["PTDR_II", i - 1, a4]) +
                       (PTDR_IN_n[a4] * la["PTDR_IN", i - 1, a4]) +
                       (NTDR_T_IIpsi_tau * la["NTDR_T_IIpsi", i - 1, a4]) +
                       (NTDR_T_INpsi_tau * la["NTDR_T_INpsi", i - 1, a4]) +
                       (PTDR_T_INpsi_tau * la["PTDR_T_INpsi", i - 1, a4]) +
                       (PTDR_T_IIpsi_tau * la["PTDR_T_IIpsi", i - 1, a4]) +
                      -(PTDR_lambda[i - 1, a4] * PTDR_p[a4] * PTDR_x * la["PTDR_R", i - 1, a4]) +
                      -(PTDR_r[a4] * la["PTDR_R", i - 1, a4]) +
                      -(PTDS_lambda[i - 1, a4] * PTDS_p[a4] * PTDS_x * la["PTDR_R", i - 1, a4]) +
                       # Accumulate
                       (surv3 * la["PTDR_R", i - 1, a3]) +
                       (PTDR_II_n[a3] * la["PTDR_II", i - 1, a3]) +
                       (PTDR_IN_n[a3] * la["PTDR_IN", i - 1, a3]) +
                       (NTDR_T_IIpsi_tau * la["NTDR_T_IIpsi", i - 1, a3]) +
                       (NTDR_T_INpsi_tau * la["NTDR_T_INpsi", i - 1, a3]) +
                       (PTDR_T_INpsi_tau * la["PTDR_T_INpsi", i - 1, a3]) +
                       (PTDR_T_IIpsi_tau * la["PTDR_T_IIpsi", i - 1, a3]) +
                      -(PTDR_lambda[i - 1, a3] * PTDR_p[a3] * PTDR_x * la["PTDR_R", i - 1, a3]) +
                      -(PTDR_r[a3] * la["PTDR_R", i - 1, a3]) +
                      -(PTDS_lambda[i - 1, a3] * PTDS_p[a3] * PTDS_x * la["PTDR_R", i - 1, a3])

## Misdiagnosed and Treated - Infectious
la["PTDR_mdt", i, a4] <-  (surv4 * la["PTDR_mdt", i - 1, a4]) +
                          # Misdiagnosed and treated
                          (PTDR_II_mdt[a4] * la["PTDR_II", i - 1, a4]) +
                          # Exit outwards
                          -(la["PTDR_mdt", i - 1, a4] * PTDR_mdt_exit) +
                          -(PTDR_mdt_u[a4] * la["PTDR_mdt", i - 1, a4]) +
                          (surv3 * la["PTDR_mdt", i - 1, a3]) +
                          # Misdiagnosed and treated
                          (PTDR_II_mdt[a3] * la["PTDR_II", i - 1, a3]) +
                          # Exit outwards
                          -(la["PTDR_mdt", i - 1, a3] * PTDR_mdt_exit) +
                          -(PTDR_mdt_u[a3] * la["PTDR_mdt", i - 1, a3])

## Misdiagnosed and Treated - Non-Infectious
la["PTDR_nid", i, a4] <-  (surv4 * la["PTDR_nid", i - 1, a4]) +
                          # Misdiagnosed and treated
                          (PTDR_IN_nid[a4] * la["PTDR_IN", i - 1, a4]) +
                          # Exit outwards
                          -(la["PTDR_nid", i - 1, a4] * PTDR_nid_exit) +
                          -(PTDR_nid_u[a4] * la["PTDR_nid", i - 1, a4]) +
                          (surv3 * la["PTDR_mdt", i - 1, a3]) +
                          # Misdiagnosed and treated
                          (PTDR_IN_nid[a3] * la["PTDR_IN", i - 1, a3]) +
                          # Exit outwards
                          -(la["PTDR_nid", i - 1, a3] * PTDR_nid_exit) +
                          -(PTDR_nid_u[a3] * la["PTDR_nid", i - 1, a3])

      # Check if scaling needed in 1950
      if (current_year == 1950) {
        scaling_table <- la[, current_step, ]
        age_props <- matrix(0, nrow(scaling_table), ncol(scaling_table))
        scol_tot <- colSums(scaling_table)
        for (scol in 1:ncol(scaling_table)) {
          if (sum(scaling_table[, scol]) > 0) {
            age_props[, scol] <- scaling_table[, scol]/scol_tot[scol]
          } else {
            age_props[1, scol] <- 1
          }
        }
        # Multiply new compartment proportions
        for (age in 1:length(total_population[51, ])) {
          la[, current_step, age] <- (age_props[, age] * total_population[51, age])
        }
      }

      # Calculate SDR
      # Check if current year = 2012, if so calculate 2011 prevalence rate.
      if (current_year == (sdr_start_year + 1)) {
        # Calculate 2011 AllTB prevalence rate
        sdr_start_yr_steps <- ((sdr_start_year - year1) * (1/dt)) + seq(1:(1/dt))
        bp_prevcmp <- c("NTDS_II", "PTDS_II", "NTDR_II", "PTDR_II")
        sdr_start_yr_prev <- sum(colSums(la[bp_prevcmp, sdr_start_yr_steps, ], dim = 1))/sum(la[, sdr_start_yr_steps, ]) * 1e+05
        rm(sdr_start_yr_steps, bp_prevcmp)
      }

      # Calculate SDR of current step
      if (current_year >= (sdr_start_year + 1)) {
        # Calculate prevalence rate of current step
        bp_prevcmp <- c("NTDS_II", "PTDS_II", "NTDR_II", "PTDR_II")
        bp_prev_current_step <- sum(la[bp_prevcmp, current_step, ])/sum(la[, current_step, ]) * 1e+05
        sdr[current_step] <- sdr_base * bp_prev_current_step / sdr_start_yr_prev
        rm(bp_prev_current_step, bp_prevcmp)
      }

      # Source and run the transit equations
      rm(a1, a0)
      a1 <- 1:(max_age)
      a0 <- 1:(max_age)
      a4 <- 100
      a3 <- 99
      a5 <- 1:max_age


# DSTB deaths
ta["DSTB_deaths", i - 1, a5] <- (NTDS_II_u[a5] * la["NTDS_II", i - 1, a5]) +
                            (NTDS_IN_u[a5] * la["NTDS_IN", i - 1, a5]) +
                            (NTDS_T_u[a5] * la["NTDS_T", i - 1, a5]) +
                            (PTDS_II_u[a5] * la["PTDS_II", i - 1, a5]) +
                            (PTDS_IN_u[a5] * la["PTDS_IN", i - 1, a5]) +
                            (PTDS_T_u[a5] * la["PTDS_T", i - 1, a5])

# Backup calculation - DSTB deaths
dsma[current_step - 1, a5]  <-  (NTDS_II_u[a5] * la["NTDS_II", i - 1, a5]) +
                            (NTDS_IN_u[a5] * la["NTDS_IN", i - 1, a5]) +
                            (NTDS_T_u[a5] * la["NTDS_T", i - 1, a5]) +
                            (PTDS_II_u[a5] * la["PTDS_II", i - 1, a5]) +
                            (PTDS_IN_u[a5] * la["PTDS_IN", i - 1, a5]) +
                            (PTDS_T_u[a5] * la["PTDS_T", i - 1, a5])

# DRTB Deaths
ta["DRTB_deaths", i - 1, a5] <- (NTDR_II_u[a5] * la["NTDR_II", i - 1, a5]) +
                            (NTDS_IN_u[a5] * la["NTDR_IN", i - 1, a5]) +
                            (NTDR_T_u[a5] * la["NTDR_T_IIpsi", i - 1, a5]) +
                            (NTDR_II_u[a5] * la["NTDR_T_IIphi", i - 1, a5]) +
                            (NTDR_T_u[a5] * la["NTDR_T_INpsi", i - 1, a5]) +
                            (NTDR_IN_u[a5] * la["NTDR_T_INphi", i - 1, a5]) +
                            (PTDR_II_u[a5] * la["PTDR_II", i - 1, a5]) +
                            (PTDS_IN_u[a5] * la["PTDR_IN", i - 1, a5]) +
                            (PTDR_T_u[a5] * la["PTDR_T_IIpsi", i - 1, a5]) +
                            (PTDR_II_u[a5] * la["PTDR_T_IIphi", i - 1, a5]) +
                            (PTDR_T_u[a5] * la["PTDR_T_INpsi", i - 1, a5]) +
                            (PTDR_IN_u[a5] * la["PTDR_T_INphi", i - 1, a5]) +
                            (NTDR_mdt_u[a5] * la["NTDR_mdt", i - 1, a5])
                            (PTDR_mdt_u[a5] * la["PTDR_mdt", i - 1, a5])
                            (NTDR_nid_u[a5] * la["NTDR_nid", i - 1, a5])
                            (PTDR_nid_u[a5] * la["PTDR_nid", i - 1, a5])

# DSTB Incidence
ta["DSTB_inc", i - 1, a5] <-  ## Deconstructed
                          (NTDS_p[a5] * NTDS_f[a5] * NTDS_lambda[i - 1, a5] * la["S", i - 1, a5]) +
                          (NTDS_p[a5] * NTDS_f[a5] * NTDS_lambda[i - 1, a5] * NTDS_x * la["NTDR_L", i - 1, a5]) +
                          (NTDS_p[a5] * NTDS_f[a5] * NTDS_lambda[i - 1, a5] * NTDS_x * la["NTDS_L", i - 1, a5]) +
                          ## End deconstruct
                          ## Reactivation from Latent
                          (NTDS_v[a5] * NTDS_f[a5] * la["NTDS_L", i - 1, a5]) +
                          ## Reinfection of Resolved - DS
                          (NTDS_lambda[i - 1, a5] * NTDS_p[a5] * NTDS_f[a5] * NTDS_x * la["NTDS_R", i - 1, a5]) +
                          ## Reactivation from Resolved
                          (NTDS_f[a5] * NTDS_r[a5] * la["NTDS_R", i - 1, a5]) +
                          ## Reinfection of Resolved - DR
                          (NTDS_lambda[i - 1, a5] * NTDS_p[a5] * NTDS_f[a5] * NTDS_x * la["NTDR_R", i - 1, a5]) +
                          ## Deconstructed (NTDS)
                          (NTDS_p[a5] * (1 - NTDS_f[a5]) * NTDS_lambda[i - 1, a5] * la["S", i - 1, a5]) +
                          (NTDS_p[a5] * (1 - NTDS_f[a5]) * NTDS_lambda[i - 1, a5] * NTDS_x * la["NTDR_L", i - 1, a5]) +
                          (NTDS_p[a5] * (1 - NTDS_f[a5]) * NTDS_lambda[i - 1, a5] * NTDS_x * la["NTDS_L", i - 1, a5]) +
                          ## End deconstruct
                          ## Reactivation from Latent
                          (NTDS_v[a5] * (1 - NTDS_f[a5]) * la["NTDS_L", i - 1, a5]) +
                          ## Reinfection of Resolved - DS
                          (NTDS_lambda[i - 1, a5] * NTDS_p[a5] * (1 - NTDS_f[a5]) * NTDS_x * la["NTDS_R", i - 1,a5]) +
                          ## Reactivation from Resolved
                          ((1 - NTDS_f[a5]) * NTDS_r[a5] * la["NTDS_R", i-1,a5]) +
                          ## Reinfection of Resolved - DR
                          (NTDS_lambda[i - 1, a5] * NTDS_p[a5] * (1 - NTDS_f[a5]) * NTDS_x * la["NTDR_R", i - 1, a5]) +
                          ## Reinfection of Resolved - DS
                          (PTDS_lambda[i - 1, a5] * PTDS_p[a5] * PTDS_f[a5] * PTDS_x * la["PTDS_R", i - 1, a5]) +
                          ## Reactivation from Resolved
                          (PTDS_f[a5] * PTDS_r[a5] * la["PTDS_R", i - 1, a5]) +
                          ## Reinfection of Resolved - DR
                          (PTDS_lambda[i - 1, a5] * PTDS_p[a5] * PTDS_f[a5] * PTDS_x * la["PTDR_R", i - 1, a5]) +
                          ## Reinfection of Resolved - DS
                          (PTDS_lambda[i - 1, a5] * PTDS_p[a5] * (1 - PTDS_f[a5]) * PTDS_x * la["PTDS_R", i - 1, a5]) +
                          ## Reactivation from Resolved
                          ((1 - PTDS_f[a5]) * PTDS_r[a5] * la["PTDS_R", i - 1, a5]) +
                          ## Reinfection of Resolved - DR
                          (PTDS_lambda[i - 1, a5] * PTDS_p[a5] * (1 - PTDS_f[a5]) * PTDS_x * la["PTDR_R", i - 1, a5])

dsia[current_step - 1, a5]  <-  # Backup calculation
                            # New Infections of Susceptible and Latent
                            (NTDS_p[a5] * NTDS_lambda[i - 1, a5] * (la["S", i - 1, a5] + (NTDS_x * la["NTDR_L", i - 1, a5]) + (NTDS_x * la["NTDS_L", i - 1, a5]))) +
                            (NTDS_lambda[i - 1, a5] * NTDS_p[a5] * NTDS_x * (la["NTDS_R", i - 1, a5] + la["NTDR_R", i - 1, a5])) +
                            (NTDS_v[a5] * la["NTDS_L", i - 1, a5]) +
                            (NTDS_r[a5] * la["NTDS_R", i - 1, a5]) +
                            # PTDS
                            (PTDS_r[a5] * la["PTDS_R", i - 1, a5]) +
                            (PTDS_lambda[i - 1, a5] * PTDS_p[a5] * PTDS_x * (la["PTDR_R", i - 1, a5] + la["PTDS_R", i - 1, a5]))


dria[current_step - 1, a5] <-       # New infections of S/NTDS_L/NTDR_L
                                (NTDR_lambda[i - 1, a5] * NTDR_p[a5] * (la["S", i - 1, a5] + (NTDR_x * la["NTDR_L", i - 1, a5]) + (NTDR_x * la["NTDS_L", i - 1, a5]))) +
                                (NTDR_v[a5] * la["NTDR_L", i - 1, a5]) +
                                (NTDR_r[a5] * la["NTDR_R", i - 1, a5]) +
                                (NTDR_lambda[i - 1, a5] * NTDR_p[a5] * NTDR_x * (la["NTDR_R", i - 1, a5] + la["NTDS_R", i - 1, a5])) +
                                ## PT
                                (PTDR_v[a5] * la["PTDR_L", i - 1, a5]) + #ZERO COMPARTMENT
                                (PTDR_r[a5] * la["PTDR_R", i - 1, a5]) +
                                (PTDR_lambda[i - 1, a5] * PTDR_p[a5] * PTDR_x * (la["PTDS_R", i - 1, a5] + la["PTDR_R", i - 1, a5])) +
                                (NT_xi * PTDR_p[a5] * la["NTDS_T", i - 1, a5] * PTDR_f[a5] * nt_mdt_loss) +
                                (PT_xi * PTDR_p[a5] * la["PTDS_T", i - 1, a5] * PTDR_f[a5] * pt_mdt_loss) +
                                (NT_xi * PTDR_p[a5] * la["NTDS_T", i - 1, a5] * (1 - PTDR_f[a5]) * nt_nid_loss) +
                                (PT_xi * PTDR_p[a5] * la["PTDS_T", i - 1, a5] * (1 - PTDR_f[a5]) * pt_nid_loss)

# DRTB Incident Cases per timestep
ta["DRTB_inc", i - 1, a5]   <-  ## Deconstructed compound term - new infections in susceptble, NTDR-Latent and NTDS-Latent
                                (NTDR_p[a5] * NTDR_f[a5] * NTDR_lambda[i - 1, a5] * la["S", i - 1, a5]) +
                                (NTDR_p[a5] * NTDR_f[a5] * NTDR_lambda[i - 1, a5] * NTDR_x * la["NTDR_L", i - 1, a5]) +
                                (NTDR_p[a5] * NTDR_f[a5] * NTDR_lambda[i - 1, a5] * NTDR_x * la["NTDS_L", i - 1, a5]) +
                                ## End deconstruct
                                ## Reactivation from Latent
                                (NTDR_v[a5] * NTDR_f[a5] * la["NTDR_L", i - 1, a5]) +
                                ## Reinfection of Resolved - DR
                                (NTDR_lambda[i - 1, a5] * NTDR_p[a5] * NTDR_f[a5] * NTDR_x * la["NTDR_R", i - 1, a5]) +
                                ## Reactivation from Resolved
                                (NTDR_f[a5] * NTDR_r[a5] * la["NTDR_R", i - 1, a5]) +
                                ## Reinfection of Resolved - DS
                                (NTDR_lambda[i - 1, a5] * NTDR_p[a5] * NTDR_f[a5] * NTDR_x * la["NTDS_R", i - 1, a5]) +
                                # Deconstructed
                                (NTDR_p[a5] * (1 - NTDR_f[a5]) * NTDR_lambda[i - 1, a5] * la["S", i - 1, a5]) +
                                (NTDR_p[a5] * (1 - NTDR_f[a5]) * NTDR_lambda[i - 1, a5] * NTDR_x * la["NTDR_L", i - 1, a5]) +
                                (NTDR_p[a5] * (1 - NTDR_f[a5]) * NTDR_lambda[i - 1, a5] * NTDR_x * la["NTDS_L", i - 1, a5]) +
                                # End deconstruct
                                ## Reactivation from Latent
                                (NTDR_v[a5] * (1 - NTDR_f[a5]) * la["NTDR_L", i - 1, a5]) +
                                ## Reinfection of Resolved - DR
                                (NTDR_lambda[i - 1, a5] * NTDR_p[a5] * (1 - NTDR_f[a5]) * NTDR_x * la["NTDR_R", i - 1, a5]) +
                                ## Reactivation from Resolved
                                ((1 - NTDR_f[a5]) * NTDR_r[a5] * la["NTDR_R", i - 1, a5]) +
                                ## Reinfection of Resolved - DS
                                (NTDR_lambda[i - 1, a5] * NTDR_p[a5] * (1 - NTDR_f[a5]) * NTDR_x * la["NTDS_R", i - 1, a5]) +
                                ## Reactivation from Latent
                                (PTDR_v[a5] * PTDR_f[a5] * la["PTDR_L", i - 1, a5]) +
                                ## Reinfection of Resolved - DR
                                (PTDR_lambda[i - 1, a5] * PTDR_p[a5] * PTDR_f[a5] * PTDR_x * la["PTDR_R", i - 1, a5]) +
                                ## Reactivation from Resolved
                                (PTDR_f[a5] * (PTDR_r[a5] ) * la["PTDR_R", i - 1, a5]) +
                                ## Reinfection of Resolved - DS
                                (PTDR_lambda[i - 1, a5] * PTDR_p[a5] * PTDR_f[a5] * PTDR_x * la["PTDS_R", i - 1, a5]) +
                                ## Reactivation from Latent
                                (PTDR_v[a5] * (1 - PTDR_f[a5]) * la["PTDR_L", i - 1, a5]) +
                                ## Reinfection of Resolved - DS
                                (PTDR_lambda[i - 1, a5] * PTDR_p[a5] * (1 - PTDR_f[a5]) * PTDR_x * la["PTDS_R", i - 1, a5]) +
                                ## Reactivation from Resolved
                                ((1 - PTDR_f[a5]) * PTDR_r[a5] * la["PTDR_R", i - 1, a5]) +
                                ## Reinfection of Resolved - DR
                                (PTDR_lambda[i - 1, a5] * PTDR_p[a5] * (1 - PTDR_f[a5]) * PTDR_x * la["PTDR_R", i - 1, a5]) +
                                ## MDR Acquisition
                                (NT_xi * PTDR_p[a5] * la["NTDS_T", i - 1, a5] * PTDR_f[a5] * nt_mdt_loss) +
                                (PT_xi * PTDR_p[a5] * la["PTDS_T", i - 1, a5] * PTDR_f[a5] * pt_mdt_loss) +
                                (NT_xi * PTDR_p[a5] * la["NTDS_T", i - 1, a5] * (1 - PTDR_f[a5]) * nt_nid_loss) +
                                (PT_xi * PTDR_p[a5] * la["PTDS_T", i - 1, a5] * (1 - PTDR_f[a5]) * pt_nid_loss)

# ALL TB incidence backup calc
tbia[current_step - 1, a5]  <- dsia[current_step - 1, a5] + dria[current_step - 1, a5]

# ALL TB incidence
ta["TB_inc", i - 1, a5]     <- ta["DSTB_inc", i - 1, a5] + ta["DRTB_inc", i - 1, a5]

### Incidence by treatment history and bacteriologic status ###
## NTDR incidence
ta["DR_nt_inc", i - 1, a5]  <-  # New transmission by infection of susceptible and latent pools
                            (NTDR_p[a5] * NTDR_lambda[i - 1, a5] * la["S", i - 1, a5]) +
                            (NTDR_p[a5] * NTDR_lambda[i - 1, a5] * NTDR_x * la["NTDR_L", i - 1, a5]) +
                            (NTDR_p[a5] * NTDR_lambda[i - 1, a5] * NTDR_x * la["NTDS_L", i - 1, a5]) +
                            ## Reinfection of Resolved - DR
                            (NTDR_lambda[i - 1, a5] * NTDR_p[a5] * NTDR_x * la["NTDR_R", i - 1, a5]) +
                            ## Reinfection of Resolved - DS
                            (NTDR_lambda[i - 1, a5] * NTDR_p[a5] * NTDR_x * la["NTDS_R", i - 1, a5]) +
                            ## Reactivation from Latent
                            (NTDR_v[a5] * la["NTDR_L", i - 1, a5]) +
                            ## Reactivation from Resolved
                            (NTDR_r[a5] * la["NTDR_R", i - 1, a5])

# NTDR Bact+ Incidence
ta["DR_nt_inc_f", i - 1, a5] <- (ta["DR_nt_inc", i - 1, a5] * NTDR_f[a5]) +
                            (NTDR_omega * la["NTDR_IN", i - 1, a5])

# PTDR Incidence
ta["DR_pt_inc", i - 1, a5] <-   # Reinfection of Resolved - DR
                            (PTDR_lambda[i - 1, a5] * PTDR_p[a5] * PTDR_x * la["PTDR_R", i - 1, a5]) +
                            ## Reinfection of Resolved - DS
                            (PTDR_lambda[i - 1, a5] * PTDR_p[a5] * PTDR_x * la["PTDS_R", i - 1, a5]) +
                            ## Reactivation from Resolved
                            ((PTDR_r[a5]) * la["PTDR_R", i - 1, a5]) +
                            ## Reactivation from Latent
                            (PTDR_v[a5] * la["PTDR_L", i - 1, a5])
                            ## MDR Acquisition - not included.

# PTDR Bact+ Incidence
ta["DR_pt_inc_f", i - 1, a5] <- (ta["DR_pt_inc", i - 1, a5] * PTDS_f[a5]) +
                            (PTDR_omega * la["PTDR_IN", i - 1, a5]) +
                            # Bacteriologically positive MDR acquisition
                            (NT_xi * PTDR_p[a5] * la["NTDS_T", i - 1, a5] * PTDR_f[a5] * nt_mdt_loss) +
                            (PT_xi * PTDR_p[a5] * la["PTDS_T", i - 1, a5] * PTDR_f[a5] * pt_mdt_loss)

# NTDS Incidence
ta["DS_nt_inc", i - 1, a5] <- # Transmission by infection of suceptible and latent pools
                            (NTDS_p[a5] * NTDS_lambda[i - 1, a5] * la["S", i - 1, a5]) +
                            (NTDS_p[a5] * NTDS_lambda[i - 1, a5] * NTDS_x * la["NTDR_L", i - 1, a5]) +
                            (NTDS_p[a5] * NTDS_lambda[i - 1, a5] * NTDS_x * la["NTDS_L", i - 1, a5]) +
                            ## Reactivation from Latent
                            (NTDS_v[a5] * la["NTDS_L", i - 1, a5]) +
                            ## Reinfection of Resolved - DS
                            (NTDS_lambda[i - 1, a5] * NTDS_p[a5] * NTDS_x * la["NTDS_R", i - 1, a5]) +
                            ## Reactivation from Resolved
                            (NTDS_r[a5]  * la["NTDS_R", i - 1, a5]) +
                            ## Reinfection of Resolved - DR
                            (NTDS_lambda[i - 1, a5] * NTDS_p[a5] * NTDS_x * la["NTDR_R", i - 1, a5])

# NTDS Bact+ Incidence
ta["DS_nt_inc_f", i - 1, a5] <- (ta["DS_nt_inc", i - 1, a5] * NTDS_f[a5]) +
                            (NTDS_omega  * la["NTDS_IN", i - 1, a5])

# PTDS Incidence
ta["DS_pt_inc", i - 1, a5] <- ## Reinfection of Resolved - DS
                          (PTDS_lambda[i - 1, a5] * PTDS_p[a5] * PTDS_x * la["PTDS_R", i - 1, a5]) +
                          ## Reactivation from Resolved
                          (PTDS_r[a5] * la["PTDS_R", i - 1, a5]) +
                          ## Reinfection of Resolved - DR
                          (PTDS_lambda[i - 1, a5] * PTDS_p[a5] * PTDS_x * la["PTDR_R", i - 1, a5])

# PTDS Bact+ Incidence
ta["DS_pt_inc_f", i - 1, a5] <- (ta["DS_pt_inc", i - 1, a5] * PTDS_f[a5]) +
                            (PTDS_omega * la["PTDS_IN", i - 1, a5])

# AllTB Never Treated Incidence
ta["All_nt_inc", i - 1, a5] <- ta["DR_nt_inc", i - 1, a5] + ta["DS_nt_inc", i - 1, a5]

# AllTB Previously Treated Incidence
ta["All_pt_inc", i - 1, a5] <- ta["DR_pt_inc", i - 1, a5] + ta["DS_pt_inc", i - 1, a5]

# AllTB Never Treated Bact+ Incidence
ta["All_nt_inc_f", i - 1, a5] <- ta["DR_nt_inc_f", i - 1, a5] + ta["DS_nt_inc_f", i - 1, a5]

# AllTB Previously Treated Bact+ Incidence
ta["All_pt_inc_f", i - 1, a5] <- ta["DR_pt_inc_f", i - 1, a5] + ta["DS_pt_inc_f", i - 1, a5]

############################ DRTB Treatment Initiations [Notifications] #################################

##########################
########## TZ ############

# 'Fitting factor', 'ff' for diagnostics only.
# ff = 1 i.e. no effect during normal model functioning.
ff <- 1

##########################
##########################

# NTDR Treatment Initiations to any destination - infectious
ta["DR_nt_intx_f", i - 1, a5]  <- (NTDR_II_kappa[a5] * NTDR_II_psi * la["NTDR_II", i - 1, a5]) +
                                  (NTDR_II_kappa[a5] * NTDR_II_phi * la["NTDR_II", i - 1, a5]) +
                                  (NTDR_II_mdt[a5] * la["NTDR_II", i - 1, a5])

# PTDR Treatment Initiations to any destination - infectious
ta["DR_pt_intx_f", i - 1, a5]  <- (PTDR_II_kappa[a5] * PTDR_II_psi * la["PTDR_II", i - 1, a5]) +
                                  (PTDR_II_kappa[a5] * PTDR_II_phi * la["PTDR_II", i - 1, a5]) +
                                  (PTDR_II_mdt[a5] * la["PTDR_II", i - 1, a5]) +
                                  (ff * (la["NTDR_mdt", i - 1, a5] * NTDR_mdt_exit * (1 - nt_mdt_loss) * PTDR_II_psi)) +
                                  (ff * (la["NTDR_mdt", i - 1, a5] * NTDR_mdt_exit * (1 - nt_mdt_loss) * PTDR_II_phi)) +
                                  (ff * (la["PTDR_mdt", i - 1, a5] * PTDR_mdt_exit * (1 - pt_mdt_loss) * PTDR_II_psi)) +
                                  (ff * (la["PTDR_mdt", i - 1, a5] * PTDR_mdt_exit * (1 - pt_mdt_loss) * PTDR_II_phi)) +
                                  (ff * (NT_xi * PTDR_p[a5] * PTDR_f[a5] * la["NTDS_T", i - 1, a5] * (1 - nt_mdt_loss) * PTDR_II_psi)) +
                                  (ff * (PT_xi * PTDR_p[a5] * PTDR_f[a5] * la["PTDS_T", i - 1, a5] * (1 - pt_mdt_loss) * PTDR_II_psi)) +
                                  (ff * (NT_xi * PTDR_p[a5] * PTDR_f[a5] * la["NTDS_T", i - 1, a5] * (1 - nt_mdt_loss) * PTDR_II_phi)) +
                                  (ff * (PT_xi * PTDR_p[a5] * PTDR_f[a5] * la["PTDS_T", i - 1, a5] * (1 - pt_mdt_loss) * PTDR_II_phi))

# NTDS Treatment initiations to any destination - infectious
ta["DS_nt_intx_f", i - 1, a5]  <- (NTDS_II_kappa[a5] * la["NTDS_II", i - 1, a5])

# PTDS Treatment initiations to any destination - infectious
ta["DS_pt_intx_f", i - 1, a5]  <- (PTDS_II_kappa[a5] * la["PTDS_II", i - 1, a5])

# All-Never-Treated TB Treatment initiations to any destination - infectious
ta["All_nt_intx_f", i - 1, a5] <- (ta["DR_nt_intx_f", i - 1, a5] + ta["DS_nt_intx_f", i - 1, a5])

# All-Previously-Treated TB Treatment initiations to any destination - infectious
ta["All_pt_intx_f", i - 1, a5] <- (ta["DR_pt_intx_f", i - 1, a5] + ta["DS_pt_intx_f", i - 1, a5])

# Treatment Initations by _type of regimen_
## Treatment Initiations of DSTB Regimen
ta["DSTB_initRx", i - 1, a5] <- (NTDS_II_kappa[a5] * la["NTDS_II", i - 1, a5]) +
                                (NTDS_IN_kappa[a5] * la["NTDS_IN", i - 1, a5]) +
                                (PTDS_II_kappa[a5] * la["PTDS_II", i - 1, a5]) +
                                (PTDS_IN_kappa[a5] * la["PTDS_IN", i - 1, a5]) +
                                (NTDR_II_mdt[a5] * la["NTDR_II", i - 1, a5]) +
                                (NTDR_IN_nid[a5] * la["NTDR_IN", i - 1, a5]) +
                                (PTDR_II_mdt[a5] * la["PTDR_II", i - 1, a5]) +
                                (PTDR_IN_nid[a5] * la["PTDR_IN", i - 1, a5])

## Treatment Initiations of DRTB Regimen
ta["DRTB_initRx", i - 1, a5] <- (NTDR_II_kappa[a5] * NTDR_II_psi * la["NTDR_II", i - 1, a5]) +
                                (NTDR_II_kappa[a5] * NTDR_II_phi * la["NTDR_II", i - 1, a5]) +
                                (NTDR_IN_kappa[a5] * NTDR_IN_psi * la["NTDR_IN", i - 1, a5]) +
                                (NTDR_IN_kappa[a5] * NTDR_IN_phi * la["NTDR_IN", i - 1, a5]) +
                                (PTDR_II_kappa[a5] * PTDR_II_psi * la["PTDR_II", i - 1, a5]) +
                                (PTDR_II_kappa[a5] * PTDR_II_phi * la["PTDR_II", i - 1, a5]) +
                                (PTDR_IN_kappa[a5] * PTDR_IN_psi * la["PTDR_IN", i - 1, a5]) +
                                (PTDR_IN_kappa[a5] * PTDR_IN_phi * la["PTDR_IN", i - 1, a5]) +
                                (la["NTDR_mdt", i - 1, a5] * NTDR_mdt_exit * (1 - nt_mdt_loss) * PTDR_II_psi) +
                                (la["PTDR_mdt", i - 1, a5] * PTDR_mdt_exit * (1 - pt_mdt_loss) * PTDR_II_psi) +
                                (la["NTDR_mdt", i - 1, a5] * NTDR_mdt_exit * (1 - nt_mdt_loss) * PTDR_II_phi) +
                                (la["PTDR_mdt", i - 1, a5] * PTDR_mdt_exit * (1 - pt_mdt_loss) * PTDR_II_phi) +
                                (la["NTDR_nid", i - 1, a5] * NTDR_nid_exit * (1 - nt_nid_loss) * PTDR_IN_psi) +
                                (la["PTDR_nid", i - 1, a5] * PTDR_nid_exit * (1 - pt_nid_loss) * PTDR_IN_psi) +
                                (la["NTDR_nid", i - 1, a5] * NTDR_nid_exit * (1 - nt_mdt_loss) * PTDR_IN_phi) +
                                (la["PTDR_nid", i - 1, a5] * PTDR_nid_exit * (1 - pt_mdt_loss) * PTDR_IN_phi) +
                                (NT_xi * PTDR_p[a5] * PTDR_f[a5] * la["NTDS_T", i - 1, a5] * (1 - nt_mdt_loss) * PTDR_II_psi) +
                                (PT_xi * PTDR_p[a5] * PTDR_f[a5] * la["PTDS_T", i - 1, a5] * (1 - pt_mdt_loss) * PTDR_II_psi) +
                                (NT_xi * PTDR_p[a5] * PTDR_f[a5] * la["NTDS_T", i - 1, a5] * (1 - nt_mdt_loss) * PTDR_II_phi) +
                                (PT_xi * PTDR_p[a5] * PTDR_f[a5] * la["PTDS_T", i - 1, a5] * (1 - pt_mdt_loss) * PTDR_II_phi) +
                                (NT_xi * PTDR_p[a5] * (1 - PTDR_f[a5]) * la["NTDS_T", i - 1, a5] * (1 - nt_nid_loss) * PTDR_IN_psi) +
                                (PT_xi * PTDR_p[a5] * (1 - PTDR_f[a5]) * la["PTDS_T", i - 1, a5] * (1 - pt_nid_loss) * PTDR_IN_psi) +
                                (NT_xi * PTDR_p[a5] * (1 - PTDR_f[a5]) * la["NTDS_T", i - 1, a5] * (1 - nt_nid_loss) * PTDR_IN_phi) +
                                (PT_xi * PTDR_p[a5] * (1 - PTDR_f[a5]) * la["PTDS_T", i - 1, a5] * (1 - pt_nid_loss) * PTDR_IN_phi)

## Laboratory Confirmed Treatment Initiations of DRTB Regimen
ta["DRTB_initRxLab", i - 1, a5] <-(NTDR_II_kappa[a5] * NTDR_II_psi * la["NTDR_II", i - 1, a5] * nt_dst_p[current_year - year1  + 1]) +
                                  (NTDR_II_kappa[a5] * NTDR_II_phi * la["NTDR_II", i - 1, a5] * nt_dst_p[current_year - year1  + 1]) +
                                  (PTDR_II_kappa[a5] * PTDR_II_psi * la["PTDR_II", i - 1, a5] * pt_dst_p[current_year - year1  + 1]) +
                                  (PTDR_II_kappa[a5] * PTDR_II_phi * la["PTDR_II", i - 1, a5] * pt_dst_p[current_year - year1  + 1]) +
                                  (la["NTDR_mdt", i - 1, a5] * NTDR_mdt_exit * nt_dst_p[current_year - year1  + 1] * PTDR_II_psi) +
                                  (la["PTDR_mdt", i - 1, a5] * PTDR_mdt_exit * pt_dst_p[current_year - year1  + 1] * PTDR_II_psi) +
                                  (la["NTDR_mdt", i - 1, a5] * NTDR_mdt_exit * nt_dst_p[current_year - year1  + 1] * PTDR_II_phi) +
                                  (la["PTDR_mdt", i - 1, a5] * PTDR_mdt_exit * pt_dst_p[current_year - year1  + 1] * PTDR_II_phi) +
                                  (NT_xi * PTDR_p[a5] * PTDR_f[a5] * la["NTDS_T", i - 1, a5] * nt_dst_p[current_year - year1  + 1]) +
                                  (PT_xi * PTDR_p[a5] * PTDR_f[a5] * la["PTDS_T", i - 1, a5] * pt_dst_p[current_year - year1  + 1]) +
                                  (NT_xi * PTDR_p[a5] * PTDR_f[a5] * la["NTDS_T", i - 1, a5] * nt_dst_p[current_year - year1  + 1] * PTDR_II_phi) +
                                  (PT_xi * PTDR_p[a5] * PTDR_f[a5] * la["PTDS_T", i - 1, a5] * pt_dst_p[current_year - year1  + 1] * PTDR_II_phi)

# DSTB person-time on treatment - ** IN PERSON-MONTHS **
ta["DSTB_onRx", i - 1, a5] <- (12 * dt) *  (la["NTDS_T", i - 1, a5] +
                                        la["PTDS_T", i - 1, a5] +
                                        la["NTDR_mdt", i - 1, a5] +
                                        la["PTDR_mdt", i - 1, a5] +
                                        la["NTDR_nid", i - 1, a5] +
                                        la["PTDR_nid", i - 1, a5])

# DRTB person-time on treatment ** IN PERSON-MONTHS**
ta["DRTB_onRx", i - 1, a5] <- (12 * dt) *  (la["NTDR_T_IIphi", i - 1, a5] +
                                        la["NTDR_T_IIpsi", i - 1, a5] +
                                        la["NTDR_T_INphi", i - 1, a5] +
                                        la["NTDR_T_INpsi", i - 1, a5] +
                                        la["PTDR_T_IIphi", i - 1, a5] +
                                        la["PTDR_T_IIpsi", i - 1, a5] +
                                        la["PTDR_T_INphi", i - 1, a5] +
                                        la["PTDR_T_INpsi", i - 1, a5])

# DSTB Diagnostic Cost
ca["ds_dx", i - 1, a5] <- (ta["DSTB_initRx", i - 1, a5] * sdr[i - 1] * ds_dx_cost) * 1e+03

# DRTB Diagnostic Cost (Excl DST)
ca["dr_dx", i - 1, a5] <- (ta["DRTB_initRx", i - 1, a5] * sdr[i - 1] * dr_dx_cost) * 1000

# DST Cost
ca["dst", i - 1, a5]   <- (ta["DRTB_initRxLab", i - 1, a5] * dst_cost) * 1000

# DSTB Treatment Cost
ca["ds_tx", i - 1, a5] <- (ta["DSTB_onRx", i - 1, a5] * ds_tx_cost) * 1000

# DRTB Treatment Cost
ca["dr_tx", i - 1, a5] <- (ta["DRTB_onRx", i - 1, a5] * dr_tx_cost) * 1000

# Total treatment cost including fractional inflation for programme cost
ca["tbrx", i - 1, a5]  <- #
                        (ca["ds_dx", i - 1, a5] +
                        ca["dr_dx", i - 1, a5] +
                        ca["dst", i - 1, a5] +
                        ca["ds_tx", i - 1, a5] +
                        ca["dr_tx", i - 1, a5]) * (1 + prog_cost)

      # Recalculate matrices for infection parameters
      # Popsize
      psizematrix[i, 1] <- sum(la[, i, 1:6])
      psizematrix[i, 2] <- sum(la[, i, 7:20])
      psizematrix[i, 3] <- sum(la[, i, 21:65])
      psizematrix[i, 4] <- sum(la[, i, 66:max_age])

      ## Total Infectious Cases by contact matrix age classes - DSTB
      DS_Imatrix[i, 1]  <- sum(la["PTDS_II", i, 1:6], la["NTDS_II", i, 1:6])
      DS_Imatrix[i, 2]  <- sum(la["PTDS_II", i, 7:20], la["NTDS_II", i, 7:20])
      DS_Imatrix[i, 3]  <- sum(la["PTDS_II", i, 21:65], la["NTDS_II", i, 21:65])
      DS_Imatrix[i, 4]  <- sum(la["PTDS_II", i, 66:max_age], la["NTDS_II", i, 66:max_age])

      ## Total Infectious Cases by contact matrix age classes - DRTB
      DR_Imatrix[i, 1]  <- sum(la["PTDR_II", i, 1:6], la["NTDR_II", i, 1:6], la["PTDR_T_IIphi", i, 1:6], la["NTDR_T_IIphi", i, 1:6], la["NTDR_mdt", i, 1:6], la["PTDR_mdt", i, 1:6])
      DR_Imatrix[i, 2]  <- sum(la["PTDR_II", i, 7:20], la["NTDR_II", i, 7:20], la["PTDR_T_IIphi", i, 7:20], la["NTDR_T_IIphi", i, 7:20], la["NTDR_mdt", i, 7:20], la["PTDR_mdt", i, 7:20])
      DR_Imatrix[i, 3]  <- sum(la["PTDR_II", i, 21:65], la["NTDR_II", i, 21:65], la["PTDR_T_IIphi", i, 21:65], la["NTDR_T_IIphi", i, 21:65], la["NTDR_mdt", i, 21:65], la["PTDR_mdt", i, 21:65])
      DR_Imatrix[i, 4]  <- sum(la["PTDR_II", i, 66:max_age], la["NTDR_II", i, 66:max_age], la["PTDR_T_IIphi", i, 66:max_age], la["NTDR_T_IIphi", i, 66:max_age], la["NTDR_mdt", i, 66:max_age], la["PTDR_mdt", i, 66:max_age])

    } else {

      # Start Step-Not-1-of-year Calculations

      # Transmission term calculation - lambda
      NTDS_lambda_raw <- colSums(-(myneta[1:4, 1:max_age]) * z * ((DS_Imatrix[i - 1, 1:4])/(psizematrix[i - 1, 1:4])))
      NTDS_lambda[i - 1, 1:max_age] <- t(DS_neta * (1 - exp(NTDS_lambda_raw)))
      PTDS_lambda <- NTDS_lambda

      NTDR_lambda_raw <- colSums(-(myneta[1:4, 1:max_age]) * z * ((DR_Imatrix[i - 1, 1:4])/(psizematrix[i - 1, 1:4])))
      NTDR_lambda[i - 1, 1:max_age] <- t(DR_neta * (1 - exp(NTDR_lambda_raw)))
      PTDR_lambda <- NTDR_lambda

      # Source and run the epi equations
      rm(a1, a0)
      a1 <- 1:max_age
      a0 <- 1:max_age


rm(surv)
surv <- 1 - (u[a0])

# Susceptibles
la["S", i, a1] <-   (surv * la["S", i - 1, a0]) +
                   -(NTDS_lambda[i - 1, a0] * la["S", i - 1, a0]) +
                   -(NTDR_lambda[i - 1, a0] * la["S", i - 1, a0])

rm(surv)
surv <- 1 - (u[a0])

## Never Treated Drug Sensitive TB
## NTDS-Latent
la["NTDS_L", i, a1] <- (surv * la["NTDS_L", i - 1, a0]) +
                       ((1 - NTDS_p[a0]) * NTDS_lambda[i - 1, a0] * la["S", i - 1, a0]) +
                       ((1 - NTDS_p[a0]) * NTDS_lambda[i - 1, a0] * NTDS_x * la["NTDR_L", i - 1, a0]) +
                      -(NTDS_v[a0]  * la["NTDS_L", i - 1, a0]) +
                      -(NTDS_x * NTDS_lambda[i - 1, a0] * NTDS_p[a0] * la["NTDS_L", i - 1, a0]) +
                      -(NTDR_x * NTDR_lambda[i - 1, a0] * la["NTDS_L", i - 1, a0])

## Infectious NTDS
la["NTDS_II", i, a1] <- (surv * la["NTDS_II", i - 1, a0]) +
                        ## Conversion of Non-infectious to Infectious
                        (NTDS_omega  * la["NTDS_IN", i - 1, a0]) +
                        ## Deconstructed
                        (NTDS_p[a0] * NTDS_f[a0] * NTDS_lambda[i - 1, a0] * la["S", i - 1, a0]) +
                        (NTDS_p[a0] * NTDS_f[a0] * NTDS_lambda[i - 1, a0] * NTDS_x * la["NTDR_L", i - 1, a0]) +
                        (NTDS_p[a0] * NTDS_f[a0] * NTDS_lambda[i - 1, a0] * NTDS_x * la["NTDS_L", i - 1, a0]) +
                        ## End deconstruct
                        ## Reactivation from Latent
                        (NTDS_v[a0]  * NTDS_f[a0] * la["NTDS_L", i - 1, a0]) +
                        ## Reinfection of Resolved - DS
                        (NTDS_lambda[i - 1, a0] * NTDS_p[a0] * NTDS_f[a0] * NTDS_x * la["NTDS_R", i - 1, a0]) +
                        ## Reactivation from Resolved
                        (NTDS_f[a0] * NTDS_r[a0]  * la["NTDS_R", i - 1, a0]) +
                        ## Reinfection of Resolved - DR
                        (NTDS_lambda[i - 1, a0] * NTDS_p[a0] * NTDS_f[a0] * NTDS_x * la["NTDR_R", i - 1, a0]) +
                        ## Natural Cure
                       -(NTDS_II_n[a0]  * la["NTDS_II", i - 1, a0]) +
                        ## Case detection
                       -(NTDS_II_kappa[a0] * la["NTDS_II", i - 1, a0]) +
                        ## TB Death
                       -(NTDS_II_u[a0] * la["NTDS_II", i - 1, a0])

## Non-Infectious NTDS
la["NTDS_IN", i,a1] <- (surv * la["NTDS_IN", i-1,a0]) +
                       ## Deconstructed
                       (NTDS_p[a0] * (1 - NTDS_f[a0]) * NTDS_lambda[i - 1, a0] * la["S", i - 1, a0]) +
                       (NTDS_p[a0] * (1 - NTDS_f[a0]) * NTDS_lambda[i - 1, a0] * NTDS_x * la["NTDR_L", i - 1, a0]) +
                       (NTDS_p[a0] * (1 - NTDS_f[a0]) * NTDS_lambda[i - 1, a0] * NTDS_x * la["NTDS_L", i - 1, a0]) +
                       ## End deconstruct
                       ## Reactivation from Latent
                       (NTDS_v[a0]  * (1 - NTDS_f[a0]) * la["NTDS_L", i - 1, a0]) +
                       ## Reinfection of Resolved - DS
                       (NTDS_lambda[i - 1, a0] * NTDS_p[a0] * (1 - NTDS_f[a0]) * NTDS_x * la["NTDS_R", i - 1,a0]) +
                       ## Reactivation from Resolved
                       ((1 - NTDS_f[a0]) * NTDS_r[a0]  * la["NTDS_R", i-1,a0]) +
                       ## Reinfection of Resolved - DR
                       (NTDS_lambda[i - 1, a0] * NTDS_p[a0] * (1 - NTDS_f[a0]) * NTDS_x * la["NTDR_R", i - 1, a0]) +
                       ## Natural Cure
                      -(NTDS_IN_n[a0]  * la["NTDS_IN", i - 1, a0]) +
                       ## Case Detection
                      -(NTDS_IN_kappa[a0] * la["NTDS_IN", i - 1, a0]) +
                       ## TB Death
                      -(NTDS_IN_u[a0] * la["NTDS_IN", i - 1, a0]) +
                       ## Conversion from Non-infectious to Infectious
                      -(NTDS_omega  * la["NTDS_IN", i - 1, a0])

## NTDS on treatment
la["NTDS_T", i,a1]  <- (surv * la["NTDS_T", i-1,a0]) +
                       ## Treatment initiation from Infectious
                       (NTDS_II_kappa[a0]*la["NTDS_II", i-1,a0]) +
                       ## Treatment initiation from Non-infectious
                       (NTDS_IN_kappa[a0]*la["NTDS_IN", i-1,a0]) +
                       ## Resistance acquisition
                      -(NT_xi * la["NTDS_T", i-1,a0]) +
                       ## Treatment failures
                      -(NTDS_phi * la["NTDS_T", i-1,a0]) +
                       ## Treatment successes
                      -(NTDS_psi * la["NTDS_T", i-1,a0]) +
                       ## Death on treatment
                      -(NTDS_T_u[a0] * la["NTDS_T", i-1,a0])

# NTDS Resolved
la["NTDS_R", i, a1] <- (surv * la["NTDS_R", i - 1, a0]) +
                       ## Natural cure Infectious
                       (NTDS_II_n[a0]  * la["NTDS_II", i - 1, a0]) +
                       ## Natural cure Non-infectious
                       (NTDS_IN_n[a0]  * la["NTDS_IN", i - 1, a0]) +
                       ## Reinfection of Resolved - DS
                      -(NTDS_lambda[i - 1, a0] * NTDS_p[a0] * NTDS_x * la["NTDS_R", i - 1, a0]) +
                       ## Reactivation of Resolved
                      -(NTDS_r[a0]  * la["NTDS_R", i - 1, a0]) +
                       ## Reinfection of Resolved - DR
                      -(NTDR_lambda[i - 1, a0] * NTDR_p[a0] * NTDR_x * la["NTDS_R", i - 1, a0])

rm(surv)
surv <- 1 - (u[a0] )

## Never Treated Drug Resistant

## NTDR-Latent
la["NTDR_L", i, a1] <- (surv * la["NTDR_L", i - 1, a0]) +
                       ## New Infections from Susceptible
                       ((1 - NTDR_p[a0]) * NTDR_lambda[i - 1, a0] * la["S", i - 1, a0]) +
                       ## New Infections from NTDS-Latent
                       ((1 - NTDR_p[a0]) * NTDR_lambda[i - 1, a0] * NTDR_x * la["NTDS_L", i - 1, a0]) +
                       ## Reactivation from Latent
                      -(NTDR_v[a0] * la["NTDR_L", i - 1, a0]) +
                       ## DRTB Reinfection of NTDR-Latent
                      -(NTDR_x * NTDR_lambda[i - 1, a0] * NTDR_p[a0] * la["NTDR_L", i - 1, a0]) +
                       ## DSTB Reinfection of NTDR-Latent
                      -(NTDR_x * NTDS_lambda[i - 1, a0] * la["NTDR_L", i - 1, a0])

## Infectious NTDR
la["NTDR_II", i,a1] <- (surv * la["NTDR_II", i - 1, a0]) +
                       ## Conversion of non-infectious to infectious
                       (NTDR_omega * la["NTDR_IN", i - 1, a0]) +
                       ## Deconstructed compound term - new infections in susceptble, NTDR-Latent and NTDS-Latent
                       (NTDR_p[a0] * NTDR_f[a0] * NTDR_lambda[i - 1, a0] * la["S", i - 1, a0]) +
                       (NTDR_p[a0] * NTDR_f[a0] * NTDR_lambda[i - 1, a0] * NTDR_x * la["NTDR_L", i - 1, a0]) +
                       (NTDR_p[a0] * NTDR_f[a0] * NTDR_lambda[i - 1, a0] * NTDR_x * la["NTDS_L", i - 1, a0]) +
                       ## End deconstruct
                       ## Reactivation from Latent
                       (NTDR_v[a0] * NTDR_f[a0] * la["NTDR_L", i - 1, a0]) +
                       ## Reinfection of Resolved - DR
                       (NTDR_lambda[i - 1, a0] * NTDR_p[a0] * NTDR_f[a0] * NTDR_x * la["NTDR_R", i - 1, a0]) +
                       ## Reactivation from Resolved
                       (NTDR_f[a0] * NTDR_r[a0] * la["NTDR_R", i - 1, a0]) +
                       ## Reinfection of Resolved - DS
                       (NTDR_lambda[i - 1, a0] * NTDR_p[a0] * NTDR_f[a0] * NTDR_x * la["NTDS_R", i - 1, a0]) +
                       ## Natural cure
                      -(NTDR_II_n[a0] * la["NTDR_II", i - 1, a0]) +
                       ## Case detection
                      -(NTDR_II_kappa[a0] * la["NTDR_II", i - 1, a0]) +
                       ## Misdiagnosis and Treatment
                      -(NTDR_II_mdt[a0] * la["NTDR_II", i - 1, a0]) +
                       ## TB death
                      -(NTDR_II_u[a0] * la["NTDR_II", i - 1, a0])

## Non-Infectious NTDR
la["NTDR_IN", i,a1] <- (surv * la["NTDR_IN", i - 1, a0]) +
                       # Deconstructed
                       (NTDR_p[a0] * (1 - NTDR_f[a0]) * NTDR_lambda[i - 1, a0] * la["S", i - 1, a0]) +
                       (NTDR_p[a0] * (1 - NTDR_f[a0]) * NTDR_lambda[i - 1, a0] * NTDR_x * la["NTDR_L", i - 1, a0]) +
                       (NTDR_p[a0] * (1 - NTDR_f[a0]) * NTDR_lambda[i - 1, a0] * NTDR_x * la["NTDS_L", i - 1, a0]) +
                       # End deconstruct
                       ## Reactivation from Latent
                       (NTDR_v[a0] * (1 - NTDR_f[a0]) * la["NTDR_L", i - 1, a0]) +
                       ## Reinfection of Resolved - DR
                       (NTDR_lambda[i - 1, a0] * NTDR_p[a0] * (1 - NTDR_f[a0]) * NTDR_x * la["NTDR_R", i - 1, a0]) +
                       ## Reactivation from Resolved
                       ((1 - NTDR_f[a0]) * NTDR_r[a0] * la["NTDR_R", i - 1, a0]) +
                       ## Reinfection of Resolved - DS
                       (NTDR_lambda[i - 1, a0] * NTDR_p[a0] * (1 - NTDR_f[a0]) * NTDR_x * la["NTDS_R", i - 1, a0]) +
                       ## Natural cure
                      -(NTDR_IN_n[a0] * la["NTDR_IN", i - 1, a0]) +
                       ## Case detection
                      -(NTDR_IN_kappa[a0] * la["NTDR_IN", i - 1, a0]) +
                       ## Misdiagnosis and Treatment
                      -(NTDR_IN_nid[a0] * la["NTDR_IN", i - 1, a0]) +
                       ## TB death
                      -(NTDR_IN_u[a0] * la["NTDR_IN", i - 1, a0]) +
                       ## Conversion from non-infectious to infectious
                      -(NTDR_omega * la["NTDR_IN", i - 1, a0])

## Treatment - NTDR
## Treatment Succeeding - Originating from Infectious DR-TB
la["NTDR_T_IIpsi", i, a1] <- (surv * la["NTDR_T_IIpsi", i - 1, a0]) +
                             (NTDR_II_kappa[a0] * NTDR_II_psi * la["NTDR_II", i - 1, a0]) +
                            -(NTDR_T_IIpsi_tau * la["NTDR_T_IIpsi", i - 1, a0]) +
                            -(NTDR_T_u[a0] * la["NTDR_T_IIpsi", i - 1, a0])

## Treatment Failing - Originating from Infectious DR-TB
la["NTDR_T_IIphi", i, a1] <- (surv * la["NTDR_T_IIphi", i - 1, a0]) +
                             (NTDR_II_kappa[a0] * NTDR_II_phi * la["NTDR_II", i - 1, a0]) +
                            -(NTDR_T_IIphi_tau * la["NTDR_T_IIphi", i - 1, a0]) +
                            -(NTDR_II_u[a0] * la["NTDR_T_IIphi", i - 1, a0])

## Treatment Succeeding - Originating from Non-Infectious DR-TB
la["NTDR_T_INpsi", i, a1] <- (surv * la["NTDR_T_INpsi", i - 1, a0]) +
                             (NTDR_IN_kappa[a0] * NTDR_IN_psi * la["NTDR_IN", i - 1, a0]) +
                            -(NTDR_T_INpsi_tau * la["NTDR_T_INpsi", i - 1, a0]) +
                            -(NTDR_T_u[a0] * la["NTDR_T_INpsi", i - 1, a0])

## Treatment Failing - Originating from Non-Infectious DR-TB
la["NTDR_T_INphi", i, a1] <- (surv * la["NTDR_T_INphi", i - 1, a0]) +
                             (NTDR_IN_kappa[a0] * NTDR_IN_phi * la["NTDR_IN", i - 1, a0]) +
                            -(NTDR_T_INphi_tau * la["NTDR_T_INphi", i - 1, a0]) +
                            -(NTDR_IN_u[a0] * la["NTDR_T_INphi", i - 1, a0])

## Resolved - DR TB
la["NTDR_R", i, a1] <- (surv * la["NTDR_R", i - 1, a0]) +
                       ## Natural cure - Infectious
                       (NTDR_II_n[a0] * la["NTDR_II", i - 1, a0]) +
                       ## Natural cure - Non-infectious
                       (NTDR_IN_n[a0] * la["NTDR_IN", i - 1, a0]) +
                       ## Reinfection of Resolved - to DR
                      -(NTDR_lambda[i - 1, a0] * NTDR_p[a0] * NTDR_x * la["NTDR_R", i - 1, a0]) +
                       ## Reactivation of Resolved
                      -(NTDR_r[a0] * la["NTDR_R", i - 1, a0]) +
                       ## Reinfection of Resolved - to DS
                      -(NTDS_lambda[i - 1, a0] * NTDS_p[a0] * NTDS_x * la["NTDR_R", i - 1, a0])

## Misdiagnosed and Treated - Infectious
la["NTDR_mdt", i, a1] <-  (surv * la["NTDR_mdt", i - 1, a0]) +
                          # Misdiagnosed and treated
                          (NTDR_II_mdt[a0] * la["NTDR_II", i - 1, a0]) +
                          # Exit outwards
                          -(la["NTDR_mdt", i - 1, a0] * NTDR_mdt_exit) +
                          -(NTDR_mdt_u[a0] * la["NTDR_mdt", i - 1, a0])

## Misdiagnosed and Treated - Npn-Infectious
la["NTDR_nid", i, a1] <-   (surv * la["NTDR_nid", i - 1, a0]) +
                           # Misdiagnosed and treated
                          (NTDR_IN_nid[a0] * la["NTDR_IN", i - 1, a0]) +
                           # Exit outwards
                          -(la["NTDR_nid", i - 1, a0] * NTDR_nid_exit) +
                          -(NTDR_nid_u[a0] * la["NTDR_nid", i - 1, a0])

rm(surv)
surv <- 1 - (u[a0] )

## Previously Treated Drug Sensitive TB

## Infectious PTDS
la["PTDS_II", i, a1] <- (surv * la["PTDS_II", i - 1, a0]) +
                        ## Conversion of non-infectious to infectious
                        (PTDS_omega * la["PTDS_IN", i - 1, a0]) +
                        ## Reinfection of Resolved - DS
                        (PTDS_lambda[i - 1, a0] * PTDS_p[a0] * PTDS_f[a0] * PTDS_x * la["PTDS_R", i - 1, a0]) +
                        ## Reactivation from Resolved
                        (PTDS_f[a0] * PTDS_r[a0] * la["PTDS_R", i - 1, a0]) +
                        ## Reinfection of Resolved - DR
                        (PTDS_lambda[i - 1, a0] * PTDS_p[a0] * PTDS_f[a0] * PTDS_x * la["PTDR_R", i - 1, a0]) +
                        ## Natural cure
                       -(PTDS_II_n[a0] * la["PTDS_II", i - 1, a0]) +
                        ## Case detection
                       -(PTDS_II_kappa[a0] * la["PTDS_II", i - 1, a0]) +
                        ## TB death
                       -(PTDS_II_u[a0] * la["PTDS_II", i - 1, a0])

## Non-Infectious PTDS
la["PTDS_IN", i, a1] <- (surv * la["PTDS_IN", i - 1, a0]) +
                        ## NTDS Treatment Failures
                        (NTDS_phi * la["NTDS_T", i - 1, a0]) +
                        ## PTDS Treatment Failures
                        (PTDS_phi * la["PTDS_T", i - 1, a0]) +
                        ## Reinfection of Resolved - DS
                        (PTDS_lambda[i - 1, a0] * PTDS_p[a0] * (1 - PTDS_f[a0]) * PTDS_x * la["PTDS_R", i - 1, a0]) +
                        ## Reactivation from Resolved
                        ((1 - PTDS_f[a0]) * PTDS_r[a0] * la["PTDS_R", i - 1, a0]) +
                        ## Reinfection of Resolved - DR
                        (PTDS_lambda[i - 1, a0] * PTDS_p[a0] * (1 - PTDS_f[a0]) * PTDS_x * la["PTDR_R", i - 1, a0]) +
                        ## Natural cure
                       -(PTDS_IN_n[a0] * la["PTDS_IN", i - 1, a0]) +
                        ## Case detection
                       -(PTDS_IN_kappa[a0] * la["PTDS_IN", i - 1, a0]) +
                        ## TB death
                       -(PTDS_IN_u[a0] * la["PTDS_IN", i - 1, a0]) +
                        ## Conversion from non-infectious to infectious
                       -(PTDS_omega * la["PTDS_IN", i - 1, a0])

## PTDS on treatment
la["PTDS_T", i, a1] <-  (surv * la["PTDS_T", i - 1, a0]) +
                        ## Treatment initiation from Infectious
                        (PTDS_II_kappa[a0] * la["PTDS_II", i - 1, a0]) +
                        ## Treatment initiation from Non-infectious
                        (PTDS_IN_kappa[a0] * la["PTDS_IN", i - 1, a0]) +
                        ## Resistance acquisition
                       -(PT_xi * la["PTDS_T", i - 1, a0]) +
                        ## Treatment failures
                       -(PTDS_phi * la["PTDS_T", i - 1, a0]) +
                        ## Treatment successes
                       -(PTDS_psi * la["PTDS_T", i - 1, a0]) +
                        ## Death on treatment
                       -(PTDS_T_u[a0] * la["PTDS_T", i - 1, a0])

## PTDS Resolved
la["PTDS_R", i, a1] <- (surv * la["PTDS_R", i - 1, a0]) +
                       ## Natural cure - Infectious
                       (PTDS_II_n[a0] * la["PTDS_II", i - 1, a0]) +
                       ## Natural cure - Non-infectious
                       (PTDS_IN_n[a0] * la["PTDS_IN", i - 1, a0]) +
                       ## NTDS Treatment Success
                       (NTDS_psi * la["NTDS_T", i - 1, a0]) +
                       ## PTDS Treatment Success
                       (PTDS_psi * la["PTDS_T", i - 1, a0]) +
                       ## Reinfection of Resolved - to DS
                      -(PTDS_lambda[i - 1, a0] * PTDS_p[a0] * PTDS_x * la["PTDS_R", i - 1, a0]) +
                       ## Reactivation of Resolved
                      -(PTDS_r[a0] * la["PTDS_R", i - 1, a0]) +
                       ## Reinfection of Resolved - to DR
                      -(PTDR_lambda[i - 1, a0] * PTDR_p[a0] * PTDR_x * la["PTDS_R", i - 1, a0])

rm(surv)
surv <- 1 - (u[a0] )

## Previously Treated Drug Resistant
## PTDR-Latent
la["PTDR_L", i, a1] <- (surv * la["PTDR_L", i - 1, a0]) +
                       ((PT_xi * (1 - PTDR_p[a0])) * (la["PTDS_T", i - 1, a0])) +
                       ((NT_xi * (1 - PTDR_p[a0])) * (la["NTDS_T", i - 1, a0])) +
                      -(PTDR_v[a0] * la["PTDR_L", i - 1, a0])

## Infectious PTDR
la["PTDR_II", i, a1] <- (surv * la["PTDR_II", i - 1, a0]) +
                        ## Conversion of non-infectious to infectious
                        (PTDR_omega * la["PTDR_IN", i - 1, a0]) +
                        ## Reactivation from Latent
                        (PTDR_v[a0] * PTDR_f[a0] * la["PTDR_L", i - 1, a0]) +
                        ## Reinfection of Resolved - DR
                        (PTDR_lambda[i - 1, a0] * PTDR_p[a0] * PTDR_f[a0] * PTDR_x * la["PTDR_R", i - 1, a0]) +
                        ## Reactivation from Resolved
                        (PTDR_f[a0] * (PTDR_r[a0] ) * la["PTDR_R", i - 1, a0]) +
                        ## Reinfection of Resolved - DS
                        (PTDR_lambda[i - 1, a0] * PTDR_p[a0] * PTDR_f[a0] * PTDR_x * la["PTDS_R", i - 1, a0]) +
                        ## MDR Acquisition
                        (NT_xi * PTDR_p[a0] * PTDR_f[a0] * la["NTDS_T", i - 1, a0] * nt_mdt_loss) +
                        (PT_xi * PTDR_p[a0] * PTDR_f[a0] * la["PTDS_T", i - 1, a0] * pt_mdt_loss) +
                        ## PTDR Treatment FAILURES - Infectious
                        (PTDR_T_IIphi_tau * la["PTDR_T_IIphi", i - 1, a0]) +
                        ## NTDR Treatment FAILURES - Infectious
                        (NTDR_T_IIphi_tau * la["NTDR_T_IIphi", i - 1, a0]) +
                        # Misdiagnosed - re-entry from Infectious MDT compartment
                        (la["NTDR_mdt", i - 1, a0] * NTDR_mdt_exit * nt_mdt_loss) +
                        (la["PTDR_mdt", i - 1, a0] * PTDR_mdt_exit * pt_mdt_loss) +
                        ## Natural cure
                       -(PTDR_II_n[a0] * la["PTDR_II", i - 1, a0]) +
                        ## Case detection
                       -(PTDR_II_kappa[a0] * la["PTDR_II", i - 1, a0]) +
                        ## Misdiagnosis and treatment
                       -(PTDR_II_mdt[a0] * la["PTDR_II", i - 1, a0]) +
                        ## TB death
                       -(PTDR_II_u[a0] * la["PTDR_II", i - 1, a0])

## Non-Infectious PTDR
la["PTDR_IN", i, a1] <- (surv * la["PTDR_IN", i - 1, a0]) +
                        ## Reactivation from Latent
                        (PTDR_v[a0] * (1 - PTDR_f[a0]) * la["PTDR_L", i - 1, a0]) +
                        ## Reinfection of Resolved - DS
                        (PTDR_lambda[i - 1, a0] * PTDR_p[a0] * (1 - PTDR_f[a0]) * PTDR_x * la["PTDS_R", i - 1, a0]) +
                        ## Reactivation from Resolved
                        ((1 - PTDR_f[a0]) * PTDR_r[a0] * la["PTDR_R", i - 1, a0]) +
                        ## Reinfection of Resolved - DR
                        (PTDR_lambda[i - 1, a0] * PTDR_p[a0] * (1 - PTDR_f[a0]) * PTDR_x * la["PTDR_R", i - 1, a0]) +
                        ## MDR Acquisition
                        (PT_xi * PTDR_p[a0] * (1 - PTDR_f[a0]) * la["PTDS_T", i - 1, a0] * pt_nid_loss) +
                        (NT_xi * PTDR_p[a0] * (1 - PTDR_f[a0]) * la["NTDS_T", i - 1, a0] * nt_nid_loss) +
                        ## PTDR Treatment FAILURES - Non-Infectious
                        (PTDR_T_INphi_tau * la["PTDR_T_INphi", i - 1, a0]) +
                        ## NTDR Treatment FAILURES - Non-Infectious
                        (NTDR_T_INphi_tau * la["NTDR_T_INphi", i - 1, a0]) +
                        # Misdiagnosed - re-entry from Non-Infectious NID compartment
                        (la["NTDR_nid", i - 1, a0] * NTDR_nid_exit * nt_nid_loss) +
                        (la["PTDR_nid", i - 1, a0] * PTDR_nid_exit * pt_nid_loss) +
                        ## Natural cure
                       -(PTDR_IN_n[a0] *la["PTDR_IN", i - 1, a0]) +
                        ## Case detection
                       -(PTDR_IN_kappa[a0] * la["PTDR_IN", i - 1, a0]) +
                        ## Misdiagnosis and treatment
                       -(PTDR_IN_nid[a0] * la["PTDR_IN", i - 1, a0]) +
                        ## TB death
                       -(PTDR_IN_u[a0] * la["PTDR_IN", i - 1, a0]) +
                        ## Conversion from non-infectious to infectious
                       -(PTDR_omega * la["PTDR_IN", i - 1, a0])

## PTDR on treatment
## Treatment Succeeding - Originating from Infectious DR-TB
la["PTDR_T_IIpsi", i, a1] <-  (surv * la["PTDR_T_IIpsi", i - 1, a0]) +
                              (PTDR_II_kappa[a0] * PTDR_II_psi * la["PTDR_II", i - 1, a0]) +
                              # Misdiagnosed - re-entry from MDR compartment
                              (la["NTDR_mdt", i - 1, a0] * NTDR_mdt_exit * (1 - nt_mdt_loss) * PTDR_II_psi) +
                              (la["PTDR_mdt", i - 1, a0] * PTDR_mdt_exit * (1 - pt_mdt_loss) * PTDR_II_psi) +
                              (NT_xi * PTDR_p[a0] * PTDR_f[a0] * la["NTDS_T", i - 1, a0] * (1 - nt_mdt_loss) * PTDR_II_psi) +
                              (PT_xi * PTDR_p[a0] * PTDR_f[a0] * la["PTDS_T", i - 1, a0] * (1 - pt_mdt_loss) * PTDR_II_psi) +
                             -(PTDR_T_IIpsi_tau * la["PTDR_T_IIpsi", i - 1, a0]) +
                             -(PTDR_T_u[a0] * la["PTDR_T_IIpsi", i - 1, a0])

## Treatment Failing - Originating from Infectious DR-TB
la["PTDR_T_IIphi", i, a1] <-  (surv * la["PTDR_T_IIphi", i - 1, a0]) +
                              (PTDR_II_kappa[a0] * PTDR_II_phi * la["PTDR_II", i - 1, a0]) +
                              # Misdiagnosed - re-entry from MDR compartment
                              (la["NTDR_mdt", i - 1, a0] * NTDR_mdt_exit * (1 - nt_mdt_loss) * PTDR_II_phi) +
                              (la["PTDR_mdt", i - 1, a0] * PTDR_mdt_exit * (1 - pt_mdt_loss) * PTDR_II_phi) +
                              (NT_xi * PTDR_p[a0] * PTDR_f[a0] * la["NTDS_T", i - 1, a0] * (1 - nt_mdt_loss) * PTDR_II_phi) +
                              (PT_xi * PTDR_p[a0] * PTDR_f[a0] * la["PTDS_T", i - 1, a0] * (1 - pt_mdt_loss) * PTDR_II_phi) +
                             -(PTDR_T_IIphi_tau * la["PTDR_T_IIphi", i - 1, a0]) +
                             -(PTDR_II_u[a0] * la["PTDR_T_IIphi", i - 1, a0])

## Treatment Succeeding - Originating from Non-Infectious DR-TB
la["PTDR_T_INpsi", i, a1] <-  (surv * la["PTDR_T_INpsi", i - 1, a0]) +
                              (PTDR_IN_kappa[a0] * PTDR_IN_psi * la["PTDR_IN", i - 1, a0]) +
                              # Misdiagnosed - re-entry from MDR compartment
                              (la["NTDR_nid", i - 1, a0] * NTDR_nid_exit * (1 - nt_nid_loss) * PTDR_IN_psi) +
                              (la["PTDR_nid", i - 1, a0] * PTDR_nid_exit * (1 - pt_nid_loss) * PTDR_IN_psi) +
                              (NT_xi * PTDR_p[a0] * (1 - PTDR_f[a0]) * la["NTDS_T", i - 1, a0] * (1 - nt_nid_loss) * PTDR_IN_psi) +
                              (PT_xi * PTDR_p[a0] * (1 - PTDR_f[a0]) * la["PTDS_T", i - 1, a0] * (1 - pt_nid_loss) * PTDR_IN_psi) +
                             -(PTDR_T_INpsi_tau * la["PTDR_T_INpsi", i - 1, a0]) +
                             -(PTDR_T_u[a0] * la["PTDR_T_INpsi", i - 1, a0])

## Treatment Failing - Originating from Non-Infectious DR-TB
la["PTDR_T_INphi", i, a1] <-  (surv * la["PTDR_T_INphi", i - 1, a0]) +
                              (PTDR_IN_kappa[a0] * PTDR_IN_phi * la["PTDR_IN", i - 1, a0]) +
                              # Misdiagnosed - re-entry from MDR compartment
                              (la["NTDR_nid", i - 1, a0] * NTDR_nid_exit * (1 - nt_nid_loss) * PTDR_IN_phi) +
                              (la["PTDR_nid", i - 1, a0] * PTDR_nid_exit * (1 - pt_nid_loss) * PTDR_IN_phi) +
                              (NT_xi * PTDR_p[a0] * (1 - PTDR_f[a0]) * la["NTDS_T", i - 1, a0] * (1 - nt_nid_loss) * PTDR_IN_phi) +
                              (PT_xi * PTDR_p[a0] * (1 - PTDR_f[a0]) * la["PTDS_T", i - 1, a0] * (1 - pt_nid_loss) * PTDR_IN_phi) +
                             -(PTDR_T_INphi_tau * la["PTDR_T_INphi", i - 1, a0]) +
                             -(PTDR_II_u[a0] * la["PTDR_T_INphi", i - 1, a0])

## PTDR Resolved
la["PTDR_R", i, a1] <- (surv * la["PTDR_R", i - 1, a0]) +
                       ## Natural cure - Infectious
                       (PTDR_II_n[a0] * la["PTDR_II", i - 1, a0]) +
                       ## Natural cure - Non-infectious
                       (PTDR_IN_n[a0] * la["PTDR_IN", i - 1, a0]) +
                       ## NTDR Treatment Success - Infectious
                       (NTDR_T_IIpsi_tau * la["NTDR_T_IIpsi", i - 1, a0]) +
                       ## NTDR Treatment Success - Non-infectious
                       (NTDR_T_INpsi_tau * la["NTDR_T_INpsi", i - 1, a0]) +
                       ## PTDR Treatment Success - Infectious
                       (PTDR_T_INpsi_tau * la["PTDR_T_INpsi", i - 1, a0]) +
                       ## PTDR Treatment Success - Non-infectious
                       (PTDR_T_IIpsi_tau * la["PTDR_T_IIpsi", i - 1, a0]) +
                       ## Reinfection of Resolved - to DS
                      -(PTDS_lambda[i - 1, a0] * PTDS_p[a0] * PTDS_x * la["PTDR_R", i - 1, a0]) +
                       ## Reactivation of Resolved
                      -(PTDR_r[a0] * la["PTDR_R", i - 1, a0]) +
                      ## Reinfection of Resolved - to DR
                      -(PTDR_lambda[i - 1, a0] * PTDR_p[a0] * PTDR_x * la["PTDR_R", i - 1, a0])

## Misdiagnosed and Treated - Infectious
la["PTDR_mdt", i, a1] <-  (surv * la["PTDR_mdt", i - 1, a0]) +
                          # Misdiagnosed and treated
                          (PTDR_II_mdt[a0] * la["PTDR_II", i - 1, a0]) +
                          # Exit outwards
                          -(la["PTDR_mdt", i - 1, a0] * PTDR_mdt_exit) +
                          -(PTDR_mdt_u[a0] * la["PTDR_mdt", i - 1, a0])

## Misdiagnosed and Treated - Npn-Infectious
la["PTDR_nid", i, a1] <-   (surv * la["PTDR_nid", i - 1, a0]) +
                           # Misdiagnosed and treated
                          (PTDR_IN_nid[a0] * la["PTDR_IN", i - 1, a0]) +
                           # Exit outwards
                          -(la["PTDR_nid", i - 1, a0] * PTDR_nid_exit) +
                          -(PTDR_nid_u[a0] * la["PTDR_nid", i - 1, a0])

      # Calculate SDR of current step
      if (current_year >= (sdr_start_year + 1)) {
        # Calculate prevalence rate of current step
        bp_prevcmp <- c("NTDS_II", "PTDS_II", "NTDR_II", "PTDR_II")
        bp_prev_current_step <- sum(la[bp_prevcmp, current_step, ])/sum(la[, current_step, ]) * 1e+05
        sdr[current_step] <- sdr_base * bp_prev_current_step / sdr_start_yr_prev
        rm(bp_prev_current_step, bp_prevcmp)
      }

      # Source and run the transit Equations
      rm(a1, a0)
      a1 <- 1:max_age
      a0 <- 1:max_age
      a5 <- 1:max_age


# DSTB deaths
ta["DSTB_deaths", i - 1, a5] <- (NTDS_II_u[a5] * la["NTDS_II", i - 1, a5]) +
                            (NTDS_IN_u[a5] * la["NTDS_IN", i - 1, a5]) +
                            (NTDS_T_u[a5] * la["NTDS_T", i - 1, a5]) +
                            (PTDS_II_u[a5] * la["PTDS_II", i - 1, a5]) +
                            (PTDS_IN_u[a5] * la["PTDS_IN", i - 1, a5]) +
                            (PTDS_T_u[a5] * la["PTDS_T", i - 1, a5])

# Backup calculation - DSTB deaths
dsma[current_step - 1, a5]  <-  (NTDS_II_u[a5] * la["NTDS_II", i - 1, a5]) +
                            (NTDS_IN_u[a5] * la["NTDS_IN", i - 1, a5]) +
                            (NTDS_T_u[a5] * la["NTDS_T", i - 1, a5]) +
                            (PTDS_II_u[a5] * la["PTDS_II", i - 1, a5]) +
                            (PTDS_IN_u[a5] * la["PTDS_IN", i - 1, a5]) +
                            (PTDS_T_u[a5] * la["PTDS_T", i - 1, a5])

# DRTB Deaths
ta["DRTB_deaths", i - 1, a5] <- (NTDR_II_u[a5] * la["NTDR_II", i - 1, a5]) +
                            (NTDS_IN_u[a5] * la["NTDR_IN", i - 1, a5]) +
                            (NTDR_T_u[a5] * la["NTDR_T_IIpsi", i - 1, a5]) +
                            (NTDR_II_u[a5] * la["NTDR_T_IIphi", i - 1, a5]) +
                            (NTDR_T_u[a5] * la["NTDR_T_INpsi", i - 1, a5]) +
                            (NTDR_IN_u[a5] * la["NTDR_T_INphi", i - 1, a5]) +
                            (PTDR_II_u[a5] * la["PTDR_II", i - 1, a5]) +
                            (PTDS_IN_u[a5] * la["PTDR_IN", i - 1, a5]) +
                            (PTDR_T_u[a5] * la["PTDR_T_IIpsi", i - 1, a5]) +
                            (PTDR_II_u[a5] * la["PTDR_T_IIphi", i - 1, a5]) +
                            (PTDR_T_u[a5] * la["PTDR_T_INpsi", i - 1, a5]) +
                            (PTDR_IN_u[a5] * la["PTDR_T_INphi", i - 1, a5]) +
                            (NTDR_mdt_u[a5] * la["NTDR_mdt", i - 1, a5])
                            (PTDR_mdt_u[a5] * la["PTDR_mdt", i - 1, a5])
                            (NTDR_nid_u[a5] * la["NTDR_nid", i - 1, a5])
                            (PTDR_nid_u[a5] * la["PTDR_nid", i - 1, a5])

# DSTB Incidence
ta["DSTB_inc", i - 1, a5] <-  ## Deconstructed
                          (NTDS_p[a5] * NTDS_f[a5] * NTDS_lambda[i - 1, a5] * la["S", i - 1, a5]) +
                          (NTDS_p[a5] * NTDS_f[a5] * NTDS_lambda[i - 1, a5] * NTDS_x * la["NTDR_L", i - 1, a5]) +
                          (NTDS_p[a5] * NTDS_f[a5] * NTDS_lambda[i - 1, a5] * NTDS_x * la["NTDS_L", i - 1, a5]) +
                          ## End deconstruct
                          ## Reactivation from Latent
                          (NTDS_v[a5] * NTDS_f[a5] * la["NTDS_L", i - 1, a5]) +
                          ## Reinfection of Resolved - DS
                          (NTDS_lambda[i - 1, a5] * NTDS_p[a5] * NTDS_f[a5] * NTDS_x * la["NTDS_R", i - 1, a5]) +
                          ## Reactivation from Resolved
                          (NTDS_f[a5] * NTDS_r[a5] * la["NTDS_R", i - 1, a5]) +
                          ## Reinfection of Resolved - DR
                          (NTDS_lambda[i - 1, a5] * NTDS_p[a5] * NTDS_f[a5] * NTDS_x * la["NTDR_R", i - 1, a5]) +
                          ## Deconstructed (NTDS)
                          (NTDS_p[a5] * (1 - NTDS_f[a5]) * NTDS_lambda[i - 1, a5] * la["S", i - 1, a5]) +
                          (NTDS_p[a5] * (1 - NTDS_f[a5]) * NTDS_lambda[i - 1, a5] * NTDS_x * la["NTDR_L", i - 1, a5]) +
                          (NTDS_p[a5] * (1 - NTDS_f[a5]) * NTDS_lambda[i - 1, a5] * NTDS_x * la["NTDS_L", i - 1, a5]) +
                          ## End deconstruct
                          ## Reactivation from Latent
                          (NTDS_v[a5] * (1 - NTDS_f[a5]) * la["NTDS_L", i - 1, a5]) +
                          ## Reinfection of Resolved - DS
                          (NTDS_lambda[i - 1, a5] * NTDS_p[a5] * (1 - NTDS_f[a5]) * NTDS_x * la["NTDS_R", i - 1,a5]) +
                          ## Reactivation from Resolved
                          ((1 - NTDS_f[a5]) * NTDS_r[a5] * la["NTDS_R", i-1,a5]) +
                          ## Reinfection of Resolved - DR
                          (NTDS_lambda[i - 1, a5] * NTDS_p[a5] * (1 - NTDS_f[a5]) * NTDS_x * la["NTDR_R", i - 1, a5]) +
                          ## Reinfection of Resolved - DS
                          (PTDS_lambda[i - 1, a5] * PTDS_p[a5] * PTDS_f[a5] * PTDS_x * la["PTDS_R", i - 1, a5]) +
                          ## Reactivation from Resolved
                          (PTDS_f[a5] * PTDS_r[a5] * la["PTDS_R", i - 1, a5]) +
                          ## Reinfection of Resolved - DR
                          (PTDS_lambda[i - 1, a5] * PTDS_p[a5] * PTDS_f[a5] * PTDS_x * la["PTDR_R", i - 1, a5]) +
                          ## Reinfection of Resolved - DS
                          (PTDS_lambda[i - 1, a5] * PTDS_p[a5] * (1 - PTDS_f[a5]) * PTDS_x * la["PTDS_R", i - 1, a5]) +
                          ## Reactivation from Resolved
                          ((1 - PTDS_f[a5]) * PTDS_r[a5] * la["PTDS_R", i - 1, a5]) +
                          ## Reinfection of Resolved - DR
                          (PTDS_lambda[i - 1, a5] * PTDS_p[a5] * (1 - PTDS_f[a5]) * PTDS_x * la["PTDR_R", i - 1, a5])

dsia[current_step - 1, a5]  <-  # Backup calculation
                            # New Infections of Susceptible and Latent
                            (NTDS_p[a5] * NTDS_lambda[i - 1, a5] * (la["S", i - 1, a5] + (NTDS_x * la["NTDR_L", i - 1, a5]) + (NTDS_x * la["NTDS_L", i - 1, a5]))) +
                            (NTDS_lambda[i - 1, a5] * NTDS_p[a5] * NTDS_x * (la["NTDS_R", i - 1, a5] + la["NTDR_R", i - 1, a5])) +
                            (NTDS_v[a5] * la["NTDS_L", i - 1, a5]) +
                            (NTDS_r[a5] * la["NTDS_R", i - 1, a5]) +
                            # PTDS
                            (PTDS_r[a5] * la["PTDS_R", i - 1, a5]) +
                            (PTDS_lambda[i - 1, a5] * PTDS_p[a5] * PTDS_x * (la["PTDR_R", i - 1, a5] + la["PTDS_R", i - 1, a5]))


dria[current_step - 1, a5] <-       # New infections of S/NTDS_L/NTDR_L
                                (NTDR_lambda[i - 1, a5] * NTDR_p[a5] * (la["S", i - 1, a5] + (NTDR_x * la["NTDR_L", i - 1, a5]) + (NTDR_x * la["NTDS_L", i - 1, a5]))) +
                                (NTDR_v[a5] * la["NTDR_L", i - 1, a5]) +
                                (NTDR_r[a5] * la["NTDR_R", i - 1, a5]) +
                                (NTDR_lambda[i - 1, a5] * NTDR_p[a5] * NTDR_x * (la["NTDR_R", i - 1, a5] + la["NTDS_R", i - 1, a5])) +
                                ## PT
                                (PTDR_v[a5] * la["PTDR_L", i - 1, a5]) + #ZERO COMPARTMENT
                                (PTDR_r[a5] * la["PTDR_R", i - 1, a5]) +
                                (PTDR_lambda[i - 1, a5] * PTDR_p[a5] * PTDR_x * (la["PTDS_R", i - 1, a5] + la["PTDR_R", i - 1, a5])) +
                                (NT_xi * PTDR_p[a5] * la["NTDS_T", i - 1, a5] * PTDR_f[a5] * nt_mdt_loss) +
                                (PT_xi * PTDR_p[a5] * la["PTDS_T", i - 1, a5] * PTDR_f[a5] * pt_mdt_loss) +
                                (NT_xi * PTDR_p[a5] * la["NTDS_T", i - 1, a5] * (1 - PTDR_f[a5]) * nt_nid_loss) +
                                (PT_xi * PTDR_p[a5] * la["PTDS_T", i - 1, a5] * (1 - PTDR_f[a5]) * pt_nid_loss)

# DRTB Incident Cases per timestep
ta["DRTB_inc", i - 1, a5]   <-  ## Deconstructed compound term - new infections in susceptble, NTDR-Latent and NTDS-Latent
                                (NTDR_p[a5] * NTDR_f[a5] * NTDR_lambda[i - 1, a5] * la["S", i - 1, a5]) +
                                (NTDR_p[a5] * NTDR_f[a5] * NTDR_lambda[i - 1, a5] * NTDR_x * la["NTDR_L", i - 1, a5]) +
                                (NTDR_p[a5] * NTDR_f[a5] * NTDR_lambda[i - 1, a5] * NTDR_x * la["NTDS_L", i - 1, a5]) +
                                ## End deconstruct
                                ## Reactivation from Latent
                                (NTDR_v[a5] * NTDR_f[a5] * la["NTDR_L", i - 1, a5]) +
                                ## Reinfection of Resolved - DR
                                (NTDR_lambda[i - 1, a5] * NTDR_p[a5] * NTDR_f[a5] * NTDR_x * la["NTDR_R", i - 1, a5]) +
                                ## Reactivation from Resolved
                                (NTDR_f[a5] * NTDR_r[a5] * la["NTDR_R", i - 1, a5]) +
                                ## Reinfection of Resolved - DS
                                (NTDR_lambda[i - 1, a5] * NTDR_p[a5] * NTDR_f[a5] * NTDR_x * la["NTDS_R", i - 1, a5]) +
                                # Deconstructed
                                (NTDR_p[a5] * (1 - NTDR_f[a5]) * NTDR_lambda[i - 1, a5] * la["S", i - 1, a5]) +
                                (NTDR_p[a5] * (1 - NTDR_f[a5]) * NTDR_lambda[i - 1, a5] * NTDR_x * la["NTDR_L", i - 1, a5]) +
                                (NTDR_p[a5] * (1 - NTDR_f[a5]) * NTDR_lambda[i - 1, a5] * NTDR_x * la["NTDS_L", i - 1, a5]) +
                                # End deconstruct
                                ## Reactivation from Latent
                                (NTDR_v[a5] * (1 - NTDR_f[a5]) * la["NTDR_L", i - 1, a5]) +
                                ## Reinfection of Resolved - DR
                                (NTDR_lambda[i - 1, a5] * NTDR_p[a5] * (1 - NTDR_f[a5]) * NTDR_x * la["NTDR_R", i - 1, a5]) +
                                ## Reactivation from Resolved
                                ((1 - NTDR_f[a5]) * NTDR_r[a5] * la["NTDR_R", i - 1, a5]) +
                                ## Reinfection of Resolved - DS
                                (NTDR_lambda[i - 1, a5] * NTDR_p[a5] * (1 - NTDR_f[a5]) * NTDR_x * la["NTDS_R", i - 1, a5]) +
                                ## Reactivation from Latent
                                (PTDR_v[a5] * PTDR_f[a5] * la["PTDR_L", i - 1, a5]) +
                                ## Reinfection of Resolved - DR
                                (PTDR_lambda[i - 1, a5] * PTDR_p[a5] * PTDR_f[a5] * PTDR_x * la["PTDR_R", i - 1, a5]) +
                                ## Reactivation from Resolved
                                (PTDR_f[a5] * (PTDR_r[a5] ) * la["PTDR_R", i - 1, a5]) +
                                ## Reinfection of Resolved - DS
                                (PTDR_lambda[i - 1, a5] * PTDR_p[a5] * PTDR_f[a5] * PTDR_x * la["PTDS_R", i - 1, a5]) +
                                ## Reactivation from Latent
                                (PTDR_v[a5] * (1 - PTDR_f[a5]) * la["PTDR_L", i - 1, a5]) +
                                ## Reinfection of Resolved - DS
                                (PTDR_lambda[i - 1, a5] * PTDR_p[a5] * (1 - PTDR_f[a5]) * PTDR_x * la["PTDS_R", i - 1, a5]) +
                                ## Reactivation from Resolved
                                ((1 - PTDR_f[a5]) * PTDR_r[a5] * la["PTDR_R", i - 1, a5]) +
                                ## Reinfection of Resolved - DR
                                (PTDR_lambda[i - 1, a5] * PTDR_p[a5] * (1 - PTDR_f[a5]) * PTDR_x * la["PTDR_R", i - 1, a5]) +
                                ## MDR Acquisition
                                (NT_xi * PTDR_p[a5] * la["NTDS_T", i - 1, a5] * PTDR_f[a5] * nt_mdt_loss) +
                                (PT_xi * PTDR_p[a5] * la["PTDS_T", i - 1, a5] * PTDR_f[a5] * pt_mdt_loss) +
                                (NT_xi * PTDR_p[a5] * la["NTDS_T", i - 1, a5] * (1 - PTDR_f[a5]) * nt_nid_loss) +
                                (PT_xi * PTDR_p[a5] * la["PTDS_T", i - 1, a5] * (1 - PTDR_f[a5]) * pt_nid_loss)

# ALL TB incidence backup calc
tbia[current_step - 1, a5]  <- dsia[current_step - 1, a5] + dria[current_step - 1, a5]

# ALL TB incidence
ta["TB_inc", i - 1, a5]     <- ta["DSTB_inc", i - 1, a5] + ta["DRTB_inc", i - 1, a5]

### Incidence by treatment history and bacteriologic status ###
## NTDR incidence
ta["DR_nt_inc", i - 1, a5]  <-  # New transmission by infection of susceptible and latent pools
                            (NTDR_p[a5] * NTDR_lambda[i - 1, a5] * la["S", i - 1, a5]) +
                            (NTDR_p[a5] * NTDR_lambda[i - 1, a5] * NTDR_x * la["NTDR_L", i - 1, a5]) +
                            (NTDR_p[a5] * NTDR_lambda[i - 1, a5] * NTDR_x * la["NTDS_L", i - 1, a5]) +
                            ## Reinfection of Resolved - DR
                            (NTDR_lambda[i - 1, a5] * NTDR_p[a5] * NTDR_x * la["NTDR_R", i - 1, a5]) +
                            ## Reinfection of Resolved - DS
                            (NTDR_lambda[i - 1, a5] * NTDR_p[a5] * NTDR_x * la["NTDS_R", i - 1, a5]) +
                            ## Reactivation from Latent
                            (NTDR_v[a5] * la["NTDR_L", i - 1, a5]) +
                            ## Reactivation from Resolved
                            (NTDR_r[a5] * la["NTDR_R", i - 1, a5])

# NTDR Bact+ Incidence
ta["DR_nt_inc_f", i - 1, a5] <- (ta["DR_nt_inc", i - 1, a5] * NTDR_f[a5]) +
                            (NTDR_omega * la["NTDR_IN", i - 1, a5])

# PTDR Incidence
ta["DR_pt_inc", i - 1, a5] <-   # Reinfection of Resolved - DR
                            (PTDR_lambda[i - 1, a5] * PTDR_p[a5] * PTDR_x * la["PTDR_R", i - 1, a5]) +
                            ## Reinfection of Resolved - DS
                            (PTDR_lambda[i - 1, a5] * PTDR_p[a5] * PTDR_x * la["PTDS_R", i - 1, a5]) +
                            ## Reactivation from Resolved
                            ((PTDR_r[a5]) * la["PTDR_R", i - 1, a5]) +
                            ## Reactivation from Latent
                            (PTDR_v[a5] * la["PTDR_L", i - 1, a5])
                            ## MDR Acquisition - not included.

# PTDR Bact+ Incidence
ta["DR_pt_inc_f", i - 1, a5] <- (ta["DR_pt_inc", i - 1, a5] * PTDS_f[a5]) +
                            (PTDR_omega * la["PTDR_IN", i - 1, a5]) +
                            # Bacteriologically positive MDR acquisition
                            (NT_xi * PTDR_p[a5] * la["NTDS_T", i - 1, a5] * PTDR_f[a5] * nt_mdt_loss) +
                            (PT_xi * PTDR_p[a5] * la["PTDS_T", i - 1, a5] * PTDR_f[a5] * pt_mdt_loss)

# NTDS Incidence
ta["DS_nt_inc", i - 1, a5] <- # Transmission by infection of suceptible and latent pools
                            (NTDS_p[a5] * NTDS_lambda[i - 1, a5] * la["S", i - 1, a5]) +
                            (NTDS_p[a5] * NTDS_lambda[i - 1, a5] * NTDS_x * la["NTDR_L", i - 1, a5]) +
                            (NTDS_p[a5] * NTDS_lambda[i - 1, a5] * NTDS_x * la["NTDS_L", i - 1, a5]) +
                            ## Reactivation from Latent
                            (NTDS_v[a5] * la["NTDS_L", i - 1, a5]) +
                            ## Reinfection of Resolved - DS
                            (NTDS_lambda[i - 1, a5] * NTDS_p[a5] * NTDS_x * la["NTDS_R", i - 1, a5]) +
                            ## Reactivation from Resolved
                            (NTDS_r[a5]  * la["NTDS_R", i - 1, a5]) +
                            ## Reinfection of Resolved - DR
                            (NTDS_lambda[i - 1, a5] * NTDS_p[a5] * NTDS_x * la["NTDR_R", i - 1, a5])

# NTDS Bact+ Incidence
ta["DS_nt_inc_f", i - 1, a5] <- (ta["DS_nt_inc", i - 1, a5] * NTDS_f[a5]) +
                            (NTDS_omega  * la["NTDS_IN", i - 1, a5])

# PTDS Incidence
ta["DS_pt_inc", i - 1, a5] <- ## Reinfection of Resolved - DS
                          (PTDS_lambda[i - 1, a5] * PTDS_p[a5] * PTDS_x * la["PTDS_R", i - 1, a5]) +
                          ## Reactivation from Resolved
                          (PTDS_r[a5] * la["PTDS_R", i - 1, a5]) +
                          ## Reinfection of Resolved - DR
                          (PTDS_lambda[i - 1, a5] * PTDS_p[a5] * PTDS_x * la["PTDR_R", i - 1, a5])

# PTDS Bact+ Incidence
ta["DS_pt_inc_f", i - 1, a5] <- (ta["DS_pt_inc", i - 1, a5] * PTDS_f[a5]) +
                            (PTDS_omega * la["PTDS_IN", i - 1, a5])

# AllTB Never Treated Incidence
ta["All_nt_inc", i - 1, a5] <- ta["DR_nt_inc", i - 1, a5] + ta["DS_nt_inc", i - 1, a5]

# AllTB Previously Treated Incidence
ta["All_pt_inc", i - 1, a5] <- ta["DR_pt_inc", i - 1, a5] + ta["DS_pt_inc", i - 1, a5]

# AllTB Never Treated Bact+ Incidence
ta["All_nt_inc_f", i - 1, a5] <- ta["DR_nt_inc_f", i - 1, a5] + ta["DS_nt_inc_f", i - 1, a5]

# AllTB Previously Treated Bact+ Incidence
ta["All_pt_inc_f", i - 1, a5] <- ta["DR_pt_inc_f", i - 1, a5] + ta["DS_pt_inc_f", i - 1, a5]

############################ DRTB Treatment Initiations [Notifications] #################################

##########################
########## TZ ############

# 'Fitting factor', 'ff' for diagnostics only.
# ff = 1 i.e. no effect during normal model functioning.
ff <- 1

##########################
##########################

# NTDR Treatment Initiations to any destination - infectious
ta["DR_nt_intx_f", i - 1, a5]  <- (NTDR_II_kappa[a5] * NTDR_II_psi * la["NTDR_II", i - 1, a5]) +
                                  (NTDR_II_kappa[a5] * NTDR_II_phi * la["NTDR_II", i - 1, a5]) +
                                  (NTDR_II_mdt[a5] * la["NTDR_II", i - 1, a5])

# PTDR Treatment Initiations to any destination - infectious
ta["DR_pt_intx_f", i - 1, a5]  <- (PTDR_II_kappa[a5] * PTDR_II_psi * la["PTDR_II", i - 1, a5]) +
                                  (PTDR_II_kappa[a5] * PTDR_II_phi * la["PTDR_II", i - 1, a5]) +
                                  (PTDR_II_mdt[a5] * la["PTDR_II", i - 1, a5]) +
                                  (ff * (la["NTDR_mdt", i - 1, a5] * NTDR_mdt_exit * (1 - nt_mdt_loss) * PTDR_II_psi)) +
                                  (ff * (la["NTDR_mdt", i - 1, a5] * NTDR_mdt_exit * (1 - nt_mdt_loss) * PTDR_II_phi)) +
                                  (ff * (la["PTDR_mdt", i - 1, a5] * PTDR_mdt_exit * (1 - pt_mdt_loss) * PTDR_II_psi)) +
                                  (ff * (la["PTDR_mdt", i - 1, a5] * PTDR_mdt_exit * (1 - pt_mdt_loss) * PTDR_II_phi)) +
                                  (ff * (NT_xi * PTDR_p[a5] * PTDR_f[a5] * la["NTDS_T", i - 1, a5] * (1 - nt_mdt_loss) * PTDR_II_psi)) +
                                  (ff * (PT_xi * PTDR_p[a5] * PTDR_f[a5] * la["PTDS_T", i - 1, a5] * (1 - pt_mdt_loss) * PTDR_II_psi)) +
                                  (ff * (NT_xi * PTDR_p[a5] * PTDR_f[a5] * la["NTDS_T", i - 1, a5] * (1 - nt_mdt_loss) * PTDR_II_phi)) +
                                  (ff * (PT_xi * PTDR_p[a5] * PTDR_f[a5] * la["PTDS_T", i - 1, a5] * (1 - pt_mdt_loss) * PTDR_II_phi))

# NTDS Treatment initiations to any destination - infectious
ta["DS_nt_intx_f", i - 1, a5]  <- (NTDS_II_kappa[a5] * la["NTDS_II", i - 1, a5])

# PTDS Treatment initiations to any destination - infectious
ta["DS_pt_intx_f", i - 1, a5]  <- (PTDS_II_kappa[a5] * la["PTDS_II", i - 1, a5])

# All-Never-Treated TB Treatment initiations to any destination - infectious
ta["All_nt_intx_f", i - 1, a5] <- (ta["DR_nt_intx_f", i - 1, a5] + ta["DS_nt_intx_f", i - 1, a5])

# All-Previously-Treated TB Treatment initiations to any destination - infectious
ta["All_pt_intx_f", i - 1, a5] <- (ta["DR_pt_intx_f", i - 1, a5] + ta["DS_pt_intx_f", i - 1, a5])

# Treatment Initations by _type of regimen_
## Treatment Initiations of DSTB Regimen
ta["DSTB_initRx", i - 1, a5] <- (NTDS_II_kappa[a5] * la["NTDS_II", i - 1, a5]) +
                                (NTDS_IN_kappa[a5] * la["NTDS_IN", i - 1, a5]) +
                                (PTDS_II_kappa[a5] * la["PTDS_II", i - 1, a5]) +
                                (PTDS_IN_kappa[a5] * la["PTDS_IN", i - 1, a5]) +
                                (NTDR_II_mdt[a5] * la["NTDR_II", i - 1, a5]) +
                                (NTDR_IN_nid[a5] * la["NTDR_IN", i - 1, a5]) +
                                (PTDR_II_mdt[a5] * la["PTDR_II", i - 1, a5]) +
                                (PTDR_IN_nid[a5] * la["PTDR_IN", i - 1, a5])

## Treatment Initiations of DRTB Regimen
ta["DRTB_initRx", i - 1, a5] <- (NTDR_II_kappa[a5] * NTDR_II_psi * la["NTDR_II", i - 1, a5]) +
                                (NTDR_II_kappa[a5] * NTDR_II_phi * la["NTDR_II", i - 1, a5]) +
                                (NTDR_IN_kappa[a5] * NTDR_IN_psi * la["NTDR_IN", i - 1, a5]) +
                                (NTDR_IN_kappa[a5] * NTDR_IN_phi * la["NTDR_IN", i - 1, a5]) +
                                (PTDR_II_kappa[a5] * PTDR_II_psi * la["PTDR_II", i - 1, a5]) +
                                (PTDR_II_kappa[a5] * PTDR_II_phi * la["PTDR_II", i - 1, a5]) +
                                (PTDR_IN_kappa[a5] * PTDR_IN_psi * la["PTDR_IN", i - 1, a5]) +
                                (PTDR_IN_kappa[a5] * PTDR_IN_phi * la["PTDR_IN", i - 1, a5]) +
                                (la["NTDR_mdt", i - 1, a5] * NTDR_mdt_exit * (1 - nt_mdt_loss) * PTDR_II_psi) +
                                (la["PTDR_mdt", i - 1, a5] * PTDR_mdt_exit * (1 - pt_mdt_loss) * PTDR_II_psi) +
                                (la["NTDR_mdt", i - 1, a5] * NTDR_mdt_exit * (1 - nt_mdt_loss) * PTDR_II_phi) +
                                (la["PTDR_mdt", i - 1, a5] * PTDR_mdt_exit * (1 - pt_mdt_loss) * PTDR_II_phi) +
                                (la["NTDR_nid", i - 1, a5] * NTDR_nid_exit * (1 - nt_nid_loss) * PTDR_IN_psi) +
                                (la["PTDR_nid", i - 1, a5] * PTDR_nid_exit * (1 - pt_nid_loss) * PTDR_IN_psi) +
                                (la["NTDR_nid", i - 1, a5] * NTDR_nid_exit * (1 - nt_mdt_loss) * PTDR_IN_phi) +
                                (la["PTDR_nid", i - 1, a5] * PTDR_nid_exit * (1 - pt_mdt_loss) * PTDR_IN_phi) +
                                (NT_xi * PTDR_p[a5] * PTDR_f[a5] * la["NTDS_T", i - 1, a5] * (1 - nt_mdt_loss) * PTDR_II_psi) +
                                (PT_xi * PTDR_p[a5] * PTDR_f[a5] * la["PTDS_T", i - 1, a5] * (1 - pt_mdt_loss) * PTDR_II_psi) +
                                (NT_xi * PTDR_p[a5] * PTDR_f[a5] * la["NTDS_T", i - 1, a5] * (1 - nt_mdt_loss) * PTDR_II_phi) +
                                (PT_xi * PTDR_p[a5] * PTDR_f[a5] * la["PTDS_T", i - 1, a5] * (1 - pt_mdt_loss) * PTDR_II_phi) +
                                (NT_xi * PTDR_p[a5] * (1 - PTDR_f[a5]) * la["NTDS_T", i - 1, a5] * (1 - nt_nid_loss) * PTDR_IN_psi) +
                                (PT_xi * PTDR_p[a5] * (1 - PTDR_f[a5]) * la["PTDS_T", i - 1, a5] * (1 - pt_nid_loss) * PTDR_IN_psi) +
                                (NT_xi * PTDR_p[a5] * (1 - PTDR_f[a5]) * la["NTDS_T", i - 1, a5] * (1 - nt_nid_loss) * PTDR_IN_phi) +
                                (PT_xi * PTDR_p[a5] * (1 - PTDR_f[a5]) * la["PTDS_T", i - 1, a5] * (1 - pt_nid_loss) * PTDR_IN_phi)

## Laboratory Confirmed Treatment Initiations of DRTB Regimen
ta["DRTB_initRxLab", i - 1, a5] <-(NTDR_II_kappa[a5] * NTDR_II_psi * la["NTDR_II", i - 1, a5] * nt_dst_p[current_year - year1  + 1]) +
                                  (NTDR_II_kappa[a5] * NTDR_II_phi * la["NTDR_II", i - 1, a5] * nt_dst_p[current_year - year1  + 1]) +
                                  (PTDR_II_kappa[a5] * PTDR_II_psi * la["PTDR_II", i - 1, a5] * pt_dst_p[current_year - year1  + 1]) +
                                  (PTDR_II_kappa[a5] * PTDR_II_phi * la["PTDR_II", i - 1, a5] * pt_dst_p[current_year - year1  + 1]) +
                                  (la["NTDR_mdt", i - 1, a5] * NTDR_mdt_exit * nt_dst_p[current_year - year1  + 1] * PTDR_II_psi) +
                                  (la["PTDR_mdt", i - 1, a5] * PTDR_mdt_exit * pt_dst_p[current_year - year1  + 1] * PTDR_II_psi) +
                                  (la["NTDR_mdt", i - 1, a5] * NTDR_mdt_exit * nt_dst_p[current_year - year1  + 1] * PTDR_II_phi) +
                                  (la["PTDR_mdt", i - 1, a5] * PTDR_mdt_exit * pt_dst_p[current_year - year1  + 1] * PTDR_II_phi) +
                                  (NT_xi * PTDR_p[a5] * PTDR_f[a5] * la["NTDS_T", i - 1, a5] * nt_dst_p[current_year - year1  + 1]) +
                                  (PT_xi * PTDR_p[a5] * PTDR_f[a5] * la["PTDS_T", i - 1, a5] * pt_dst_p[current_year - year1  + 1]) +
                                  (NT_xi * PTDR_p[a5] * PTDR_f[a5] * la["NTDS_T", i - 1, a5] * nt_dst_p[current_year - year1  + 1] * PTDR_II_phi) +
                                  (PT_xi * PTDR_p[a5] * PTDR_f[a5] * la["PTDS_T", i - 1, a5] * pt_dst_p[current_year - year1  + 1] * PTDR_II_phi)

# DSTB person-time on treatment - ** IN PERSON-MONTHS **
ta["DSTB_onRx", i - 1, a5] <- (12 * dt) *  (la["NTDS_T", i - 1, a5] +
                                        la["PTDS_T", i - 1, a5] +
                                        la["NTDR_mdt", i - 1, a5] +
                                        la["PTDR_mdt", i - 1, a5] +
                                        la["NTDR_nid", i - 1, a5] +
                                        la["PTDR_nid", i - 1, a5])

# DRTB person-time on treatment ** IN PERSON-MONTHS**
ta["DRTB_onRx", i - 1, a5] <- (12 * dt) *  (la["NTDR_T_IIphi", i - 1, a5] +
                                        la["NTDR_T_IIpsi", i - 1, a5] +
                                        la["NTDR_T_INphi", i - 1, a5] +
                                        la["NTDR_T_INpsi", i - 1, a5] +
                                        la["PTDR_T_IIphi", i - 1, a5] +
                                        la["PTDR_T_IIpsi", i - 1, a5] +
                                        la["PTDR_T_INphi", i - 1, a5] +
                                        la["PTDR_T_INpsi", i - 1, a5])

# DSTB Diagnostic Cost
ca["ds_dx", i - 1, a5] <- (ta["DSTB_initRx", i - 1, a5] * sdr[i - 1] * ds_dx_cost) * 1e+03

# DRTB Diagnostic Cost (Excl DST)
ca["dr_dx", i - 1, a5] <- (ta["DRTB_initRx", i - 1, a5] * sdr[i - 1] * dr_dx_cost) * 1000

# DST Cost
ca["dst", i - 1, a5]   <- (ta["DRTB_initRxLab", i - 1, a5] * dst_cost) * 1000

# DSTB Treatment Cost
ca["ds_tx", i - 1, a5] <- (ta["DSTB_onRx", i - 1, a5] * ds_tx_cost) * 1000

# DRTB Treatment Cost
ca["dr_tx", i - 1, a5] <- (ta["DRTB_onRx", i - 1, a5] * dr_tx_cost) * 1000

# Total treatment cost including fractional inflation for programme cost
ca["tbrx", i - 1, a5]  <- #
                        (ca["ds_dx", i - 1, a5] +
                        ca["dr_dx", i - 1, a5] +
                        ca["dst", i - 1, a5] +
                        ca["ds_tx", i - 1, a5] +
                        ca["dr_tx", i - 1, a5]) * (1 + prog_cost)

      # Recalculate matrices for infection parameters
      # Popsize
      psizematrix[i, 1] <- sum(la[, i, 1:6])
      psizematrix[i, 2] <- sum(la[, i, 7:20])
      psizematrix[i, 3] <- sum(la[, i, 21:65])
      psizematrix[i, 4] <- sum(la[, i, 66:max_age])

      ## Total Infectious Cases by contact matrix age classes - DSTB
      DS_Imatrix[i, 1]  <- sum(la["PTDS_II", i, 1:6], la["NTDS_II", i, 1:6])
      DS_Imatrix[i, 2]  <- sum(la["PTDS_II", i, 7:20], la["NTDS_II", i, 7:20])
      DS_Imatrix[i, 3]  <- sum(la["PTDS_II", i, 21:65], la["NTDS_II", i, 21:65])
      DS_Imatrix[i, 4]  <- sum(la["PTDS_II", i, 66:max_age], la["NTDS_II", i, 66:max_age])

      ## Total Infectious Cases by contact matrix age classes - DRTB
      DR_Imatrix[i, 1]  <- sum(la["PTDR_II", i, 1:6], la["NTDR_II", i, 1:6], la["PTDR_T_IIphi", i, 1:6], la["NTDR_T_IIphi", i, 1:6], la["NTDR_mdt", i, 1:6], la["PTDR_mdt", i, 1:6])
      DR_Imatrix[i, 2]  <- sum(la["PTDR_II", i, 7:20], la["NTDR_II", i, 7:20], la["PTDR_T_IIphi", i, 7:20], la["NTDR_T_IIphi", i, 7:20], la["NTDR_mdt", i, 7:20], la["PTDR_mdt", i, 7:20])
      DR_Imatrix[i, 3]  <- sum(la["PTDR_II", i, 21:65], la["NTDR_II", i, 21:65], la["PTDR_T_IIphi", i, 21:65], la["NTDR_T_IIphi", i, 21:65], la["NTDR_mdt", i, 21:65], la["PTDR_mdt", i, 21:65])
      DR_Imatrix[i, 4]  <- sum(la["PTDR_II", i, 66:max_age], la["NTDR_II", i, 66:max_age], la["PTDR_T_IIphi", i, 66:max_age], la["NTDR_T_IIphi", i, 66:max_age], la["NTDR_mdt", i, 66:max_age], la["PTDR_mdt", i, 66:max_age])

    }

    # Store yearly cdr/kappa arrays
    # hcdr[current_year - year1 +1, ]       <- NTDS.CDR.II
    hkappa[current_year - year1 + 1, ]      <- NTDS_II_kappa
    hntdr_kappa[current_year - year1 + 1, ] <- NTDR_II_kappa
    hptdr_kappa[current_year - year1 + 1, ] <- PTDR_II_kappa
    hcdrscaling[current_year - year1 + 1, ] <- c(DS_CDRscale, DS_CDRscaleO, DS_CDRscaleE)
    hn[current_year - year1 + 1, ]          <- NTDS_II_n
    hui[current_year - year1 + 1, ]         <- NTDS_II_u

  }


  # Common core calculations

  # Mode specific calculations and output
  if (mode == 0) {
# Local calculations

# 1. Housekeeping
# 1.1 New holding environment for results
calc                                 <- new.env()
# 1.2 Working time horizons
calc$wth                             <- 1950:yearend
# 1.3 Indices
# Indices: Susceptible, Recovered and Latent [i.e. "Non-Prevalent"] Compartments
calc$nonprev                         <- match(c("S", "NTDS_R", "PTDS_R", "NTDR_R", "PTDR_R", "NTDR_L", "NTDS_L", "PTDR_L"), unlist(dimnames(la)[1]))
# Indices: DR compartments
calc$dr_comps                        <- grep("DR", unlist(dimnames(la)[1]))
# Indices: DS compartments
calc$ds_comps                        <- grep("DS", unlist(dimnames(la)[1]))

# 2. Demographics - annual age-wise population
calc$ann_agewise                     <- annual_mean(colSums(la, dims = 1), dt, year1, yearend)[-c(1:50), ]

# 3. DSTB
# 3.1 DSTB Prevalence
# Annual average DSTB prevalence, Infectious TB only, raw
calc$ann_ds_prev_ii                  <- (annual_mean(colSums(la[c("NTDS_II", "PTDS_II"), , ], dims = 1), dt, year1, yearend))[-c(1:50), ]
# Annual average DSTB prevalence by age **group**, Infectious TB only, population normalised
calc$ann_ds_prev_norm_ii             <- as.data.frame(cbind(calc$wth, psadj(calc$ann_ds_prev_ii, 1:100), psadj(calc$ann_ds_prev_ii, 16:100), psadj(calc$ann_ds_prev_ii, 1:15), psadj(calc$ann_ds_prev_ii, 16:30), psadj(calc$ann_ds_prev_ii, 31:45), psadj(calc$ann_ds_prev_ii, 46:60), psadj(calc$ann_ds_prev_ii, 61:100)))

# 3.2 DSTB Incidence
# Annual DSTB incidence by age, raw
calc$ann_dstb_inc                    <- annual_sum(dsia[201:800, ], dt, 1950, yearend)
# Annual DSTB incidence by age **group**, population normalised
calc$ann_dstb_inc_norm               <- as.data.frame(cbind(calc$wth, psadj(calc$ann_dstb_inc, 1:100), psadj(calc$ann_dstb_inc, 16:100), psadj(calc$ann_dstb_inc, 1:15), psadj(calc$ann_dstb_inc, 16:65), psadj(calc$ann_dstb_inc, 66:100)))

# 3.3 DSTB Mortality
# Annual DSTB mortality by age, raw
calc$ann_dstb_mort                   <- annual_sum(dsma[201:800,], dt, 1950, yearend)
# Annual DSTB mortality by age **group**, population normalised
calc$ann_dstb_mort_norm              <- as.data.frame(cbind(calc$wth, psadj(calc$ann_dstb_mort, 1:100),psadj(calc$ann_dstb_mort, 16:100), psadj(calc$ann_dstb_mort, 1:15), psadj(calc$ann_dstb_mort, 16:65), psadj(calc$ann_dstb_mort, 66:100)))

# 4. DRTB
# 4.1 DRTB Incidence
# Annual DRTB incidence by age, raw
calc$ann_drtb_inc                    <- annual_sum(ta["DRTB_inc", 201:800, ], dt, 1950, yearend)
# Annual DRTB incidence by age **group**, population-normalised
calc$ann_drtb_inc_norm               <- as.data.frame(cbind(calc$wth, psadj(calc$ann_drtb_inc, 1:100), psadj(calc$ann_drtb_inc, 16:100), psadj(calc$ann_drtb_inc, 1:15), psadj(calc$ann_drtb_inc, 16:65), psadj(calc$ann_drtb_inc, 66:100)))

# 4.2 DRTB Mortality
calc$ann_drtb_mort <- annual_sum(ta["DRTB_deaths", 201:800, ], dt, 1950, yearend)
# Annual drtb mortality by age **group**, population normalised
calc$ann_drtb_mort_norm <- as.data.frame(cbind(calc$wth, psadj(calc$ann_drtb_mort, 1:100), psadj(calc$ann_drtb_mort, 16:100), psadj(calc$ann_drtb_mort, 1:15), psadj(calc$ann_drtb_mort, 16:65), psadj(calc$ann_drtb_mort, 66:100)))

# 4.3 DRTB Notifications
# DRTB Treatment Initiations - Including Misdiagnoses converted to correct treatment
calc$ann_drtb_notif               <- annual_sum(1e03 * rowSums(ta["DRTB_initRx", 201:800, ]), dt, 1950, yearend)
# DRTB Treatment Initiations - Only transitions from Prevalent Compartments to Treatment Compartments
calc$ann_drtb_notif_lab           <- annual_sum(1e03 * rowSums(ta["DRTB_initRxLab", 201:800, ]), dt, 1950, yearend)

# 4.4 DRTB Proportional Incidence
tot_nt_inc_f                                            <- annual_sum(rowSums(ta["All_nt_inc_f", , ]), dt, year1, yearend)
tot_pt_inc_f                                            <- annual_sum(rowSums(ta["All_pt_inc_f", , ]), dt, year1, yearend)
mdr_nt_inc_f                                            <- annual_sum(rowSums(ta["DR_nt_inc_f", , ]), dt, year1, yearend)
mdr_pt_inc_f                                            <- annual_sum(rowSums(ta["DR_pt_inc_f", , ]), dt, year1, yearend)

calc$mdr_nt_iprop_f                                      <- mdr_nt_inc_f[,2] / tot_nt_inc_f[,2]
calc$mdr_pt_iprop_f                                      <- mdr_pt_inc_f[,2] / tot_pt_inc_f[,2]

if (any(is.nan(calc$mdr_nt_iprop_f))) {
  calc$mdr_nt_iprop_f[which(is.nan(calc$mdr_nt_iprop_f))] <- 0
}

if (any(is.nan(calc$mdr_pt_iprop_f))) {
  calc$mdr_pt_iprop_f[which(is.nan(calc$mdr_pt_iprop_f))] <- 0
}

# MDR Incident Proportion by Treatment History
calc$mdr_nt_iprop_f <- as.data.frame(cbind(mdr_nt_inc_f[, 1], calc$mdr_nt_iprop_f)[-c(1:50), ])
calc$mdr_pt_iprop_f <- as.data.frame(cbind(mdr_pt_inc_f[, 1], calc$mdr_pt_iprop_f)[-c(1:50), ])

# MDR Bacteriologically Positive Incidence by Treatment History
calc$mdr_nt_inc_f <- as.data.frame(cbind(calc$wth, (mdr_nt_inc_f[-c(1:50), 2] / rowSums(calc$ann_agewise[, -1]) * 1e05)))
calc$mdr_pt_inc_f <- as.data.frame(cbind(calc$wth, (mdr_pt_inc_f[-c(1:50), 2] / rowSums(calc$ann_agewise[, -1]) * 1e05)))

# 4.5 DRTB Proportional Notifications
tot_nt_intx_f                                            <- annual_sum(rowSums(ta["All_nt_intx_f", , ]), dt, year1, yearend)
tot_pt_intx_f                                            <- annual_sum(rowSums(ta["All_pt_intx_f", , ]), dt, year1, yearend)
mdr_nt_intx_f                                            <- annual_sum(rowSums(ta["DR_nt_intx_f", , ]), dt, year1, yearend)
mdr_pt_intx_f                                            <- annual_sum(rowSums(ta["DR_pt_intx_f", , ]), dt, year1, yearend)

calc$mdr_nt_nprop_f                                      <- mdr_nt_intx_f[,2] / tot_nt_intx_f[,2]
calc$mdr_pt_nprop_f                                      <- mdr_pt_intx_f[,2] / tot_pt_intx_f[,2]

if (any(is.nan(calc$mdr_nt_nprop_f))) {
  calc$mdr_nt_nprop_f[which(is.nan(calc$mdr_nt_nprop_f))] <- 0
}

if (any(is.nan(calc$mdr_pt_nprop_f))) {
  calc$mdr_pt_nprop_f[which(is.nan(calc$mdr_pt_nprop_f))] <- 0
}

# MDR Notification Proportion by Treatment History
calc$mdr_nt_nprop_f <- as.data.frame(cbind(mdr_nt_intx_f[, 1], calc$mdr_nt_nprop_f)[-c(1:50), ])
calc$mdr_pt_nprop_f <- as.data.frame(cbind(mdr_pt_intx_f[, 1], calc$mdr_pt_nprop_f)[-c(1:50), ])

calc$mdr_nt_intx_f <- as.data.frame(cbind(calc$wth, (mdr_nt_intx_f[-c(1:50), 2] / rowSums(calc$ann_agewise[, -1]) * 1e05)))
calc$mdr_pt_intx_f <- as.data.frame(cbind(calc$wth, (mdr_pt_intx_f[-c(1:50), 2] / rowSums(calc$ann_agewise[, -1]) * 1e05)))

# 5. Health Economics - Costing
# 5.1 DSTB
# Person-time on treatment, by year in 3 month units, DSTB by age
calc$ann_dstb_rx_dtu                 <- (annual_sum(ta["DSTB_onRx", , ], dt, year1, yearend))[-c(1:50), ]
# Person-time on treatmentn, by year in 1 month units, DSTB by age
calc$ann_dstb_rx_mo                  <- cbind(calc$ann_dstb_rx_dtu[, 1], (calc$ann_dstb_rx_dtu[, -1] * 3))
# Treatment initiations per year - DSTB by age
calc$ann_dstb_rx_init                <- (annual_sum(ta["DSTB_initRx", , ], dt, year1, yearend))[-c(1:50), ]

# 5.2 DRTB
# Person-time on treatment, by year in 3 month units, DRTB by age
calc$ann_drtb_rx_dtu                 <- annual_sum(ta["DRTB_onRx", , ], dt, year1, yearend)[c(1:50), ]
# Person-time on treatmentn, by year in 1 month units, DRTB
calc$ann_drtb_rx_mo                  <- cbind(calc$ann_drtb_rx_dtu[, 1], (calc$ann_drtb_rx_dtu[, -1] * 3))
# Treatment initiations per year - DRTB
calc$ann_drtb_rx_init                <- (annual_sum(ta["DRTB_initRx", , ], dt, year1, yearend))[-c(1:50), ]

# 5.3 Total Cost
# calc$ann_total_cost <- annual_sum(rowSums(ta["total_cost", , ]), dt, year1, yearend)[-c(1:50), ]
calc$cost_array                     <- t(rowSums(ca[, , ], dims = 2))
calc$ann_cost_array                 <- annual_sum(calc$cost_array, dt, year1, yearend)[-c(1:50), ]

# 6 ALL TB
# 6.1 ALL TB Prevalence
# Annual all TB prevalence, by age, infectious TB only
calc$ann_alltb_prev_ii                  <- (annual_mean(colSums(la[c("NTDS_II", "PTDS_II", "NTDR_II", "PTDR_II"), , ], dims = 1), dt, year1, yearend))[-c(1:50), ]
# Annual all TB prevalence by age **group**, Infectious TB only, population normalised
calc$ann_alltb_prev_norm_ii             <- as.data.frame(cbind(calc$wth, psadj(calc$ann_alltb_prev_ii, 1:100), psadj(calc$ann_alltb_prev_ii, 16:100), psadj(calc$ann_alltb_prev_ii, 1:15), psadj(calc$ann_alltb_prev_ii, 16:30), psadj(calc$ann_alltb_prev_ii, 31:45), psadj(calc$ann_alltb_prev_ii, 46:60), psadj(calc$ann_alltb_prev_ii, 61:100)))
# Annual All TB incidence, by age
calc$ann_alltb_inc <- cbind(calc$ann_dstb_inc[, 1], (calc$ann_dstb_inc[, -1] + calc$ann_drtb_inc[, -1]))
# Annual all TB incidence by age **group**, population-normalised
calc$ann_alltb_inc_norm <- as.data.frame(cbind(calc$wth, psadj(calc$ann_alltb_inc, 1:100), psadj(calc$ann_alltb_inc, 16:100), psadj(calc$ann_alltb_inc, 1:15), psadj(calc$ann_alltb_inc, 16:65), psadj(calc$ann_alltb_inc, 66:100)))

# 6.3 ALL TB Mortality Annual All TB mortality, by age
calc$ann_alltb_mort <- cbind(calc$ann_dstb_mort[, 1], (calc$ann_dstb_mort[, -1] + calc$ann_drtb_mort[, -1]))
# Annual All TB mortality, by age **group**, population-normalised
calc$ann_alltb_mort_norm <- as.data.frame(cbind(calc$wth, psadj(calc$ann_alltb_mort, 1:100), psadj(calc$ann_alltb_mort, 16:100), psadj(calc$ann_alltb_mort, 1:15), psadj(calc$ann_alltb_mort, 16:65), psadj(calc$ann_alltb_mort, 66:100)))

# 6.4 ALL TB Notifications
# Raw allTB notifications
calc$alltb_notif <- ta["DSTB_initRx", , ] + ta["DRTB_initRx", , ]
# Annual all TB notifications by Age
calc$ann_alltb_notif <- annual_sum(calc$alltb_notif[201:800,], dt, 1950, yearend)
# Annual all TB notifications by age **group**, population normalised
calc$ann_alltb_notif_norm <- as.data.frame(cbind(calc$wth, psadj(calc$ann_alltb_notif, 1:100), psadj(calc$ann_alltb_notif, 16:100), psadj(calc$ann_alltb_notif, 1:15), psadj(calc$ann_alltb_notif, 16:65), psadj(calc$ann_alltb_notif, 66:100)))

# 7 Miscellaneous
# 7.1 Scaled Mortality terms
calc$scaled_ui                       <- ui
calc$scaled_uni                      <- uni
calc$scaled_ut                       <- ut

# 7.2 Annual age-wise background mortality
calc$ann_bgmort                      <- annual_mean(ta["bgmort", , ], dt, year1, yearend)[-c(1:50), ]

# 7.3 Average births per year
calc$ann_births                      <- annual_mean(colSums(la[, , 1], dims = 1), dt, year1, yearend)[-c(1:50), ]

# 7.4 Annual compartment-wise population
compartments                         <- t(rowSums(la, dim = 2))
calc$ann_compartments                <- (annual_mean(compartments, dt, year1, yearend))[-c(1:50), ]
rm(compartments)

# Table colnames

### Age class names
calc$prevnames                       <- c("Year", "All", "AllAdults", "A0-14", "A15-29", "A30-44", "A45-59", "A60+")
calc$incnames                        <- c("Year", "All", "AllAdults", "A0-14", "A15-64", "A65+")

### Column names
colnames(calc$ann_ds_prev_norm_ii)     <- calc$prevnames
colnames(calc$ann_alltb_prev_norm_ii)  <- calc$prevnames
colnames(calc$ann_dstb_inc_norm)       <- calc$incnames
colnames(calc$ann_dstb_mort_norm)      <- calc$incnames
colnames(calc$ann_drtb_mort_norm)      <- calc$incnames
colnames(calc$ann_alltb_mort_norm)     <- calc$incnames
colnames(calc$ann_alltb_notif_norm)    <- calc$incnames
colnames(calc$ann_alltb_inc_norm)      <- calc$incnames
colnames(calc$ann_drtb_inc_norm)       <- calc$incnames
colnames(calc$ann_drtb_notif)          <- c("Year", "DRTB Notifications")
colnames(calc$ann_drtb_notif_lab)      <- c("Year", "DRTB Notifications (Lab)")
colnames(calc$mdr_nt_iprop_f)          <- c("Year", "NTDR_inc_prop")
colnames(calc$mdr_pt_iprop_f)          <- c("Year", "PTDR_inc_prop")
colnames(calc$mdr_nt_nprop_f)          <- c("Year", "NTDR_not_prop")
colnames(calc$mdr_pt_nprop_f)          <- c("Year", "PTDR_not_prop")
colnames(calc$mdr_nt_intx_f)           <- c("Year", "NTDR Notifications")
colnames(calc$mdr_pt_intx_f)           <- c("Year", "PTDR Notifications")
colnames(calc$mdr_nt_inc_f)            <- c("Year", "NTDR Incidence")
colnames(calc$mdr_pt_inc_f)            <- c("Year", "PTDR Incidence")
colnames(calc$ann_cost_array)          <- c("Year", "DS Diagnosis", "DR Diagnosis", "DST", "DS Treatment", "DR Treatment", "Total")

### Detailed outputs

output_pop                           <- list2env(list(
  la = la,
  ta = ta,
  DS_Imatrix = DS_Imatrix,
  DR_Imatrix = DR_Imatrix,
  psizematrix = psizematrix,
  dsia = dsia,
  dria = dria,
  tbia = tbia),
  hash = T
)

output_internals                     <- list2env(list(
  NTDS_lambda = NTDS_lambda,
  PTDS_lambda = PTDS_lambda,
  NTDR_lambda = NTDR_lambda,
  PTDR_lambda = PTDR_lambda,
  NTDS_lambda_raw = NTDS_lambda_raw,
  NTDR_lambda_raw = NTDR_lambda_raw,
  DS_CDR = DS_CDR,
  NTDS_II_kappa = NTDS_II_kappa,
  NTDS_IN_kappa = NTDS_IN_kappa,
  PTDS_II_kappa = PTDS_II_kappa,
  PTDS_IN_kappa = PTDS_IN_kappa,
  NTDR_II_kappa = NTDR_II_kappa,
  NTDR_IN_kappa = NTDR_IN_kappa,
  PTDR_II_kappa = PTDR_II_kappa,
  PTDR_IN_kappa = PTDR_IN_kappa,
  scaling_table = scaling_table,
  age_props = age_props,
  NTDS_II_u = NTDS_II_u,
  PTDS_II_u = PTDS_II_u,
  NTDS_IN_u = NTDS_IN_u,
  PTDS_IN_u = PTDS_IN_u,
  NTDS_II_n = NTDS_II_n,
  PTDS_II_n = PTDS_II_n,
  NTDS_IN_n = NTDS_IN_n,
  PTDS_IN_n = PTDS_IN_n,
  hcdr = hcdr,
  hkappa = hkappa,
  death_rate = death_rate,
  nt_dst_p = nt_dst_p,
  pt_dst_p = pt_dst_p,
  NTDR_II_mdt = NTDR_II_mdt,
  PTDR_II_mdt = PTDR_II_mdt,
  NTDR_IN_nid = NTDR_IN_nid,
  PTDR_IN_nid = PTDR_IN_nid,
  NTDR_II_kap_mdt = NTDR_II_kap_mdt,
  NTDR_II_kap_main = NTDR_II_kap_main,
  pt_mdt_loss = pt_mdt_loss,
  nt_mdt_loss = nt_mdt_loss,
  pt_nid_loss = pt_nid_loss,
  nt_nid_loss = nt_nid_loss,
  ptii_tx_corr = ptii_tx_corr,
  ptin_tx_corr = ptin_tx_corr,
  ntii_tx_corr = ntii_tx_corr,
  ntin_tx_corr = ntin_tx_corr,
  sdr_start_yr_prev = sdr_start_yr_prev,
  sdr = sdr),
  hash = T
)

output_input                         <- list2env(list(
  year1 = year1,
  yearend = yearend,
  dt = dt,
  para_static = para_static,
  para_variable = para_variable),
  hash = T
)

output                               <- list(
  pop = output_pop,
  internals = output_internals,
  inputs = output_input,
  calc = calc
)
  } else if (mode == 1) {
    # cluster mode
fo                              <- list()
tmp                             <- list()

# Total tmp$populations of interest (large array slices)
tmp$y2000                       <- (la[, 401, ] + la[, 402, ] + la[, 403, ] + la[, 404, ])/4
tmp$y2010                       <- (la[, 441, ] + la[, 442, ] + la[, 443, ] + la[, 444, ])/4
tmp$y2015                       <- (la[, 461, ] + la[, 462, ] + la[, 463, ] + la[, 464, ])/4
tmp$y2016                       <- (la[, 465, ] + la[, 466, ] + la[, 467, ] + la[, 468, ])/4
tmp$y2013                       <- (la[, 453, ] + la[, 454, ] + la[, 455, ] + la[, 456, ])/4
tmp$y2007                       <- (la[, 429, ] + la[, 430, ] + la[, 431, ] + la[, 432, ])/4

tmp$pop_total2000               <- sum(tmp$y2000)
tmp$pop_total2010               <- sum(tmp$y2010)
tmp$pop_total2016               <- sum(tmp$y2016)
tmp$pop_total2000_allad         <- sum(tmp$y2000[, 16:100])
tmp$pop_total2010_allad         <- sum(tmp$y2010[, 16:100])
tmp$pop_total2015               <- sum(tmp$y2015)
tmp$pop_total2000_014           <- sum(tmp$y2000[, 1:15])
tmp$pop_total2010_014           <- sum(tmp$y2010[, 1:15])
tmp$pop_total2015_014           <- sum(tmp$y2015[, 1:15])
tmp$pop_total2000_1564          <- sum(tmp$y2000[, 16:65])
tmp$pop_total2010_1564          <- sum(tmp$y2010[, 16:65])
tmp$pop_total2015_1564          <- sum(tmp$y2015[, 16:65])
tmp$pop_total2000_65            <- sum(tmp$y2000[, 66:100])
tmp$pop_total2010_65            <- sum(tmp$y2010[, 66:100])
tmp$pop_total2015_65            <- sum(tmp$y2015[, 66:100])

tmp$prevcmp                     <- c("NTDS_II", "PTDS_II", "NTDR_II", "PTDR_II")

# Bacteriologically positive prevalent array
tmp$bprevarray2000              <- colMeans(colSums(la[tmp$prevcmp, 401:404, ]))
tmp$bprevarray2010              <- colMeans(colSums(la[tmp$prevcmp, 441:444, ]))


# All TB Prevalence 2000, normalised all ages, normalised
fo$alltb_previi_2000_norm_allad <- 1e+05 * sum(tmp$bprevarray2000[16:100])/tmp$pop_total2000_allad
# 2010, normalised all ages, normalised
fo$alltb_previi_2010_norm_allad <- 1e+05 * sum(tmp$bprevarray2010[16:100])/tmp$pop_total2010_allad
# 2000, normalised 15-29 ages, normalised
fo$alltb_previi_2000_norm_1529  <- 1e+05 * sum(tmp$bprevarray2000[16:30])/sum(tmp$y2000[, 16:30])
# 2010, normalised 15-29 ages, normalised
fo$alltb_previi_2010_norm_1529  <- 1e+05 * sum(tmp$bprevarray2010[16:30])/sum(tmp$y2010[, 16:30])
# 2000, normalised 30-44 ages, normalised
fo$alltb_previi_2000_norm_3044  <- 1e+05 * sum(tmp$bprevarray2000[31:45])/sum(tmp$y2000[, 31:45])
# 2010, normalised 30-44 ages, normalised
fo$alltb_previi_2010_norm_3044  <- 1e+05 * sum(tmp$bprevarray2010[31:45])/sum(tmp$y2010[, 31:45])
# 2000, normalised 45-59 ages, normalised
fo$alltb_previi_2000_norm_4559  <- 1e+05 * sum(tmp$bprevarray2000[46:60])/sum(tmp$y2000[, 46:60])
# 2010, normalised 45-59 ages, normalised
fo$alltb_previi_2010_norm_4559  <- 1e+05 * sum(tmp$bprevarray2010[46:60])/sum(tmp$y2010[, 46:60])
# 2000, normalised 60+ ages, normalised
fo$alltb_previi_2000_norm_60    <- 1e+05 * sum(tmp$bprevarray2000[61:100])/sum(tmp$y2000[, 61:100])
# 2010, normalised 60+ ages, normalised
fo$alltb_previi_2010_norm_60    <- 1e+05 * sum(tmp$bprevarray2010[61:100])/sum(tmp$y2010[, 61:100])

# All TB Incidence
# 2000, incidence normalised, all ages
fo$alltb_inc_2010_norm_all      <- 1e+05 * sum(ta[c("DSTB_inc", "DRTB_inc"), 441:444, ])/tmp$pop_total2010
# 2010, incidence normalised, all ages
fo$alltb_inc_2016_norm_all      <- 1e+05 * sum(ta[c("DSTB_inc", "DRTB_inc"), 465:468, ])/tmp$pop_total2016


# All TB mortality rate 2000 mortality normalised, all ages
fo$alltb_mort_2010_all          <- 1e+05 * sum(ta[c("DSTB_deaths", "DRTB_deaths"), 441:444, ])/tmp$pop_total2010
# 2000 mortality normalised, 0-14
fo$alltb_mort_2010_014          <- 1e+05 * sum(ta[c("DSTB_deaths", "DRTB_deaths"), 441:444, 1:15])/tmp$pop_total2010_014
# 2000 mortality normalised, 0-14
fo$alltb_mort_2010_1564         <- 1e+05 * sum(ta[c("DSTB_deaths", "DRTB_deaths"), 441:444, 16:65])/tmp$pop_total2010_1564
# 2000 mortality normalised, 0-14
fo$alltb_mort_2010_65           <- 1e+05 * sum(ta[c("DSTB_deaths", "DRTB_deaths"), 441:444, 66:100])/tmp$pop_total2010_65

# All TB notification rate 2000 notification normalised, all ages
fo$alltb_notif_2015_all         <- 1e+05 * sum(ta[c("DSTB_initRx", "DRTB_initRx"), 461:464, ])/tmp$pop_total2015
# 2000 notification normalised, 0-14
fo$alltb_notif_2015_014         <- 1e+05 * sum(ta[c("DSTB_initRx", "DRTB_initRx"), 461:464, 1:15])/tmp$pop_total2015_014
# 2000 notification normalised, 0-14
fo$alltb_notif_2015_1564        <- 1e+05 * sum(ta[c("DSTB_initRx", "DRTB_initRx"), 461:464, 16:65])/tmp$pop_total2015_1564
# 2000 notification normalised, 0-14
fo$alltb_notif_2015_65          <- 1e+05 * sum(ta[c("DSTB_initRx", "DRTB_initRx"), 461:464, 66:100])/tmp$pop_total2015_65

# # DSTB Notifications - 2017 (lif = lfu inflation factor)
# tmp$dx_2017                   <- sum(ta[c("DSTB_initRx", 481:484, )]) * lif
# # DSTB total treatment person-time - 2017
# tmp$ptot_2017                 <- sum(la[c("NTDS_T", "PTDS_T", 481:484, )]) * 12 * dt


# DRTB Calculations
# TOtal incidence 2016
fo$drtb_inc_2016                <- 1e+05 * sum(ta["DRTB_inc", 465:468, ])/tmp$pop_total2016

# Notification proportion - 2007
fo$drtb_nt_fnot_prop_2007       <- 100 * sum(ta["DR_nt_intx_f", 429:432, ])/sum(ta["All_nt_intx_f", 429:432, ])
fo$drtb_pt_fnot_prop_2007       <- 100 * sum(ta["DR_pt_intx_f", 429:432, ])/sum(ta["All_pt_intx_f", 429:432, ])

# Notification proportion - 2013
fo$drtb_nt_fnot_prop_2013       <- 100 * sum(ta["DR_nt_intx_f", 453:456, ])/sum(ta["All_nt_intx_f", 453:456, ])
fo$drtb_pt_fnot_prop_2013       <- 100 * sum(ta["DR_pt_intx_f", 453:456, ])/sum(ta["All_pt_intx_f", 453:456, ])

# Lab confirmed treatment Initiation
fo$drtb_not_lab_2013            <- 1e+03 * sum(ta["DRTB_initRxLab", 453:456, ])

# Health Economics - 2017 total cost (inc programme cost)
fo$total_tbrx_cost_2017         <- sum(ca["tbrx", 469:472, ])

rm(tmp)

# Step finder
# Given a particular calendar year, dt and starting year, isolate the steps for that calendar year
# sf                            <- function(year, year1 = year1, dt = dt) {
# stepspa                       <- 1/dt
# t1step                        <- (year - year1) * stepspa
# tastep                        <- c(rep(t1step, stepspa))
# tastep                        <- tastep + c(1:stepspa)
# return(tastep)
# }
  }

  if (mode == 0) {
    # Local mode - return detailed output
    return(output)
  } else if (mode == 1) {
    # Cluster mode - return only fitting values
    return(fo)
  }

}
