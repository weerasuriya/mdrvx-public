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

# Vaccine compartments
p_compartments <- c(p_compartments, paste0("v_", p_compartments))

# Initialise the "Large Array" (la) - the axes of this 3D array are: compartment name, timestep, age
la <- array(0, dim = c(length(p_compartments), steps, age_cls), dimnames = list(p_compartments, c(as.character(1:steps)), age_nms))

# Background mortality matrix
bg_mort <- array(0, dim = c(length(p_compartments), steps, age_cls), dimnames = list(p_compartments, c(as.character(1:steps)), age_nms))

# List of Transits (i.e. flows)
transit <- c(
  "DSTB_deaths",
  "DRTB_deaths",
  "DSTB_inc",
  "DRTB_inc",
  "DSTB_onRx",
  "DRTB_onRx",
  "bgmort",
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
  "total_cost",
  "v_DSTB_deaths",
  "v_DRTB_deaths",
  "v_DSTB_inc",
  "v_DRTB_inc",
  "v_DSTB_onRx",
  "v_DRTB_onRx",
  "v_bgmort",
  "v_TB_inc",
  "v_DSTB_initRx",
  "v_DRTB_initRx",
  "v_DS_nt_intx",
  "v_DS_nt_intx_f",
  "v_DS_pt_intx",
  "v_DS_pt_intx_f",
  "v_DR_nt_intx",
  "v_DR_nt_intx_f",
  "v_DR_pt_intx",
  "v_DR_pt_intx_f",
  "v_All_nt_intx_f",
  "v_All_pt_intx_f",
  "v_All_nt_intx",
  "v_All_pt_intx",
  "v_DRTB_initRxLab",
  "v_DS_nt_inc",
  "v_DS_nt_inc_f",
  "v_DS_pt_inc",
  "v_DS_pt_inc_f",
  "v_DR_nt_inc",
  "v_DR_nt_inc_f",
  "v_DR_pt_inc",
  "v_DR_pt_inc_f",
  "v_All_nt_inc_f",
  "v_All_pt_inc_f",
  "v_All_nt_inc",
  "v_All_pt_inc",
  "v_total_cost",
  "dsi_new",
  "dsi_rr",
  "dri_new",
  "dri_rr",
  "ti_new",
  "ti_rr",
  "ti_new_s",
  "ti_new_lr",
  "dsi_new_s",
  "dsi_new_lr",
  "dri_new_s",
  "dri_new_lr",
  "ti",
  "v_dsi_new",
  "v_dsi_rr",
  "v_dri_new",
  "v_dri_rr",
  "v_ti_new",
  "v_ti_rr",
  "v_ti_new_s",
  "v_ti_new_lr",
  "v_dsi_new_s",
  "v_dsi_new_lr",
  "v_dri_new_s",
  "v_dri_new_lr",
  "v_ti"
)

# Transit non/vaccine compartment list
tvc  <- grep(x = transit, pattern = "v_", value = T)
tnvc <- gsub(x = tvc, pattern = "v_", replacement = "")

# Initialise the 'Transit Array' (ta) - the axes of this 3D array are:
# compartment name (i.e. flow), timestep, age
ta   <- array(0, dim = c(length(transit), steps, age_cls), dimnames = list(transit, c(as.character(1:steps)), age_nms))

cost_types <- c(
  "ds_dx",
  "dr_dx",
  "dst",
  "ds_tx",
  "dr_tx",
  "tbrx",
  "nv_ds_dx",
  "nv_dr_dx",
  "nv_dst",
  "nv_ds_tx",
  "nv_dr_tx",
  "nv_tbrx",
  "v_ds_dx",
  "v_dr_dx",
  "v_dst",
  "v_ds_tx",
  "v_dr_tx",
  "v_tbrx",
  "v_immR",
  "v_immM",
  "v_ucvxR",
  "v_ucvxM",
  "v_prog",
  "grand_total",
  "nv_grand_total",
  "v_grand_total"
)

cac   <- c("ds_dx", "dr_dx", "dst", "ds_tx", "dr_tx", "tbrx")
cvc   <- paste0("v_", cac)
cnvc  <- paste0("nv_", cac)

# Initialise the 'Cost Array' (ta) - the axes of this 3D array are:
# compartment name (i.e. cost), timestep, age
# ca    <- array(0,
#               dim = c(
#                 length(cost_types),
#                 steps,
#                 age_cls),
#               dimnames = list(
#                 cost_types,
#                 c(as.character(1:steps)),
#                 age_nms))

ca <- matrix(0, nrow = steps, ncol = length(cost_types))
colnames(ca) <- cost_types
# Cost sub array - for vaccine parameters independent of dynamic effects
if (vaccine == 1) {
  vxa_c <- c("costR", "costM", "costT", "deliv_costR", "deliv_costM", "deliv_costT")
  vxa   <- array(0,
                 dim = c(
                   length(vxa_c),
                   steps,
                   length(vx_regimens)
                 ),
                 dimnames = list(vxa_c,
                                 as.character(1:steps),
                                 as.character(paste0("USD", vx_regimens))
                 )
  )
}

# Special arrays
# DSTB incidence array
dsia           <- matrix(0, steps, age_cls)
rownames(dsia) <- paste("Step", as.character(1:steps))
colnames(dsia) <- age_nms

# DRTB incidence array
dria           <- matrix(0, steps, age_cls)
rownames(dria) <- paste("Step", as.character(1:steps))
colnames(dria) <- age_nms

# AllTB incidence array
tbia           <- matrix(0, steps, age_cls)
rownames(tbia) <- paste("Step", as.character(1:steps))
colnames(tbia) <- age_nms

# DSTB TB mortality array
dsma           <- matrix(0, steps, age_cls)
rownames(dsma) <- paste("Step", as.character(1:steps))
colnames(dsma) <- age_nms

# Vaccines - number immunised per step
immunised           <- matrix(0, steps, age_cls)
rownames(immunised) <- paste("Step", as.character(1:steps))
colnames(immunised) <- age_nms

if (vaccine == 1) {
  v_dsia <- matrix(0, steps, age_cls)
  v_dria <- matrix(0, steps, age_cls)
  v_tbia <- matrix(0, steps, age_cls)
  v_dsma <- matrix(0, steps, age_cls)
}

# Initialise Step1 Compartment Populations ####
# Assign initial TB proportions
# For each year of age (column in total_pop), multiply by the compartment distribution (init_tb, age-invariant) and assign into la[, 1, ]
init_tb   <- init_tb_prop[, prop_col]
la_cmps   <- length(p_compartments)
la[, 1, ] <- outer(init_tb, total_population[1, ])

# Contact Matrix components ####
# psize for contact matirx

# Age classes in contact matrix
infage                    <- list(c(1:6), c(7:20), c(21:65), c(66:age_cls))

# Infectious compartments used in contact matrix
# DSTB Infectious compartments
dsicmp                    <- c("PTDS_II", "NTDS_II", "v_PTDS_II", "v_NTDS_II")
# DRTB Infectious compartments
dricmp                    <- c("PTDR_II", "NTDR_II", "PTDR_T_IIphi", "NTDR_T_IIphi", "NTDR_mdt", "PTDR_mdt",
  "v_PTDR_II", "v_NTDR_II", "v_PTDR_T_IIphi", "v_NTDR_T_IIphi", "v_NTDR_mdt", "v_PTDR_mdt")

# Set up population size (denominator) Matrix
psizematrix               <- matrix(0, steps, length(infage))

# Set up Infectious matrices for DS and DRTB
DS_Imatrix                <- matrix(0, steps, length(infage))
DR_Imatrix                <- matrix(0, steps, length(infage))

# Calculate number of infectious cases by age class and enter into Infectious matrices
psizematrix[1, ]          <- vapply(infage, function(x) sum(la[, 1, x]), FUN.VALUE = numeric(1))
DS_Imatrix[1, ]           <- vapply(infage, function(x) sum(la[dsicmp, 1, x]), FUN.VALUE = numeric(1))
DR_Imatrix[1, ]           <- vapply(infage, function(x) sum(la[dricmp, 1, x]), FUN.VALUE = numeric(1))

# Force Of Infection ####
## Drug Sensitive
NTDS_lambda               <- PTDS_lambda <- matrix(0, steps, age_cls)

## Drug Resistant
NTDR_lambda               <- PTDR_lambda <- matrix(0, steps, age_cls)

# Apply# Denominator matrix the DR_ts scaling factor
DR_neta                   <- DS_neta * DR_ts

## Initialise the first step force of infection
# Calculate transmission terms - Lambda
NTDS_lambda_raw           <- colSums(-(myneta[1:4, 1:age_cls]) * z * ((DS_Imatrix[1, 1:4])/(psizematrix[1, 1:4])))
NTDS_lambda[1, 1:age_cls] <- t(DS_neta * (1 - exp(NTDS_lambda_raw)))
PTDS_lambda               <- NTDS_lambda

NTDR_lambda_raw           <- colSums(-(myneta[1:4, 1:age_cls]) * z * ((DR_Imatrix[1, 1:4])/(psizematrix[1, 1:4])))
NTDR_lambda[1, 1:age_cls] <- t(DR_neta * (1 - exp(NTDR_lambda_raw)))
PTDR_lambda               <- NTDR_lambda

if (vaccine == 1) {
  v_NTDS_lambda <- NTDS_lambda * (1 - effI)
  v_PTDS_lambda <- v_NTDS_lambda * (1 - effI)
  v_NTDR_lambda <- NTDR_lambda * (1 - effI)
  v_PTDR_lambda <- v_NTDR_lambda * (1 - effI)
}

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

# Initialise DST and NID probability vectors
nt_dst_p <- pt_dst_p <- ntii_emp_p <- ptii_emp_p <- ntin_emp_p <- ptin_emp_p <- ntii_tx_corr <- ptii_tx_corr <- ntin_tx_corr <- ptin_tx_corr <- vector(mode = "numeric", length = (yearend - year1 + 1))

# Array to hold calculated SDR
sdr                     <- vector(mode = "numeric", length = steps)
sdr_start_yr_steps      <- ((sdr_start_year - year1) * (1/dt)) + seq(1:(1/dt))
sdr[sdr_start_yr_steps] <- sdr_base

# Vaccinated vs non-vaccinated p_compartments
nonvx_cmp <- c(
  "NTDS_II",
  "NTDS_IN",
  "NTDS_T",
  "NTDR_II",
  "NTDR_IN",
  "NTDR_T_IIphi",
  "NTDR_T_IIpsi",
  "NTDR_T_INphi",
  "NTDR_T_INpsi",
  "PTDS_II",
  "PTDS_IN",
  "PTDS_T",
  "PTDR_II",
  "PTDR_IN",
  "PTDR_T_IIphi",
  "PTDR_T_IIpsi",
  "PTDR_T_INphi",
  "PTDR_T_INpsi",
  "NTDR_mdt",
  "NTDR_nid",
  "PTDR_mdt",
  "PTDR_nid"
)

# Garbage collection
rm(sdr_start_yr_steps)

# Error / warning conditions and signal handling
e_ntu       <- simpleError(message = "kappa>1", call = list(rtype = 1))
e_la        <- simpleError(message = "negative values in la[]", call = list(rtype = 3, step = 1))
e_nan       <- simpleError(message = "NA/NaN in results", call = list(rtype = 4))
e_psz       <- simpleError(message = "psizematrix <= 0", call=list(rtype = 5, extra = psizematrix))

# Vector to hold ini_proportion
ini_prop <- matrix(0, nrow = steps, ncol = age_cls)
prev_denom <- c("NTDS_II", "v_NTDS_II", "NTDS_IN", "v_NTDS_IN", "PTDS_II", "v_PTDS_II", "PTDS_IN", "v_PTDS_IN")
prev_numer <- c("NTDS_IN", "v_NTDS_IN", "PTDS_IN", "v_PTDS_IN")
