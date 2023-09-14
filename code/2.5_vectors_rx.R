# Interpolation functions for pnpr/ppm
npr_l <- list(
    x = c(2012, 2017),
    y = c(0, 0))

npr<- approxfun(npr_l, rule=2)

ppm_l <- list(
    x = c(2018, 2025),
    y = c(ppm, if (scen == 1) ppm else if (scen == 2) alt_ppm))

ppm <- approxfun(ppm_l, rule=2)

# GLF/CDR interpolation function
glf_f <- function(A = 0, K = 1, C = 1, B = 1, Q = 1, v = 1, M = 2003, tstart = 1970, tend = 2050) {
  climb = A + ((K - A)/((C + (Q * exp(-B * ((tstart:tend) - M)))^(1/v))))
  approxfun(x = c(tstart:tend), y = climb, rule = 2)
}

cdr <- glf_f(tstart = 1970, tend = 2050, A = glf_a, B = glf_b, K = glf_k, v = glf_v, Q = glf_q, C = glf_c, M = glf_m)

# Master detection and treatment table
DXTX <- data.table(
  Year = year1:yearend
)

DXTX[, c("ppm", "npr", "npu", "cdr") := .(ppm(Year), npr(Year), 1-npr(Year), cdr(Year))]
DXTX[, "tir":= (cdr(Year) * (1 - npr(Year)) / (..pnpu * (1 - ppm(Year))))]
DXTX[, 'not_fr':= (1-ppm(Year)) * (tir * ..pnpu)/npu ]
DXTX[, 'pnpr':= not_fr * npr(Year) / ppm(Year)]

# the not_fr works out to be the same as the CDR
# back calculate pnpr from this term

# Case Detection Rate
# Create age and year-wise CDR matrices
DS_CDR <- matrix(0, (yearend - year1 + 1), age_cls)

# Apply adjusted TIR to DS_CDR matrix
DS_CDR[] <- DXTX[, tir]

# Age-wise scaling of CDR DS_CDRScaled, DS_CDRScaledE, DR_CDRScaled,
# Older adult CDRScaling
DS_CDRscaleO <- (DS_CDRscale + DS_CDRscaleE)/2

# ## DSTB
DS_CDR[, 1:(chiyrs + yaduyrs)]                     <- CDR_scaler(scl = DS_CDRscale, cdr = DS_CDR[, 1:(chiyrs + yaduyrs)])
DS_CDR[, (chiyrs + yaduyrs + 1):(chiyrs + aduyrs)] <- CDR_scaler(scl = DS_CDRscaleO, cdr = DS_CDR[, (chiyrs + yaduyrs + 1):(chiyrs + aduyrs)])
DS_CDR[, (chiyrs + aduyrs + 1):(age_cls)]          <- CDR_scaler(scl = DS_CDRscaleE, cdr = DS_CDR[, (chiyrs + aduyrs + 1):(age_cls)])

# Scale down to per time step
DS_CDR <- 1 - (1 - DS_CDR)^(dt)

# Add a column of years
DS_CDR <- cbind(year1:yearend, DS_CDR)

# Case detection switch
if (rx == 0) {
  DS_CDR[, 2:101] <- 0
}

# Create matrix to hold new treatment initiation rates
NTDS_II_kap_main <- PTDS_II_kap_main <- NTDS_IN_kap_main <- PTDS_IN_kap_main <- NTDR_II_kap_main <- PTDR_II_kap_main <- NTDR_IN_kap_main <- PTDR_IN_kap_main <- matrix(0, yearend - year1 + 1, age_cls + 1)
NTDR_II_kap_mdt  <- PTDR_II_kap_mdt  <- PTDR_IN_kap_nid <- NTDR_IN_kap_nid   <- matrix(0, yearend - year1 + 1, age_cls + 1)

NTDS_II_kap_main[, 1] <- PTDS_II_kap_main[, 1]<- NTDS_IN_kap_main[, 1]  <- PTDS_IN_kap_main[, 1] <- NTDR_II_kap_main[, 1] <- PTDR_II_kap_main[, 1] <- NTDR_IN_kap_main[, 1] <- PTDR_IN_kap_main[, 1] <- year1:yearend
NTDR_II_kap_mdt[, 1]  <- PTDR_II_kap_mdt[, 1]  <- PTDR_IN_kap_nid[, 1] <- NTDR_IN_kap_nid[, 1]   <- year1:yearend

######################### Treatment Probabilities for DRTB ##################################
# Create linear scale up of empirical treatment probability until DST start year1

dx_p <- data.table(
  Year = year1:yearend,
  nt_dst_p = 0,
  ntii_emp_p = 0,
  ntin_emp_p = 0,
  ntii_tx_corr = 0,
  ntin_tx_corr = 0,
  pt_dst_p = 0,
  ptii_emp_p = 0,
  ptin_emp_p = 0,
  ptii_tx_corr = 0,
  ptin_tx_corr = 0
)
setkey(dx_p, Year)

# Linear interpolation functions for scale up of DST
nt_dst <- approxfun(
  x = c(dst_start_year, f_dst_year, dstp_scaleup$year),
  y = c(0, nt_dst_prob, {if (scen == 1) dstp_scaleup$sq else if (scen == 2) dstp_scaleup$pol}),
  rule = 2
)

pt_dst <- approxfun(
  x = c(dst_start_year, f_dst_year, dstp_scaleup$year),
  y = c(0, pt_dst_prob, {if (scen == 1) dstp_scaleup$sq else if (scen == 2) dstp_scaleup$pol}),
  rule = 2
)

dx_p[,c("nt_dst_p", "pt_dst_p") := .(nt_dst(Year), pt_dst(Year))]

# Initialise DRTB Empirical treatment probability
# Empirical treatment probability = 0 until 1 year before when MDR acquisition begins
wrap_emp <- function(emp_tx) {
  approxfun(x = c(mdrac_yr, last_fit_year), y = c(0, emp_tx), rule = 2)
}

ntin_emp <- wrap_emp(ntin_emp_tx_p)
ptin_emp <- wrap_emp(ptin_emp_tx_p)
ntii_emp <- wrap_emp(ntii_emp_tx_p)
ptii_emp <- wrap_emp(ptii_emp_tx_p)

dx_p[, c("ntin_emp_p", "ptin_emp_p", "ntii_emp_p", "ptii_emp_p") := .(ntin_emp(Year), ptin_emp(Year),((1 - nt_dst(Year)) * ntii_emp(Year)), ((1 - pt_dst(Year)) * ptii_emp(Year)) )]

dx_p[, pub_ntii_tx_corr := ntii_emp_p + nt_dst_p]
dx_p[, pub_ntin_tx_corr := ntin_emp_p]
dx_p[, pub_ptii_tx_corr := ptii_emp_p + pt_dst_p]
dx_p[, pub_ptin_tx_corr := ptin_emp_p]

dx_p[, pri_ntii_tx_corr := 0]
dx_p[, pri_ntin_tx_corr := 0]
dx_p[, pri_ptii_tx_corr := 0]
dx_p[, pri_ptin_tx_corr := 0]

# Weighted average correct treatment proportion between private and public sector
dx_p[, ntii_tx_corr := (ppm(Year) * pri_ntii_tx_corr) + ((1 - ppm(Year)) * pub_ntii_tx_corr)]
dx_p[, ptii_tx_corr := (ppm(Year) * pri_ptii_tx_corr) + ((1 - ppm(Year)) * pub_ptii_tx_corr)]
dx_p[, ptin_tx_corr := (ppm(Year) * pri_ptin_tx_corr) + ((1 - ppm(Year)) * pub_ptin_tx_corr)]
dx_p[, ntin_tx_corr := (ppm(Year) * pri_ntin_tx_corr) + ((1 - ppm(Year)) * pub_ntin_tx_corr)]

# Proportion of private sector MDR treatment initiation _among correct treatment initiations_
dx_p[, ntii_prtp := (ppm(Year) * pri_ntii_tx_corr)/ntii_tx_corr]
dx_p[, ptii_prtp := (ppm(Year) * pri_ptii_tx_corr)/ptii_tx_corr]
dx_p[, ptin_prtp := (ppm(Year) * pri_ptin_tx_corr)/ptin_tx_corr]
dx_p[, ntin_prtp := (ppm(Year) * pri_ntin_tx_corr)/ntin_tx_corr]

############################################################################################
# Calculate treatment initiation rates
NTDS_II_kap_main[, -1] <- sweep(x = DS_CDR[, -1], MARGIN = 2, FUN = "*", STATS = (NTDS_II_u + NTDS_II_n))/(1 - DS_CDR[, -1])
PTDS_II_kap_main[, -1] <- sweep(x = DS_CDR[, -1], MARGIN = 2, FUN = "*", STATS = (PTDS_II_u + PTDS_II_n))/(1 - DS_CDR[, -1])
NTDS_IN_kap_main[, -1] <- sweep(x = DS_CDR[, -1] * e, MARGIN = 2, FUN = "*", STATS = (NTDS_IN_u + NTDS_IN_n))/(1 - (DS_CDR[, -1] * e))
PTDS_IN_kap_main[, -1] <- sweep(x = DS_CDR[, -1] * e, MARGIN = 2, FUN = "*", STATS = (PTDS_IN_u + PTDS_IN_n))/(1 - (DS_CDR[, -1] * e))

NTDR_II_kap_main[NTDR_II_kap_main[, 1] >= mdrac_yr, -1] <- NTDS_II_kap_main[NTDS_II_kap_main[, 1]>= mdrac_yr, -1] * dx_p[Year>=mdrac_yr, ntii_tx_corr]
NTDR_II_kap_mdt[NTDR_II_kap_mdt[, 1] >= mdrac_yr, -1 ]  <- NTDS_II_kap_main[NTDS_II_kap_main[, 1]>= mdrac_yr, -1] * (1 - dx_p[Year>=mdrac_yr, ntii_tx_corr])

PTDR_II_kap_main[PTDR_II_kap_main[, 1] >= mdrac_yr, -1] <- PTDS_II_kap_main[PTDS_II_kap_main[, 1]>= mdrac_yr, -1] * dx_p[Year>=mdrac_yr, ptii_tx_corr]
PTDR_II_kap_mdt[PTDR_II_kap_mdt[, 1] >= mdrac_yr, -1 ]  <- PTDS_II_kap_main[PTDS_II_kap_main[, 1]>= mdrac_yr, -1] * (1 - dx_p[Year>=mdrac_yr, ptii_tx_corr])

NTDR_IN_kap_main[NTDR_IN_kap_main[, 1] >= mdrac_yr, -1] <- NTDS_IN_kap_main[NTDS_IN_kap_main[, 1]>= mdrac_yr, -1] * dx_p[Year>=mdrac_yr, ntin_tx_corr]
NTDR_IN_kap_nid[NTDR_IN_kap_nid[, 1] >= mdrac_yr, -1 ]  <- NTDS_IN_kap_main[NTDS_IN_kap_main[, 1]>= mdrac_yr, -1] * (1 - dx_p[Year>=mdrac_yr, ntin_tx_corr])

PTDR_IN_kap_main[PTDR_IN_kap_main[, 1] >= mdrac_yr, -1] <- PTDS_IN_kap_main[PTDS_IN_kap_main[, 1]>= mdrac_yr, -1] * dx_p[Year>=mdrac_yr, ptin_tx_corr]
PTDR_IN_kap_nid[PTDR_IN_kap_nid[, 1] >= mdrac_yr, -1 ]  <- PTDS_IN_kap_main[PTDS_IN_kap_main[, 1]>= mdrac_yr, -1] * (1 - dx_p[Year>=mdrac_yr, ptin_tx_corr])

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
  stop(e_ntu)
}

# Weight treatment success according to PPM and generate interpolation function
psi_pp[, ds_psi := (ds_pri_psi * ppm(Year)) + (ds_pub_psi * (1-ppm(Year)))]
ds_psi <- approxfun(psi_pp[, .(Year, ds_psi)], rule=2)
# For DRTB the proportion of private sector MDR treatment must be taken from the dx_p table
psi_pp[, dr_psi := (dr_pri_psi * dx_p[Year %in% psi_pp[, Year], .(ptii_prtp)]) + (dr_pub_psi * (1-dx_p[Year %in% psi_pp[, Year], .(ptii_prtp)]))]
dr_psi <- approxfun(psi_pp[, .(Year, dr_psi)], rule=2)

# Interpolation functions for yearwise-lookup
# dst_probabilities
if_pt_dst_p <- approxfun(dx_p[, .(x = Year, y = pt_dst_p)], rule = 2)
if_nt_dst_p <- approxfun(dx_p[, .(x = Year, y = nt_dst_p)], rule = 2)

# correct treatment probabilities
if_ptii_tx_corr <- approxfun(dx_p[, .(x = Year, y = ptii_tx_corr)], rule = 2)
if_ntii_tx_corr <- approxfun(dx_p[, .(x = Year, y = ntii_tx_corr)], rule = 2)
if_ntin_tx_corr <- approxfun(dx_p[, .(x = Year, y = ntin_tx_corr)], rule = 2)
if_ptin_tx_corr <- approxfun(dx_p[, .(x = Year, y = ptin_tx_corr)], rule = 2)

# Cost scaling
# China - R1/R2 weighting for MDR-TB treatment
mdrtb_rx[, R1_cost := R1 / 100 * mdr_r1_cost]
mdrtb_rx[, R2_cost := R2 / 100 * mdr_r2_cost]
mdrtb_rx[, SC_cost := SC / 100 * mdr_sc_cost]
mdrtb_rx[, MDR_cost_pol := (R1_cost + R2_cost + SC_cost) * 12 * ..dt]
mdrtb_rx[, MDR_cost_sq := mdr_r2_cost * 12 * ..dt]

arc_mdrtb_rx <- copy(mdrtb_rx)

ds_tx_cost <- mds_tx_cost * 12 * dt
