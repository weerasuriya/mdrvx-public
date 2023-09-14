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

#Ribbon mode
if (mode == 2) {
  library(readr)
  library(digest)
  para_hash <- digest(para_variable)
  write_csv(cbind(para_hash, calc$ann_alltb_notif_norm), append = T, "../output/trajectories/alltb_notifications.csv")
  write_csv(cbind(para_hash, calc$ann_alltb_prev_norm_ii), append = T, "../output/trajectories/alltb_prevalence.csv")
  write_csv(cbind(para_hash, calc$ann_alltb_mort_norm), append = T, "../output/trajectories/alltb_mortality.csv")
  write_csv(cbind(para_hash, calc$ann_alltb_inc_norm), append = T, "../output/trajectories/alltb_incidence.csv")
  write_csv(cbind(para_hash, calc$ann_drtb_inc_norm), append = T, "../output/trajectories/drtb_incidence.csv")
  write_csv(cbind(para_hash, calc$mdr_nt_nprop_f), append = T, "../output/trajectories/drtb_nt_prop.csv")
  write_csv(cbind(para_hash, calc$mdr_pt_nprop_f), append = T, "../output/trajectories/drtb_pt_prop.csv")
  write_csv(cbind(para_hash, calc$ann_drtb_notif_lab), append = T, "../output/trajectories/drtb_lab_init.csv")
  write_csv(cbind(para_hash, calc$ann_cost_array), append = T, "../output/trajectories/cost_array.csv")

  return(1)
}

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
