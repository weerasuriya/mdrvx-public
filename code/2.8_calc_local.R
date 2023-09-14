# Local calculations

# 1. Housekeeping
# 1.1 New holding environment for results and intermediate results
calc              <- new.env(hash = TRUE)
ir                <- new.env(hash = TRUE)
ir$premelt  <- list()

# 1.2 Working time horizons
ir$Years          <- rep(year1:yearend, each = 1/dt)
if (vaccine == 1 || vaccine == 2) {
  if (exists("extra_vx_start")) {
    ir$export_start <- extra_vx_start
  } else {
      ir$export_start <- 2018
      }
  ir$export_end   <- 2050
} else {
  ir$export_start <- 1990
  ir$export_end   <- 2050
}
ir$wth            <- ir$export_start:ir$export_end
# 1.3 Indices Indices: Susceptible, Recovered and Latent [i.e. 'Non-Prevalent']
# Compartments
ir$nonprev        <- match(c("S", "NTDS_R", "PTDS_R", "NTDR_R", "PTDR_R", "NTDR_L", "NTDS_L", "PTDR_L"), unlist(dimnames(la)[1]))
# Indices: DR compartments
ir$dr_comps       <- grep("DR", unlist(dimnames(la)[1]))
# Indices: DS compartments
ir$ds_comps       <- grep("DS", unlist(dimnames(la)[1]))

if (!is.null(modespec)) {
  aglp <- modespec$aglp
  aglo <- modespec$aglo
} else {
# 1.4 Age groups
# Prevalence tables
aglp <- list(All = c(0:max_age), A15 = c(15:max_age), A014 = c(0:14), A1529 = c(15:29), A3044 = c(30:44), A4559 = c(45:59), A60 = c(60:max_age))
# Other
  aglo <- list(All = c(0:max_age), A014 = c(0:14), A15 = c(15:max_age))
}


# 1.5 Data.table attributes
ir$reshape <- list(
  id.vars = c("Year", "para_hash", "scen"),
  melt.grp = c(names(aglo), "All")
)

# 2. Demographics - annual age-wise population
ir$ann_agewise     <- data.table(Year = ir$Years, colSums(la, dims = 1))[Year >= ir$export_start & Year <= ir$export_end,lapply(.SD,mean) , keyby = Year]

ir$premelt$ann_agegrps        <- data.table(Year = ir$wth,
                              A014 = ir$ann_agewise[, rowSums(.SD), .SDcol = as.character(0:14)],
                              A1564 = ir$ann_agewise[, rowSums(.SD), .SDcol = as.character(15:64)],
                              A65 = ir$ann_agewise[, rowSums(.SD), .SDcol = as.character(64:max_age)],
                              All = ir$ann_agewise[, rowSums(.SD), .SDcol = as.character(0:max_age)],
                              AllAdults = ir$ann_agewise[, rowSums(.SD), .SDcol = as.character(15:max_age)],
                              A15 = ir$ann_agewise[, rowSums(.SD), .SDcol = as.character(15:max_age)],
                              A1529 = ir$ann_agewise[, rowSums(.SD), .SDcol = as.character(15:29)],
                              A3044 = ir$ann_agewise[, rowSums(.SD), .SDcol = as.character(30:44)],
                              A4559 = ir$ann_agewise[, rowSums(.SD), .SDcol = as.character(45:59)],
                              A60 = ir$ann_agewise[, rowSums(.SD), .SDcol = as.character(60:max_age)])[, c("para_hash", "scen") := .(para_hash, scen)]
setkey(ir$premelt$ann_agegrps, Year)

# 3. DSTB
# 3.1 DSTB Prevalence
# Annual average DSTB prevalence, Infectious TB only, raw
ir$ann_dstb_prev_raw <- data.table(Year = ir$Years, Prevalence = colSums(la[c("NTDS_II", "PTDS_II", "v_NTDS_II", "v_PTDS_II"), , ], dims = 1))[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, mean), keyby = Year]
# Annual average DSTB prevalence by age **group**, Infectious TB only, population normalised
ir$premelt$ann_dstb_prev     <- setDT(lapply(X = aglp, function(x) psadj(ir$ann_dstb_prev_raw, x)))[, Year := ir$wth][, c("para_hash", "scen") := .(para_hash, scen)]
setcolorder(ir$premelt$ann_dstb_prev, c("Year"))

# 3.2 DSTB Incidence
# Annual DSTB incidence by age, raw
ir$ann_dstb_inc_raw   <- data.table(Year = ir$Years, dsia)[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, sum), keyby = Year]
# Annual DSTB incidence by age **group**, population normalised
ir$premelt$ann_dstb_inc       <- setDT(lapply(X = aglo, function(x) psadj(ir$ann_dstb_inc_raw, x)))[, Year := ir$wth][, c("para_hash", "scen") := .(para_hash, scen)]
setcolorder(ir$premelt$ann_dstb_inc, c("Year"))

# 3.3 DSTB Mortality
# Annual DSTB mortality by age, raw
ir$ann_dstb_mort_raw  <- data.table(Year = ir$Years, dsma)[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, sum), keyby = Year]
# Annual DSTB mortality by age **group**, population normalised
ir$premelt$ann_dstb_mort      <- setDT(lapply(X = aglo, function(x) psadj(ir$ann_dstb_mort_raw, x)))[, Year := ir$wth][, c("para_hash", "scen") := .(para_hash, scen)]
setcolorder(ir$premelt$ann_dstb_mort, c("Year"))

# 3.4 DSTB Notifications
# Annual
ir$ann_dstb_notif_raw  <- data.table(Year = ir$Years, ta["DSTB_initRx", , ])[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, sum), keyby = Year]

# By age
ir$premelt$ann_dstb_notif     <- setDT(lapply(X = aglo, function(x) psadj(ir$ann_dstb_notif_raw, x)))[, Year := ir$wth][, c("para_hash", "scen") := .(para_hash, scen)]
setcolorder(ir$premelt$ann_dstb_notif, c("Year"))

# 4. DRTB
# 4.1 DRTB Incidence
# Annual DRTB incidence by age, raw
ir$ann_drtb_inc_raw   <- data.table(Year = ir$Years, dria)[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, sum), keyby = Year]
# Annual DRTB incidence by age **group**, population-normalised
ir$premelt$ann_drtb_inc       <- setDT(lapply(X = aglo, function(x) psadj(ir$ann_drtb_inc_raw, x)))[, Year := ir$wth][, c("para_hash", "scen") := .(para_hash, scen)]
setcolorder(ir$premelt$ann_drtb_inc, c("Year"))

# 4.2 DRTB Mortality
ir$ann_drtb_mort_raw  <- data.table(Year = ir$Years, ta["DRTB_deaths", , ])[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, sum), keyby = Year]
# Annual drtb mortality by age **group**, population normalised
ir$premelt$ann_drtb_mort      <- setDT(lapply(X = aglo, function(x) psadj(ir$ann_drtb_mort_raw, x)))[, Year := ir$wth][, c("para_hash", "scen") := .(para_hash, scen)]
setcolorder(ir$premelt$ann_drtb_mort, c("Year"))
setnames(ir$ann_drtb_mort_raw, c("Year", as.character(0:99)))

# 4.3 DRTB Notifications
# DRTB Treatment Initiations - Including Misdiagnoses converted to correct treatment
ir$premelt$ann_drtb_notif     <- data.table(Year = ir$Years, 1e03 * rowSums(ta["DRTB_initRx", , ]))[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, sum), keyby = Year][, c("para_hash", "scen") := .(para_hash, scen)]

# Annual
ir$ann_drtb_anotif_raw  <- data.table(Year = ir$Years, ta["DRTB_initRx", , ])[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, sum), keyby = Year]
# By age
ir$premelt$ann_drtb_anotif     <- setDT(lapply(X = aglo, function(x) psadj(ir$ann_drtb_anotif_raw, x)))[, Year := ir$wth][, c("para_hash", "scen") := .(para_hash, scen)]
setcolorder(ir$premelt$ann_drtb_anotif, c("Year"))

# DRTB Treatment Initiations - Scaled to CCDC Only
setnames(ir$premelt$ann_drtb_notif, old = "V2", new = "Total")
ir$premelt$ann_drtb_notif[, CCDC := Total * ..chr]
ir$premelt$ann_drtb_notif[, Hospital := Total - CCDC]

# DRTB Treatment Initiations - Only transitions from Prevalent Compartments to Treatment Compartments
ir$premelt$ann_drtb_notif_lab <- data.table(Year = ir$Years, 1e03 * rowSums(ta["DRTB_initRxLab", , ]))[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, sum), keyby = Year][, c("para_hash", "scen") := .(para_hash, scen)]

# 4.4 DRTB Proportional Incidence
tot_nt_inc_f            <- data.table(Year = ir$Years, NT = rowSums(ta["All_nt_inc_f", , ]))[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, sum), keyby = Year]
tot_pt_inc_f            <- data.table(Year = ir$Years, PT = rowSums(ta["All_pt_inc_f", , ]))[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, sum), keyby = Year]
mdr_nt_inc_f            <- data.table(Year = ir$Years, NT = rowSums(ta["DR_nt_inc_f", , ]))[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, sum), keyby = Year]
mdr_pt_inc_f            <- data.table(Year = ir$Years, PT = rowSums(ta["DR_pt_inc_f", , ]))[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, sum), keyby = Year]

ir$premelt$mdr_iprop_f        <- data.table(Year = ir$wth, (mdr_nt_inc_f[,2] / tot_nt_inc_f[,2]), (mdr_pt_inc_f[,2] / tot_pt_inc_f[,2]))[, c("para_hash", "scen") := .(para_hash, scen)]
ir$premelt$mdr_iprop_f[is.na(NT), NT := 0]
ir$premelt$mdr_iprop_f[is.na(PT), PT := 0]

# 4.5 DRTB Proportional Notifications
tot_nt_intx_f           <- data.table(Year = ir$Years, NT = rowSums(ta["All_nt_intx_f", , ]))[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, sum), keyby = Year]
tot_pt_intx_f           <- data.table(Year = ir$Years, PT = rowSums(ta["All_pt_intx_f", , ]))[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, sum), keyby = Year]
mdr_nt_intx_f           <- data.table(Year = ir$Years, NT = rowSums(ta["DR_nt_intx_f", , ]))[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, sum), keyby = Year]
mdr_pt_intx_f           <- data.table(Year = ir$Years, PT = rowSums(ta["DR_pt_intx_f", , ]))[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, sum), keyby = Year]

ir$premelt$mdr_nprop_f        <- data.table(Year = ir$wth, (mdr_nt_intx_f[,2] / tot_nt_intx_f[,2]), (mdr_pt_intx_f[,2] / tot_pt_intx_f[,2]))[, c("para_hash", "scen") := .(para_hash, scen)]
ir$premelt$mdr_nprop_f[is.na(NT), NT := 0]
ir$premelt$mdr_nprop_f[is.na(PT), PT := 0]

ir$premelt$disintx <- cbind(tot_nt_intx_f, tot_pt_intx_f[, 2], M = mdr_nt_intx_f[, 2], M = mdr_pt_intx_f[, 2])[, c("para_hash", "scen") := .(para_hash, scen)]

# 4.6 DRTB Prevalence
# Annual average DSTB prevalence, Infectious TB only, raw
ir$ann_drtb_prev_raw      <- data.table(Year = ir$Years, Prevalence = colSums(la[c("NTDR_II", "PTDR_II", "v_NTDR_II", "v_PTDR_II"), , ], dims = 1))[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, mean), keyby = Year]
# Annual average DSTB prevalence by age **group**, Infectious TB only, population normalised
ir$premelt$ann_drtb_prev        <- setDT(lapply(X = aglp, function(x) psadj(ir$ann_drtb_prev_raw, x)))[, Year := ir$wth][, c("para_hash", "scen") := .(para_hash, scen)]
setcolorder(ir$premelt$ann_drtb_prev, c("Year"))

# 5. Health Economics - Costing
# 5.1 Time on Treatment
# Person-time on treatment, by year in 3 month units and 1 month units
ir$ann_onrx           <- data.table(Year = ir$Years,
                            DSTBonRx = rowSums(ta["DSTB_onRx", , ]),
                            DRTBonRx = rowSums(ta["DRTB_onRx", , ]))[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, sum), keyby = Year]

ir$ann_onrx[, DSTBonRx_mo := (12 * ..dt) * DSTBonRx]
ir$ann_onrx[, DRTBonRx_mo := (12 * ..dt) * DRTBonRx]

# Treatment initiations per year
ir$premelt$ann_initRx       <- data.table(Year = ir$Years,
                              DSTB_initRx = rowSums(ta["DSTB_initRx", , ]),
                              DRTB_initRx = rowSums(ta["DRTB_initRx", , ]))[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, sum), keyby = Year][, c("para_hash", "scen") := .(para_hash, scen)]

# 5.3 Total Cost
# ir$premelt$ann_total_cost <- annual_sum(rowSums(ta["total_cost", , ]), dt, year1, yearend)[-c(1:50), ]
ir$cost_array         <- data.table(Year = ir$Years, ca[, c(cac, cnvc, cvc)])
ir$cost_array[, ccdc_tbrx := ((ds_dx + ds_tx + ((dr_dx + dr_tx + dst) * ..chr))) * (1 + prog_cost)]

ir$premelt$ann_cost_array   <- ir$cost_array[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, sum), keyby = Year][, c("para_hash", "scen") := .(para_hash, scen)]

# 5.4 DALY weight compartments
# Full "any active TB disease mode"
# dalycmps            <- c("NTDS_II", "NTDS_IN", "NTDS_T", "NTDR_II", "NTDR_IN", "NTDR_T_IIphi", "NTDR_T_IIpsi",
# "NTDR_T_INphi", "NTDR_T_INpsi", "PTDS_II", "PTDS_IN", "PTDS_T", "PTDR_II", "PTDR_IN",
# "PTDR_T_IIphi", "PTDR_T_IIpsi", "PTDR_T_INphi", "PTDR_T_INpsi", "NTDR_mdt", "NTDR_nid",
# "PTDR_mdt", "PTDR_nid", "v_NTDS_II", "v_NTDS_IN", "v_NTDS_T", "v_NTDR_II", "v_NTDR_IN",
# "v_NTDR_T_IIphi", "v_NTDR_T_IIpsi", "v_NTDR_T_INphi", "v_NTDR_T_INpsi", "v_PTDS_II",
# "v_PTDS_IN", "v_PTDS_T", "v_PTDR_II", "v_PTDR_IN", "v_PTDR_T_IIphi", "v_PTDR_T_IIpsi",
# "v_PTDR_T_INphi", "v_PTDR_T_INpsi", "v_NTDR_mdt", "v_NTDR_nid", "v_PTDR_mdt",
# "v_PTDR_nid")
# Untreated active TB only mode
dalycmps                       <- c("NTDS_II", "NTDS_IN", "NTDR_II", "NTDR_IN", "PTDS_II", "PTDS_IN", "PTDR_II",
                              "PTDR_IN", "v_NTDS_II", "v_NTDS_IN", "v_NTDR_II", "v_NTDR_IN", "v_PTDS_II", "v_PTDS_IN",
                              "v_PTDR_II", "v_PTDR_IN")

drtb_dalycmps                  <- c("NTDR_II", "NTDR_IN", "PTDR_II", "PTDR_IN", "v_NTDR_II", "v_NTDR_IN",
                              "v_PTDR_II", "v_PTDR_IN")

dstb_dalycmps                  <- c("NTDS_II", "NTDS_IN", "PTDS_II", "PTDS_IN", "v_NTDS_II", "v_NTDS_IN",
                              "v_PTDS_II", "v_PTDS_IN")

ir$premelt$ann_dalycmp               <- data.table(Year = ir$Years, colSums(la[dalycmps, , ], dims = 1) * 1e+03)[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, mean), keyby = Year][, c("para_hash", "scen") := .(para_hash, scen)]

ir$premelt$ann_dstb_dalycmp          <- data.table(Year = ir$Years, colSums(la[dstb_dalycmps, , ], dims = 1) * 1e+03)[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, mean), keyby = Year][, c("para_hash", "scen") := .(para_hash, scen)]

ir$premelt$ann_drtb_dalycmp          <- data.table(Year = ir$Years, colSums(la[drtb_dalycmps, , ], dims = 1) * 1e+03)[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, mean), keyby = Year][, c("para_hash", "scen") := .(para_hash, scen)]

# NB DALY table is **already multiplied up by 1000**
comment(ir$premelt$ann_dalycmp)      <- "Values are multiplied by 1000 within mainloop itself"
comment(ir$premelt$ann_dstb_dalycmp) <- "Values are multiplied by 1000 within mainloop itself"
comment(ir$premelt$ann_drtb_dalycmp) <- "Values are multiplied by 1000 within mainloop itself"

# 6 ALL TB
# 6.1 ALL TB Prevalence
# Annual all TB prevalence, by age, infectious TB only
ir$ann_alltb_prev_agewise <- data.table(Year = ir$Years, colSums(la[c("NTDS_II", "PTDS_II", "NTDR_II", "PTDR_II", "v_NTDS_II", "v_PTDS_II", "v_NTDR_II", "v_PTDR_II"), , ], dims = 1))[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, mean), keyby = Year]
# Annual all TB prevalence by age **group**, Infectious TB only, population adjusted
ir$premelt$ann_alltb_prev       <- setDT(lapply(X = aglp, function(x) psadj(ir$ann_alltb_prev_agewise, x)))[, Year := ir$wth][, c("para_hash", "scen") := .(para_hash, scen)]
setcolorder(ir$premelt$ann_alltb_prev, c("Year"))

# 6.2 ALL TB Incidence, by age groups
# Annual All TB incidence, by age groups, population adjusted
ir$premelt$ann_alltb_inc     <- data.table(Year = ir$wth, para_hash = para_hash, scen = scen, ir$premelt$ann_dstb_inc[, -c("Year", "para_hash", "scen")] + ir$premelt$ann_drtb_inc[, -c("Year", "para_hash", "scen")])

# 6.3 ALL TB Mortality Annual All TB mortality, by age groups, population adjusted
ir$premelt$ann_alltb_mort    <- data.table(Year = ir$wth, para_hash = para_hash, scen = scen, ir$premelt$ann_dstb_mort[, -c("Year", "para_hash", "scen")] + ir$premelt$ann_drtb_mort[, -c("Year", "para_hash", "scen")])

# 6.4 ALL TB Notifications
# Raw allTB notifications
ir$ann_alltb_notif_raw <- data.table(Year = ir$Years, (ta["DSTB_initRx", , ] + ta["DRTB_initRx", , ]))[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, sum), keyby = Year]
# Annual all TB notifications by age **group**, population normalised
ir$premelt$ann_alltb_notif   <- setDT(lapply(X = aglo, function(x) psadj(ir$ann_alltb_notif_raw, x)))[, Year := ir$wth][, c("para_hash", "scen") := .(para_hash, scen)]
setcolorder(ir$premelt$ann_alltb_notif, c("Year"))

# 6.6 LTBI %
# 6.6 LTBI %
ir$ann_ltbi_raw        <- data.table(Year = ir$Years, (la["NTDS_L", , ] + la["NTDR_L", , ] + la["v_NTDS_L", , ] + la["v_NTDR_L", , ]))[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, mean), keyby = Year]
ir$premelt$ann_ltbi          <- setDT(lapply(X = aglo, function(x) psadj(ir$ann_ltbi_raw, x, denom = 100)))[, Year := ir$wth][, c("para_hash", "scen") := .(para_hash, scen)]

# 6.6.1 LTBI by DRS
ir$ann_ltbi_raw_ds        <- data.table(Year = ir$Years, (la["NTDS_L", , ]))[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, mean), keyby = Year]
ir$premelt$ann_ltbi_ds          <- setDT(lapply(X = aglo, function(x) psadj(ir$ann_ltbi_raw_ds, x, denom = 100)))[, Year := ir$wth][, c("para_hash", "scen") := .(para_hash, scen)]
ir$ann_ltbi_raw_dr        <- data.table(Year = ir$Years, (la["NTDR_L", , ]))[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, mean), keyby = Year]
ir$premelt$ann_ltbi_dr          <- setDT(lapply(X = aglo, function(x) psadj(ir$ann_ltbi_raw_dr, x, denom = 100)))[, Year := ir$wth][, c("para_hash", "scen") := .(para_hash, scen)]

# 6.7 RESOLVED
ir$ann_resolved_raw        <- data.table(Year = ir$Years, (la["NTDS_R", , ] + la["NTDR_R", , ] + la["PTDS_R", , ] + la["PTDR_R", , ]))[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, mean), keyby = Year]
ir$premelt$ann_resolved          <- setDT(lapply(X = aglo, function(x) psadj(ir$ann_resolved_raw, x, denom = 100)))[, Year := ir$wth][, c("para_hash", "scen") := .(para_hash, scen)]

ir$ann_resolved_raw_ds        <- data.table(Year = ir$Years, (la["NTDS_R", , ] + la["PTDS_R", , ]))[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, mean), keyby = Year]
ir$premelt$ann_resolved_ds          <- setDT(lapply(X = aglo, function(x) psadj(ir$ann_resolved_raw_ds, x, denom = 100)))[, Year := ir$wth][, c("para_hash", "scen") := .(para_hash, scen)]
ir$ann_resolved_raw_dr        <- data.table(Year = ir$Years, (la["NTDR_R", , ] + la["PTDR_R", , ]))[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, mean), keyby = Year]
ir$premelt$ann_resolved_dr          <- setDT(lapply(X = aglo, function(x) psadj(ir$ann_resolved_raw_dr, x, denom = 100)))[, Year := ir$wth][, c("para_hash", "scen") := .(para_hash, scen)]

# 6.8 SUSCEPTIBLE
ir$ann_susceptible_raw        <- data.table(Year = ir$Years, (la["S", , ]))[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, mean), keyby = Year]
ir$premelt$ann_susceptible         <- setDT(lapply(X = aglo, function(x) psadj(ir$ann_susceptible_raw, x, denom = 100)))[, Year := ir$wth][, c("para_hash", "scen") := .(para_hash, scen)]

# 7 Miscellaneous
# 7.1 Scaled Mortality terms
ir$scaled_ui           <- ui
ir$scaled_uni          <- uni
ir$scaled_ut           <- ut

# 7.2 Annual age-wise background mortality - Not calculated
# ir$premelt$ann_bgmort           <- data.table(Year = ir$Years, ta["bgmort", , ])[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, sum), keyby = Year]

# 7.3 Average births per year
ir$ann_births       <- data.table(Year = ir$Years, Births = colSums(la[, , 1], dims = 1))[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, mean), keyby = Year]

# 7.4 Annual compartment-wise population
compartments        <- t(rowSums(la, dim = 2))

ir$ann_compartments <- data.table(Year = ir$Years, t(rowSums(la, dim = 2)))[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, mean), keyby = Year]
rm(compartments)

# 7.5 Combined table for source of Incidence
ir$premelt$ann_inc_source <- data.table(
  Year = ir$Years,
  Total_New = rowSums(ta["ti_new", , ])/rowSums((ta["ti_rr", , ] + ta["ti_new", , ])),
  Total_New_LR = rowSums(ta["ti_new_lr", , ])/rowSums((ta["ti_rr", , ] + ta["ti_new", , ])),
  Total_New_S = rowSums(ta["ti_new_s", , ])/rowSums((ta["ti_rr", , ] + ta["ti_new", , ])),
  DR_New = rowSums(ta["dri_new", , ])/rowSums((ta["dri_rr", , ] + ta["dri_new", , ])),
  DS_New = rowSums(ta["dsi_new", , ])/rowSums((ta["dsi_rr", , ] + ta["dsi_new", , ])),
  DS_New_LR = rowSums(ta["dsi_new_lr", , ])/rowSums((ta["dsi_rr", , ] + ta["dsi_new", , ])),
  DS_New_S = rowSums(ta["dsi_new_s", , ])/rowSums((ta["dsi_rr", , ] + ta["dsi_new", , ])),
  DR_New_LR = rowSums(ta["dri_new_lr", , ])/rowSums((ta["dri_rr", , ] + ta["dri_new", , ])),
  DR_New_S = rowSums(ta["dri_new_s", , ])/rowSums((ta["dri_rr", , ] + ta["dri_new", , ]))
)[, Total_Relapse := 1 - Total_New, ][, DR_Relapse := 1 - DR_New][, DS_Relapse := 1 - DS_New][Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, mean), keyby = Year][, c("para_hash", "scen") := .(para_hash, scen)]

# 7.6 IN/II ratio
# 7.6 Combined IS Table - Experimental new method
try({
age_grouper <- stepfun(x = c(0,15,65, 100), y = c(0, A014 = 1, 2, 3, 4))

ir$melt_pop <- melt(ir$premelt$ann_agegrps, id.vars = c("Year", "para_hash", "scen"), variable.name = "MGrp", value.name = "MPop", variable.factor = F)[MGrp %chin% c('A014', 'A1564', 'A65', 'All')]
#ir$melt_pop[, MGrp := factor(MGrp, levels = c(1,2,3,4), labels = c('A014', 'A1564', 'A65', 'All'))]

psadj2 <- function(dtta, id = NULL) {
    idv = c("Year", id)
    mddta <- melt(dtta, id.vars = idv, variable.name = 'MGrp', value.name = 'MVal', variable.factor = F)
    mall <- mddta[, .(MGrp = factor(4, label = 'All'), MVal = sum(MVal)), by = c("Year", id)]
    mddta[, MGrp := factor(age_grouper(MGrp), levels = c(1,2,3), labels = c('A014', 'A1564', 'A65'))]
    opdt <- rbind(mall, mddta[, .(MVal = sum(MVal)), by = .(Year, MGrp, MGrp2)])
    #[, MGrp := factor(MGrp, levels = c(1,2,3,4), labels = c('A014', 'A1564', 'A65', 'All'))]
    setkey(opdt, Year, MGrp)
    opdt
}

ir$mcis <- rbind(
data.table(Year = ir$Years, ta["ti_new", , ], MGrp2 = "ti_new"),
data.table(Year = ir$Years, ta["ti_new_lr", , ], MGrp2 = "ti_new_lr"),
data.table(Year = ir$Years, ta["ti_new_s", , ], MGrp2 = "ti_new_s"),
data.table(Year = ir$Years, ta["dri_new", , ], MGrp2 = "dri_new"),
data.table(Year = ir$Years, ta["dsi_new", , ], MGrp2 = "dsi_new"),
data.table(Year = ir$Years, ta["dsi_new_lr", , ], MGrp2 = "dsi_new_lr"),
data.table(Year = ir$Years, ta["dsi_new_s", , ], MGrp2 = "dsi_new_s"),
data.table(Year = ir$Years, ta["dri_new_lr", , ], MGrp2 = "dri_new_lr"),
data.table(Year = ir$Years, ta["dri_new_s", , ], MGrp2 = "dri_new_s"),
data.table(Year = ir$Years, ta["dri_rr", , ], MGrp2 = "dri_rr"),
data.table(Year = ir$Years, ta["dsi_rr", , ], MGrp2 = "dsi_rr"),
data.table(Year = ir$Years, ta["ti_rr", , ], MGrp2 = "ti_rr"),
data.table(Year = ir$Years, ta["ti", , ], MGrp2 = "ti")
)

ir$mcis_a <- psadj2(ir$mcis, id = 'MGrp2')
ir$comis <- merge(ir$melt_pop, ir$mcis_a, by = c('Year', 'MGrp'))[, .(Year, para_hash, scen, MGrp, MGrp2, MVal = MVal/MPop*1E+05)]

})

# 7.6 IN/II ratio
dr.ii <- c("NTDR_II", "PTDR_II", "v_NTDR_II", "v_PTDR_II")
ds.ii <- c("NTDS_II", "PTDS_II", "v_NTDS_II", "v_PTDS_II")
dr.in <- c("NTDR_IN", "PTDR_IN", "v_NTDR_IN", "v_PTDR_IN")
ds.in <- c("NTDS_IN", "PTDS_IN", "v_NTDS_IN", "v_PTDS_IN")
ntds.ii <- c("NTDS_II", "v_NTDS_II")
ptds.ii <- c("PTDS_II", "v_PTDS_II")
ntds.in <- c("NTDS_IN", "v_NTDS_IN")
ptds.in <- c("PTDS_IN", "v_PTDS_IN")

ir$ini_ratio <- data.table(
  Year = ir$Years,
  alltb = rowSums(colSums(la[c(ds.ii, ds.in, dr.ii, dr.in), , ], dims = 1)),
  all.ii = rowSums(colSums(la[c(ds.ii, dr.ii), , ], dims = 1)),
  all.in = rowSums(colSums(la[c(ds.in, dr.in), , ], dims = 1)),
  ds.ii = rowSums(colSums(la[c(ds.ii), , ], dims = 1)),
  dr.ii = rowSums(colSums(la[c(dr.ii), , ], dims = 1)),
  ds.in = rowSums(colSums(la[c(ds.in), , ], dims = 1)),
  dr.in = rowSums(colSums(la[c(dr.in), , ], dims = 1)),
  ntds.ii = rowSums(colSums(la[c(ntds.ii), , ], dims = 1)),
  ptds.ii = rowSums(colSums(la[c(ptds.ii), , ], dims = 1)),
  ntds.in = rowSums(colSums(la[c(ntds.in), , ], dims = 1)),
  ptds.in = rowSums(colSums(la[c(ptds.in), , ], dims = 1))
)[Year >= 1900 & Year <= 2099, lapply(.SD, mean), keyby = Year]

ir$ini_ratio[, c("all_in.ii", "ds_in.ii", "dr_in.ii") := .(all.in/all.ii, ds.in/ds.ii, dr.in/dr.ii)]
ir$ini_ratio[, c("para_hash", "scen") := .(para_hash, scen)]
ir$ini_prop <- ini_prop

# Clear any NA
for (j in seq_len(ncol(ir$premelt$ann_inc_source))) {
  set(ir$premelt$ann_inc_source, which(is.na(ir$premelt$ann_inc_source[[j]])), j, 0)
}

#Melt into calc list
ir$premelt_names <- names(ir$premelt)
ir$melt_idv <- c("Year", "para_hash", "scen")

calc <- lapply(ir$premelt_names, function(x) {
  melt(ir$premelt[[x]], id.vars = ir$melt_idv, variable.name = 'MGrp', value.name = 'MVal')
  }
)

calc <- setNames(calc, ir$premelt_names)

calc$comis <- ir$comis

# 8 Vaccines
if (vaccine == 1 || vaccine == 2) {
  vxl = list()
  vxl$difflist = list()
  vxl$point_estimates <- data.table(
    para_hash = para_hash,
    alltb_inc_rate_2030 = ir$premelt$ann_alltb_inc[Year == 2030, All],
    alltb_inc_rate_2035 = ir$premelt$ann_alltb_inc[Year == 2035, All],
    alltb_inc_rate_2050 = ir$premelt$ann_alltb_inc[Year == 2050, All],
    alltb_mort_rate_2030 = ir$premelt$ann_alltb_mort[Year == 2030, All],
    alltb_mort_rate_2035 = ir$premelt$ann_alltb_mort[Year == 2035, All],
    alltb_mort_rate_2050 = ir$premelt$ann_alltb_mort[Year == 2050, All],
    drtb_inc_rate_2030 = ir$premelt$ann_drtb_inc[Year == 2030, All],
    drtb_inc_rate_2035 = ir$premelt$ann_drtb_inc[Year == 2035, All],
    drtb_inc_rate_2050 = ir$premelt$ann_drtb_inc[Year == 2050, All],
    drtb_mort_rate_2030 = ir$premelt$ann_drtb_mort[Year == 2030, All],
    drtb_mort_rate_2035 = ir$premelt$ann_drtb_mort[Year == 2035, All],
    drtb_mort_rate_2050 = ir$premelt$ann_drtb_mort[Year == 2050, All]
  )
  vxl$point_estimates <- melt(vxl$point_estimates, id.vars = "para_hash", variable.name = "MGrp", value.name = "MVal")
  vxl$point_estimates[, scen := scen]
  setkey(vxl$point_estimates, scen, para_hash, MGrp)

  vxl$difflist$ann_alltb_inc_raw = rbindlist(list(ir$ann_drtb_inc_raw, ir$ann_dstb_inc_raw))[, lapply(.SD, sum), by = Year][, c("para_hash", "scen") := .(para_hash, scen)]
  vxl$difflist$ann_alltb_mort_raw = rbindlist(list(ir$ann_drtb_mort_raw, ir$ann_dstb_mort_raw))[, lapply(.SD, sum), by = Year][, c("para_hash", "scen") := .(para_hash, scen)]
  vxl$difflist$ann_drtb_inc_raw = ir$ann_drtb_inc_raw[, c("para_hash", "scen") := .(para_hash, scen)]
  vxl$difflist$ann_drtb_mort_raw = ir$ann_drtb_mort_raw[, c("para_hash", "scen") := .(para_hash, scen)]
  vxl$difflist$ann_dstb_inc_raw = ir$ann_dstb_inc_raw[, c("para_hash", "scen") := .(para_hash, scen)]
  vxl$difflist$ann_dstb_mort_raw = ir$ann_dstb_mort_raw[, c("para_hash", "scen") := .(para_hash, scen)]
  vxl$difflist$ann_dalycmp = ir$premelt$ann_dalycmp
  vxl$difflist$ann_dstb_dalycmp = ir$premelt$ann_dstb_dalycmp
  vxl$difflist$ann_drtb_dalycmp = ir$premelt$ann_drtb_dalycmp
  vxl$difflist$ann_initRx = ir$premelt$ann_initRx

  vxl$para_hash = para_hash

  if (cost ==1 ) {
    vxl$difflist$ann_cost_array = ir$premelt$ann_cost_array
  }

  if (vaccine == 1) {
    require(R.utils)
    vxl$ann_imms <- data.table(
      Year = ir$Years,
      Routine = ca[, "v_immR"],
      Mass = ca[, "v_immM"],
      Total = ca[, "v_immR"] + ca[, "v_immM"]
    )[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, sum), keyby = Year][, c("para_hash", "scen") := .(para_hash, scen)]

    if (cost == 1) {
      vxa["costR", , ] <- ca[, "v_immR"] %o% (vxdelivR + vx_regimens)
      vxa["costM", , ] <- ca[, "v_immM"] %o% (vx_cmp_prog_cost + vxdelivM + vx_regimens)
      vxa["deliv_costR", , ] <- ca[, "v_immR"] %o% (vxdelivR + rep.int(x = 0, times = length(vx_regimens)))
      vxa["deliv_costM", , ] <- ca[, "v_immM"] %o% (vx_cmp_prog_cost + vxdelivM + rep.int(x = 0, times = length(vx_regimens)))
      vxa["costT", , ] <- vxa["costM", , ] + vxa["costR", , ]
      vxa["deliv_costT", , ] <- vxa["deliv_costM", , ] + vxa["deliv_costR", , ]
      vxl$vxa          <- data.table(Year = ir$Years, wrap.array(vxa, map = list(c(2), c(1, 3)), sep = "_"))[Year >= ir$export_start & Year <= ir$export_end, lapply(.SD, sum), keyby = Year][, c("para_hash", "scen") := .(para_hash, scen)]
    }
    ir$premelt$ann_vxcost <- vxl$vxa
    ir$premelt$ann_ims    <- vxl$ann_imms
    calc$ann_ims          <- melt(ir$premelt$ann_ims, id.vars = ir$melt_idv, variable.name = 'MGrp', value.name = 'MVal')
    calc$ann_vxcost       <- melt(ir$premelt$ann_vxcost, id.vars = ir$melt_idv, variable.name = 'MGrp', value.name = 'MVal')
  }

  vxl$mdifflist <- lapply(names(vxl$difflist),
         function(x) {
           melt(vxl$difflist[[x]], id.vars = c("Year", "para_hash", "scen"), variable.name = "MGrp", value.name = "MVal")
         })
  vxl$mdifflist <- setNames(vxl$mdifflist, names(vxl$difflist))
  # Key output data tables
  lapply(vxl$mdifflist, setkeyv, c("Year", "scen", "para_hash", "MGrp"))
}

#Ribbon mode
if (mode == 4) {
  setwd(opp)
  setDTthreads(1)

if (vaccine == 1) {
    # Construct filename stem
    fns<- sprintf("vx%s__d%s__e%s__aM%s__aR%s__wM%s__cM%s__cR%s__MI%s",
                  para_vax$vxtype,
                  para_vax$duration,
                  para_vax$effD,
                  para_vax$ageM,
                  para_vax$ageR,
                  para_vax$widthM,
                  para_vax$coverageM,
                  para_vax$coverageR,
                  para_vax$mass_interval)

    vxcl <- list(   "vxtype" = para_vax$vxtype,
                    "duration" = para_vax$duration,
                    "effD" = para_vax$effD,
                    "ageM" = para_vax$ageM,
                    "ageR" = para_vax$ageR,
                    "widthM" = para_vax$widthM,
                    "coverageM" = para_vax$coverageM,
                    "coverageR" = para_vax$coverageR,
                    "mass_interval" = para_vax$mass_interval,
                    "scen" = scen,
                    "para_hash" = para_hash)

    lapply(ls(calc), function(x) fwrite(calc[[x]][, names(vxcl) := vxcl], sprintf("%s__%s__%s_%s.csv", fns, x, scen, para_hash) ))
    if (!is.null(modespec)) {
      fwrite(ir$ann_alltb_prev_agewise[, names(vxcl) := vxcl], sprintf("%s__%s__%s_%s.csv", fns, "prevagewise", scen, para_hash), na = "NA")
    }
} else {
    # Construct filename stem
    fns<- "vx0__d0__e0__aM0__aR0__wM0__cM0__cR0__MI0"

    vxcl <- list(   "vxtype" = NA,
                    "duration" = NA,
                    "effD" = NA,
                    "ageM" = NA,
                    "ageR" = NA,
                    "widthM" = NA,
                    "coverageM" = NA,
                    "coverageR" = NA,
                    "mass_interval" = NA,
                    "scen" = scen,
                    "para_hash" = para_hash)

    lapply(ls(calc), function(x) fwrite(calc[[x]][, names(vxcl) := vxcl], sprintf("%s__%s__%s_%s.csv", fns, x, scen, para_hash), na = "NA"))
    if (!is.null(modespec)) {
      fwrite(ir$ann_alltb_prev_agewise[, names(vxcl) := vxcl], sprintf("%s__%s__%s_%s.csv", fns, "prevagewise", scen, para_hash), na = "NA")
}
  }

  return(calc)
} else if (mode == 2) {
  setDTthreads(1)
  lapply(ls(calc), function(x) {
    y   <- as.name(x)
    fwrite(eval(substitute(ir$premelt$y)), file.path("../output/trajectories", paste0(x, ".csv")))
  })
  return(1)
} else if (mode == 3) {
  return(vxl)
}

output_pop                           <- list(
                                              la = la,
                                              ta = ta,
                                              ca = ca,
                                              DS_Imatrix = DS_Imatrix,
                                              DR_Imatrix = DR_Imatrix,
                                              psizematrix = psizematrix,
                                              dsia = dsia,
                                              dria = dria,
                                              tbia = tbia)

if (vaccine == 1) output_pop$vxa <- vxa

output_internals                     <- list(
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
  comp_prop_table = comp_prop_table,
  NTDS_II_u = NTDS_II_u,
  PTDS_II_u = PTDS_II_u,
  NTDS_IN_u = NTDS_IN_u,
  PTDS_IN_u = PTDS_IN_u,
  NTDS_II_n = NTDS_II_n,
  PTDS_II_n = PTDS_II_n,
  NTDS_IN_n = NTDS_IN_n,
  PTDS_IN_n = PTDS_IN_n,
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
   prev_ratio_sdr_sy = prev_ratio_sdr_sy,
  sdr = sdr,
  dx_p = dx_p,
  mdrtb_rx = mdrtb_rx,
  arc_mdrtb_rx = arc_mdrtb_rx,
  DXTX = DXTX,
  if_mdrtb_xr = if_mdrtb_xr,
  if_mdrtb_cost = if_mdrtb_cost,
  cdr = cdr,
  ds_psi = ds_psi,
  dr_psi = dr_psi)

if (vaccine != 0) output_internals$vxl <- vxl

output_input                         <- list(
  year1 = year1,
  yearend = yearend,
  dt = dt,
  para_static = para_static,
  para_variable = para_variable)

output                               <- list(
  pop = output_pop,
  internals = output_internals,
  inputs = output_input,
  calc = calc,
  ir = ir
)
