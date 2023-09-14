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
