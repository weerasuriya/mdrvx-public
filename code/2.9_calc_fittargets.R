fo                              <- list()
tmp                             <- list()

# Total tmp$populations of interest (large array slices)
tmp$y2000                       <- (la[, 401, ] + la[, 402, ] + la[, 403, ] + la[, 404, ])/4
tmp$y2010                       <- (la[, 441, ] + la[, 442, ] + la[, 443, ] + la[, 444, ])/4
tmp$y2015                       <- (la[, 461, ] + la[, 462, ] + la[, 463, ] + la[, 464, ])/4
tmp$y2017                       <- (la[, 469, ] + la[, 470, ] + la[, 471, ] + la[, 472, ])/4
tmp$y2013                       <- (la[, 453, ] + la[, 454, ] + la[, 455, ] + la[, 456, ])/4
tmp$y2007                       <- (la[, 429, ] + la[, 430, ] + la[, 431, ] + la[, 432, ])/4

tmp$pop_total2000               <- sum(tmp$y2000)
tmp$pop_total2010               <- sum(tmp$y2010)
tmp$pop_total2017               <- sum(tmp$y2017)
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
fo$alltb_prev_allad_2000 <- 1e+05 * sum(tmp$bprevarray2000[16:100])/tmp$pop_total2000_allad
# 2010, normalised all ages, normalised
fo$alltb_prev_allad_2010 <- 1e+05 * sum(tmp$bprevarray2010[16:100])/tmp$pop_total2010_allad
# 2000, normalised 15-29 ages, normalised
fo$alltb_prev_1529_2000  <- 1e+05 * sum(tmp$bprevarray2000[16:30])/sum(tmp$y2000[, 16:30])
# 2010, normalised 15-29 ages, normalised
fo$alltb_prev_1529_2010  <- 1e+05 * sum(tmp$bprevarray2010[16:30])/sum(tmp$y2010[, 16:30])
# 2000, normalised 30-44 ages, normalised
fo$alltb_prev_3044_2000  <- 1e+05 * sum(tmp$bprevarray2000[31:45])/sum(tmp$y2000[, 31:45])
# 2010, normalised 30-44 ages, normalised
fo$alltb_prev_3044_2010  <- 1e+05 * sum(tmp$bprevarray2010[31:45])/sum(tmp$y2010[, 31:45])
# 2000, normalised 45-59 ages, normalised
fo$alltb_prev_4559_2000  <- 1e+05 * sum(tmp$bprevarray2000[46:60])/sum(tmp$y2000[, 46:60])
# 2010, normalised 45-59 ages, normalised
fo$alltb_prev_4559_2010  <- 1e+05 * sum(tmp$bprevarray2010[46:60])/sum(tmp$y2010[, 46:60])
# 2000, normalised 60+ ages, normalised
fo$alltb_prev_60p_2000    <- 1e+05 * sum(tmp$bprevarray2000[61:100])/sum(tmp$y2000[, 61:100])
# 2010, normalised 60+ ages, normalised
fo$alltb_prev_60p_2010    <- 1e+05 * sum(tmp$bprevarray2010[61:100])/sum(tmp$y2010[, 61:100])

# All TB Incidence
# 2000, incidence normalised, all ages
fo$alltb_inc_all_2000      <- 1e+05 * sum(ta[c("DSTB_inc", "DRTB_inc"), 401:404, ])/tmp$pop_total2000
# 2010, incidence normalised, all ages
fo$alltb_inc_all_2017      <- 1e+05 * sum(ta[c("DSTB_inc", "DRTB_inc"), 469:472, ])/tmp$pop_total2017


# All TB mortality rate 2000 mortality normalised, all ages
fo$alltb_mort_all_2010          <- 1e+05 * sum(ta[c("DSTB_deaths", "DRTB_deaths"), 441:444, ])/tmp$pop_total2010
# 2000 mortality normalised, 0-14
fo$alltb_mort_014_2010          <- 1e+05 * sum(ta[c("DSTB_deaths", "DRTB_deaths"), 441:444, 1:15])/tmp$pop_total2010_014
# 2000 mortality normalised, 0-14
fo$alltb_mort_1564_2010         <- 1e+05 * sum(ta[c("DSTB_deaths", "DRTB_deaths"), 441:444, 16:65])/tmp$pop_total2010_1564
# 2000 mortality normalised, 0-14
fo$alltb_mort_65p_2010           <- 1e+05 * sum(ta[c("DSTB_deaths", "DRTB_deaths"), 441:444, 66:100])/tmp$pop_total2010_65

# All TB notification rate 2000 notification normalised, all ages
fo$alltb_notif_all_2015         <- 1e+05 * sum(ta[c("DSTB_initRx", "DRTB_initRx"), 461:464, ])/tmp$pop_total2015
# 2000 notification normalised, 0-14
fo$alltb_notif_014_2015         <- 1e+05 * sum(ta[c("DSTB_initRx", "DRTB_initRx"), 461:464, 1:15])/tmp$pop_total2015_014
# 2000 notification normalised, 0-14
fo$alltb_notif_1564_2015        <- 1e+05 * sum(ta[c("DSTB_initRx", "DRTB_initRx"), 461:464, 16:65])/tmp$pop_total2015_1564
# 2000 notification normalised, 0-14
fo$alltb_notif_65p_2015          <- 1e+05 * sum(ta[c("DSTB_initRx", "DRTB_initRx"), 461:464, 66:100])/tmp$pop_total2015_65

# DRTB Calculations
# TOtal incidence 2016
fo$drtb_inc_2017                <- 1e+05 * sum(ta["DRTB_inc", 469:472, ])/tmp$pop_total2017

# Notification proportion - 2007
fo$drtb_nt_notp_2007       <- 100 * sum(ta["DR_nt_intx_f", 429:432, ])/sum(ta["All_nt_intx_f", 429:432, ])
fo$drtb_pt_notp_2007       <- 100 * sum(ta["DR_pt_intx_f", 429:432, ])/sum(ta["All_pt_intx_f", 429:432, ])

# Notification proportion - 2013
fo$drtb_nt_notp_2013       <- 100 * sum(ta["DR_nt_intx_f", 453:456, ])/sum(ta["All_nt_intx_f", 453:456, ])
fo$drtb_pt_notp_2013       <- 100 * sum(ta["DR_pt_intx_f", 453:456, ])/sum(ta["All_pt_intx_f", 453:456, ])

# DRTB Treatment initiations - 2017
fo$dr_treated2017               <- 1e+03 * sum(ta["DRTB_initRx", 469:472, ]) * chr

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
