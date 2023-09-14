#>>>
mainloop <- function(data, mode = 0, para_static, para_variable, para_cost, startcon, vaccine = 0, para_vax = NULL, transmission = 1, mdrac = 1, rx = 1, cost = 1, scen = 1, debug = FALSE, modespec = NULL) {
  time.start <- Sys.time()
  require(here)
  require(data.table)
  require(digest)
  #  init_dir <- setwd(here("code"))
  int_env <- environment()

  if (!exists("para_hash", envir = list2env(para_variable))) para_hash <- digest(para_variable)

  list2env(data, int_env, hash = TRUE)
  list2env(para_static, int_env, hash = TRUE)
  list2env(para_variable, int_env, hash = TRUE)
  list2env(para_cost, int_env, hash = TRUE)
  list2env(startcon, int_env, hash = TRUE)
  if (vaccine == 1) {
    list2env(para_vax, int_env, hash = TRUE)
  }
  # Set up age classes - children, adults, younger adults, elderly
  chiyrs    <- 15
  aduyrs    <- 50
  yaduyrs   <- 40
  eldyrs    <- (age_cls - chiyrs - aduyrs)

  # Global MDR switch. If mdrac==0, then disable MDR acquisition in this run.
  if (mdrac == 0) {
    if (exists("xi")) {
      rm(xi)
    }
    assign("NT_xi", 0, envir = int_env)
    assign("PT_xi", 0, envir = int_env)
    assign("v_NT_xi", 0, envir = int_env)
    assign("v_PT_xi", 0, envir = int_env)
    mdrac_yr <- year1
  } else {
    mdrac_yr <- mdrac
  }

  if (transmission == 0) {z <- 0}

  #<<<
  # Call helper functions
  source("1.1_library.R", local = TRUE)

  # Create major matrices, including the largearray la[] and loading the base line
  # population into very first time step.
  source("2.1_initialise.R", local = TRUE)

  # Assign the imported parameters to variable names for use in loop
  source("2.2_name.R", local = TRUE)

  # Source natural history (nh) outcome parameter vector calculations
  source("2.3_vectors_nh.R", local = TRUE)

  # Source TB and treatment associated associated mortality parameter vector
  # calculations
  source("2.4_vectors_mort.R", local = TRUE)

  # Source Case Detection Rate and Treatment Success matrices
  source("2.5_vectors_rx.R", local = TRUE)

  # Source vaccine related settings
  source("2.5.2_vectors_vx.R", local = TRUE)
  #>>>

  if (scen == 1) {
    mdrtb_rx <- mdrtb_rx[, .(Year, xr_sq, MDR_cost_sq)]
    setnames(mdrtb_rx, new = "xr", old = "xr_sq")
    setnames(mdrtb_rx, new = "MDR_cost", old = "MDR_cost_sq")
  } else if (scen ==2) {
    mdrtb_rx <- mdrtb_rx[, .(Year, xr_pol, MDR_cost_pol)]
    setnames(mdrtb_rx, new = 'xr', old = 'xr_pol')
    setnames(mdrtb_rx, new = 'MDR_cost', old = 'MDR_cost_pol')
  }

  if_mdrtb_xr   <- approxfun(mdrtb_rx[, .(x = Year, y = xr)], rule = 2)
  if_mdrtb_cost <- approxfun(mdrtb_rx[, .(x = Year, y = MDR_cost)], rule = 2)

  ## Main iterator
  for (i in 2:steps) {
      if (debug == TRUE) {
        if (any(la<0)) stop(e_la)
     }

    current_step <- i
    # Set the substep and current_year variables for this iteration
    current_year <- timekeeper(i, dt, year1)[2]
    substep      <- timekeeper(i, dt, year1)[1]
    YrIndex      <- current_year - year1 + 1

    # # Deaths - select age-specific background death rates for current year
    if (bgu == 1) {
      if (substep == 1) {
        u <- death_rate[YrIndex - 1, 2:(age_cls + 1)]
      } else if (substep != 1) {
        u <- death_rate[YrIndex, 2:(age_cls + 1)]
      }
    } else if (bgu == 0) {
      u <- rep(0, 100)
    }

    # If MDR enabled, then start resistance acquisition in this year
    if (mdrac != 0) {
      if (current_year < mdrac) {
        assign("NT_xi", 0, envir = int_env)
        assign("PT_xi", 0, envir = int_env)
        assign("v_NT_xi", 0, envir = int_env)
        assign("v_PT_xi", 0, envir = int_env)
      } else if (current_year >= mdrac) {
        assign("NT_xi", xi_init, envir = int_env)
        assign("PT_xi", xi_init, envir = int_env)
        assign("v_NT_xi", xi_init, envir = int_env)
        assign("v_PT_xi", xi_init, envir = int_env)
      }
    }

    # Set treatment initiation rate - DSTB
    NTDS_II_kappa <- NTDS_II_kap_main[YrIndex, 2:(age_cls + 1)]
    PTDS_II_kappa <- PTDS_II_kap_main[YrIndex, 2:(age_cls + 1)]
    NTDS_IN_kappa <- NTDS_IN_kap_main[YrIndex, 2:(age_cls + 1)]
    PTDS_IN_kappa <- PTDS_IN_kap_main[YrIndex, 2:(age_cls + 1)]

    # Set treatment initiation rate - DRTB
    NTDR_II_kappa <- NTDR_II_kap_main[YrIndex, 2:(age_cls + 1)]
    PTDR_II_kappa <- PTDR_II_kap_main[YrIndex, 2:(age_cls + 1)]
    NTDR_IN_kappa <- NTDR_IN_kap_main[YrIndex, 2:(age_cls + 1)]
    PTDR_IN_kappa <- PTDR_IN_kap_main[YrIndex, 2:(age_cls + 1)]

    # Set misdiagnosed treatment initiation rate
    NTDR_II_mdt <- NTDR_II_kap_mdt[YrIndex, 2:(age_cls + 1)]
    PTDR_II_mdt <- PTDR_II_kap_mdt[YrIndex, 2:(age_cls + 1)]
    NTDR_IN_nid <- NTDR_IN_kap_nid[YrIndex, 2:(age_cls + 1)]
    PTDR_IN_nid <- PTDR_IN_kap_nid[YrIndex, 2:(age_cls + 1)]

    if (any(
      any(NTDS_II_n + NTDS_II_kappa + NTDS_II_u + u > 1) |
      any(NTDS_IN_n + NTDS_IN_kappa + NTDS_IN_u + u > 1) |
      any(PTDS_II_n + PTDS_II_kappa + PTDS_II_u + u > 1) |
      any(PTDS_IN_n + PTDS_IN_kappa + PTDS_IN_u + u > 1) |
      any(NTDR_II_n + NTDR_II_mdt + NTDR_II_kappa + NTDR_II_u + u > 1) |
      any(NTDR_IN_n + NTDR_IN_nid + NTDR_IN_kappa + NTDR_IN_u + u > 1) |
      any(PTDR_II_n + PTDR_II_mdt + PTDR_II_kappa + PTDR_II_u + u > 1) |
      any(PTDR_IN_n + PTDR_IN_nid + PTDR_IN_kappa + PTDR_IN_u + u > 1)
   )) {
         stop("Prev Outflow Error")
      }

    # Select treatment success/failure rates
    # DSTB
    NTDS_psi    <- ds_psi(current_year) * (1 - exp(-(1/0.5) * dt)) #Exit timing constant
    NTDS_phi    <- (1 - ds_psi(current_year)) * (1 - exp(-(1/0.5) * dt))

    PTDS_psi    <- ds_psi(current_year)* (1 - exp(-(1/0.5) * dt))
    PTDS_phi    <- (1 - ds_psi(current_year)) * (1 - exp(-(1/0.5) * dt))

    # DRTB
    PTDR_IN_psi <- PTDR_II_psi <- NTDR_IN_psi <- NTDR_II_psi <- dr_psi(current_year)
    PTDR_IN_phi <- PTDR_II_phi <- NTDR_IN_phi <- NTDR_II_phi <- 1 - NTDR_IN_psi

    NTDR_mdt_exit <- PTDR_mdt_exit <- (1 - exp(-(1/0.5) * dt))
    NTDR_nid_exit <- PTDR_nid_exit <- (1 - exp(-(1/0.5) * dt))

    # DRTB Treatment compartment exit risk
    NTDR_T_IIpsi_tau <- NTDR_T_INpsi_tau <- PTDR_T_IIpsi_tau <- PTDR_T_INpsi_tau <- if_mdrtb_xr(current_year)
    NTDR_T_IIphi_tau <- NTDR_T_INphi_tau <- PTDR_T_IIphi_tau <- PTDR_T_INphi_tau <- if_mdrtb_xr(current_year)

    # Duplicate for vaccine stratum
    if (vaccine == 1) {
      v_NTDS_II_kappa <- NTDS_II_kappa
      v_PTDS_II_kappa <- PTDS_II_kappa
      v_NTDS_IN_kappa <- NTDS_IN_kappa
      v_PTDS_IN_kappa <- PTDS_IN_kappa

      v_NTDR_II_kappa <- NTDR_II_kappa
      v_PTDR_II_kappa <- PTDR_II_kappa
      v_NTDR_IN_kappa <- NTDR_IN_kappa
      v_PTDR_IN_kappa <- PTDR_IN_kappa

      v_NTDR_II_mdt <- NTDR_II_mdt
      v_PTDR_II_mdt <- PTDR_II_mdt
      v_NTDR_IN_nid <- NTDR_IN_nid
      v_PTDR_IN_nid <- PTDR_IN_nid

      v_NTDS_psi    <- NTDS_psi
      v_NTDS_phi    <- NTDS_phi

      v_PTDS_psi    <- PTDS_psi
      v_PTDS_phi    <- PTDS_phi

      v_NTDR_II_psi <- NTDR_II_psi
      v_NTDR_II_phi <- NTDR_II_phi

      v_NTDR_IN_psi <- NTDR_IN_psi
      v_NTDR_IN_phi <- NTDR_IN_phi

      v_PTDR_II_psi <- PTDR_II_psi
      v_PTDR_II_phi <- PTDR_II_phi

      v_PTDR_IN_psi <- PTDR_IN_psi
      v_PTDR_IN_phi <- PTDR_IN_phi

      v_NTDR_mdt_exit <- NTDR_mdt_exit
      v_PTDR_mdt_exit <- PTDR_mdt_exit
      v_NTDR_nid_exit <- NTDR_nid_exit
      v_PTDR_nid_exit <- PTDR_nid_exit

      v_NTDR_T_IIpsi_tau <- v_NTDR_T_INpsi_tau <- v_PTDR_T_IIpsi_tau <- v_PTDR_T_INpsi_tau <- if_mdrtb_xr(current_year)
      v_NTDR_T_IIphi_tau <- v_NTDR_T_INphi_tau <- v_PTDR_T_IIphi_tau <- v_PTDR_T_INphi_tau <- if_mdrtb_xr(current_year)
    }

    # Set up "lost to prevalent pool" parameter for MDR misdiagnosed compartments
    if (current_year < dst_start_year) {
      nt_mdt_loss  <- pt_mdt_loss  <- nt_nid_loss  <- pt_nid_loss  <- 1
    } else if (current_year >= dst_start_year) {
      nt_mdt_loss <- 1 - if_ntii_tx_corr(current_year)
      pt_mdt_loss <- 1 - if_ptii_tx_corr(current_year)
      nt_nid_loss <- 1 - if_ntin_tx_corr(current_year)
      pt_nid_loss <- 1 - if_ptin_tx_corr(current_year)
    }

    # Set up dst probs
    pt_dst_p <- if_pt_dst_p(current_year)
    nt_dst_p <- if_nt_dst_p(current_year)

    # Equation age range initialisation
    a1    <- NULL
    a0    <- NULL
    surv  <- NULL

    # Substep 1 Specific Code

    if (substep == 1) {
      # Birth model - add absolute number of births for this year into population
      # vector
      if (fert == 1) {
        la["S", i, 1] <- births[YrIndex, 2]
      } else if (fert == 0) {
        la["S", i, 1] <- 0
      }
      rm(a1, a0)
      a1 <- 2:(age_cls)
      a0 <- 1:(age_cls - 1)
    } else {
      rm(a1, a0)
      a1 <- 1:age_cls
      a0 <- 1:age_cls
    }

    # Calculate transmission terms - Lambda
    # Psizematrix check
    if (any(psizematrix[i - 1,] <= 0)) stop(e_psz)

    NTDS_lambda_raw <- -(myneta * z * DS_Imatrix[i - 1, ]/psizematrix[i - 1,])
    NTDS_lambda_raw <- colSums(NTDS_lambda_raw)
    PTDS_lambda[i, 1:age_cls] <- NTDS_lambda[i, 1:age_cls] <- t(DS_neta * (1 - exp(NTDS_lambda_raw)))

    NTDR_lambda_raw <- -(myneta * z * DR_Imatrix[i - 1, ]/psizematrix[i - 1,])
    NTDR_lambda_raw <- colSums(NTDR_lambda_raw)
    PTDR_lambda[i, 1:age_cls] <- NTDR_lambda[i, 1:age_cls] <- t(DR_neta * (1 - exp(NTDR_lambda_raw)))

    # Vaccine efficacy against infection - disabled
    if (vaccine == 1) {
      v_NTDS_lambda <- NTDS_lambda# * (1 - effI)
      v_PTDS_lambda <- PTDS_lambda# * (1 - effI)
      v_NTDR_lambda <- NTDR_lambda# * (1 - effI)
      v_PTDR_lambda <- PTDR_lambda# * (1 - effI)
    }

    # Source and run the epi Equations

    # Proportion of non-infectious among infectious DSTB

    ini_prop[current_step, ] <- colSums(la[prev_numer, (current_step-1), ])/colSums(la[prev_denom, (current_step-1), ])
    if (is.nan(ini_prop[current_step, 1])) ini_prop[current_step, 1] <- 0

    #<<<
    source("2.6.0_eq_susceptible.R", local = TRUE)
    source("2.6.1_eq_ntds.R", local = TRUE)
    source("2.6.2_eq_ntdr.R", local = TRUE)
    source("2.6.3_eq_ptds.R", local = TRUE)
    source("2.6.4_eq_ptdr.R", local = TRUE)
    #>>>

    if (vaccine == 1) {
      #<<<
      source("2.6.0_eq_susceptible_vx.R", local = TRUE)
      source("2.6.1_eq_ntds_vx.R", local = TRUE)
      source("2.6.2_eq_ntdr_vx.R", local = TRUE)
      source("2.6.3_eq_ptds_vx.R", local = TRUE)
      source("2.6.4_eq_ptdr_vx.R", local = TRUE)
      #>>>
    }

    # Reset population in 1950 using new compartment proportions
    if (substep == 1 && current_year == 1950) {
      # Generate table of proportion-in-each-compartment by age
      comp_prop_table      <- prop.table(la[, current_step, ], margin = 2)
      la[, current_step, ] <- comp_prop_table %*% diag(total_population[51, ])
      }

    # Calculate SDR
    # Check if current year = 2012, if so calculate 2011 prevalence rate.
   #  if (current_year == (sdr_start_year + 1)) {
   #    # Calculate 2011 AllTB prevalence rate
   #    sdr_start_yr_steps <- ((sdr_start_year - year1) * (1/dt)) + seq(1:(1/dt))
   #    bp_prevcmp <- c("NTDS_II", "PTDS_II", "NTDR_II", "PTDR_II", "v_NTDS_II", "v_PTDS_II", "v_NTDR_II", "v_PTDR_II")
   #    sdr_start_yr_prev <- sum(colSums(la[bp_prevcmp, sdr_start_yr_steps, ], dim = 1))/sum(la[, sdr_start_yr_steps, ]) * 1e+05
   #    rm(sdr_start_yr_steps, bp_prevcmp)
   #  }

   #  # Calculate SDR of current step
   #  if (current_year >= (sdr_start_year + 1)) {
   #    # Calculate prevalence rate of current step
   #    bp_prevcmp <- c("NTDS_II", "PTDS_II", "NTDR_II", "PTDR_II", "v_NTDS_II", "v_PTDS_II", "v_NTDR_II", "v_PTDR_II")
   #    bp_prev_current_step <- sum(la[bp_prevcmp, current_step, ])/sum(la[, current_step, ]) * 1e+05
   #    sdr[current_step] <- sdr_base * bp_prev_current_step / sdr_start_yr_prev
   #    rm(bp_prev_current_step, bp_prevcmp)
   #  }

    ## Vaccination ##

    if (vaccine == 1) {
      if (current_year >= vx_start_year) {
        # Immunisation and waning in vaccinated compartments
        # Susceptibles
        tmp_S       <- la["S", current_step,] - (la["S", current_step,] * vaxmatS[current_step, ]) + ((1 - vaxmatS[current_step, ]) * wanes[current_step, ] * la["v_S", current_step,])
        tmp_v_S     <- la["v_S", current_step,] + (la["S", current_step,] * vaxmatS[current_step, ]) - ((1 - vaxmatS[current_step, ]) * wanes[current_step, ] * la["v_S", current_step,])

        # Latents and Resolved
        tmp_NTDS_L   <- la["NTDS_L", current_step, ] - (la["NTDS_L", current_step, ] * vaxmatL[current_step, ]) + ((1 - vaxmatL[current_step, ]) * wanes[current_step, ] * (la["v_NTDS_L", current_step, ]))
        tmp_v_NTDS_L <- la["v_NTDS_L", current_step, ] + (la["NTDS_L", current_step, ] * vaxmatL[current_step, ]) - ((1 - vaxmatL[current_step, ]) * wanes[current_step, ] * (la["v_NTDS_L", current_step, ]))
        tmp_NTDR_L   <- la["NTDR_L", current_step, ] - (la["NTDR_L", current_step, ] * vaxmatL[current_step, ]) + ((1 - vaxmatL[current_step, ]) * wanes[current_step, ] * (la["v_NTDR_L", current_step, ]))
        tmp_v_NTDR_L <- la["v_NTDR_L", current_step, ] + (la["NTDR_L", current_step, ] * vaxmatL[current_step, ]) - ((1 - vaxmatL[current_step, ]) * wanes[current_step, ] * (la["v_NTDR_L", current_step, ]))
        tmp_NTDR_R   <- la["NTDR_R", current_step, ] - (la["NTDR_R", current_step, ] * vaxmatR[current_step, ]) + ((1 - vaxmatR[current_step, ]) * wanes[current_step, ] * (la["v_NTDR_R", current_step, ]))
        tmp_v_NTDR_R <- la["v_NTDR_R", current_step, ] + (la["NTDR_R", current_step, ] * vaxmatR[current_step, ]) - ((1 - vaxmatR[current_step, ]) * wanes[current_step, ] * (la["v_NTDR_R", current_step, ]))
        tmp_NTDS_R   <- la["NTDS_R", current_step, ] - (la["NTDS_R", current_step, ] * vaxmatR[current_step, ]) + ((1 - vaxmatR[current_step, ]) * wanes[current_step, ] * (la["v_NTDS_R", current_step, ]))
        tmp_v_NTDS_R <- la["v_NTDS_R", current_step, ] + (la["NTDS_R", current_step, ] * vaxmatR[current_step, ]) - ((1 - vaxmatR[current_step, ]) * wanes[current_step, ] * (la["v_NTDS_R", current_step, ]))
        tmp_PTDS_R   <- la["PTDS_R", current_step, ] - (la["PTDS_R", current_step, ] * vaxmatR[current_step, ]) + ((1 - vaxmatR[current_step, ]) * wanes[current_step, ] * (la["v_PTDS_R", current_step, ]))
        tmp_v_PTDS_R <- la["v_PTDS_R", current_step, ] + (la["PTDS_R", current_step, ] * vaxmatR[current_step, ]) - ((1 - vaxmatR[current_step, ]) * wanes[current_step, ] * (la["v_PTDS_R", current_step, ]))
        tmp_PTDR_R   <- la["PTDR_R", current_step, ] - (la["PTDR_R", current_step, ] * vaxmatR[current_step, ]) + ((1 - vaxmatR[current_step, ]) * wanes[current_step, ] * (la["v_PTDR_R", current_step, ]))
        tmp_v_PTDR_R <- la["v_PTDR_R", current_step, ] + (la["PTDR_R", current_step, ] * vaxmatR[current_step, ]) - ((1 - vaxmatR[current_step, ]) * wanes[current_step, ] * (la["v_PTDR_R", current_step, ]))

        # Number vaccinated in this time step - NB order matters - this must happen before the swap.
        immunised[current_step, ] <- (la["S", current_step,] * vaxmatI[current_step, ]) + (la["v_S", current_step,] * vaxmatI[current_step, ]) +
        (la["NTDS_L", current_step, ] * vaxmatI[current_step, ]) + (la["v_NTDS_L", current_step, ] * vaxmatI[current_step, ]) +
        (la["NTDR_L", current_step, ] * vaxmatI[current_step, ]) + (la["v_NTDR_L", current_step, ] * vaxmatI[current_step, ]) +
        (la["NTDR_R", current_step, ] * vaxmatI[current_step, ]) + (la["v_NTDR_R", current_step, ] * vaxmatI[current_step, ]) +
        (la["NTDS_R", current_step, ] * vaxmatI[current_step, ]) + (la["v_NTDS_R", current_step, ] * vaxmatI[current_step, ]) +
        (la["PTDS_R", current_step, ] * vaxmatI[current_step, ]) + (la["v_PTDS_R", current_step, ] * vaxmatI[current_step, ]) +
        (la["PTDR_R", current_step, ] * vaxmatI[current_step, ]) + (la["v_PTDR_R", current_step, ] * vaxmatI[current_step, ])

        #Swap into correct compartments
        la["S", current_step, ]        <- tmp_S
        la["v_S", current_step, ]      <- tmp_v_S
        la["NTDS_L", current_step, ]   <- tmp_NTDS_L
        la["v_NTDS_L", current_step, ] <- tmp_v_NTDS_L
        la["NTDR_L", current_step, ]   <- tmp_NTDR_L
        la["v_NTDR_L", current_step, ] <- tmp_v_NTDR_L
        la["NTDR_R", current_step, ]   <- tmp_NTDR_R
        la["v_NTDR_R", current_step, ] <- tmp_v_NTDR_R
        la["NTDS_R", current_step, ]   <- tmp_NTDS_R
        la["v_NTDS_R", current_step, ] <- tmp_v_NTDS_R
        la["PTDS_R", current_step, ]   <- tmp_PTDS_R
        la["v_PTDS_R", current_step, ] <- tmp_v_PTDS_R
        la["PTDR_R", current_step, ]   <- tmp_PTDR_R
        la["v_PTDR_R", current_step, ] <- tmp_v_PTDR_R

        # Waning between non-vaccinated compartments
        for (cmp in nonvx_cmp) {
          v_cmp                     <- paste0("v_", cmp)
          tmp_cmp                   <- la[cmp, current_step, ] + (la[v_cmp, current_step, ] * wanes[current_step, ])
          tmp_v_cmp                 <- la[v_cmp, current_step, ] - (la[v_cmp, current_step, ] * wanes[current_step, ])
          la[cmp, current_step, ]   <- tmp_cmp
          la[v_cmp, current_step, ] <- tmp_v_cmp
          rm(tmp_cmp, tmp_v_cmp)
        }
        rm(tmp_S, tmp_v_S, tmp_NTDS_L, tmp_v_NTDS_L, tmp_NTDR_L, tmp_v_NTDR_L, tmp_NTDR_R, tmp_v_NTDR_R, tmp_NTDS_R, tmp_v_NTDS_R, tmp_PTDS_R, tmp_v_PTDS_R, tmp_PTDR_R, tmp_v_PTDR_R)
      }
    }

    # Source and run the transit equations
    a5 <- 1:age_cls
    a1 <- 1:age_cls
    a0 <- 1:age_cls
    #<<<
    source("2.7_eq_transit.R", local = TRUE)
    #>>>

    # Recalculate matrices for force of infection calculation
    psizematrix[i, ] <- vapply(infage, vsf, cmp = p_compartments, FUN.VALUE = numeric(1))
    DS_Imatrix[i, ]  <- vapply(infage, vsf, cmp = dsicmp, FUN.VALUE = numeric(1))
    DR_Imatrix[i, ]  <- vapply(infage, vsf, cmp = dricmp, FUN.VALUE = numeric(1))

  }

  # Large array positivity check
  if (any(la < 0)) stop(e_la)

  # Post-loop cost calculations
  if (cost == 1) {

   sdr                                    <- vector(length = steps)
   sdr[1:(sf(sdr_start_year)[1] - 1)]     <- 0
   sdr[(sf(sdr_start_year)[1]):steps]     <- sdr_base
   # Prevalence ratio in SDR start year
   bp_prevcmp                             <- c("NTDS_II", "PTDS_II", "NTDR_II", "PTDR_II", "v_NTDS_II", "v_PTDS_II", "v_NTDR_II", "v_PTDR_II")
   prev_ratio_sdr_sy                      <- sum(rowSums(colSums(la[bp_prevcmp, sf(sdr_start_year), ], dims = 1))) / sum(rowSums(colSums(la[, sf(sdr_start_year), ], dims = 1))) * 1e+05
   prev_ratio <- rowSums(colSums(la[bp_prevcmp, , ], dims = 1)) / rowSums(colSums(la[, , ], dims = 1)) * 1e+05
   sdr                                    <- ifelse(test = sdr * prev_ratio / prev_ratio_sdr_sy < 1, yes = 1, no = sdr * prev_ratio / prev_ratio_sdr_sy)
   sdr[sf(sdr_start_year)]                <- sdr_base

    vec_drtx_cost <- if_mdrtb_cost(yrfinder(1:800))

    a5                     <- 1:age_cls
    ca[, "nv_ds_dx"]       <- (rowSums(ta["DSTB_initRx", , a5]) * sdr * ds_dx_cost) * 1000
    ca[, "nv_dr_dx"]       <- (rowSums(ta["DRTB_initRx", , a5]) * sdr * dr_dx_cost) * 1000
    ca[, "nv_dst"]         <- (rowSums(ta["DRTB_initRxLab", , a5]) * dst_cost) * 1000
    ca[, "nv_ds_tx"]       <- (rowSums(ta["DSTB_onRx", , a5]) * ds_tx_cost) * 1000
    ca[, "nv_dr_tx"]       <- (rowSums(ta["DRTB_onRx", , a5]) * vec_drtx_cost) * 1000
    ca[, "nv_tbrx"]        <- (ca[, "nv_ds_dx"] + ca[, "nv_dr_dx"] + ca[, "nv_dst"] + ca[, "nv_ds_tx"] + ca[, "nv_dr_tx"]) * (1 + prog_cost)
    ca[, "nv_grand_total"] <- ca[, "tbrx"]

    if (vaccine == 1) {
      ca[, "v_ds_dx"]      <- (rowSums(ta["v_DSTB_initRx", , a5] * sdr) * ds_dx_cost) * 1000
      ca[, "v_dr_dx"]      <- (rowSums(ta["v_DRTB_initRx", , a5] * sdr) * dr_dx_cost) * 1000
      ca[, "v_dst"]        <- (rowSums(ta["v_DRTB_initRxLab", , a5]) * dst_cost) * 1000
      ca[, "v_ds_tx"]      <- (rowSums(ta["v_DSTB_onRx", , a5]) * ds_tx_cost) * 1000
      ca[, "v_dr_tx"]      <- (rowSums(ta["v_DRTB_onRx", , a5]) * vec_drtx_cost) * 1000
      ca[, "v_tbrx"]       <- (ca[, "v_ds_dx"] + ca[, "v_dr_dx"] + ca[, "v_dst"] + ca[, "v_ds_tx"] + ca[, "v_dr_tx"]) * (1 + prog_cost)
      ca[, "v_immM"]       <- rowSums(immunised[, (ageM + 1):age_cls]) * 1000
      ca[, "v_immR"]       <- immunised[, ageR + 1] * 1000
      ca[, "v_prog"]       <- rowSums(immunised[, (ageM + 1):age_cls]) * 1000 * vx_cmp_prog_cost
    }
    ca[, cac]              <- ca[, cnvc] + ca[, cvc]
    # browser()
  }

  ### **Vaccine Flows**
  if (vaccine == 1) {
    ta[tnvc, , ]           <- ta[tvc, , ] + ta[tnvc, , ]
  }

  # Mode specific calculations and output
  # Local mode (0) or ribbon mode (2)
  if (mode == 0 | mode == 2 | mode == 3 | mode == 4) {
    #<<<
    source("2.8_calc_local.R", local = TRUE)
    #>>>
  } else if (mode == 1) {
    # cluster mode
    #<<<
    source("2.9_calc_fittargets.R", local = TRUE)
    #>>>
  }

  if (mode == 0) {
    # Local mode - return detailed output
    return(output)
  } else if (mode == 1) {
    if (any(is.nan(unlist(fo))) | any(is.na(unlist(fo)))) {stop(e_nan)}
    # Cluster mode - return only fitting values
    return(fo)
  }

}
#<<<
