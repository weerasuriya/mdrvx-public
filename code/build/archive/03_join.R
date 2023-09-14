
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

