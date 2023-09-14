
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

