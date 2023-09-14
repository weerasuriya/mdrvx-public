
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

