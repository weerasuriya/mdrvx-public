
      # Recalculate matrices for infection parameters
      # Popsize
      psizematrix[i, 1] <- sum(la[, i, 1:6])
      psizematrix[i, 2] <- sum(la[, i, 7:20])
      psizematrix[i, 3] <- sum(la[, i, 21:65])
      psizematrix[i, 4] <- sum(la[, i, 66:max_age])

      ## Total Infectious Cases by contact matrix age classes - DSTB
      DS_Imatrix[i, 1]  <- sum(la["PTDS_II", i, 1:6], la["NTDS_II", i, 1:6])
      DS_Imatrix[i, 2]  <- sum(la["PTDS_II", i, 7:20], la["NTDS_II", i, 7:20])
      DS_Imatrix[i, 3]  <- sum(la["PTDS_II", i, 21:65], la["NTDS_II", i, 21:65])
      DS_Imatrix[i, 4]  <- sum(la["PTDS_II", i, 66:max_age], la["NTDS_II", i, 66:max_age])

      ## Total Infectious Cases by contact matrix age classes - DRTB
      DR_Imatrix[i, 1]  <- sum(la["PTDR_II", i, 1:6], la["NTDR_II", i, 1:6], la["PTDR_T_IIphi", i, 1:6], la["NTDR_T_IIphi", i, 1:6], la["NTDR_mdt", i, 1:6], la["PTDR_mdt", i, 1:6])
      DR_Imatrix[i, 2]  <- sum(la["PTDR_II", i, 7:20], la["NTDR_II", i, 7:20], la["PTDR_T_IIphi", i, 7:20], la["NTDR_T_IIphi", i, 7:20], la["NTDR_mdt", i, 7:20], la["PTDR_mdt", i, 7:20])
      DR_Imatrix[i, 3]  <- sum(la["PTDR_II", i, 21:65], la["NTDR_II", i, 21:65], la["PTDR_T_IIphi", i, 21:65], la["NTDR_T_IIphi", i, 21:65], la["NTDR_mdt", i, 21:65], la["PTDR_mdt", i, 21:65])
      DR_Imatrix[i, 4]  <- sum(la["PTDR_II", i, 66:max_age], la["NTDR_II", i, 66:max_age], la["PTDR_T_IIphi", i, 66:max_age], la["NTDR_T_IIphi", i, 66:max_age], la["NTDR_mdt", i, 66:max_age], la["PTDR_mdt", i, 66:max_age])

    } else {

      # Start Step-Not-1-of-year Calculations

      # Transmission term calculation - lambda
      NTDS_lambda_raw <- colSums(-(myneta[1:4, 1:max_age]) * z * ((DS_Imatrix[i - 1, 1:4])/(psizematrix[i - 1, 1:4])))
      NTDS_lambda[i - 1, 1:max_age] <- t(DS_neta * (1 - exp(NTDS_lambda_raw)))
      PTDS_lambda <- NTDS_lambda

      NTDR_lambda_raw <- colSums(-(myneta[1:4, 1:max_age]) * z * ((DR_Imatrix[i - 1, 1:4])/(psizematrix[i - 1, 1:4])))
      NTDR_lambda[i - 1, 1:max_age] <- t(DR_neta * (1 - exp(NTDR_lambda_raw)))
      PTDR_lambda <- NTDR_lambda

      # Source and run the epi equations
      rm(a1, a0)
      a1 <- 1:max_age
      a0 <- 1:max_age

