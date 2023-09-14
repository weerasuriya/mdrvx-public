
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

    }

    # Store yearly cdr/kappa arrays
    # hcdr[current_year - year1 +1, ]       <- NTDS.CDR.II
    hkappa[current_year - year1 + 1, ]      <- NTDS_II_kappa
    hntdr_kappa[current_year - year1 + 1, ] <- NTDR_II_kappa
    hptdr_kappa[current_year - year1 + 1, ] <- PTDR_II_kappa
    hcdrscaling[current_year - year1 + 1, ] <- c(DS_CDRscale, DS_CDRscaleO, DS_CDRscaleE)
    hn[current_year - year1 + 1, ]          <- NTDS_II_n
    hui[current_year - year1 + 1, ]         <- NTDS_II_u

  }


  # Common core calculations

  # Mode specific calculations and output
  if (mode == 0) {
