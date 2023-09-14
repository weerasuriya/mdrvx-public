
# Create age-wise parameter vectors for key natural history (nh) outcome parameters
# Params: p,f,r,v,n

# Natural history 'Older adult' (age 55-65) parameter calculations
# Older adult parameters for natural cure rate, n
NTDS_II_noldadu <- (NTDS_II_n + NTDS_II_nelderly)/2
NTDS_IN_noldadu <- (NTDS_IN_n + NTDS_IN_nelderly)/2
PTDS_II_noldadu <- (PTDS_II_n + PTDS_II_nelderly)/2
PTDS_IN_noldadu <- (PTDS_IN_n + PTDS_IN_nelderly)/2

NTDR_II_noldadu <- (NTDR_II_n + NTDR_II_nelderly)/2
NTDR_IN_noldadu <- (NTDR_IN_n + NTDR_IN_nelderly)/2
PTDR_II_noldadu <- (PTDR_II_n + PTDR_II_nelderly)/2
PTDR_IN_noldadu <- (PTDR_IN_n + PTDR_IN_nelderly)/2

# Older adult parameters for reactivation Rate, v
NTDS_voldadu    <- (NTDS_vadult + NTDS_velderly)/2
NTDR_voldadu    <- (NTDR_vadult + NTDR_velderly)/2
PTDS_voldadu    <- (PTDS_vadult + PTDS_velderly)/2
PTDR_voldadu    <- (PTDR_vadult + PTDR_velderly)/2

# Older adult parameters for reactivation Rate, r
NTDS_roldadu    <- (NTDS_radult + NTDS_relderly)/2
NTDR_roldadu    <- (NTDR_radult + NTDR_relderly)/2
PTDS_roldadu    <- (PTDS_radult + PTDS_relderly)/2
PTDR_roldadu    <- (PTDR_radult + PTDR_relderly)/2

# Age-wise parameter vectors - Natural History Parameters (p,f,r,v,n)
NTDS_p = c((rep(NTDS_pchild, l = chiyrs)), (rep(NTDS_padult, l = aduyrs)), (rep(NTDS_pelderly, l = eldyrs)))
NTDS_f = c((rep(NTDS_fchild, l = chiyrs)), (rep(NTDS_fadult, l = aduyrs)), (rep(NTDS_felderly, l = eldyrs)))
NTDS_v = c((rep(NTDS_vchild, l = chiyrs)), (rep(NTDS_vadult, l = yaduyrs)), (rep(NTDS_voldadu, l = (aduyrs - yaduyrs))), (rep(NTDS_velderly, l = eldyrs)))
NTDS_r = c((rep(NTDS_rchild, l = chiyrs)), (rep(NTDS_radult, l = yaduyrs)), (rep(NTDS_roldadu, l = (aduyrs - yaduyrs))), (rep(NTDS_relderly, l = eldyrs)))

NTDS_II_n = c((rep(NTDS_II_n, l = chiyrs)), (rep(NTDS_II_n, l = yaduyrs)), (rep(NTDS_II_noldadu, l = (aduyrs - yaduyrs))), (rep(NTDS_II_nelderly, l = eldyrs)))
NTDS_IN_n = c((rep(NTDS_IN_n, l = chiyrs)), (rep(NTDS_II_n, l = yaduyrs)), (rep(NTDS_IN_noldadu, l = (aduyrs - yaduyrs))), (rep(NTDS_IN_nelderly, l = eldyrs)))

NTDR_p = c((rep(NTDR_pchild, l = chiyrs)), (rep(NTDR_padult, l = aduyrs)), (rep(NTDR_pelderly, l = eldyrs)))
NTDR_f = c((rep(NTDR_fchild, l = chiyrs)), (rep(NTDR_fadult, l = aduyrs)), (rep(NTDR_felderly, l = eldyrs)))
NTDR_v = c((rep(NTDR_vchild, l = chiyrs)), (rep(NTDR_vadult, l = yaduyrs)), (rep(NTDR_voldadu, l = (aduyrs - yaduyrs))), (rep(NTDR_velderly, l = eldyrs)))
NTDR_r = c((rep(NTDR_rchild, l = chiyrs)), (rep(NTDR_radult, l = yaduyrs)), (rep(NTDR_roldadu, l = (aduyrs - yaduyrs))), (rep(NTDR_relderly, l = eldyrs)))

NTDR_II_n = c((rep(NTDR_II_n, l = chiyrs)), (rep(NTDR_II_n, l = yaduyrs)), (rep(NTDR_II_noldadu, l = (aduyrs - yaduyrs))), (rep(NTDR_II_nelderly, l = eldyrs)))
NTDR_IN_n = c((rep(NTDR_IN_n, l = chiyrs)), (rep(NTDR_II_n, l = yaduyrs)), (rep(NTDR_IN_noldadu, l = (aduyrs - yaduyrs))), (rep(NTDR_IN_nelderly, l = eldyrs)))

PTDS_p = c((rep(PTDS_pchild, l = chiyrs)), (rep(PTDS_padult, l = aduyrs)), (rep(PTDS_pelderly, l = eldyrs)))
PTDS_f = c((rep(PTDS_fchild, l = chiyrs)), (rep(PTDS_fadult, l = aduyrs)), (rep(PTDS_felderly, l = eldyrs)))
PTDS_v = c((rep(PTDS_vchild, l = chiyrs)), (rep(PTDS_vadult, l = yaduyrs)), (rep(PTDS_voldadu, l = (aduyrs - yaduyrs))), (rep(PTDS_velderly, l = eldyrs)))
PTDS_r = c((rep(PTDS_rchild, l = chiyrs)), (rep(PTDS_radult, l = yaduyrs)), (rep(PTDS_roldadu, l = (aduyrs - yaduyrs))), (rep(PTDS_relderly, l = eldyrs)))

PTDS_II_n = c((rep(PTDS_II_n, l = chiyrs)), (rep(PTDS_II_n, l = yaduyrs)), (rep(PTDS_II_noldadu, l = (aduyrs - yaduyrs))), (rep(PTDS_II_nelderly, l = eldyrs)))
PTDS_IN_n = c((rep(PTDS_IN_n, l = chiyrs)), (rep(PTDS_II_n, l = yaduyrs)), (rep(PTDS_IN_noldadu, l = (aduyrs - yaduyrs))), (rep(PTDS_IN_nelderly, l = eldyrs)))

# PTDR_p = c((rep(PTDR_pchild, l = chiyrs)), (rep(PTDR_padult, l = aduyrs)), (rep(PTDR_pelderly, l = eldyrs)))
# Disable PTDR_L compartment
PTDR_p = c(rep(1, 100))
PTDR_f = c((rep(PTDR_fchild, l = chiyrs)), (rep(PTDR_fadult, l = aduyrs)), (rep(PTDR_felderly, l = eldyrs)))
PTDR_v = c((rep(PTDR_vchild, l = chiyrs)), (rep(PTDR_vadult, l = yaduyrs)), (rep(PTDR_voldadu, l = (aduyrs - yaduyrs))), (rep(PTDR_velderly, l = eldyrs)))
PTDR_r = c((rep(PTDR_rchild, l = chiyrs)), (rep(PTDR_radult, l = yaduyrs)), (rep(PTDR_roldadu, l = (aduyrs - yaduyrs))), (rep(PTDR_relderly, l = eldyrs)))

PTDR_II_n = c((rep(PTDR_II_n, l = chiyrs)), (rep(PTDR_II_n, l = yaduyrs)), (rep(PTDR_II_noldadu, l = (aduyrs - yaduyrs))), (rep(PTDR_II_nelderly, l = eldyrs)))
PTDR_IN_n = c((rep(PTDR_IN_n, l = chiyrs)), (rep(PTDR_II_n, l = yaduyrs)), (rep(PTDR_IN_noldadu, l = (aduyrs - yaduyrs))), (rep(PTDR_IN_nelderly, l = eldyrs)))

# PTDR_r <- c(rep(ptdr_r, max_age))

if (vaccine == 1) {
  v_NTDS_p <- NTDS_p * (1 - effD)
  v_NTDS_f <- NTDS_f
  v_NTDS_v <- NTDS_v * (1 - effD)
  v_NTDS_r <- NTDS_r * (1 - effD)

  v_NTDS_II_n <- NTDS_II_n
  v_NTDS_IN_n <- NTDS_IN_n

  v_NTDR_p <- NTDR_p * (1- effD)
  v_NTDR_f <- NTDR_f
  v_NTDR_v <- NTDR_v * (1 - effD)
  v_NTDR_r <- NTDR_r * (1 - effD)

  v_NTDR_II_n <- NTDR_II_n
  v_NTDR_IN_n <- NTDR_IN_n

  v_PTDS_p <- PTDS_p * (1 - effD)
  v_PTDS_f <- PTDS_f
  v_PTDS_v <- PTDS_v * (1 - effD)
  v_PTDS_r <- PTDS_r * (1 - effD)

  v_PTDS_II_n <- PTDS_II_n
  v_PTDS_IN_n <- PTDS_IN_n

  v_PTDR_p <- PTDR_p * (1 - effD)
  v_PTDR_f <- PTDR_f
  v_PTDR_v <- PTDR_v * (1 - effD)
  v_PTDR_r <- PTDR_r * (1 - effD)

  v_PTDR_II_n <- PTDR_II_n
  v_PTDR_IN_n <- PTDR_IN_n
}
