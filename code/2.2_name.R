# The parameters fed into the function

NTDS_II_n            <- NTDS_IN_n            <- PTDS_II_n            <- PTDS_IN_n            <- n
NTDS_II_nelderly     <- NTDS_IN_nelderly     <- PTDS_II_nelderly     <- PTDS_IN_nelderly     <- nelderly

NTDR_II_n            <- n
NTDR_II_nelderly     <- nelderly

NTDR_IN_n            <- n
NTDR_IN_nelderly     <- nelderly

PTDR_II_n            <- NTDR_II_n
PTDR_II_nelderly     <- NTDR_II_nelderly

PTDR_IN_n            <- NTDR_IN_n
PTDR_IN_nelderly     <- NTDR_IN_nelderly

NTDS_fchild          <- PTDS_fchild <- NTDR_fchild <- PTDR_fchild <- fchild

NTDS_fadult          <- PTDS_fadult <- NTDR_fadult <- PTDR_fadult <- fadult

NTDS_felderly        <- PTDS_felderly <- NTDR_felderly <- PTDR_felderly <- felderly

NTDS_rchild          <- rchild
PTDS_rchild          <- rchild
NTDR_rchild          <- rchild
PTDR_rchild          <- rchild

NTDS_radult          <- radult
PTDS_radult          <- radult
NTDR_radult          <- radult
PTDR_radult          <- radult

NTDS_relderly        <- relderly
PTDS_relderly        <- relderly
NTDR_relderly        <- relderly
PTDR_relderly        <- relderly

NTDS_pchild          <- pchild
PTDS_pchild          <- pchild
NTDR_pchild          <- pchild
PTDR_pchild          <- pchild

NTDS_padult          <- padult
PTDS_padult          <- padult
NTDR_padult          <- padult
PTDR_padult          <- padult

NTDS_pelderly        <- pelderly
PTDS_pelderly        <- pelderly
NTDR_pelderly        <- pelderly
PTDR_pelderly        <- pelderly

NTDS_vchild          <- vchild
PTDS_vchild          <- vchild
NTDR_vchild          <- vchild
PTDR_vchild          <- vchild

NTDS_vadult          <- vadult
PTDS_vadult          <- vadult
NTDR_vadult          <- vadult
PTDR_vadult          <- vadult

NTDS_velderly        <- velderly
PTDS_velderly        <- velderly
NTDR_velderly        <- velderly
PTDR_velderly        <- velderly

NTDS_omega           <- omega
PTDS_omega           <- omega
NTDR_omega           <- omega
PTDR_omega           <- omega

DS_CDRscale          <- CDRscale
DR_CDRscale          <- CDRscale

DS_CDRscaleE         <- CDRscaleE
DR_CDRscaleE         <- CDRscaleE

NTDS_x               <- x
PTDS_x               <- x
NTDR_x               <- x
PTDR_x               <- x

nt_dst_prob <- pt_dst_prob <- f_dst_prob
ntin_emp_tx_p <- ntii_emp_tx_p <- ptin_emp_tx_p <- ptii_emp_tx_p <- emp_tx_p
pt_alt_dst_prob <- nt_alt_dst_prob <- alt_dst_prob 

# Parameters copied for vaccine stratum
if (vaccine == 1) {
  v_NTDS_II_n        <- n
  v_NTDS_II_nelderly <- nelderly

  v_NTDS_IN_n        <- n
  v_NTDS_IN_nelderly <- nelderly

  v_PTDS_II_n        <- v_NTDS_II_n
  v_PTDS_II_nelderly <- v_NTDS_II_nelderly

  v_PTDS_IN_n        <- v_NTDS_IN_n
  v_PTDS_IN_nelderly <- v_NTDS_IN_nelderly

  v_NTDR_II_n        <- n
  v_NTDR_II_nelderly <- nelderly

  v_NTDR_IN_n        <- n
  v_NTDR_IN_nelderly <- nelderly

  v_PTDR_II_n        <- v_NTDR_II_n
  v_PTDR_II_nelderly <- v_NTDR_II_nelderly

  v_PTDR_IN_n        <- v_NTDR_IN_n
  v_PTDR_IN_nelderly <- v_NTDR_IN_nelderly

  v_NTDS_fchild      <- fchild
  v_PTDS_fchild      <- fchild
  v_NTDR_fchild      <- fchild
  v_PTDR_fchild      <- fchild

  v_NTDS_fadult      <- fadult
  v_PTDS_fadult      <- fadult
  v_NTDR_fadult      <- fadult
  v_PTDR_fadult      <- fadult

  v_NTDS_felderly    <- felderly
  v_PTDS_felderly    <- felderly
  v_NTDR_felderly    <- felderly
  v_PTDR_felderly    <- felderly

  v_NTDS_rchild      <- rchild
  v_PTDS_rchild      <- rchild
  v_NTDR_rchild      <- rchild
  v_PTDR_rchild      <- rchild

  v_NTDS_radult      <- radult
  v_PTDS_radult      <- radult
  v_NTDR_radult      <- radult
  v_PTDR_radult      <- radult

  v_NTDS_relderly    <- relderly
  v_PTDS_relderly    <- relderly
  v_NTDR_relderly    <- relderly
  v_PTDR_relderly    <- relderly

  v_NTDS_pchild      <- pchild
  v_PTDS_pchild      <- pchild
  v_NTDR_pchild      <- pchild
  v_PTDR_pchild      <- pchild

  v_NTDS_padult      <- padult
  v_PTDS_padult      <- padult
  v_NTDR_padult      <- padult
  v_PTDR_padult      <- padult

  v_NTDS_pelderly    <- pelderly
  v_PTDS_pelderly    <- pelderly
  v_NTDR_pelderly    <- pelderly
  v_PTDR_pelderly    <- pelderly

  v_NTDS_vchild      <- vchild
  v_PTDS_vchild      <- vchild
  v_NTDR_vchild      <- vchild
  v_PTDR_vchild      <- vchild

  v_NTDS_vadult      <- vadult
  v_PTDS_vadult      <- vadult
  v_NTDR_vadult      <- vadult
  v_PTDR_vadult      <- vadult

  v_NTDS_velderly    <- velderly
  v_PTDS_velderly    <- velderly
  v_NTDR_velderly    <- velderly
  v_PTDR_velderly    <- velderly

  v_NTDS_omega       <- omega
  v_PTDS_omega       <- omega
  v_NTDR_omega       <- omega
  v_PTDR_omega       <- omega

  v_NTDS_x           <- x
  v_PTDS_x           <- x
  v_NTDR_x           <- x
  v_PTDR_x           <- x
}
