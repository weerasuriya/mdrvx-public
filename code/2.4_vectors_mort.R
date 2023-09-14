# Apply uiscaleC
if (uiscaleC < 0 ) {
  uichild  <- ui + (uiscaleC * ui)
  unichild <- uni + (uiscaleC * uni)
  utchild  <- ut + (uiscaleC * ut)
} else if (uiscaleC >= 0) {
  uichild  <- ui + (uiscaleC * (1 - ui))
  unichild <- uni + (uiscaleC * (1 - uni))
  utchild  <- ut + (uiscaleC * (1 - ut))
}

# Apply uiscaleA
if (uiscaleA < 0 ) {
  uiadult  <- ui + (uiscaleA * ui)
  uniadult <- uni + (uiscaleA * uni)
  utadult  <- ut + (uiscaleA * ut)
} else if (uiscaleA >= 0) {
  uiadult  <- ui + (uiscaleA * (1 - ui))
  uniadult <- uni + (uiscaleA * (1 - uni))
  utadult  <- ut + (uiscaleA * (1 - ut))
}

# Apply uiscaleE
if (uiscaleE < 0 ) {
  uielderly  <- ui + (uiscaleE * ui)
  unielderly <- uni + (uiscaleE * uni)
  utelderly  <- ut + (uiscaleE * ut)
} else if (uiscaleE >= 0) {
  uielderly  <- ui + (uiscaleE * (1 - ui))
  unielderly <- uni + (uiscaleE * (1 - uni))
  utelderly  <- ut + (uiscaleE * (1 - ut))
}

# Construct age-wise mortality vectors with scaled mortality rates
NTDS_II_u  <- PTDS_II_u  <- NTDR_II_u  <- PTDR_II_u  <- c((rep(uichild, l = chiyrs)), (rep(uiadult, l = aduyrs)), (rep(uielderly, l = eldyrs)))
NTDS_IN_u  <- PTDS_IN_u  <- NTDR_IN_u  <- PTDR_IN_u  <- c((rep(unichild, l = chiyrs)), (rep(uniadult, l = aduyrs)), (rep(unielderly, l = eldyrs)))
NTDS_T_u   <- PTDS_T_u   <- c((rep(utchild, l = chiyrs)), (rep(utadult, l = aduyrs)), (rep(utelderly, l = eldyrs)))
NTDR_T_u   <- PTDR_T_u   <- rep(0, age_cls)

NTDR_mdt_u <- PTDR_mdt_u <- c((rep(uichild, l = chiyrs)), (rep(uiadult, l = aduyrs)), (rep(uielderly, l = eldyrs)))
NTDR_nid_u <- PTDR_nid_u <- c((rep(unichild, l = chiyrs)), (rep(uniadult, l = aduyrs)), (rep(unielderly, l = eldyrs)))

if (vaccine == 1) {
  v_NTDS_II_u <- NTDS_II_u
  v_NTDS_IN_u <- NTDS_IN_u
  v_NTDS_T_u <- NTDS_T_u

  v_PTDS_II_u <- PTDS_II_u
  v_PTDS_IN_u <- PTDS_IN_u
  v_PTDS_T_u <- PTDS_T_u

  v_NTDR_II_u <- NTDR_II_u
  v_NTDR_IN_u <- NTDR_IN_u
  v_NTDR_T_u <- NTDR_T_u

  v_PTDR_II_u <- PTDR_II_u
  v_PTDR_IN_u <- PTDR_IN_u
  v_PTDR_T_u <- PTDR_T_u

  v_NTDR_mdt_u <- NTDR_mdt_u
  v_PTDR_mdt_u <- PTDR_mdt_u
  v_NTDR_nid_u <- NTDR_nid_u
  v_PTDR_nid_u <- PTDR_nid_u
}
