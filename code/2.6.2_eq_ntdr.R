
# rm(surv)
surv <- 1 - u[a0]

## Never Treated Drug Resistant

## NTDR-Latent
la["NTDR_L", i, a1] <- (surv * la["NTDR_L", i - 1, a0]) +
                       ## New Infections from Susceptible
                       ((1 - NTDR_p[a0]) * NTDR_lambda[i - 1, a0] * la["S", i - 1, a0]) +
                       ## New Infections from NTDS-Latent
                       ((1 - NTDR_p[a0]) * NTDR_lambda[i - 1, a0] * NTDR_x * la["NTDS_L", i - 1, a0]) +
                       ## Reactivation from Latent
                      -(NTDR_v[a0] * la["NTDR_L", i - 1, a0]) +
                       ## DRTB Reinfection of NTDR-Latent
                      -(NTDR_x * NTDR_lambda[i - 1, a0] * NTDR_p[a0] * la["NTDR_L", i - 1, a0]) +
                       ## DSTB Reinfection of NTDR-Latent
                      -(NTDR_x * NTDS_lambda[i - 1, a0] * la["NTDR_L", i - 1, a0])

## Infectious NTDR
la["NTDR_II", i,a1] <- (surv * la["NTDR_II", i - 1, a0]) +
                       ## Conversion of non-infectious to infectious
                       (NTDR_omega * la["NTDR_IN", i - 1, a0]) +
                       ## Deconstructed compound term - new infections in susceptble, NTDR-Latent and NTDS-Latent
                       (NTDR_p[a0] * NTDR_f[a0] * NTDR_lambda[i - 1, a0] * la["S", i - 1, a0]) +
                       (NTDR_p[a0] * NTDR_f[a0] * NTDR_lambda[i - 1, a0] * NTDR_x * la["NTDR_L", i - 1, a0]) +
                       (NTDR_p[a0] * NTDR_f[a0] * NTDR_lambda[i - 1, a0] * NTDR_x * la["NTDS_L", i - 1, a0]) +
                       ## End deconstruct
                       ## Reactivation from Latent
                       (NTDR_v[a0] * NTDR_f[a0] * la["NTDR_L", i - 1, a0]) +
                       ## Reinfection of Resolved - DR
                       (NTDR_lambda[i - 1, a0] * NTDR_p[a0] * NTDR_f[a0] * NTDR_x * la["NTDR_R", i - 1, a0]) +
                       ## Reactivation from Resolved
                       (NTDR_f[a0] * NTDR_r[a0] * la["NTDR_R", i - 1, a0]) +
                       ## Reinfection of Resolved - DS
                       (NTDR_lambda[i - 1, a0] * NTDR_p[a0] * NTDR_f[a0] * NTDR_x * la["NTDS_R", i - 1, a0]) +
                       ## Natural cure
                      -(NTDR_II_n[a0] * la["NTDR_II", i - 1, a0]) +
                       ## Case detection
                      -(NTDR_II_kappa[a0] * la["NTDR_II", i - 1, a0]) +
                       ## Misdiagnosis and Treatment
                      -(NTDR_II_mdt[a0] * la["NTDR_II", i - 1, a0]) +
                       ## TB death
                      -(NTDR_II_u[a0] * la["NTDR_II", i - 1, a0])

## Non-Infectious NTDR
la["NTDR_IN", i,a1] <- (surv * la["NTDR_IN", i - 1, a0]) +
                       # Deconstructed
                       (NTDR_p[a0] * (1 - NTDR_f[a0]) * NTDR_lambda[i - 1, a0] * la["S", i - 1, a0]) +
                       (NTDR_p[a0] * (1 - NTDR_f[a0]) * NTDR_lambda[i - 1, a0] * NTDR_x * la["NTDR_L", i - 1, a0]) +
                       (NTDR_p[a0] * (1 - NTDR_f[a0]) * NTDR_lambda[i - 1, a0] * NTDR_x * la["NTDS_L", i - 1, a0]) +
                       # End deconstruct
                       ## Reactivation from Latent
                       (NTDR_v[a0] * (1 - NTDR_f[a0]) * la["NTDR_L", i - 1, a0]) +
                       ## Reinfection of Resolved - DR
                       (NTDR_lambda[i - 1, a0] * NTDR_p[a0] * (1 - NTDR_f[a0]) * NTDR_x * la["NTDR_R", i - 1, a0]) +
                       ## Reactivation from Resolved
                       ((1 - NTDR_f[a0]) * NTDR_r[a0] * la["NTDR_R", i - 1, a0]) +
                       ## Reinfection of Resolved - DS
                       (NTDR_lambda[i - 1, a0] * NTDR_p[a0] * (1 - NTDR_f[a0]) * NTDR_x * la["NTDS_R", i - 1, a0]) +
                       ## Natural cure
                      -(NTDR_IN_n[a0] * la["NTDR_IN", i - 1, a0]) +
                       ## Case detection
                      -(NTDR_IN_kappa[a0] * la["NTDR_IN", i - 1, a0]) +
                       ## Misdiagnosis and Treatment
                      -(NTDR_IN_nid[a0] * la["NTDR_IN", i - 1, a0]) +
                       ## TB death
                      -(NTDR_IN_u[a0] * la["NTDR_IN", i - 1, a0]) +
                       ## Conversion from non-infectious to infectious
                      -(NTDR_omega * la["NTDR_IN", i - 1, a0])

## Treatment - NTDR
## Treatment Succeeding - Originating from Infectious DR-TB
la["NTDR_T_IIpsi", i, a1] <- (surv * la["NTDR_T_IIpsi", i - 1, a0]) +
                             (NTDR_II_kappa[a0] * NTDR_II_psi * la["NTDR_II", i - 1, a0]) +
                            -(NTDR_T_IIpsi_tau * la["NTDR_T_IIpsi", i - 1, a0]) +
                            -(NTDR_T_u[a0] * la["NTDR_T_IIpsi", i - 1, a0])

## Treatment Failing - Originating from Infectious DR-TB
la["NTDR_T_IIphi", i, a1] <- (surv * la["NTDR_T_IIphi", i - 1, a0]) +
                             (NTDR_II_kappa[a0] * NTDR_II_phi * la["NTDR_II", i - 1, a0]) +
                            -(NTDR_T_IIphi_tau * la["NTDR_T_IIphi", i - 1, a0]) +
                            -(NTDR_II_u[a0] * la["NTDR_T_IIphi", i - 1, a0])

## Treatment Succeeding - Originating from Non-Infectious DR-TB
la["NTDR_T_INpsi", i, a1] <- (surv * la["NTDR_T_INpsi", i - 1, a0]) +
                             (NTDR_IN_kappa[a0] * NTDR_IN_psi * la["NTDR_IN", i - 1, a0]) +
                            -(NTDR_T_INpsi_tau * la["NTDR_T_INpsi", i - 1, a0]) +
                            -(NTDR_T_u[a0] * la["NTDR_T_INpsi", i - 1, a0])

## Treatment Failing - Originating from Non-Infectious DR-TB
la["NTDR_T_INphi", i, a1] <- (surv * la["NTDR_T_INphi", i - 1, a0]) +
                             (NTDR_IN_kappa[a0] * NTDR_IN_phi * la["NTDR_IN", i - 1, a0]) +
                            -(NTDR_T_INphi_tau * la["NTDR_T_INphi", i - 1, a0]) +
                            -(NTDR_IN_u[a0] * la["NTDR_T_INphi", i - 1, a0])

## Resolved - DR TB
la["NTDR_R", i, a1] <- (surv * la["NTDR_R", i - 1, a0]) +
                       ## Natural cure - Infectious
                       (NTDR_II_n[a0] * la["NTDR_II", i - 1, a0]) +
                       ## Natural cure - Non-infectious
                       (NTDR_IN_n[a0] * la["NTDR_IN", i - 1, a0]) +
                       ## Reinfection of Resolved - to DR
                      -(NTDR_lambda[i - 1, a0] * NTDR_p[a0] * NTDR_x * la["NTDR_R", i - 1, a0]) +
                       ## Reactivation of Resolved
                      -(NTDR_r[a0] * la["NTDR_R", i - 1, a0]) +
                       ## Reinfection of Resolved - to DS
                      -(NTDS_lambda[i - 1, a0] * NTDS_p[a0] * NTDS_x * la["NTDR_R", i - 1, a0])

## Misdiagnosed and Treated - Infectious
la["NTDR_mdt", i, a1] <-  (surv * la["NTDR_mdt", i - 1, a0]) +
                          # Misdiagnosed and treated
                          (NTDR_II_mdt[a0] * la["NTDR_II", i - 1, a0]) +
                          # Exit outwards
                          -(la["NTDR_mdt", i - 1, a0] * NTDR_mdt_exit) +
                          -(NTDR_mdt_u[a0] * la["NTDR_mdt", i - 1, a0])

## Misdiagnosed and Treated - Npn-Infectious
la["NTDR_nid", i, a1] <-   (surv * la["NTDR_nid", i - 1, a0]) +
                           # Misdiagnosed and treated
                          (NTDR_IN_nid[a0] * la["NTDR_IN", i - 1, a0]) +
                           # Exit outwards
                          -(la["NTDR_nid", i - 1, a0] * NTDR_nid_exit) +
                          -(NTDR_nid_u[a0] * la["NTDR_nid", i - 1, a0])
