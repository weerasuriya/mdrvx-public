
#rm(surv)
surv <- 1 - u[a0]

## Never Treated Drug Resistant

## v_NTDR-Latent
la["v_NTDR_L", i, a1] <- (surv * la["v_NTDR_L", i - 1, a0]) +
                       ## New Infections from Susceptible
                       ((1 - v_NTDR_p[a0]) * v_NTDR_lambda[i - 1, a0] * la["v_S", i - 1, a0]) +
                       ## New Infections from v_NTDS-Latent
                       ((1 - v_NTDR_p[a0]) * v_NTDR_lambda[i - 1, a0] * v_NTDR_x * la["v_NTDS_L", i - 1, a0]) +
                       ## Reactivation from Latent
                      -(v_NTDR_v[a0] * la["v_NTDR_L", i - 1, a0]) +
                       ## DRTB Reinfection of v_NTDR-Latent
                      -(v_NTDR_x * v_NTDR_lambda[i - 1, a0] * v_NTDR_p[a0] * la["v_NTDR_L", i - 1, a0]) +
                       ## DSTB Reinfection of v_NTDR-Latent
                      -(v_NTDR_x * v_NTDS_lambda[i - 1, a0] * la["v_NTDR_L", i - 1, a0])

## Infectious v_NTDR
la["v_NTDR_II", i,a1] <- (surv * la["v_NTDR_II", i - 1, a0]) +
                       ## Conversion of non-infectious to infectious
                       (v_NTDR_omega * la["v_NTDR_IN", i - 1, a0]) +
                       ## Deconstructed compound term - new infections in susceptble, v_NTDR-Latent and v_NTDS-Latent
                       (v_NTDR_p[a0] * v_NTDR_f[a0] * v_NTDR_lambda[i - 1, a0] * la["v_S", i - 1, a0]) +
                       (v_NTDR_p[a0] * v_NTDR_f[a0] * v_NTDR_lambda[i - 1, a0] * v_NTDR_x * la["v_NTDR_L", i - 1, a0]) +
                       (v_NTDR_p[a0] * v_NTDR_f[a0] * v_NTDR_lambda[i - 1, a0] * v_NTDR_x * la["v_NTDS_L", i - 1, a0]) +
                       ## End deconstruct
                       ## Reactivation from Latent
                       (v_NTDR_v[a0] * v_NTDR_f[a0] * la["v_NTDR_L", i - 1, a0]) +
                       ## Reinfection of Resolved - DR
                       (v_NTDR_lambda[i - 1, a0] * v_NTDR_p[a0] * v_NTDR_f[a0] * v_NTDR_x * la["v_NTDR_R", i - 1, a0]) +
                       ## Reactivation from Resolved
                       (v_NTDR_f[a0] * v_NTDR_r[a0] * la["v_NTDR_R", i - 1, a0]) +
                       ## Reinfection of Resolved - DS
                       (v_NTDR_lambda[i - 1, a0] * v_NTDR_p[a0] * v_NTDR_f[a0] * v_NTDR_x * la["v_NTDS_R", i - 1, a0]) +
                       ## Natural cure
                      -(v_NTDR_II_n[a0] * la["v_NTDR_II", i - 1, a0]) +
                       ## Case detection
                      -(v_NTDR_II_kappa[a0] * la["v_NTDR_II", i - 1, a0]) +
                       ## Misdiagnosis and Treatment
                      -(v_NTDR_II_mdt[a0] * la["v_NTDR_II", i - 1, a0]) +
                       ## TB death
                      -(v_NTDR_II_u[a0] * la["v_NTDR_II", i - 1, a0])

## Non-Infectious v_NTDR
la["v_NTDR_IN", i,a1] <- (surv * la["v_NTDR_IN", i - 1, a0]) +
                       # Deconstructed
                       (v_NTDR_p[a0] * (1 - v_NTDR_f[a0]) * v_NTDR_lambda[i - 1, a0] * la["v_S", i - 1, a0]) +
                       (v_NTDR_p[a0] * (1 - v_NTDR_f[a0]) * v_NTDR_lambda[i - 1, a0] * v_NTDR_x * la["v_NTDR_L", i - 1, a0]) +
                       (v_NTDR_p[a0] * (1 - v_NTDR_f[a0]) * v_NTDR_lambda[i - 1, a0] * v_NTDR_x * la["v_NTDS_L", i - 1, a0]) +
                       # End deconstruct
                       ## Reactivation from Latent
                       (v_NTDR_v[a0] * (1 - v_NTDR_f[a0]) * la["v_NTDR_L", i - 1, a0]) +
                       ## Reinfection of Resolved - DR
                       (v_NTDR_lambda[i - 1, a0] * v_NTDR_p[a0] * (1 - v_NTDR_f[a0]) * v_NTDR_x * la["v_NTDR_R", i - 1, a0]) +
                       ## Reactivation from Resolved
                       ((1 - v_NTDR_f[a0]) * v_NTDR_r[a0] * la["v_NTDR_R", i - 1, a0]) +
                       ## Reinfection of Resolved - DS
                       (v_NTDR_lambda[i - 1, a0] * v_NTDR_p[a0] * (1 - v_NTDR_f[a0]) * v_NTDR_x * la["v_NTDS_R", i - 1, a0]) +
                       ## Natural cure
                      -(v_NTDR_IN_n[a0] * la["v_NTDR_IN", i - 1, a0]) +
                       ## Case detection
                      -(v_NTDR_IN_kappa[a0] * la["v_NTDR_IN", i - 1, a0]) +
                       ## Misdiagnosis and Treatment
                      -(v_NTDR_IN_nid[a0] * la["v_NTDR_IN", i - 1, a0]) +
                       ## TB death
                      -(v_NTDR_IN_u[a0] * la["v_NTDR_IN", i - 1, a0]) +
                       ## Conversion from non-infectious to infectious
                      -(v_NTDR_omega * la["v_NTDR_IN", i - 1, a0])

## Treatment - v_NTDR
## Treatment Succeeding - Originating from Infectious DR-TB
la["v_NTDR_T_IIpsi", i, a1] <- (surv * la["v_NTDR_T_IIpsi", i - 1, a0]) +
                             (v_NTDR_II_kappa[a0] * v_NTDR_II_psi * la["v_NTDR_II", i - 1, a0]) +
                            -(v_NTDR_T_IIpsi_tau * la["v_NTDR_T_IIpsi", i - 1, a0]) +
                            -(v_NTDR_T_u[a0] * la["v_NTDR_T_IIpsi", i - 1, a0])

## Treatment Failing - Originating from Infectious DR-TB
la["v_NTDR_T_IIphi", i, a1] <- (surv * la["v_NTDR_T_IIphi", i - 1, a0]) +
                             (v_NTDR_II_kappa[a0] * v_NTDR_II_phi * la["v_NTDR_II", i - 1, a0]) +
                            -(v_NTDR_T_IIphi_tau * la["v_NTDR_T_IIphi", i - 1, a0]) +
                            -(v_NTDR_II_u[a0] * la["v_NTDR_T_IIphi", i - 1, a0])

## Treatment Succeeding - Originating from Non-Infectious DR-TB
la["v_NTDR_T_INpsi", i, a1] <- (surv * la["v_NTDR_T_INpsi", i - 1, a0]) +
                             (v_NTDR_IN_kappa[a0] * v_NTDR_IN_psi * la["v_NTDR_IN", i - 1, a0]) +
                            -(v_NTDR_T_INpsi_tau * la["v_NTDR_T_INpsi", i - 1, a0]) +
                            -(v_NTDR_T_u[a0] * la["v_NTDR_T_INpsi", i - 1, a0])

## Treatment Failing - Originating from Non-Infectious DR-TB
la["v_NTDR_T_INphi", i, a1] <- (surv * la["v_NTDR_T_INphi", i - 1, a0]) +
                             (v_NTDR_IN_kappa[a0] * v_NTDR_IN_phi * la["v_NTDR_IN", i - 1, a0]) +
                            -(v_NTDR_T_INphi_tau * la["v_NTDR_T_INphi", i - 1, a0]) +
                            -(v_NTDR_IN_u[a0] * la["v_NTDR_T_INphi", i - 1, a0])

## Resolved - DR TB
la["v_NTDR_R", i, a1] <- (surv * la["v_NTDR_R", i - 1, a0]) +
                       ## Natural cure - Infectious
                       (v_NTDR_II_n[a0] * la["v_NTDR_II", i - 1, a0]) +
                       ## Natural cure - Non-infectious
                       (v_NTDR_IN_n[a0] * la["v_NTDR_IN", i - 1, a0]) +
                       ## Reinfection of Resolved - to DR
                      -(v_NTDR_lambda[i - 1, a0] * v_NTDR_p[a0] * v_NTDR_x * la["v_NTDR_R", i - 1, a0]) +
                       ## Reactivation of Resolved
                      -(v_NTDR_r[a0] * la["v_NTDR_R", i - 1, a0]) +
                       ## Reinfection of Resolved - to DS
                      -(v_NTDS_lambda[i - 1, a0] * v_NTDS_p[a0] * v_NTDS_x * la["v_NTDR_R", i - 1, a0])

## Misdiagnosed and Treated - Infectious
la["v_NTDR_mdt", i, a1] <-  (surv * la["v_NTDR_mdt", i - 1, a0]) +
                          # Misdiagnosed and treated
                          (v_NTDR_II_mdt[a0] * la["v_NTDR_II", i - 1, a0]) +
                          # Exit outwards
                          -(la["v_NTDR_mdt", i - 1, a0] * v_NTDR_mdt_exit) +
                          -(v_NTDR_mdt_u[a0] * la["v_NTDR_mdt", i - 1, a0])

## Misdiagnosed and Treated - Npn-Infectious
la["v_NTDR_nid", i, a1] <-   (surv * la["v_NTDR_nid", i - 1, a0]) +
                           # Misdiagnosed and treated
                          (v_NTDR_IN_nid[a0] * la["v_NTDR_IN", i - 1, a0]) +
                           # Exit outwards
                          -(la["v_NTDR_nid", i - 1, a0] * v_NTDR_nid_exit) +
                          -(v_NTDR_nid_u[a0] * la["v_NTDR_nid", i - 1, a0])
