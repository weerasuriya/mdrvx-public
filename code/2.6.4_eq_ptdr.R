
# rm(surv)
surv <- 1 - u[a0]

## Previously Treated Drug Resistant
## PTDR-Latent
la["PTDR_L", i, a1] <- (surv * la["PTDR_L", i - 1, a0]) +
                       ((PT_xi * (1 - PTDR_p[a0])) * (la["PTDS_T", i - 1, a0])) +
                       ((NT_xi * (1 - PTDR_p[a0])) * (la["NTDS_T", i - 1, a0])) +
                      -(PTDR_v[a0] * la["PTDR_L", i - 1, a0])

## Infectious PTDR
la["PTDR_II", i, a1] <- (surv * la["PTDR_II", i - 1, a0]) +
                        ## Conversion of non-infectious to infectious
                        (PTDR_omega * la["PTDR_IN", i - 1, a0]) +
                        ## Reactivation from Latent
                        (PTDR_v[a0] * PTDR_f[a0] * la["PTDR_L", i - 1, a0]) +
                        ## Reinfection of Resolved - DR
                        (PTDR_lambda[i - 1, a0] * PTDR_p[a0] * PTDR_f[a0] * PTDR_x * la["PTDR_R", i - 1, a0]) +
                        ## Reactivation from Resolved
                        (PTDR_f[a0] * (PTDR_r[a0] ) * la["PTDR_R", i - 1, a0]) +
                        ## Reinfection of Resolved - DS
                        (PTDR_lambda[i - 1, a0] * PTDR_p[a0] * PTDR_f[a0] * PTDR_x * la["PTDS_R", i - 1, a0]) +
                        ## MDR Acquisition
                        (NT_xi * PTDR_p[a0] * PTDR_f[a0] * la["NTDS_T", i - 1, a0] * nt_mdt_loss) +
                        (PT_xi * PTDR_p[a0] * PTDR_f[a0] * la["PTDS_T", i - 1, a0] * pt_mdt_loss) +
                        ## PTDR Treatment FAILURES - Infectious
                        (PTDR_T_IIphi_tau * la["PTDR_T_IIphi", i - 1, a0]) +
                        ## NTDR Treatment FAILURES - Infectious
                        (NTDR_T_IIphi_tau * la["NTDR_T_IIphi", i - 1, a0]) +
                        # Misdiagnosed - re-entry from Infectious MDT compartment
                        (la["NTDR_mdt", i - 1, a0] * NTDR_mdt_exit * nt_mdt_loss) +
                        (la["PTDR_mdt", i - 1, a0] * PTDR_mdt_exit * pt_mdt_loss) +
                        ## Natural cure
                       -(PTDR_II_n[a0] * la["PTDR_II", i - 1, a0]) +
                        ## Case detection
                       -(PTDR_II_kappa[a0] * la["PTDR_II", i - 1, a0]) +
                        ## Misdiagnosis and treatment
                       -(PTDR_II_mdt[a0] * la["PTDR_II", i - 1, a0]) +
                        ## TB death
                       -(PTDR_II_u[a0] * la["PTDR_II", i - 1, a0])

## Non-Infectious PTDR
la["PTDR_IN", i, a1] <- (surv * la["PTDR_IN", i - 1, a0]) +
                        ## Reactivation from Latent
                        (PTDR_v[a0] * (1 - PTDR_f[a0]) * la["PTDR_L", i - 1, a0]) +
                        ## Reinfection of Resolved - DS
                        (PTDR_lambda[i - 1, a0] * PTDR_p[a0] * (1 - PTDR_f[a0]) * PTDR_x * la["PTDS_R", i - 1, a0]) +
                        ## Reactivation from Resolved
                        ((1 - PTDR_f[a0]) * PTDR_r[a0] * la["PTDR_R", i - 1, a0]) +
                        ## Reinfection of Resolved - DR
                        (PTDR_lambda[i - 1, a0] * PTDR_p[a0] * (1 - PTDR_f[a0]) * PTDR_x * la["PTDR_R", i - 1, a0]) +
                        ## MDR Acquisition
                        (PT_xi * PTDR_p[a0] * (1 - PTDR_f[a0]) * la["PTDS_T", i - 1, a0] * pt_nid_loss) +
                        (NT_xi * PTDR_p[a0] * (1 - PTDR_f[a0]) * la["NTDS_T", i - 1, a0] * nt_nid_loss) +
                        ## PTDR Treatment FAILURES - Non-Infectious
                        (PTDR_T_INphi_tau * la["PTDR_T_INphi", i - 1, a0]) +
                        ## NTDR Treatment FAILURES - Non-Infectious
                        (NTDR_T_INphi_tau * la["NTDR_T_INphi", i - 1, a0]) +
                        # Misdiagnosed - re-entry from Non-Infectious NID compartment
                        (la["NTDR_nid", i - 1, a0] * NTDR_nid_exit * nt_nid_loss) +
                        (la["PTDR_nid", i - 1, a0] * PTDR_nid_exit * pt_nid_loss) +
                        ## Natural cure
                       -(PTDR_IN_n[a0] *la["PTDR_IN", i - 1, a0]) +
                        ## Case detection
                       -(PTDR_IN_kappa[a0] * la["PTDR_IN", i - 1, a0]) +
                        ## Misdiagnosis and treatment
                       -(PTDR_IN_nid[a0] * la["PTDR_IN", i - 1, a0]) +
                        ## TB death
                       -(PTDR_IN_u[a0] * la["PTDR_IN", i - 1, a0]) +
                        ## Conversion from non-infectious to infectious
                       -(PTDR_omega * la["PTDR_IN", i - 1, a0])

## PTDR on treatment
## Treatment Succeeding - Originating from Infectious DR-TB
la["PTDR_T_IIpsi", i, a1] <-  (surv * la["PTDR_T_IIpsi", i - 1, a0]) +
                              (PTDR_II_kappa[a0] * PTDR_II_psi * la["PTDR_II", i - 1, a0]) +
                              # Misdiagnosed - re-entry from MDR compartment
                              (la["NTDR_mdt", i - 1, a0] * NTDR_mdt_exit * (1 - nt_mdt_loss) * PTDR_II_psi) +
                              (la["PTDR_mdt", i - 1, a0] * PTDR_mdt_exit * (1 - pt_mdt_loss) * PTDR_II_psi) +
                              (NT_xi * PTDR_p[a0] * PTDR_f[a0] * la["NTDS_T", i - 1, a0] * (1 - nt_mdt_loss) * PTDR_II_psi) +
                              (PT_xi * PTDR_p[a0] * PTDR_f[a0] * la["PTDS_T", i - 1, a0] * (1 - pt_mdt_loss) * PTDR_II_psi) +
                             -(PTDR_T_IIpsi_tau * la["PTDR_T_IIpsi", i - 1, a0]) +
                             -(PTDR_T_u[a0] * la["PTDR_T_IIpsi", i - 1, a0])

## Treatment Failing - Originating from Infectious DR-TB
la["PTDR_T_IIphi", i, a1] <-  (surv * la["PTDR_T_IIphi", i - 1, a0]) +
                              (PTDR_II_kappa[a0] * PTDR_II_phi * la["PTDR_II", i - 1, a0]) +
                              # Misdiagnosed - re-entry from MDR compartment
                              (la["NTDR_mdt", i - 1, a0] * NTDR_mdt_exit * (1 - nt_mdt_loss) * PTDR_II_phi) +
                              (la["PTDR_mdt", i - 1, a0] * PTDR_mdt_exit * (1 - pt_mdt_loss) * PTDR_II_phi) +
                              (NT_xi * PTDR_p[a0] * PTDR_f[a0] * la["NTDS_T", i - 1, a0] * (1 - nt_mdt_loss) * PTDR_II_phi) +
                              (PT_xi * PTDR_p[a0] * PTDR_f[a0] * la["PTDS_T", i - 1, a0] * (1 - pt_mdt_loss) * PTDR_II_phi) +
                             -(PTDR_T_IIphi_tau * la["PTDR_T_IIphi", i - 1, a0]) +
                             -(PTDR_II_u[a0] * la["PTDR_T_IIphi", i - 1, a0])

## Treatment Succeeding - Originating from Non-Infectious DR-TB
la["PTDR_T_INpsi", i, a1] <-  (surv * la["PTDR_T_INpsi", i - 1, a0]) +
                              (PTDR_IN_kappa[a0] * PTDR_IN_psi * la["PTDR_IN", i - 1, a0]) +
                              # Misdiagnosed - re-entry from MDR compartment
                              (la["NTDR_nid", i - 1, a0] * NTDR_nid_exit * (1 - nt_nid_loss) * PTDR_IN_psi) +
                              (la["PTDR_nid", i - 1, a0] * PTDR_nid_exit * (1 - pt_nid_loss) * PTDR_IN_psi) +
                              (NT_xi * PTDR_p[a0] * (1 - PTDR_f[a0]) * la["NTDS_T", i - 1, a0] * (1 - nt_nid_loss) * PTDR_IN_psi) +
                              (PT_xi * PTDR_p[a0] * (1 - PTDR_f[a0]) * la["PTDS_T", i - 1, a0] * (1 - pt_nid_loss) * PTDR_IN_psi) +
                             -(PTDR_T_INpsi_tau * la["PTDR_T_INpsi", i - 1, a0]) +
                             -(PTDR_T_u[a0] * la["PTDR_T_INpsi", i - 1, a0])

## Treatment Failing - Originating from Non-Infectious DR-TB
la["PTDR_T_INphi", i, a1] <-  (surv * la["PTDR_T_INphi", i - 1, a0]) +
                              (PTDR_IN_kappa[a0] * PTDR_IN_phi * la["PTDR_IN", i - 1, a0]) +
                              # Misdiagnosed - re-entry from MDR compartment
                              (la["NTDR_nid", i - 1, a0] * NTDR_nid_exit * (1 - nt_nid_loss) * PTDR_IN_phi) +
                              (la["PTDR_nid", i - 1, a0] * PTDR_nid_exit * (1 - pt_nid_loss) * PTDR_IN_phi) +
                              (NT_xi * PTDR_p[a0] * (1 - PTDR_f[a0]) * la["NTDS_T", i - 1, a0] * (1 - nt_nid_loss) * PTDR_IN_phi) +
                              (PT_xi * PTDR_p[a0] * (1 - PTDR_f[a0]) * la["PTDS_T", i - 1, a0] * (1 - pt_nid_loss) * PTDR_IN_phi) +
                             -(PTDR_T_INphi_tau * la["PTDR_T_INphi", i - 1, a0]) +
                             -(PTDR_IN_u[a0] * la["PTDR_T_INphi", i - 1, a0])

## PTDR Resolved
la["PTDR_R", i, a1] <- (surv * la["PTDR_R", i - 1, a0]) +
                       ## Natural cure - Infectious
                       (PTDR_II_n[a0] * la["PTDR_II", i - 1, a0]) +
                       ## Natural cure - Non-infectious
                       (PTDR_IN_n[a0] * la["PTDR_IN", i - 1, a0]) +
                       ## NTDR Treatment Success - Infectious
                       (NTDR_T_IIpsi_tau * la["NTDR_T_IIpsi", i - 1, a0]) +
                       ## NTDR Treatment Success - Non-infectious
                       (NTDR_T_INpsi_tau * la["NTDR_T_INpsi", i - 1, a0]) +
                       ## PTDR Treatment Success - Infectious
                       (PTDR_T_INpsi_tau * la["PTDR_T_INpsi", i - 1, a0]) +
                       ## PTDR Treatment Success - Non-infectious
                       (PTDR_T_IIpsi_tau * la["PTDR_T_IIpsi", i - 1, a0]) +
                       ## Reinfection of Resolved - to DS
                      -(PTDS_lambda[i - 1, a0] * PTDS_p[a0] * PTDS_x * la["PTDR_R", i - 1, a0]) +
                       ## Reactivation of Resolved
                      -(PTDR_r[a0] * la["PTDR_R", i - 1, a0]) +
                      ## Reinfection of Resolved - to DR
                      -(PTDR_lambda[i - 1, a0] * PTDR_p[a0] * PTDR_x * la["PTDR_R", i - 1, a0])

## Misdiagnosed and Treated - Infectious
la["PTDR_mdt", i, a1] <-  (surv * la["PTDR_mdt", i - 1, a0]) +
                          # Misdiagnosed and treated
                          (PTDR_II_mdt[a0] * la["PTDR_II", i - 1, a0]) +
                          # Exit outwards
                          -(la["PTDR_mdt", i - 1, a0] * PTDR_mdt_exit) +
                          -(PTDR_mdt_u[a0] * la["PTDR_mdt", i - 1, a0])

## Misdiagnosed and Treated - Npn-Infectious
la["PTDR_nid", i, a1] <-   (surv * la["PTDR_nid", i - 1, a0]) +
                           # Misdiagnosed and treated
                          (PTDR_IN_nid[a0] * la["PTDR_IN", i - 1, a0]) +
                           # Exit outwards
                          -(la["PTDR_nid", i - 1, a0] * PTDR_nid_exit) +
                          -(PTDR_nid_u[a0] * la["PTDR_nid", i - 1, a0])
