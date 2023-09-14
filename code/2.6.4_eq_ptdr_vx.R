
#rm(surv)
surv <- 1 - u[a0]

## Previously Treated Drug Resistant
## v_PTDR-Latent
la["v_PTDR_L", i, a1] <- (surv * la["v_PTDR_L", i - 1, a0]) +
                       ((v_PT_xi * (1 - v_PTDR_p[a0])) * (la["v_PTDS_T", i - 1, a0])) +
                       ((v_NT_xi * (1 - v_PTDR_p[a0])) * (la["v_NTDS_T", i - 1, a0])) +
                      -(v_PTDR_v[a0] * la["v_PTDR_L", i - 1, a0])

## Infectious v_PTDR
la["v_PTDR_II", i, a1] <- (surv * la["v_PTDR_II", i - 1, a0]) +
                        ## Conversion of non-infectious to infectious
                        (v_PTDR_omega * la["v_PTDR_IN", i - 1, a0]) +
                        ## Reactivation from Latent
                        (v_PTDR_v[a0] * v_PTDR_f[a0] * la["v_PTDR_L", i - 1, a0]) +
                        ## Reinfection of Resolved - DR
                        (v_PTDR_lambda[i - 1, a0] * v_PTDR_p[a0] * v_PTDR_f[a0] * v_PTDR_x * la["v_PTDR_R", i - 1, a0]) +
                        ## Reactivation from Resolved
                        (v_PTDR_f[a0] * (v_PTDR_r[a0] ) * la["v_PTDR_R", i - 1, a0]) +
                        ## Reinfection of Resolved - DS
                        (v_PTDR_lambda[i - 1, a0] * v_PTDR_p[a0] * v_PTDR_f[a0] * v_PTDR_x * la["v_PTDS_R", i - 1, a0]) +
                        ## MDR Acquisition
                        (v_NT_xi * v_PTDR_p[a0] * v_PTDR_f[a0] * la["v_NTDS_T", i - 1, a0] * nt_mdt_loss) +
                        (v_PT_xi * v_PTDR_p[a0] * v_PTDR_f[a0] * la["v_PTDS_T", i - 1, a0] * pt_mdt_loss) +
                        ## v_PTDR Treatment FAILURES - Infectious
                        (v_PTDR_T_IIphi_tau * la["v_PTDR_T_IIphi", i - 1, a0]) +
                        ## v_NTDR Treatment FAILURES - Infectious
                        (v_NTDR_T_IIphi_tau * la["v_NTDR_T_IIphi", i - 1, a0]) +
                        # Misdiagnosed - re-entry from Infectious MDT compartment
                        (la["v_NTDR_mdt", i - 1, a0] * v_NTDR_mdt_exit * nt_mdt_loss) +
                        (la["v_PTDR_mdt", i - 1, a0] * v_PTDR_mdt_exit * pt_mdt_loss) +
                        ## Natural cure
                       -(v_PTDR_II_n[a0] * la["v_PTDR_II", i - 1, a0]) +
                        ## Case detection
                       -(v_PTDR_II_kappa[a0] * la["v_PTDR_II", i - 1, a0]) +
                        ## Misdiagnosis and treatment
                       -(v_PTDR_II_mdt[a0] * la["v_PTDR_II", i - 1, a0]) +
                        ## TB death
                       -(v_PTDR_II_u[a0] * la["v_PTDR_II", i - 1, a0])

## Non-Infectious v_PTDR
la["v_PTDR_IN", i, a1] <- (surv * la["v_PTDR_IN", i - 1, a0]) +
                        ## Reactivation from Latent
                        (v_PTDR_v[a0] * (1 - v_PTDR_f[a0]) * la["v_PTDR_L", i - 1, a0]) +
                        ## Reinfection of Resolved - DS
                        (v_PTDR_lambda[i - 1, a0] * v_PTDR_p[a0] * (1 - v_PTDR_f[a0]) * v_PTDR_x * la["v_PTDS_R", i - 1, a0]) +
                        ## Reactivation from Resolved
                        ((1 - v_PTDR_f[a0]) * v_PTDR_r[a0] * la["v_PTDR_R", i - 1, a0]) +
                        ## Reinfection of Resolved - DR
                        (v_PTDR_lambda[i - 1, a0] * v_PTDR_p[a0] * (1 - v_PTDR_f[a0]) * v_PTDR_x * la["v_PTDR_R", i - 1, a0]) +
                        ## MDR Acquisition
                        (v_PT_xi * v_PTDR_p[a0] * (1 - v_PTDR_f[a0]) * la["v_PTDS_T", i - 1, a0] * pt_nid_loss) +
                        (v_NT_xi * v_PTDR_p[a0] * (1 - v_PTDR_f[a0]) * la["v_NTDS_T", i - 1, a0] * nt_nid_loss) +
                        ## v_PTDR Treatment FAILURES - Non-Infectious
                        (v_PTDR_T_INphi_tau * la["v_PTDR_T_INphi", i - 1, a0]) +
                        ## v_NTDR Treatment FAILURES - Non-Infectious
                        (v_NTDR_T_INphi_tau * la["v_NTDR_T_INphi", i - 1, a0]) +
                        # Misdiagnosed - re-entry from Non-Infectious NID compartment
                        (la["v_NTDR_nid", i - 1, a0] * v_NTDR_nid_exit * nt_nid_loss) +
                        (la["v_PTDR_nid", i - 1, a0] * v_PTDR_nid_exit * pt_nid_loss) +
                        ## Natural cure
                       -(v_PTDR_IN_n[a0] *la["v_PTDR_IN", i - 1, a0]) +
                        ## Case detection
                       -(v_PTDR_IN_kappa[a0] * la["v_PTDR_IN", i - 1, a0]) +
                        ## Misdiagnosis and treatment
                       -(v_PTDR_IN_nid[a0] * la["v_PTDR_IN", i - 1, a0]) +
                        ## TB death
                       -(v_PTDR_IN_u[a0] * la["v_PTDR_IN", i - 1, a0]) +
                        ## Conversion from non-infectious to infectious
                       -(v_PTDR_omega * la["v_PTDR_IN", i - 1, a0])

## v_PTDR on treatment
## Treatment Succeeding - Originating from Infectious DR-TB
la["v_PTDR_T_IIpsi", i, a1] <-  (surv * la["v_PTDR_T_IIpsi", i - 1, a0]) +
                              (v_PTDR_II_kappa[a0] * v_PTDR_II_psi * la["v_PTDR_II", i - 1, a0]) +
                              # Misdiagnosed - re-entry from MDR compartment
                              (la["v_NTDR_mdt", i - 1, a0] * v_NTDR_mdt_exit * (1 - nt_mdt_loss) * v_PTDR_II_psi) +
                              (la["v_PTDR_mdt", i - 1, a0] * v_PTDR_mdt_exit * (1 - pt_mdt_loss) * v_PTDR_II_psi) +
                              (v_NT_xi * v_PTDR_p[a0] * v_PTDR_f[a0] * la["v_NTDS_T", i - 1, a0] * (1 - nt_mdt_loss) * v_PTDR_II_psi) +
                              (v_PT_xi * v_PTDR_p[a0] * v_PTDR_f[a0] * la["v_PTDS_T", i - 1, a0] * (1 - pt_mdt_loss) * v_PTDR_II_psi) +
                             -(v_PTDR_T_IIpsi_tau * la["v_PTDR_T_IIpsi", i - 1, a0]) +
                             -(v_PTDR_T_u[a0] * la["v_PTDR_T_IIpsi", i - 1, a0])

## Treatment Failing - Originating from Infectious DR-TB
la["v_PTDR_T_IIphi", i, a1] <-  (surv * la["v_PTDR_T_IIphi", i - 1, a0]) +
                              (v_PTDR_II_kappa[a0] * v_PTDR_II_phi * la["v_PTDR_II", i - 1, a0]) +
                              # Misdiagnosed - re-entry from MDR compartment
                              (la["v_NTDR_mdt", i - 1, a0] * v_NTDR_mdt_exit * (1 - nt_mdt_loss) * v_PTDR_II_phi) +
                              (la["v_PTDR_mdt", i - 1, a0] * v_PTDR_mdt_exit * (1 - pt_mdt_loss) * v_PTDR_II_phi) +
                              (v_NT_xi * v_PTDR_p[a0] * v_PTDR_f[a0] * la["v_NTDS_T", i - 1, a0] * (1 - nt_mdt_loss) * v_PTDR_II_phi) +
                              (v_PT_xi * v_PTDR_p[a0] * v_PTDR_f[a0] * la["v_PTDS_T", i - 1, a0] * (1 - pt_mdt_loss) * v_PTDR_II_phi) +
                             -(v_PTDR_T_IIphi_tau * la["v_PTDR_T_IIphi", i - 1, a0]) +
                             -(v_PTDR_II_u[a0] * la["v_PTDR_T_IIphi", i - 1, a0])

## Treatment Succeeding - Originating from Non-Infectious DR-TB
la["v_PTDR_T_INpsi", i, a1] <-  (surv * la["v_PTDR_T_INpsi", i - 1, a0]) +
                              (v_PTDR_IN_kappa[a0] * v_PTDR_IN_psi * la["v_PTDR_IN", i - 1, a0]) +
                              # Misdiagnosed - re-entry from MDR compartment
                              (la["v_NTDR_nid", i - 1, a0] * v_NTDR_nid_exit * (1 - nt_nid_loss) * v_PTDR_IN_psi) +
                              (la["v_PTDR_nid", i - 1, a0] * v_PTDR_nid_exit * (1 - pt_nid_loss) * v_PTDR_IN_psi) +
                              (v_NT_xi * v_PTDR_p[a0] * (1 - v_PTDR_f[a0]) * la["v_NTDS_T", i - 1, a0] * (1 - nt_nid_loss) * v_PTDR_IN_psi) +
                              (v_PT_xi * v_PTDR_p[a0] * (1 - v_PTDR_f[a0]) * la["v_PTDS_T", i - 1, a0] * (1 - pt_nid_loss) * v_PTDR_IN_psi) +
                             -(v_PTDR_T_INpsi_tau * la["v_PTDR_T_INpsi", i - 1, a0]) +
                             -(v_PTDR_T_u[a0] * la["v_PTDR_T_INpsi", i - 1, a0])

## Treatment Failing - Originating from Non-Infectious DR-TB
la["v_PTDR_T_INphi", i, a1] <-  (surv * la["v_PTDR_T_INphi", i - 1, a0]) +
                              (v_PTDR_IN_kappa[a0] * v_PTDR_IN_phi * la["v_PTDR_IN", i - 1, a0]) +
                              # Misdiagnosed - re-entry from MDR compartment
                              (la["v_NTDR_nid", i - 1, a0] * v_NTDR_nid_exit * (1 - nt_nid_loss) * v_PTDR_IN_phi) +
                              (la["v_PTDR_nid", i - 1, a0] * v_PTDR_nid_exit * (1 - pt_nid_loss) * v_PTDR_IN_phi) +
                              (v_NT_xi * v_PTDR_p[a0] * (1 - v_PTDR_f[a0]) * la["v_NTDS_T", i - 1, a0] * (1 - nt_nid_loss) * v_PTDR_IN_phi) +
                              (v_PT_xi * v_PTDR_p[a0] * (1 - v_PTDR_f[a0]) * la["v_PTDS_T", i - 1, a0] * (1 - pt_nid_loss) * v_PTDR_IN_phi) +
                             -(v_PTDR_T_INphi_tau * la["v_PTDR_T_INphi", i - 1, a0]) +
                             -(v_PTDR_IN_u[a0] * la["v_PTDR_T_INphi", i - 1, a0])

## v_PTDR Resolved
la["v_PTDR_R", i, a1] <- (surv * la["v_PTDR_R", i - 1, a0]) +
                       ## Natural cure - Infectious
                       (v_PTDR_II_n[a0] * la["v_PTDR_II", i - 1, a0]) +
                       ## Natural cure - Non-infectious
                       (v_PTDR_IN_n[a0] * la["v_PTDR_IN", i - 1, a0]) +
                       ## v_NTDR Treatment Success - Infectious
                       (v_NTDR_T_IIpsi_tau * la["v_NTDR_T_IIpsi", i - 1, a0]) +
                       ## v_NTDR Treatment Success - Non-infectious
                       (v_NTDR_T_INpsi_tau * la["v_NTDR_T_INpsi", i - 1, a0]) +
                       ## v_PTDR Treatment Success - Infectious
                       (v_PTDR_T_INpsi_tau * la["v_PTDR_T_INpsi", i - 1, a0]) +
                       ## v_PTDR Treatment Success - Non-infectious
                       (v_PTDR_T_IIpsi_tau * la["v_PTDR_T_IIpsi", i - 1, a0]) +
                       ## Reinfection of Resolved - to DS
                      -(v_PTDS_lambda[i - 1, a0] * v_PTDS_p[a0] * v_PTDS_x * la["v_PTDR_R", i - 1, a0]) +
                       ## Reactivation of Resolved
                      -(v_PTDR_r[a0] * la["v_PTDR_R", i - 1, a0]) +
                      ## Reinfection of Resolved - to DR
                      -(v_PTDR_lambda[i - 1, a0] * v_PTDR_p[a0] * v_PTDR_x * la["v_PTDR_R", i - 1, a0])

## Misdiagnosed and Treated - Infectious
la["v_PTDR_mdt", i, a1] <-  (surv * la["v_PTDR_mdt", i - 1, a0]) +
                          # Misdiagnosed and treated
                          (v_PTDR_II_mdt[a0] * la["v_PTDR_II", i - 1, a0]) +
                          # Exit outwards
                          -(la["v_PTDR_mdt", i - 1, a0] * v_PTDR_mdt_exit) +
                          -(v_PTDR_mdt_u[a0] * la["v_PTDR_mdt", i - 1, a0])

## Misdiagnosed and Treated - Npn-Infectious
la["v_PTDR_nid", i, a1] <-   (surv * la["v_PTDR_nid", i - 1, a0]) +
                           # Misdiagnosed and treated
                          (v_PTDR_IN_nid[a0] * la["v_PTDR_IN", i - 1, a0]) +
                           # Exit outwards
                          -(la["v_PTDR_nid", i - 1, a0] * v_PTDR_nid_exit) +
                          -(v_PTDR_nid_u[a0] * la["v_PTDR_nid", i - 1, a0])
