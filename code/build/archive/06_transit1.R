
# DSTB deaths
ta["DSTB_deaths", i - 1, a5] <- (NTDS_II_u[a5] * la["NTDS_II", i - 1, a5]) +
                            (NTDS_IN_u[a5] * la["NTDS_IN", i - 1, a5]) +
                            (NTDS_T_u[a5] * la["NTDS_T", i - 1, a5]) +
                            (PTDS_II_u[a5] * la["PTDS_II", i - 1, a5]) +
                            (PTDS_IN_u[a5] * la["PTDS_IN", i - 1, a5]) +
                            (PTDS_T_u[a5] * la["PTDS_T", i - 1, a5])

# Backup calculation - DSTB deaths
dsma[current_step - 1, a5]  <-  (NTDS_II_u[a5] * la["NTDS_II", i - 1, a5]) +
                            (NTDS_IN_u[a5] * la["NTDS_IN", i - 1, a5]) +
                            (NTDS_T_u[a5] * la["NTDS_T", i - 1, a5]) +
                            (PTDS_II_u[a5] * la["PTDS_II", i - 1, a5]) +
                            (PTDS_IN_u[a5] * la["PTDS_IN", i - 1, a5]) +
                            (PTDS_T_u[a5] * la["PTDS_T", i - 1, a5])

# DRTB Deaths
ta["DRTB_deaths", i - 1, a5] <- (NTDR_II_u[a5] * la["NTDR_II", i - 1, a5]) +
                            (NTDS_IN_u[a5] * la["NTDR_IN", i - 1, a5]) +
                            (NTDR_T_u[a5] * la["NTDR_T_IIpsi", i - 1, a5]) +
                            (NTDR_II_u[a5] * la["NTDR_T_IIphi", i - 1, a5]) +
                            (NTDR_T_u[a5] * la["NTDR_T_INpsi", i - 1, a5]) +
                            (NTDR_IN_u[a5] * la["NTDR_T_INphi", i - 1, a5]) +
                            (PTDR_II_u[a5] * la["PTDR_II", i - 1, a5]) +
                            (PTDS_IN_u[a5] * la["PTDR_IN", i - 1, a5]) +
                            (PTDR_T_u[a5] * la["PTDR_T_IIpsi", i - 1, a5]) +
                            (PTDR_II_u[a5] * la["PTDR_T_IIphi", i - 1, a5]) +
                            (PTDR_T_u[a5] * la["PTDR_T_INpsi", i - 1, a5]) +
                            (PTDR_IN_u[a5] * la["PTDR_T_INphi", i - 1, a5]) +
                            (NTDR_mdt_u[a5] * la["NTDR_mdt", i - 1, a5])
                            (PTDR_mdt_u[a5] * la["PTDR_mdt", i - 1, a5])
                            (NTDR_nid_u[a5] * la["NTDR_nid", i - 1, a5])
                            (PTDR_nid_u[a5] * la["PTDR_nid", i - 1, a5])

# DSTB Incidence
ta["DSTB_inc", i - 1, a5] <-  ## Deconstructed
                          (NTDS_p[a5] * NTDS_f[a5] * NTDS_lambda[i - 1, a5] * la["S", i - 1, a5]) +
                          (NTDS_p[a5] * NTDS_f[a5] * NTDS_lambda[i - 1, a5] * NTDS_x * la["NTDR_L", i - 1, a5]) +
                          (NTDS_p[a5] * NTDS_f[a5] * NTDS_lambda[i - 1, a5] * NTDS_x * la["NTDS_L", i - 1, a5]) +
                          ## End deconstruct
                          ## Reactivation from Latent
                          (NTDS_v[a5] * NTDS_f[a5] * la["NTDS_L", i - 1, a5]) +
                          ## Reinfection of Resolved - DS
                          (NTDS_lambda[i - 1, a5] * NTDS_p[a5] * NTDS_f[a5] * NTDS_x * la["NTDS_R", i - 1, a5]) +
                          ## Reactivation from Resolved
                          (NTDS_f[a5] * NTDS_r[a5] * la["NTDS_R", i - 1, a5]) +
                          ## Reinfection of Resolved - DR
                          (NTDS_lambda[i - 1, a5] * NTDS_p[a5] * NTDS_f[a5] * NTDS_x * la["NTDR_R", i - 1, a5]) +
                          ## Deconstructed (NTDS)
                          (NTDS_p[a5] * (1 - NTDS_f[a5]) * NTDS_lambda[i - 1, a5] * la["S", i - 1, a5]) +
                          (NTDS_p[a5] * (1 - NTDS_f[a5]) * NTDS_lambda[i - 1, a5] * NTDS_x * la["NTDR_L", i - 1, a5]) +
                          (NTDS_p[a5] * (1 - NTDS_f[a5]) * NTDS_lambda[i - 1, a5] * NTDS_x * la["NTDS_L", i - 1, a5]) +
                          ## End deconstruct
                          ## Reactivation from Latent
                          (NTDS_v[a5] * (1 - NTDS_f[a5]) * la["NTDS_L", i - 1, a5]) +
                          ## Reinfection of Resolved - DS
                          (NTDS_lambda[i - 1, a5] * NTDS_p[a5] * (1 - NTDS_f[a5]) * NTDS_x * la["NTDS_R", i - 1,a5]) +
                          ## Reactivation from Resolved
                          ((1 - NTDS_f[a5]) * NTDS_r[a5] * la["NTDS_R", i-1,a5]) +
                          ## Reinfection of Resolved - DR
                          (NTDS_lambda[i - 1, a5] * NTDS_p[a5] * (1 - NTDS_f[a5]) * NTDS_x * la["NTDR_R", i - 1, a5]) +
                          ## Reinfection of Resolved - DS
                          (PTDS_lambda[i - 1, a5] * PTDS_p[a5] * PTDS_f[a5] * PTDS_x * la["PTDS_R", i - 1, a5]) +
                          ## Reactivation from Resolved
                          (PTDS_f[a5] * PTDS_r[a5] * la["PTDS_R", i - 1, a5]) +
                          ## Reinfection of Resolved - DR
                          (PTDS_lambda[i - 1, a5] * PTDS_p[a5] * PTDS_f[a5] * PTDS_x * la["PTDR_R", i - 1, a5]) +
                          ## Reinfection of Resolved - DS
                          (PTDS_lambda[i - 1, a5] * PTDS_p[a5] * (1 - PTDS_f[a5]) * PTDS_x * la["PTDS_R", i - 1, a5]) +
                          ## Reactivation from Resolved
                          ((1 - PTDS_f[a5]) * PTDS_r[a5] * la["PTDS_R", i - 1, a5]) +
                          ## Reinfection of Resolved - DR
                          (PTDS_lambda[i - 1, a5] * PTDS_p[a5] * (1 - PTDS_f[a5]) * PTDS_x * la["PTDR_R", i - 1, a5])

dsia[current_step - 1, a5]  <-  # Backup calculation
                            # New Infections of Susceptible and Latent
                            (NTDS_p[a5] * NTDS_lambda[i - 1, a5] * (la["S", i - 1, a5] + (NTDS_x * la["NTDR_L", i - 1, a5]) + (NTDS_x * la["NTDS_L", i - 1, a5]))) +
                            (NTDS_lambda[i - 1, a5] * NTDS_p[a5] * NTDS_x * (la["NTDS_R", i - 1, a5] + la["NTDR_R", i - 1, a5])) +
                            (NTDS_v[a5] * la["NTDS_L", i - 1, a5]) +
                            (NTDS_r[a5] * la["NTDS_R", i - 1, a5]) +
                            # PTDS
                            (PTDS_r[a5] * la["PTDS_R", i - 1, a5]) +
                            (PTDS_lambda[i - 1, a5] * PTDS_p[a5] * PTDS_x * (la["PTDR_R", i - 1, a5] + la["PTDS_R", i - 1, a5]))


dria[current_step - 1, a5] <-       # New infections of S/NTDS_L/NTDR_L
                                (NTDR_lambda[i - 1, a5] * NTDR_p[a5] * (la["S", i - 1, a5] + (NTDR_x * la["NTDR_L", i - 1, a5]) + (NTDR_x * la["NTDS_L", i - 1, a5]))) +
                                (NTDR_v[a5] * la["NTDR_L", i - 1, a5]) +
                                (NTDR_r[a5] * la["NTDR_R", i - 1, a5]) +
                                (NTDR_lambda[i - 1, a5] * NTDR_p[a5] * NTDR_x * (la["NTDR_R", i - 1, a5] + la["NTDS_R", i - 1, a5])) +
                                ## PT
                                (PTDR_v[a5] * la["PTDR_L", i - 1, a5]) + #ZERO COMPARTMENT
                                (PTDR_r[a5] * la["PTDR_R", i - 1, a5]) +
                                (PTDR_lambda[i - 1, a5] * PTDR_p[a5] * PTDR_x * (la["PTDS_R", i - 1, a5] + la["PTDR_R", i - 1, a5])) +
                                (NT_xi * PTDR_p[a5] * la["NTDS_T", i - 1, a5] * PTDR_f[a5] * nt_mdt_loss) +
                                (PT_xi * PTDR_p[a5] * la["PTDS_T", i - 1, a5] * PTDR_f[a5] * pt_mdt_loss) +
                                (NT_xi * PTDR_p[a5] * la["NTDS_T", i - 1, a5] * (1 - PTDR_f[a5]) * nt_nid_loss) +
                                (PT_xi * PTDR_p[a5] * la["PTDS_T", i - 1, a5] * (1 - PTDR_f[a5]) * pt_nid_loss)

# DRTB Incident Cases per timestep
ta["DRTB_inc", i - 1, a5]   <-  ## Deconstructed compound term - new infections in susceptble, NTDR-Latent and NTDS-Latent
                                (NTDR_p[a5] * NTDR_f[a5] * NTDR_lambda[i - 1, a5] * la["S", i - 1, a5]) +
                                (NTDR_p[a5] * NTDR_f[a5] * NTDR_lambda[i - 1, a5] * NTDR_x * la["NTDR_L", i - 1, a5]) +
                                (NTDR_p[a5] * NTDR_f[a5] * NTDR_lambda[i - 1, a5] * NTDR_x * la["NTDS_L", i - 1, a5]) +
                                ## End deconstruct
                                ## Reactivation from Latent
                                (NTDR_v[a5] * NTDR_f[a5] * la["NTDR_L", i - 1, a5]) +
                                ## Reinfection of Resolved - DR
                                (NTDR_lambda[i - 1, a5] * NTDR_p[a5] * NTDR_f[a5] * NTDR_x * la["NTDR_R", i - 1, a5]) +
                                ## Reactivation from Resolved
                                (NTDR_f[a5] * NTDR_r[a5] * la["NTDR_R", i - 1, a5]) +
                                ## Reinfection of Resolved - DS
                                (NTDR_lambda[i - 1, a5] * NTDR_p[a5] * NTDR_f[a5] * NTDR_x * la["NTDS_R", i - 1, a5]) +
                                # Deconstructed
                                (NTDR_p[a5] * (1 - NTDR_f[a5]) * NTDR_lambda[i - 1, a5] * la["S", i - 1, a5]) +
                                (NTDR_p[a5] * (1 - NTDR_f[a5]) * NTDR_lambda[i - 1, a5] * NTDR_x * la["NTDR_L", i - 1, a5]) +
                                (NTDR_p[a5] * (1 - NTDR_f[a5]) * NTDR_lambda[i - 1, a5] * NTDR_x * la["NTDS_L", i - 1, a5]) +
                                # End deconstruct
                                ## Reactivation from Latent
                                (NTDR_v[a5] * (1 - NTDR_f[a5]) * la["NTDR_L", i - 1, a5]) +
                                ## Reinfection of Resolved - DR
                                (NTDR_lambda[i - 1, a5] * NTDR_p[a5] * (1 - NTDR_f[a5]) * NTDR_x * la["NTDR_R", i - 1, a5]) +
                                ## Reactivation from Resolved
                                ((1 - NTDR_f[a5]) * NTDR_r[a5] * la["NTDR_R", i - 1, a5]) +
                                ## Reinfection of Resolved - DS
                                (NTDR_lambda[i - 1, a5] * NTDR_p[a5] * (1 - NTDR_f[a5]) * NTDR_x * la["NTDS_R", i - 1, a5]) +
                                ## Reactivation from Latent
                                (PTDR_v[a5] * PTDR_f[a5] * la["PTDR_L", i - 1, a5]) +
                                ## Reinfection of Resolved - DR
                                (PTDR_lambda[i - 1, a5] * PTDR_p[a5] * PTDR_f[a5] * PTDR_x * la["PTDR_R", i - 1, a5]) +
                                ## Reactivation from Resolved
                                (PTDR_f[a5] * (PTDR_r[a5] ) * la["PTDR_R", i - 1, a5]) +
                                ## Reinfection of Resolved - DS
                                (PTDR_lambda[i - 1, a5] * PTDR_p[a5] * PTDR_f[a5] * PTDR_x * la["PTDS_R", i - 1, a5]) +
                                ## Reactivation from Latent
                                (PTDR_v[a5] * (1 - PTDR_f[a5]) * la["PTDR_L", i - 1, a5]) +
                                ## Reinfection of Resolved - DS
                                (PTDR_lambda[i - 1, a5] * PTDR_p[a5] * (1 - PTDR_f[a5]) * PTDR_x * la["PTDS_R", i - 1, a5]) +
                                ## Reactivation from Resolved
                                ((1 - PTDR_f[a5]) * PTDR_r[a5] * la["PTDR_R", i - 1, a5]) +
                                ## Reinfection of Resolved - DR
                                (PTDR_lambda[i - 1, a5] * PTDR_p[a5] * (1 - PTDR_f[a5]) * PTDR_x * la["PTDR_R", i - 1, a5]) +
                                ## MDR Acquisition
                                (NT_xi * PTDR_p[a5] * la["NTDS_T", i - 1, a5] * PTDR_f[a5] * nt_mdt_loss) +
                                (PT_xi * PTDR_p[a5] * la["PTDS_T", i - 1, a5] * PTDR_f[a5] * pt_mdt_loss) +
                                (NT_xi * PTDR_p[a5] * la["NTDS_T", i - 1, a5] * (1 - PTDR_f[a5]) * nt_nid_loss) +
                                (PT_xi * PTDR_p[a5] * la["PTDS_T", i - 1, a5] * (1 - PTDR_f[a5]) * pt_nid_loss)

# ALL TB incidence backup calc
tbia[current_step - 1, a5]  <- dsia[current_step - 1, a5] + dria[current_step - 1, a5]

# ALL TB incidence
ta["TB_inc", i - 1, a5]     <- ta["DSTB_inc", i - 1, a5] + ta["DRTB_inc", i - 1, a5]

### Incidence by treatment history and bacteriologic status ###
## NTDR incidence
ta["DR_nt_inc", i - 1, a5]  <-  # New transmission by infection of susceptible and latent pools
                            (NTDR_p[a5] * NTDR_lambda[i - 1, a5] * la["S", i - 1, a5]) +
                            (NTDR_p[a5] * NTDR_lambda[i - 1, a5] * NTDR_x * la["NTDR_L", i - 1, a5]) +
                            (NTDR_p[a5] * NTDR_lambda[i - 1, a5] * NTDR_x * la["NTDS_L", i - 1, a5]) +
                            ## Reinfection of Resolved - DR
                            (NTDR_lambda[i - 1, a5] * NTDR_p[a5] * NTDR_x * la["NTDR_R", i - 1, a5]) +
                            ## Reinfection of Resolved - DS
                            (NTDR_lambda[i - 1, a5] * NTDR_p[a5] * NTDR_x * la["NTDS_R", i - 1, a5]) +
                            ## Reactivation from Latent
                            (NTDR_v[a5] * la["NTDR_L", i - 1, a5]) +
                            ## Reactivation from Resolved
                            (NTDR_r[a5] * la["NTDR_R", i - 1, a5])

# NTDR Bact+ Incidence
ta["DR_nt_inc_f", i - 1, a5] <- (ta["DR_nt_inc", i - 1, a5] * NTDR_f[a5]) +
                            (NTDR_omega * la["NTDR_IN", i - 1, a5])

# PTDR Incidence
ta["DR_pt_inc", i - 1, a5] <-   # Reinfection of Resolved - DR
                            (PTDR_lambda[i - 1, a5] * PTDR_p[a5] * PTDR_x * la["PTDR_R", i - 1, a5]) +
                            ## Reinfection of Resolved - DS
                            (PTDR_lambda[i - 1, a5] * PTDR_p[a5] * PTDR_x * la["PTDS_R", i - 1, a5]) +
                            ## Reactivation from Resolved
                            ((PTDR_r[a5]) * la["PTDR_R", i - 1, a5]) +
                            ## Reactivation from Latent
                            (PTDR_v[a5] * la["PTDR_L", i - 1, a5])
                            ## MDR Acquisition - not included.

# PTDR Bact+ Incidence
ta["DR_pt_inc_f", i - 1, a5] <- (ta["DR_pt_inc", i - 1, a5] * PTDS_f[a5]) +
                            (PTDR_omega * la["PTDR_IN", i - 1, a5]) +
                            # Bacteriologically positive MDR acquisition
                            (NT_xi * PTDR_p[a5] * la["NTDS_T", i - 1, a5] * PTDR_f[a5] * nt_mdt_loss) +
                            (PT_xi * PTDR_p[a5] * la["PTDS_T", i - 1, a5] * PTDR_f[a5] * pt_mdt_loss)

# NTDS Incidence
ta["DS_nt_inc", i - 1, a5] <- # Transmission by infection of suceptible and latent pools
                            (NTDS_p[a5] * NTDS_lambda[i - 1, a5] * la["S", i - 1, a5]) +
                            (NTDS_p[a5] * NTDS_lambda[i - 1, a5] * NTDS_x * la["NTDR_L", i - 1, a5]) +
                            (NTDS_p[a5] * NTDS_lambda[i - 1, a5] * NTDS_x * la["NTDS_L", i - 1, a5]) +
                            ## Reactivation from Latent
                            (NTDS_v[a5] * la["NTDS_L", i - 1, a5]) +
                            ## Reinfection of Resolved - DS
                            (NTDS_lambda[i - 1, a5] * NTDS_p[a5] * NTDS_x * la["NTDS_R", i - 1, a5]) +
                            ## Reactivation from Resolved
                            (NTDS_r[a5]  * la["NTDS_R", i - 1, a5]) +
                            ## Reinfection of Resolved - DR
                            (NTDS_lambda[i - 1, a5] * NTDS_p[a5] * NTDS_x * la["NTDR_R", i - 1, a5])

# NTDS Bact+ Incidence
ta["DS_nt_inc_f", i - 1, a5] <- (ta["DS_nt_inc", i - 1, a5] * NTDS_f[a5]) +
                            (NTDS_omega  * la["NTDS_IN", i - 1, a5])

# PTDS Incidence
ta["DS_pt_inc", i - 1, a5] <- ## Reinfection of Resolved - DS
                          (PTDS_lambda[i - 1, a5] * PTDS_p[a5] * PTDS_x * la["PTDS_R", i - 1, a5]) +
                          ## Reactivation from Resolved
                          (PTDS_r[a5] * la["PTDS_R", i - 1, a5]) +
                          ## Reinfection of Resolved - DR
                          (PTDS_lambda[i - 1, a5] * PTDS_p[a5] * PTDS_x * la["PTDR_R", i - 1, a5])

# PTDS Bact+ Incidence
ta["DS_pt_inc_f", i - 1, a5] <- (ta["DS_pt_inc", i - 1, a5] * PTDS_f[a5]) +
                            (PTDS_omega * la["PTDS_IN", i - 1, a5])

# AllTB Never Treated Incidence
ta["All_nt_inc", i - 1, a5] <- ta["DR_nt_inc", i - 1, a5] + ta["DS_nt_inc", i - 1, a5]

# AllTB Previously Treated Incidence
ta["All_pt_inc", i - 1, a5] <- ta["DR_pt_inc", i - 1, a5] + ta["DS_pt_inc", i - 1, a5]

# AllTB Never Treated Bact+ Incidence
ta["All_nt_inc_f", i - 1, a5] <- ta["DR_nt_inc_f", i - 1, a5] + ta["DS_nt_inc_f", i - 1, a5]

# AllTB Previously Treated Bact+ Incidence
ta["All_pt_inc_f", i - 1, a5] <- ta["DR_pt_inc_f", i - 1, a5] + ta["DS_pt_inc_f", i - 1, a5]

############################ DRTB Treatment Initiations [Notifications] #################################

##########################
########## TZ ############

# 'Fitting factor', 'ff' for diagnostics only.
# ff = 1 i.e. no effect during normal model functioning.
ff <- 1

##########################
##########################

# NTDR Treatment Initiations to any destination - infectious
ta["DR_nt_intx_f", i - 1, a5]  <- (NTDR_II_kappa[a5] * NTDR_II_psi * la["NTDR_II", i - 1, a5]) +
                                  (NTDR_II_kappa[a5] * NTDR_II_phi * la["NTDR_II", i - 1, a5]) +
                                  (NTDR_II_mdt[a5] * la["NTDR_II", i - 1, a5])

# PTDR Treatment Initiations to any destination - infectious
ta["DR_pt_intx_f", i - 1, a5]  <- (PTDR_II_kappa[a5] * PTDR_II_psi * la["PTDR_II", i - 1, a5]) +
                                  (PTDR_II_kappa[a5] * PTDR_II_phi * la["PTDR_II", i - 1, a5]) +
                                  (PTDR_II_mdt[a5] * la["PTDR_II", i - 1, a5]) +
                                  (ff * (la["NTDR_mdt", i - 1, a5] * NTDR_mdt_exit * (1 - nt_mdt_loss) * PTDR_II_psi)) +
                                  (ff * (la["NTDR_mdt", i - 1, a5] * NTDR_mdt_exit * (1 - nt_mdt_loss) * PTDR_II_phi)) +
                                  (ff * (la["PTDR_mdt", i - 1, a5] * PTDR_mdt_exit * (1 - pt_mdt_loss) * PTDR_II_psi)) +
                                  (ff * (la["PTDR_mdt", i - 1, a5] * PTDR_mdt_exit * (1 - pt_mdt_loss) * PTDR_II_phi)) +
                                  (ff * (NT_xi * PTDR_p[a5] * PTDR_f[a5] * la["NTDS_T", i - 1, a5] * (1 - nt_mdt_loss) * PTDR_II_psi)) +
                                  (ff * (PT_xi * PTDR_p[a5] * PTDR_f[a5] * la["PTDS_T", i - 1, a5] * (1 - pt_mdt_loss) * PTDR_II_psi)) +
                                  (ff * (NT_xi * PTDR_p[a5] * PTDR_f[a5] * la["NTDS_T", i - 1, a5] * (1 - nt_mdt_loss) * PTDR_II_phi)) +
                                  (ff * (PT_xi * PTDR_p[a5] * PTDR_f[a5] * la["PTDS_T", i - 1, a5] * (1 - pt_mdt_loss) * PTDR_II_phi))

# NTDS Treatment initiations to any destination - infectious
ta["DS_nt_intx_f", i - 1, a5]  <- (NTDS_II_kappa[a5] * la["NTDS_II", i - 1, a5])

# PTDS Treatment initiations to any destination - infectious
ta["DS_pt_intx_f", i - 1, a5]  <- (PTDS_II_kappa[a5] * la["PTDS_II", i - 1, a5])

# All-Never-Treated TB Treatment initiations to any destination - infectious
ta["All_nt_intx_f", i - 1, a5] <- (ta["DR_nt_intx_f", i - 1, a5] + ta["DS_nt_intx_f", i - 1, a5])

# All-Previously-Treated TB Treatment initiations to any destination - infectious
ta["All_pt_intx_f", i - 1, a5] <- (ta["DR_pt_intx_f", i - 1, a5] + ta["DS_pt_intx_f", i - 1, a5])

# Treatment Initations by _type of regimen_
## Treatment Initiations of DSTB Regimen
ta["DSTB_initRx", i - 1, a5] <- (NTDS_II_kappa[a5] * la["NTDS_II", i - 1, a5]) +
                                (NTDS_IN_kappa[a5] * la["NTDS_IN", i - 1, a5]) +
                                (PTDS_II_kappa[a5] * la["PTDS_II", i - 1, a5]) +
                                (PTDS_IN_kappa[a5] * la["PTDS_IN", i - 1, a5]) +
                                (NTDR_II_mdt[a5] * la["NTDR_II", i - 1, a5]) +
                                (NTDR_IN_nid[a5] * la["NTDR_IN", i - 1, a5]) +
                                (PTDR_II_mdt[a5] * la["PTDR_II", i - 1, a5]) +
                                (PTDR_IN_nid[a5] * la["PTDR_IN", i - 1, a5])

## Treatment Initiations of DRTB Regimen
ta["DRTB_initRx", i - 1, a5] <- (NTDR_II_kappa[a5] * NTDR_II_psi * la["NTDR_II", i - 1, a5]) +
                                (NTDR_II_kappa[a5] * NTDR_II_phi * la["NTDR_II", i - 1, a5]) +
                                (NTDR_IN_kappa[a5] * NTDR_IN_psi * la["NTDR_IN", i - 1, a5]) +
                                (NTDR_IN_kappa[a5] * NTDR_IN_phi * la["NTDR_IN", i - 1, a5]) +
                                (PTDR_II_kappa[a5] * PTDR_II_psi * la["PTDR_II", i - 1, a5]) +
                                (PTDR_II_kappa[a5] * PTDR_II_phi * la["PTDR_II", i - 1, a5]) +
                                (PTDR_IN_kappa[a5] * PTDR_IN_psi * la["PTDR_IN", i - 1, a5]) +
                                (PTDR_IN_kappa[a5] * PTDR_IN_phi * la["PTDR_IN", i - 1, a5]) +
                                (la["NTDR_mdt", i - 1, a5] * NTDR_mdt_exit * (1 - nt_mdt_loss) * PTDR_II_psi) +
                                (la["PTDR_mdt", i - 1, a5] * PTDR_mdt_exit * (1 - pt_mdt_loss) * PTDR_II_psi) +
                                (la["NTDR_mdt", i - 1, a5] * NTDR_mdt_exit * (1 - nt_mdt_loss) * PTDR_II_phi) +
                                (la["PTDR_mdt", i - 1, a5] * PTDR_mdt_exit * (1 - pt_mdt_loss) * PTDR_II_phi) +
                                (la["NTDR_nid", i - 1, a5] * NTDR_nid_exit * (1 - nt_nid_loss) * PTDR_IN_psi) +
                                (la["PTDR_nid", i - 1, a5] * PTDR_nid_exit * (1 - pt_nid_loss) * PTDR_IN_psi) +
                                (la["NTDR_nid", i - 1, a5] * NTDR_nid_exit * (1 - nt_mdt_loss) * PTDR_IN_phi) +
                                (la["PTDR_nid", i - 1, a5] * PTDR_nid_exit * (1 - pt_mdt_loss) * PTDR_IN_phi) +
                                (NT_xi * PTDR_p[a5] * PTDR_f[a5] * la["NTDS_T", i - 1, a5] * (1 - nt_mdt_loss) * PTDR_II_psi) +
                                (PT_xi * PTDR_p[a5] * PTDR_f[a5] * la["PTDS_T", i - 1, a5] * (1 - pt_mdt_loss) * PTDR_II_psi) +
                                (NT_xi * PTDR_p[a5] * PTDR_f[a5] * la["NTDS_T", i - 1, a5] * (1 - nt_mdt_loss) * PTDR_II_phi) +
                                (PT_xi * PTDR_p[a5] * PTDR_f[a5] * la["PTDS_T", i - 1, a5] * (1 - pt_mdt_loss) * PTDR_II_phi) +
                                (NT_xi * PTDR_p[a5] * (1 - PTDR_f[a5]) * la["NTDS_T", i - 1, a5] * (1 - nt_nid_loss) * PTDR_IN_psi) +
                                (PT_xi * PTDR_p[a5] * (1 - PTDR_f[a5]) * la["PTDS_T", i - 1, a5] * (1 - pt_nid_loss) * PTDR_IN_psi) +
                                (NT_xi * PTDR_p[a5] * (1 - PTDR_f[a5]) * la["NTDS_T", i - 1, a5] * (1 - nt_nid_loss) * PTDR_IN_phi) +
                                (PT_xi * PTDR_p[a5] * (1 - PTDR_f[a5]) * la["PTDS_T", i - 1, a5] * (1 - pt_nid_loss) * PTDR_IN_phi)

## Laboratory Confirmed Treatment Initiations of DRTB Regimen
ta["DRTB_initRxLab", i - 1, a5] <-(NTDR_II_kappa[a5] * NTDR_II_psi * la["NTDR_II", i - 1, a5] * nt_dst_p[current_year - year1  + 1]) +
                                  (NTDR_II_kappa[a5] * NTDR_II_phi * la["NTDR_II", i - 1, a5] * nt_dst_p[current_year - year1  + 1]) +
                                  (PTDR_II_kappa[a5] * PTDR_II_psi * la["PTDR_II", i - 1, a5] * pt_dst_p[current_year - year1  + 1]) +
                                  (PTDR_II_kappa[a5] * PTDR_II_phi * la["PTDR_II", i - 1, a5] * pt_dst_p[current_year - year1  + 1]) +
                                  (la["NTDR_mdt", i - 1, a5] * NTDR_mdt_exit * nt_dst_p[current_year - year1  + 1] * PTDR_II_psi) +
                                  (la["PTDR_mdt", i - 1, a5] * PTDR_mdt_exit * pt_dst_p[current_year - year1  + 1] * PTDR_II_psi) +
                                  (la["NTDR_mdt", i - 1, a5] * NTDR_mdt_exit * nt_dst_p[current_year - year1  + 1] * PTDR_II_phi) +
                                  (la["PTDR_mdt", i - 1, a5] * PTDR_mdt_exit * pt_dst_p[current_year - year1  + 1] * PTDR_II_phi) +
                                  (NT_xi * PTDR_p[a5] * PTDR_f[a5] * la["NTDS_T", i - 1, a5] * nt_dst_p[current_year - year1  + 1]) +
                                  (PT_xi * PTDR_p[a5] * PTDR_f[a5] * la["PTDS_T", i - 1, a5] * pt_dst_p[current_year - year1  + 1]) +
                                  (NT_xi * PTDR_p[a5] * PTDR_f[a5] * la["NTDS_T", i - 1, a5] * nt_dst_p[current_year - year1  + 1] * PTDR_II_phi) +
                                  (PT_xi * PTDR_p[a5] * PTDR_f[a5] * la["PTDS_T", i - 1, a5] * pt_dst_p[current_year - year1  + 1] * PTDR_II_phi)

# DSTB person-time on treatment - ** IN PERSON-MONTHS **
ta["DSTB_onRx", i - 1, a5] <- (12 * dt) *  (la["NTDS_T", i - 1, a5] +
                                        la["PTDS_T", i - 1, a5] +
                                        la["NTDR_mdt", i - 1, a5] +
                                        la["PTDR_mdt", i - 1, a5] +
                                        la["NTDR_nid", i - 1, a5] +
                                        la["PTDR_nid", i - 1, a5])

# DRTB person-time on treatment ** IN PERSON-MONTHS**
ta["DRTB_onRx", i - 1, a5] <- (12 * dt) *  (la["NTDR_T_IIphi", i - 1, a5] +
                                        la["NTDR_T_IIpsi", i - 1, a5] +
                                        la["NTDR_T_INphi", i - 1, a5] +
                                        la["NTDR_T_INpsi", i - 1, a5] +
                                        la["PTDR_T_IIphi", i - 1, a5] +
                                        la["PTDR_T_IIpsi", i - 1, a5] +
                                        la["PTDR_T_INphi", i - 1, a5] +
                                        la["PTDR_T_INpsi", i - 1, a5])

# DSTB Diagnostic Cost
ca["ds_dx", i - 1, a5] <- (ta["DSTB_initRx", i - 1, a5] * sdr[i - 1] * ds_dx_cost) * 1e+03

# DRTB Diagnostic Cost (Excl DST)
ca["dr_dx", i - 1, a5] <- (ta["DRTB_initRx", i - 1, a5] * sdr[i - 1] * dr_dx_cost) * 1000

# DST Cost
ca["dst", i - 1, a5]   <- (ta["DRTB_initRxLab", i - 1, a5] * dst_cost) * 1000

# DSTB Treatment Cost
ca["ds_tx", i - 1, a5] <- (ta["DSTB_onRx", i - 1, a5] * ds_tx_cost) * 1000

# DRTB Treatment Cost
ca["dr_tx", i - 1, a5] <- (ta["DRTB_onRx", i - 1, a5] * dr_tx_cost) * 1000

# Total treatment cost including fractional inflation for programme cost
ca["tbrx", i - 1, a5]  <- #
                        (ca["ds_dx", i - 1, a5] +
                        ca["dr_dx", i - 1, a5] +
                        ca["dst", i - 1, a5] +
                        ca["ds_tx", i - 1, a5] +
                        ca["dr_tx", i - 1, a5]) * (1 + prog_cost)
