
rm(surv)
surv <- 1 - (u[a0])

# Susceptibles
la["S", i, a1] <-   (surv * la["S", i - 1, a0]) +
                   -(NTDS_lambda[i - 1, a0] * la["S", i - 1, a0]) +
                   -(NTDR_lambda[i - 1, a0] * la["S", i - 1, a0])

rm(surv)
surv <- 1 - (u[a0])

## Never Treated Drug Sensitive TB
## NTDS-Latent
la["NTDS_L", i, a1] <- (surv * la["NTDS_L", i - 1, a0]) +
                       ((1 - NTDS_p[a0]) * NTDS_lambda[i - 1, a0] * la["S", i - 1, a0]) +
                       ((1 - NTDS_p[a0]) * NTDS_lambda[i - 1, a0] * NTDS_x * la["NTDR_L", i - 1, a0]) +
                      -(NTDS_v[a0]  * la["NTDS_L", i - 1, a0]) +
                      -(NTDS_x * NTDS_lambda[i - 1, a0] * NTDS_p[a0] * la["NTDS_L", i - 1, a0]) +
                      -(NTDR_x * NTDR_lambda[i - 1, a0] * la["NTDS_L", i - 1, a0])

## Infectious NTDS
la["NTDS_II", i, a1] <- (surv * la["NTDS_II", i - 1, a0]) +
                        ## Conversion of Non-infectious to Infectious
                        (NTDS_omega  * la["NTDS_IN", i - 1, a0]) +
                        ## Deconstructed
                        (NTDS_p[a0] * NTDS_f[a0] * NTDS_lambda[i - 1, a0] * la["S", i - 1, a0]) +
                        (NTDS_p[a0] * NTDS_f[a0] * NTDS_lambda[i - 1, a0] * NTDS_x * la["NTDR_L", i - 1, a0]) +
                        (NTDS_p[a0] * NTDS_f[a0] * NTDS_lambda[i - 1, a0] * NTDS_x * la["NTDS_L", i - 1, a0]) +
                        ## End deconstruct
                        ## Reactivation from Latent
                        (NTDS_v[a0]  * NTDS_f[a0] * la["NTDS_L", i - 1, a0]) +
                        ## Reinfection of Resolved - DS
                        (NTDS_lambda[i - 1, a0] * NTDS_p[a0] * NTDS_f[a0] * NTDS_x * la["NTDS_R", i - 1, a0]) +
                        ## Reactivation from Resolved
                        (NTDS_f[a0] * NTDS_r[a0]  * la["NTDS_R", i - 1, a0]) +
                        ## Reinfection of Resolved - DR
                        (NTDS_lambda[i - 1, a0] * NTDS_p[a0] * NTDS_f[a0] * NTDS_x * la["NTDR_R", i - 1, a0]) +
                        ## Natural Cure
                       -(NTDS_II_n[a0]  * la["NTDS_II", i - 1, a0]) +
                        ## Case detection
                       -(NTDS_II_kappa[a0] * la["NTDS_II", i - 1, a0]) +
                        ## TB Death
                       -(NTDS_II_u[a0] * la["NTDS_II", i - 1, a0])

## Non-Infectious NTDS
la["NTDS_IN", i,a1] <- (surv * la["NTDS_IN", i-1,a0]) +
                       ## Deconstructed
                       (NTDS_p[a0] * (1 - NTDS_f[a0]) * NTDS_lambda[i - 1, a0] * la["S", i - 1, a0]) +
                       (NTDS_p[a0] * (1 - NTDS_f[a0]) * NTDS_lambda[i - 1, a0] * NTDS_x * la["NTDR_L", i - 1, a0]) +
                       (NTDS_p[a0] * (1 - NTDS_f[a0]) * NTDS_lambda[i - 1, a0] * NTDS_x * la["NTDS_L", i - 1, a0]) +
                       ## End deconstruct
                       ## Reactivation from Latent
                       (NTDS_v[a0]  * (1 - NTDS_f[a0]) * la["NTDS_L", i - 1, a0]) +
                       ## Reinfection of Resolved - DS
                       (NTDS_lambda[i - 1, a0] * NTDS_p[a0] * (1 - NTDS_f[a0]) * NTDS_x * la["NTDS_R", i - 1,a0]) +
                       ## Reactivation from Resolved
                       ((1 - NTDS_f[a0]) * NTDS_r[a0]  * la["NTDS_R", i-1,a0]) +
                       ## Reinfection of Resolved - DR
                       (NTDS_lambda[i - 1, a0] * NTDS_p[a0] * (1 - NTDS_f[a0]) * NTDS_x * la["NTDR_R", i - 1, a0]) +
                       ## Natural Cure
                      -(NTDS_IN_n[a0]  * la["NTDS_IN", i - 1, a0]) +
                       ## Case Detection
                      -(NTDS_IN_kappa[a0] * la["NTDS_IN", i - 1, a0]) +
                       ## TB Death
                      -(NTDS_IN_u[a0] * la["NTDS_IN", i - 1, a0]) +
                       ## Conversion from Non-infectious to Infectious
                      -(NTDS_omega  * la["NTDS_IN", i - 1, a0])

## NTDS on treatment
la["NTDS_T", i,a1]  <- (surv * la["NTDS_T", i-1,a0]) +
                       ## Treatment initiation from Infectious
                       (NTDS_II_kappa[a0]*la["NTDS_II", i-1,a0]) +
                       ## Treatment initiation from Non-infectious
                       (NTDS_IN_kappa[a0]*la["NTDS_IN", i-1,a0]) +
                       ## Resistance acquisition
                      -(NT_xi * la["NTDS_T", i-1,a0]) +
                       ## Treatment failures
                      -(NTDS_phi * la["NTDS_T", i-1,a0]) +
                       ## Treatment successes
                      -(NTDS_psi * la["NTDS_T", i-1,a0]) +
                       ## Death on treatment
                      -(NTDS_T_u[a0] * la["NTDS_T", i-1,a0])

# NTDS Resolved
la["NTDS_R", i, a1] <- (surv * la["NTDS_R", i - 1, a0]) +
                       ## Natural cure Infectious
                       (NTDS_II_n[a0]  * la["NTDS_II", i - 1, a0]) +
                       ## Natural cure Non-infectious
                       (NTDS_IN_n[a0]  * la["NTDS_IN", i - 1, a0]) +
                       ## Reinfection of Resolved - DS
                      -(NTDS_lambda[i - 1, a0] * NTDS_p[a0] * NTDS_x * la["NTDS_R", i - 1, a0]) +
                       ## Reactivation of Resolved
                      -(NTDS_r[a0]  * la["NTDS_R", i - 1, a0]) +
                       ## Reinfection of Resolved - DR
                      -(NTDR_lambda[i - 1, a0] * NTDR_p[a0] * NTDR_x * la["NTDS_R", i - 1, a0])

rm(surv4, surv3)
surv4 <- 1 - (u[a4] )
surv3 <- 1 - (u[a3] )

# Susceptibles
la["S", i, a4] <-      (surv4 * la["S", i - 1, a4]) +
                      -(NTDS_lambda[i - 1, a4] * la["S", i - 1, a4]) +
                      -(NTDR_lambda[i - 1, a4] * la["S", i - 1, a4]) +
# Accumulate
                       (surv3 * la["S", i - 1, a3]) +
                      -(NTDS_lambda[i - 1, a3] * la["S", i - 1, a3]) +
                      -(NTDR_lambda[i - 1, a3] * la["S", i - 1, a3])

# Never Treated Drug Sensitive TB
# NTDS-latent
la["NTDS_L", i, a4] <- (surv4 * la["NTDS_L", i - 1, a4]) +
                       ((1 - NTDS_p[a4]) * NTDS_lambda[i - 1, a4] * la["S", i - 1, a4]) +
                       ((1 - NTDS_p[a4]) * NTDS_lambda[i - 1, a4] * NTDS_x * la["NTDR_L", i - 1, a4]) +
                      -(NTDS_v[a4] * la["NTDS_L", i - 1, a4]) +
                      -(NTDS_x * NTDS_lambda[i - 1, a4] * NTDS_p[a4] * la["NTDS_L", i - 1, a4]) +
                      -(NTDR_x * NTDR_lambda[i - 1, a4] * la["NTDS_L", i - 1, a4]) +
                      #Accumulate
                       (surv3 * la["NTDS_L", i - 1, a3]) +
                       ((1 - NTDS_p[a3]) * NTDS_lambda[i - 1, a3] * la["S", i - 1, a3]) +
                       ((1 - NTDS_p[a3]) * NTDS_lambda[i - 1, a3] * NTDS_x * la["NTDR_L", i - 1, a3]) +
                      -(NTDS_v[a3] * la["NTDS_L", i - 1, a3]) +
                      -(NTDS_x * NTDS_lambda[i - 1, a3] * NTDS_p[a3] * la["NTDS_L", i - 1, a3]) +
                      -(NTDR_x * NTDR_lambda[i - 1, a3] * la["NTDS_L", i - 1, a3])

### Infectious DS-TB
la["NTDS_II", i, a4] <- (surv4 * la["NTDS_II", i - 1, a4]) +
                        (NTDS_omega * la["NTDS_IN", i - 1, a4]) +
                        (NTDS_p[a4] * NTDS_f[a4] * NTDS_lambda[i - 1, a4] * la["S", i - 1, a4]) +
                        (NTDS_p[a4] * NTDS_f[a4] * NTDS_lambda[i - 1, a4] * NTDS_x * la["NTDR_L", i - 1, a4]) +
                        (NTDS_p[a4] * NTDS_f[a4] * NTDS_lambda[i - 1, a4] * NTDS_x * la["NTDS_L", i - 1, a4]) +
                        (NTDS_v[a4] * NTDS_f[a4] * la["NTDS_L", i - 1, a4]) +
                        (NTDS_lambda[i - 1, a4] * NTDS_p[a4] * NTDS_f[a4] * NTDS_x * la["NTDS_R", i - 1, a4]) +
                        (NTDS_f[a4] * NTDS_r[a4] * la["NTDS_R", i - 1, a4]) +
                        (NTDS_lambda[i - 1, a4] * NTDS_p[a4] * NTDS_f[a4] * NTDS_x * la["NTDR_R", i - 1, a4]) +
                       -(NTDS_II_n[a4] * la["NTDS_II", i - 1, a4]) +
                       -(NTDS_II_kappa[a4] * la["NTDS_II", i - 1, a4]) +
                       -(NTDS_II_u[a4] * la["NTDS_II", i - 1, a4]) +
                       # Accumulate
                        (surv3 * la["NTDS_II", i - 1, a3]) +
                       (NTDS_omega * la["NTDS_IN", i - 1, a3]) +
                        (NTDS_p[a3] * NTDS_f[a3] * NTDS_lambda[i - 1, a3] * la["S", i - 1, a3]) +
                        (NTDS_p[a3] * NTDS_f[a3] * NTDS_lambda[i - 1, a3] * NTDS_x * la["NTDR_L", i - 1, a3]) +
                        (NTDS_p[a3] * NTDS_f[a3] * NTDS_lambda[i - 1, a3] * NTDS_x * la["NTDS_L", i - 1, a3]) +
                        (NTDS_v[a3] * NTDS_f[a3] * la["NTDS_L", i - 1, a3]) +
                        (NTDS_lambda[i - 1, a3] * NTDS_p[a3] * NTDS_f[a3] * NTDS_x * la["NTDS_R", i - 1, a3]) +
                        (NTDS_f[a3] * NTDS_r[a3] * la["NTDS_R", i - 1, a3]) +
                        (NTDS_lambda[i - 1, a3] * NTDS_p[a3] * NTDS_f[a3] * NTDS_x * la["NTDR_R", i - 1, a3]) +
                       -(NTDS_II_n[a3] * la["NTDS_II", i - 1, a3]) +
                       -(NTDS_II_kappa[a3] * la["NTDS_II", i - 1, a3]) +
                       -(NTDS_II_u[a3] * la["NTDS_II", i - 1, a3])

## Non-infectious DS-TB
la["NTDS_IN", i,a4] <-  (surv4 * la["NTDS_IN", i-1,a4]) +
                        (NTDS_p[a4] * (1 - NTDS_f[a4]) * NTDS_lambda[i - 1, a4] * la["S", i - 1, a4]) +
                        (NTDS_p[a4] * (1 - NTDS_f[a4]) * NTDS_lambda[i - 1, a4] * NTDS_x * la["NTDR_L", i - 1, a4]) +
                        (NTDS_p[a4] * (1 - NTDS_f[a4]) * NTDS_lambda[i - 1, a4] * NTDS_x * la["NTDS_L", i - 1, a4]) +
                        (NTDS_v[a4] * (1 - NTDS_f[a4]) * la["NTDS_L", i - 1, a4]) +
                        (NTDS_lambda[i - 1, a4] * NTDS_p[a4] * (1 - NTDS_f[a4]) * NTDS_x * la["NTDS_R", i-1,a4]) +
                        ((1 - NTDS_f[a4]) * NTDS_r[a4] * la["NTDS_R", i-1,a4]) +
                        (NTDS_lambda[i - 1, a4] * NTDS_p[a4] * (1 - NTDS_f[a4]) * NTDS_x * la["NTDR_R", i - 1, a4]) +
                       -(NTDS_IN_n[a4] * la["NTDS_IN", i - 1, a4]) +
                       -(NTDS_IN_kappa[a4] * la["NTDS_IN", i - 1, a4]) +
                       -(NTDS_IN_u[a4] * la["NTDS_IN", i - 1, a4]) +
                       -(NTDS_omega * la["NTDS_IN", i - 1, a4]) +
                       # Accumulate
                        (surv3 * la["NTDS_IN", i-1,a3]) +
                        (NTDS_p[a3] * (1 - NTDS_f[a3]) * NTDS_lambda[i - 1, a3] * la["S", i - 1, a3]) +
                        (NTDS_p[a3] * (1 - NTDS_f[a3]) * NTDS_lambda[i - 1, a3] * NTDS_x * la["NTDR_L", i - 1, a3]) +
                        (NTDS_p[a3] * (1 - NTDS_f[a3]) * NTDS_lambda[i - 1, a3] * NTDS_x * la["NTDS_L", i - 1, a3]) +
                        (NTDS_v[a3] * (1 - NTDS_f[a3]) * la["NTDS_L", i - 1, a3]) +
                        (NTDS_lambda[i - 1, a3] * NTDS_p[a3] * (1 - NTDS_f[a3]) * NTDS_x * la["NTDS_R", i-1,a3]) +
                        ((1 - NTDS_f[a3]) * NTDS_r[a3] * la["NTDS_R", i-1,a3]) +
                        (NTDS_lambda[i - 1, a3] * NTDS_p[a3] * (1 - NTDS_f[a3]) * NTDS_x * la["NTDR_R", i - 1, a3]) +
                       -(NTDS_IN_n[a3] * la["NTDS_IN", i - 1, a3]) +
                       -(NTDS_IN_kappa[a3] * la["NTDS_IN", i - 1, a3]) +
                       -(NTDS_IN_u[a3] * la["NTDS_IN", i - 1, a3]) +
                       -(NTDS_omega * la["NTDS_IN", i - 1, a3])

## NTDS On Treatment
la["NTDS_T", i,a4]  <- (surv4 * la["NTDS_T", i-1,a4]) +
                      (NTDS_II_kappa[a4]*la["NTDS_II", i-1,a4]) +
                      (NTDS_IN_kappa[a4]*la["NTDS_IN", i-1,a4]) -
                      (NT_xi * la["NTDS_T", i-1,a4]) -
                      (NTDS_phi * la["NTDS_T", i-1,a4]) -
                      (NTDS_psi * la["NTDS_T", i-1,a4]) -
                      (NTDS_T_u[a4] * la["NTDS_T", i-1,a4]) +
                      #Accumulate
                      (surv3 * la["NTDS_T", i-1,a3]) +
                      (NTDS_II_kappa[a3]*la["NTDS_II", i-1,a3]) +
                      (NTDS_IN_kappa[a3]*la["NTDS_IN", i-1,a3]) -
                      (NT_xi * la["NTDS_T", i-1,a3]) -
                      (NTDS_phi * la["NTDS_T", i-1,a3]) -
                      (NTDS_psi * la["NTDS_T", i-1,a3]) -
                      (NTDS_T_u[a3] * la["NTDS_T", i-1,a3])

# NTDS Resolved
la["NTDS_R", i, a4] <-  (surv4 * la["NTDS_R", i - 1, a4]) +
                        (NTDS_II_n[a4] * la["NTDS_II", i - 1, a4]) +
                        (NTDS_IN_n[a4] * la["NTDS_IN", i - 1, a4]) +
                       -(NTDS_lambda[i - 1, a4] * NTDS_p[a4] * NTDS_x * la["NTDS_R", i - 1, a4]) +
                       -(NTDS_r[a4] * la["NTDS_R", i - 1, a4]) +
                       -(NTDR_lambda[i - 1, a4] * NTDR_p[a4] * NTDR_x * la["NTDS_R", i - 1, a4]) +
                        # Accumulate
                        (surv3 * la["NTDS_R", i - 1, a3]) +
                        (NTDS_II_n[a3] * la["NTDS_II", i - 1, a3]) +
                        (NTDS_IN_n[a3] * la["NTDS_IN", i - 1, a3]) +
                       -(NTDS_lambda[i - 1, a3] * NTDS_p[a3] * NTDS_x * la["NTDS_R", i - 1, a3]) +
                       -(NTDS_r[a3] * la["NTDS_R", i - 1, a3]) +
                       -(NTDR_lambda[i - 1, a3] * NTDR_p[a3] * NTDR_x * la["NTDS_R", i - 1, a3])

rm(surv)
surv <- 1 - (u[a0] )

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

rm(surv4, surv3)
surv4 <- 1 - (u[a4])
surv3 <- 1 - (u[a3])

## Never Treated Drug Resistant

## NTDR-Latent
la["NTDR_L", i, a4] <-  (surv4 * la["NTDR_L", i - 1, a4]) +
                        (surv3 * la["NTDR_L", i - 1, a3]) +
                        ## New Infections from Susceptible
                        ((1 - NTDR_p[a4]) * NTDR_lambda[i - 1, a4] * la["S", i - 1, a4]) +
                        ((1 - NTDR_p[a3]) * NTDR_lambda[i - 1, a3] * la["S", i - 1, a3]) +
                        ## New Infections from NTDS-Latent
                        ((1 - NTDR_p[a4]) * NTDR_lambda[i - 1, a4] * NTDR_x * la["NTDS_L", i - 1, a4]) +
                        ((1 - NTDR_p[a3]) * NTDR_lambda[i - 1, a3] * NTDR_x * la["NTDS_L", i - 1, a3]) +
                        ## Reactivation from Latent
                       -(NTDR_v[a4] * la["NTDR_L", i - 1, a4]) +
                       -(NTDR_v[a3] * la["NTDR_L", i - 1, a3]) +
                        ## DRTB Reinfection of NTDR-Latent
                       -(NTDR_x * NTDR_lambda[i - 1, a4] * NTDR_p[a4] * la["NTDR_L", i - 1, a4]) +
                       -(NTDR_x * NTDR_lambda[i - 1, a3] * NTDR_p[a3] * la["NTDR_L", i - 1, a3]) +
                        ## DSTB Reinfection of NTDR-Latent
                       -(NTDR_x * NTDS_lambda[i - 1, a4] * la["NTDR_L", i - 1, a4]) +
                       -(NTDR_x * NTDS_lambda[i - 1, a3] * la["NTDR_L", i - 1, a3])

### Infectious NTDR
la["NTDR_II", i,a4] <-  (surv4 * la["NTDR_II", i - 1, a4]) +
                        (surv3 * la["NTDR_II", i - 1, a3]) +
                        ## Conversion of non-infectious to infectious
                        (NTDR_omega * la["NTDR_IN", i - 1, a4]) +
                        (NTDR_omega * la["NTDR_IN", i - 1, a3]) +
                        ## Deconstructed
                        (NTDR_p[a4] * NTDR_f[a4] * NTDR_lambda[i - 1, a4] * la["S", i - 1, a4]) +
                        (NTDR_p[a4] * NTDR_f[a4] * NTDR_lambda[i - 1, a4] * NTDR_x * la["NTDR_L", i - 1, a4]) +
                        (NTDR_p[a4] * NTDR_f[a4] * NTDR_lambda[i - 1, a4] * NTDR_x * la["NTDS_L", i - 1, a4]) +
                        (NTDR_p[a3] * NTDR_f[a3] * NTDR_lambda[i - 1, a3] * la["S", i - 1, a3]) +
                        (NTDR_p[a3] * NTDR_f[a3] * NTDR_lambda[i - 1, a3] * NTDR_x * la["NTDR_L", i - 1, a3]) +
                        (NTDR_p[a3] * NTDR_f[a3] * NTDR_lambda[i - 1, a3] * NTDR_x * la["NTDS_L", i - 1, a3]) +
                        ## End deconstruct
                        ## Reactivation from Latent
                        (NTDR_v[a4] * NTDR_f[a4] * la["NTDR_L", i - 1, a4]) +
                        (NTDR_v[a3] * NTDR_f[a3] * la["NTDR_L", i - 1, a3]) +
                        ## Reinfection of Resolved - DR
                        (NTDR_lambda[i - 1, a4] * NTDR_p[a4] * NTDR_f[a4] * NTDR_x * la["NTDR_R", i - 1, a4]) +
                        (NTDR_lambda[i - 1, a3] * NTDR_p[a3] * NTDR_f[a3] * NTDR_x * la["NTDR_R", i - 1, a3]) +
                        ## Reactivation from Resolved
                        (NTDR_f[a4] * NTDR_r[a4] * la["NTDR_R", i - 1, a4]) +
                        (NTDR_f[a3] * NTDR_r[a3] * la["NTDR_R", i - 1, a3]) +
                        ## Reinfection of Resolved - DS
                        (NTDR_lambda[i - 1, a4] * NTDR_p[a4] * NTDR_f[a4] * NTDR_x * la["NTDS_R", i - 1, a4]) +
                        (NTDR_lambda[i - 1, a3] * NTDR_p[a3] * NTDR_f[a3] * NTDR_x * la["NTDS_R", i - 1, a3]) +
                        ## Natural cure
                       -(NTDR_II_n[a4] * la["NTDR_II", i - 1, a4]) +
                       -(NTDR_II_n[a3] * la["NTDR_II", i - 1, a3]) +
                        ## Case detection
                       -(NTDR_II_kappa[a4] * la["NTDR_II", i - 1, a4]) +
                       -(NTDR_II_kappa[a3] * la["NTDR_II", i - 1, a3]) +
                        ## Misdiagnosis and Treatment
                       -(NTDR_II_mdt[a4] * la["NTDR_II", i - 1, a4]) +
                       -(NTDR_II_mdt[a3] * la["NTDR_II", i - 1, a3]) +
                        ## TB death
                       -(NTDR_II_u[a4] * la["NTDR_II", i - 1, a4]) +
                       -(NTDR_II_u[a3] * la["NTDR_II", i - 1, a3])

### Non-Infectious NTDR
la["NTDR_IN", i,a4] <-  (surv4 * la["NTDR_IN", i - 1, a4]) +
                        (surv3 * la["NTDR_IN", i - 1, a3]) +
                        # Deconstructed
                        (NTDR_p[a4] * (1 - NTDR_f[a4]) * NTDR_lambda[i - 1, a4] * la["S", i - 1, a4]) +
                        (NTDR_p[a4] * (1 - NTDR_f[a4]) * NTDR_lambda[i - 1, a4] * NTDR_x * la["NTDR_L", i - 1, a4]) +
                        (NTDR_p[a4] * (1 - NTDR_f[a4]) * NTDR_lambda[i - 1, a4] * NTDR_x * la["NTDS_L", i - 1, a4]) +
                        (NTDR_p[a3] * (1 - NTDR_f[a3]) * NTDR_lambda[i - 1, a3] * la["S", i - 1, a3]) +
                        (NTDR_p[a3] * (1 - NTDR_f[a3]) * NTDR_lambda[i - 1, a3] * NTDR_x * la["NTDR_L", i - 1, a3]) +
                        (NTDR_p[a3] * (1 - NTDR_f[a3]) * NTDR_lambda[i - 1, a3] * NTDR_x * la["NTDS_L", i - 1, a3]) +
                        # End deconstruct
                        ## Reactivation from Latent
                        (NTDR_v[a4] * (1 - NTDR_f[a4]) * la["NTDR_L", i - 1, a4]) +
                        (NTDR_v[a3] * (1 - NTDR_f[a3]) * la["NTDR_L", i - 1, a3]) +
                        ## Reinfection of Resolved - DR
                        (NTDR_lambda[i - 1, a4] * NTDR_p[a4] * (1 - NTDR_f[a4]) * NTDR_x * la["NTDR_R", i - 1, a4]) +
                        (NTDR_lambda[i - 1, a3] * NTDR_p[a3] * (1 - NTDR_f[a3]) * NTDR_x * la["NTDR_R", i - 1, a3]) +
                        ## Reactivation from Resolved
                        ((1 - NTDR_f[a4]) * NTDR_r[a4] * la["NTDR_R", i - 1, a4]) +
                        ((1 - NTDR_f[a3]) * NTDR_r[a3] * la["NTDR_R", i - 1, a3]) +
                        ## Reinfection of Resolved - DS
                        (NTDR_lambda[i - 1, a4] * NTDR_p[a4] * (1 - NTDR_f[a4]) * NTDR_x * la["NTDS_R", i - 1, a4]) +
                        (NTDR_lambda[i - 1, a3] * NTDR_p[a3] * (1 - NTDR_f[a3]) * NTDR_x * la["NTDS_R", i - 1, a3]) +
                        ## Natural cure
                       -(NTDR_IN_n[a4] * la["NTDR_IN", i - 1, a4]) +
                       -(NTDR_IN_n[a3] * la["NTDR_IN", i - 1, a3]) +
                        ## Case detection
                       -(NTDR_IN_kappa[a4] * la["NTDR_IN", i - 1, a4]) +
                       -(NTDR_IN_kappa[a3] * la["NTDR_IN", i - 1, a3]) +
                        ## Misdiagnosis and Treatment
                       -(NTDR_IN_nid[a4] * la["NTDR_IN", i - 1, a4]) +
                       -(NTDR_IN_nid[a3] * la["NTDR_IN", i - 1, a3]) +
                        ## TB death
                       -(NTDR_IN_u[a4] * la["NTDR_IN", i - 1, a4]) +
                       -(NTDR_IN_u[a3] * la["NTDR_IN", i - 1, a3]) +
                        ## Conversion from non-infectious to infectious
                       -(NTDR_omega * la["NTDR_IN", i - 1, a4]) +
                       -(NTDR_omega * la["NTDR_IN", i - 1, a3])

## Treatment - NTDR
## Treatment Succeeding - Originating from Infectious DR-TB
la["NTDR_T_IIpsi", i, a4] <- (surv4 * la["NTDR_T_IIpsi", i - 1, a4]) +
                             (surv3 * la["NTDR_T_IIpsi", i - 1, a3]) +
                             (NTDR_II_kappa[a4] * NTDR_II_psi * la["NTDR_II", i - 1, a4]) +
                             (NTDR_II_kappa[a3] * NTDR_II_psi * la["NTDR_II", i - 1, a3]) +
                            -(NTDR_T_IIpsi_tau * la["NTDR_T_IIpsi", i - 1, a4]) +
                            -(NTDR_T_IIpsi_tau * la["NTDR_T_IIpsi", i - 1, a3]) +
                            -(NTDR_T_u[a4] * la["NTDR_T_IIpsi", i - 1, a4]) +
                            -(NTDR_T_u[a3] * la["NTDR_T_IIpsi", i - 1, a3])

## Treatment Failing - Originating from Infectious DR-TB
la["NTDR_T_IIphi", i, a4] <- (surv4 * la["NTDR_T_IIphi", i - 1, a4]) +
                             (surv3 * la["NTDR_T_IIphi", i - 1, a3]) +
                             (NTDR_II_kappa[a4] * NTDR_II_phi * la["NTDR_II", i - 1, a4]) +
                             (NTDR_II_kappa[a3] * NTDR_II_phi * la["NTDR_II", i - 1, a3]) +
                            -(NTDR_T_IIphi_tau * la["NTDR_T_IIphi", i - 1, a4]) +
                            -(NTDR_T_IIphi_tau * la["NTDR_T_IIphi", i - 1, a3]) +
                            -(NTDR_II_u[a4] * la["NTDR_T_IIphi", i - 1, a4]) +
                            -(NTDR_II_u[a3] * la["NTDR_T_IIphi", i - 1, a3])

## Treatment Succeeding - Originating from Non-Infectious DR-TB
la["NTDR_T_INpsi", i, a4] <- (surv4 * la["NTDR_T_INpsi", i - 1, a4]) +
                             (surv3 * la["NTDR_T_INpsi", i - 1, a3]) +
                             (NTDR_IN_kappa[a4] * NTDR_IN_psi * la["NTDR_IN", i - 1, a4]) +
                             (NTDR_IN_kappa[a3] * NTDR_IN_psi * la["NTDR_IN", i - 1, a3]) +
                            -(NTDR_T_INpsi_tau * la["NTDR_T_INpsi", i - 1, a4]) +
                            -(NTDR_T_INpsi_tau * la["NTDR_T_INpsi", i - 1, a3]) +
                            -(NTDR_T_u[a4] * la["NTDR_T_INpsi", i - 1, a4]) +
                            -(NTDR_T_u[a3] * la["NTDR_T_INpsi", i - 1, a3])

## Treatment Failing - Originating from Non-Infectious DR-TB
la["NTDR_T_INphi", i, a4] <- (surv4 * la["NTDR_T_INphi", i - 1, a4]) +
                             (surv3 * la["NTDR_T_INphi", i - 1, a3]) +
                             (NTDR_IN_kappa[a4] * NTDR_IN_phi * la["NTDR_IN", i - 1, a4]) +
                             (NTDR_IN_kappa[a3] * NTDR_IN_phi * la["NTDR_IN", i - 1, a3]) +
                            -(NTDR_T_INphi_tau * la["NTDR_T_INphi", i - 1, a4]) +
                            -(NTDR_T_INphi_tau * la["NTDR_T_INphi", i - 1, a3]) +
                            -(NTDR_IN_u[a4] * la["NTDR_T_INphi", i - 1, a4]) +
                            -(NTDR_IN_u[a3] * la["NTDR_T_INphi", i - 1, a3])

## Resolved - DR TB
la["NTDR_R", i, a4] <- (surv4 * la["NTDR_R", i - 1, a4]) +
                       (surv3 * la["NTDR_R", i - 1, a3]) +
                       ## Natural cure - Infectious
                       (NTDR_II_n[a4] * la["NTDR_II", i - 1, a4]) +
                       (NTDR_II_n[a3] * la["NTDR_II", i - 1, a3]) +
                       ## Natural cure - Non-infectious
                       (NTDR_IN_n[a4] * la["NTDR_IN", i - 1, a4]) +
                       (NTDR_IN_n[a3] * la["NTDR_IN", i - 1, a3]) +
                       ## Reinfection of Resolved - to DR
                      -(NTDR_lambda[i - 1, a4] * NTDR_p[a4] * NTDR_x * la["NTDR_R", i - 1, a4]) +
                      -(NTDR_lambda[i - 1, a3] * NTDR_p[a3] * NTDR_x * la["NTDR_R", i - 1, a3]) +
                       ## Reactivation of Resolved
                      -(NTDR_r[a4] * la["NTDR_R", i - 1, a4]) +
                      -(NTDR_r[a3] * la["NTDR_R", i - 1, a3]) +
                       ## Reinfection of Resolved - to DS
                      -(NTDS_lambda[i - 1, a4] * NTDS_p[a4] * NTDS_x * la["NTDR_R", i - 1, a4]) +
                      -(NTDS_lambda[i - 1, a3] * NTDS_p[a3] * NTDS_x * la["NTDR_R", i - 1, a3])

## Misdiagnosed and Treated - Infectious
la["NTDR_mdt", i, a4] <-  (surv4 * la["NTDR_mdt", i - 1, a4]) +
                          # Misdiagnosed and treated
                          (NTDR_II_mdt[a4] * la["NTDR_II", i - 1, a4]) +
                          # Exit outwards
                          -(la["NTDR_mdt", i - 1, a4] * NTDR_mdt_exit) +
                          -(NTDR_mdt_u[a4] * la["NTDR_mdt", i - 1, a4]) +
                          (surv3 * la["NTDR_mdt", i - 1, a3]) +
                          # Misdiagnosed and treated
                          (NTDR_II_mdt[a3] * la["NTDR_II", i - 1, a3]) +
                          # Exit outwards
                          -(la["NTDR_mdt", i - 1, a3] * NTDR_mdt_exit) +
                          -(NTDR_mdt_u[a3] * la["NTDR_mdt", i - 1, a3])

## Misdiagnosed and Treated - Non-Infectious
la["NTDR_nid", i, a4] <-  (surv4 * la["NTDR_nid", i - 1, a4]) +
                          # Misdiagnosed and treated
                          (NTDR_IN_nid[a4] * la["NTDR_IN", i - 1, a4]) +
                          # Exit outwards
                          -(la["NTDR_nid", i - 1, a4] * NTDR_nid_exit) +
                          -(NTDR_nid_u[a4] * la["NTDR_nid", i - 1, a4]) +
                          (surv3 * la["NTDR_mdt", i - 1, a3]) +
                          # Misdiagnosed and treated
                          (NTDR_IN_nid[a3] * la["NTDR_IN", i - 1, a3]) +
                          # Exit outwards
                          -(la["NTDR_nid", i - 1, a3] * NTDR_nid_exit) +
                          -(NTDR_nid_u[a3] * la["NTDR_nid", i - 1, a3])

rm(surv)
surv <- 1 - (u[a0] )

## Previously Treated Drug Sensitive TB

## Infectious PTDS
la["PTDS_II", i, a1] <- (surv * la["PTDS_II", i - 1, a0]) +
                        ## Conversion of non-infectious to infectious
                        (PTDS_omega * la["PTDS_IN", i - 1, a0]) +
                        ## Reinfection of Resolved - DS
                        (PTDS_lambda[i - 1, a0] * PTDS_p[a0] * PTDS_f[a0] * PTDS_x * la["PTDS_R", i - 1, a0]) +
                        ## Reactivation from Resolved
                        (PTDS_f[a0] * PTDS_r[a0] * la["PTDS_R", i - 1, a0]) +
                        ## Reinfection of Resolved - DR
                        (PTDS_lambda[i - 1, a0] * PTDS_p[a0] * PTDS_f[a0] * PTDS_x * la["PTDR_R", i - 1, a0]) +
                        ## Natural cure
                       -(PTDS_II_n[a0] * la["PTDS_II", i - 1, a0]) +
                        ## Case detection
                       -(PTDS_II_kappa[a0] * la["PTDS_II", i - 1, a0]) +
                        ## TB death
                       -(PTDS_II_u[a0] * la["PTDS_II", i - 1, a0])

## Non-Infectious PTDS
la["PTDS_IN", i, a1] <- (surv * la["PTDS_IN", i - 1, a0]) +
                        ## NTDS Treatment Failures
                        (NTDS_phi * la["NTDS_T", i - 1, a0]) +
                        ## PTDS Treatment Failures
                        (PTDS_phi * la["PTDS_T", i - 1, a0]) +
                        ## Reinfection of Resolved - DS
                        (PTDS_lambda[i - 1, a0] * PTDS_p[a0] * (1 - PTDS_f[a0]) * PTDS_x * la["PTDS_R", i - 1, a0]) +
                        ## Reactivation from Resolved
                        ((1 - PTDS_f[a0]) * PTDS_r[a0] * la["PTDS_R", i - 1, a0]) +
                        ## Reinfection of Resolved - DR
                        (PTDS_lambda[i - 1, a0] * PTDS_p[a0] * (1 - PTDS_f[a0]) * PTDS_x * la["PTDR_R", i - 1, a0]) +
                        ## Natural cure
                       -(PTDS_IN_n[a0] * la["PTDS_IN", i - 1, a0]) +
                        ## Case detection
                       -(PTDS_IN_kappa[a0] * la["PTDS_IN", i - 1, a0]) +
                        ## TB death
                       -(PTDS_IN_u[a0] * la["PTDS_IN", i - 1, a0]) +
                        ## Conversion from non-infectious to infectious
                       -(PTDS_omega * la["PTDS_IN", i - 1, a0])

## PTDS on treatment
la["PTDS_T", i, a1] <-  (surv * la["PTDS_T", i - 1, a0]) +
                        ## Treatment initiation from Infectious
                        (PTDS_II_kappa[a0] * la["PTDS_II", i - 1, a0]) +
                        ## Treatment initiation from Non-infectious
                        (PTDS_IN_kappa[a0] * la["PTDS_IN", i - 1, a0]) +
                        ## Resistance acquisition
                       -(PT_xi * la["PTDS_T", i - 1, a0]) +
                        ## Treatment failures
                       -(PTDS_phi * la["PTDS_T", i - 1, a0]) +
                        ## Treatment successes
                       -(PTDS_psi * la["PTDS_T", i - 1, a0]) +
                        ## Death on treatment
                       -(PTDS_T_u[a0] * la["PTDS_T", i - 1, a0])

## PTDS Resolved
la["PTDS_R", i, a1] <- (surv * la["PTDS_R", i - 1, a0]) +
                       ## Natural cure - Infectious
                       (PTDS_II_n[a0] * la["PTDS_II", i - 1, a0]) +
                       ## Natural cure - Non-infectious
                       (PTDS_IN_n[a0] * la["PTDS_IN", i - 1, a0]) +
                       ## NTDS Treatment Success
                       (NTDS_psi * la["NTDS_T", i - 1, a0]) +
                       ## PTDS Treatment Success
                       (PTDS_psi * la["PTDS_T", i - 1, a0]) +
                       ## Reinfection of Resolved - to DS
                      -(PTDS_lambda[i - 1, a0] * PTDS_p[a0] * PTDS_x * la["PTDS_R", i - 1, a0]) +
                       ## Reactivation of Resolved
                      -(PTDS_r[a0] * la["PTDS_R", i - 1, a0]) +
                       ## Reinfection of Resolved - to DR
                      -(PTDR_lambda[i - 1, a0] * PTDR_p[a0] * PTDR_x * la["PTDS_R", i - 1, a0])

rm(surv4, surv3)
surv4 <- 1 - (u[a4] )
surv3 <- 1 - (u[a3] )

## Previously Treated for TB : Drug Sensitive
### Infectious DS-TB
la["PTDS_II", i, a4] <- (surv4 * la["PTDS_II", i - 1, a4]) +
                        (PTDS_omega * la["PTDS_IN", i - 1, a4]) +
                        (PTDS_lambda[i - 1, a4] * PTDS_p[a4] * PTDS_f[a4] * PTDS_x * la["PTDS_R", i - 1, a4]) +
                        (PTDS_f[a4] * PTDS_r[a4] * la["PTDS_R", i - 1, a4]) +
                        (PTDS_lambda[i - 1, a4] * PTDS_p[a4] * PTDS_f[a4] * PTDS_x * la["PTDR_R", i - 1, a4]) +
                       -(PTDS_II_n[a4] * la["PTDS_II", i - 1, a4]) +
                       -(PTDS_II_kappa[a4] * la["PTDS_II", i - 1, a4]) +
                       -(PTDS_II_u[a4] * la["PTDS_II", i - 1, a4]) +
                        # Accumulate
                        (surv3 * la["PTDS_II", i - 1, a3]) +
                        (PTDS_omega * la["PTDS_IN", i - 1, a3]) +
                        (PTDS_lambda[i - 1, a3] * PTDS_p[a3] * PTDS_f[a3] * PTDS_x * la["PTDS_R", i - 1, a3]) +
                        (PTDS_f[a3] * PTDS_r[a3] * la["PTDS_R", i - 1, a3]) +
                        (PTDS_lambda[i - 1, a3] * PTDS_p[a3] * PTDS_f[a3] * PTDS_x * la["PTDR_R", i - 1, a3]) +
                       -(PTDS_II_n[a3] * la["PTDS_II", i - 1, a3]) +
                       -(PTDS_II_kappa[a3] * la["PTDS_II", i - 1, a3]) +
                       -(PTDS_II_u[a3] * la["PTDS_II", i - 1, a3])

### Non-Infectious DS-TB
la["PTDS_IN", i, a4] <- (surv4 * la["PTDS_IN", i - 1, a4]) +
                        (NTDS_phi * la["NTDS_T", i - 1, a4]) +
                        (PTDS_phi * la["PTDS_T", i - 1, a4]) +
                        (PTDS_lambda[i - 1, a4] * PTDS_p[a4] * (1 - PTDS_f[a4]) * PTDS_x * la["PTDS_R", i - 1, a4]) +
                        ((1 - PTDS_f[a4]) * PTDS_r[a4] * la["PTDS_R", i - 1, a4]) +
                        (PTDS_lambda[i - 1, a4] * PTDS_p[a4] * (1 - PTDS_f[a4]) * PTDS_x * la["PTDR_R", i - 1, a4]) +
                       -(PTDS_IN_n[a4] * la["PTDS_IN", i - 1, a4]) +
                       -(PTDS_IN_kappa[a4] * la["PTDS_IN", i - 1, a4]) +
                       -(PTDS_IN_u[a4] * la["PTDS_IN", i - 1, a4]) +
                       -(PTDS_omega * la["PTDS_IN", i - 1, a4]) +
                       # Accumulate
                        (surv3 * la["PTDS_IN", i - 1, a3]) +
                        (NTDS_phi * la["NTDS_T", i - 1, a3]) +
                        (PTDS_phi * la["PTDS_T", i - 1, a3]) +
                        (PTDS_lambda[i - 1, a3] * PTDS_p[a3] * (1 - PTDS_f[a3]) * PTDS_x * la["PTDS_R", i - 1, a3]) +
                        ((1 - PTDS_f[a3]) * PTDS_r[a3] * la["PTDS_R", i - 1, a3]) +
                        (PTDS_lambda[i - 1, a3] * PTDS_p[a3] * (1 - PTDS_f[a3]) * PTDS_x * la["PTDR_R", i - 1, a3]) +
                       -(PTDS_IN_n[a3] * la["PTDS_IN", i - 1, a3]) +
                       -(PTDS_IN_kappa[a3] * la["PTDS_IN", i - 1, a3]) +
                       -(PTDS_IN_u[a3] * la["PTDS_IN", i - 1, a3]) +
                       -(PTDS_omega * la["PTDS_IN", i - 1, a3])

### Treatment DS-TB
la["PTDS_T", i, a4] <- (surv4 * la["PTDS_T", i - 1, a4]) +
                      (PTDS_II_kappa[a4] * la["PTDS_II", i - 1, a4]) +
                      (PTDS_IN_kappa[a4] * la["PTDS_IN", i - 1, a4]) -
                      (PT_xi * la["PTDS_T", i - 1, a4]) -
                      (PTDS_phi * la["PTDS_T", i - 1, a4]) -
                      (PTDS_psi * la["PTDS_T", i - 1, a4]) -
                      (PTDS_T_u[a4] * la["PTDS_T", i - 1, a4]) +
                      #Accumulate
                      (surv3 * la["PTDS_T", i - 1, a3]) +
                      (PTDS_II_kappa[a3] * la["PTDS_II", i - 1, a3]) +
                      (PTDS_IN_kappa[a3] * la["PTDS_IN", i - 1, a3]) -
                      (PT_xi * la["PTDS_T", i - 1, a3]) -
                      (PTDS_phi * la["PTDS_T", i - 1, a3]) -
                      (PTDS_psi * la["PTDS_T", i - 1, a3]) -
                      (PTDS_T_u[a3] * la["PTDS_T", i - 1, a3])

### PTDS Resolved
la["PTDS_R", i, a4] <-  (surv4 * la["PTDS_R", i - 1, a4]) +
                        (PTDS_II_n[a4] * la["PTDS_II", i - 1, a4]) +
                        (PTDS_IN_n[a4] * la["PTDS_IN", i - 1, a4]) +
                        ((NTDS_psi * la["NTDS_T", i - 1, a4]) +
                        (PTDS_psi * la["PTDS_T", i - 1, a4])) +
                       -(PTDS_lambda[i - 1, a4] * PTDS_p[a4] * PTDS_x * la["PTDS_R", i - 1, a4]) +
                       -(PTDS_r[a4] * la["PTDS_R", i - 1, a4]) +
                       -(PTDR_lambda[i - 1, a4] * PTDR_p[a4] * PTDR_x * la["PTDS_R", i - 1, a4]) +
                        # Accumulate
                        (surv3 * la["PTDS_R", i - 1, a3]) +
                        (PTDS_II_n[a3] * la["PTDS_II", i - 1, a3]) +
                        (PTDS_IN_n[a3] * la["PTDS_IN", i - 1, a3]) +
                        (NTDS_psi * la["NTDS_T", i - 1, a3]) +
                        (PTDS_psi * la["PTDS_T", i - 1, a3]) +
                       -(PTDS_lambda[i - 1, a3] * PTDS_p[a3] * PTDS_x * la["PTDS_R", i - 1, a3]) +
                       -(PTDS_r[a3] * la["PTDS_R", i - 1, a3]) +
                       -(PTDR_lambda[i - 1, a3] * PTDR_p[a3] * PTDR_x * la["PTDS_R", i - 1, a3])

rm(surv)
surv <- 1 - (u[a0] )

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
                             -(PTDR_II_u[a0] * la["PTDR_T_INphi", i - 1, a0])

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

rm(surv4, surv3)
surv4 <- 1 - (u[a4] )
surv3 <- 1 - (u[a3] )

## Previously Treated for TB: Drug Resistant
### DR-LTBI
la["PTDR_L", i, a4] <-  (surv4 * la["PTDR_L", i - 1, a4]) +
                        (PT_xi * (1 - PTDR_p[a4]) * la["PTDS_T", i - 1, a4]) +
                        (NT_xi * (1 - PTDR_p[a4]) * la["NTDS_T", i - 1, a4]) +
                       -(PTDR_v[a4] * la["PTDR_L", i - 1, a4]) +
                        # Accumulate
                        (surv3 * la["PTDR_L", i - 1, a3]) +
                        (PT_xi * (1 - PTDR_p[a3]) * la["PTDS_T", i - 1, a3]) +
                        (NT_xi * (1 - PTDR_p[a3]) * la["NTDS_T", i - 1, a3]) +
                       -(PTDR_v[a3] * la["PTDR_L", i - 1, a3])

### Infectious DR-TB
la["PTDR_II", i, a4] <- (surv4 * la["PTDR_II", i - 1, a4]) +
                        (PTDR_omega * la["PTDR_IN", i - 1, a4]) +
                        (PTDR_v[a4] * PTDR_f[a4] * la["PTDR_L", i - 1, a4]) +
                        (PTDR_lambda[i - 1, a4] * PTDR_p[a4] * PTDR_f[a4] * PTDR_x * la["PTDR_R", i - 1, a4]) +
                        (PTDR_f[a4] * (PTDR_r[a4] ) * la["PTDR_R", i - 1, a4]) +
                        (PTDR_lambda[i - 1, a4] * PTDR_p[a4] * PTDR_f[a4] * PTDR_x * la["PTDS_R", i - 1, a4]) +
                        (NT_xi * PTDR_p[a4] * PTDR_f[a4] * la["NTDS_T", i - 1, a4] * nt_mdt_loss) +
                        (PT_xi * PTDR_p[a4] * PTDR_f[a4] * la["PTDS_T", i - 1, a4] * pt_mdt_loss) +
                        (PTDR_T_IIphi_tau * la["PTDR_T_IIphi", i - 1, a4]) +
                        (NTDR_T_IIphi_tau * la["NTDR_T_IIphi", i - 1, a4]) +
                        (la["NTDR_mdt", i - 1, a4] * NTDR_mdt_exit * nt_mdt_loss) +
                        (la["PTDR_mdt", i - 1, a4] * PTDR_mdt_exit * pt_mdt_loss) +
                       -(PTDR_II_n[a4] * la["PTDR_II", i - 1, a4]) +
                       -(PTDR_II_kappa[a4] * la["PTDR_II", i - 1, a4]) +
                       -(PTDR_II_mdt[a4] * la["PTDR_II", i - 1, a4]) +
                       -(PTDR_II_u[a4] * la["PTDR_II", i - 1, a4]) +
                        # Accumulate
                        (surv3 * la["PTDR_II", i - 1, a3]) +
                        (PTDR_omega * la["PTDR_IN", i - 1, a3]) +
                        (PTDR_v[a3] * PTDR_f[a3] * la["PTDR_L", i - 1, a3]) +
                        (PTDR_lambda[i - 1, a3] * PTDR_p[a3] * PTDR_f[a3] * PTDR_x * la["PTDR_R", i - 1, a3]) +
                        (PTDR_f[a3] * (PTDR_r[a3] ) * la["PTDR_R", i - 1, a3]) +
                        (PTDR_lambda[i - 1, a3] * PTDR_p[a3] * PTDR_f[a3] * PTDR_x * la["PTDS_R", i - 1, a3]) +
                        (NT_xi * PTDR_p[a3] * PTDR_f[a3] * la["NTDS_T", i - 1, a3] * nt_mdt_loss) +
                        (PT_xi * PTDR_p[a3] * PTDR_f[a3] * la["PTDS_T", i - 1, a3] * pt_mdt_loss) +
                        (PTDR_T_IIphi_tau * la["PTDR_T_IIphi", i - 1, a3]) +
                        (NTDR_T_IIphi_tau * la["NTDR_T_IIphi", i - 1, a3]) +
                        (la["NTDR_mdt", i - 1, a3] * NTDR_mdt_exit * nt_mdt_loss) +
                        (la["PTDR_mdt", i - 1, a3] * PTDR_mdt_exit * pt_mdt_loss) +
                       -(PTDR_II_n[a3] * la["PTDR_II", i - 1, a3]) +
                       -(PTDR_II_kappa[a3] * la["PTDR_II", i - 1, a3]) +
                       -(PTDR_II_mdt[a3] * la["PTDR_II", i - 1, a3]) +
                       -(PTDR_II_u[a3] * la["PTDR_II", i - 1, a3])

### Non-Infectious DR-TB
la["PTDR_IN", i, a4] <- (surv4 * la["PTDR_IN", i - 1, a4]) +
                        (PTDR_v[a4] * (1 - PTDR_f[a4]) * la["PTDR_L", i - 1, a4]) +
                        (PTDR_lambda[i - 1, a4] * PTDR_p[a4] * (1 - PTDR_f[a4]) * PTDR_x * la["PTDR_R", i - 1, a4]) +
                        ((1 - PTDR_f[a4]) * PTDR_r[a4] * la["PTDR_R", i - 1, a4]) +
                        (PTDR_lambda[i - 1, a4] * PTDR_p[a4] * (1 - PTDR_f[a4]) * PTDR_x * la["PTDS_R", i - 1, a4]) +
                        (PT_xi * PTDR_p[a4] * (1 - PTDR_f[a4]) * la["PTDS_T", i - 1, a4] * pt_nid_loss) +
                        (NT_xi * PTDR_p[a4] * (1 - PTDR_f[a4]) * la["NTDS_T", i - 1, a4] * nt_nid_loss) +
                        (PTDR_T_INphi_tau * la["PTDR_T_INphi", i - 1, a4]) +
                        (NTDR_T_INphi_tau * la["NTDR_T_INphi", i - 1, a4]) +
                        (la["NTDR_nid", i - 1, a4] * NTDR_nid_exit * nt_nid_loss) +
                        (la["PTDR_nid", i - 1, a4] * PTDR_nid_exit * pt_nid_loss) +
                       -(PTDR_IN_n[a4] *la["PTDR_IN", i - 1, a4]) +
                       -(PTDR_IN_kappa[a4] * la["PTDR_IN", i - 1, a4]) +
                       -(PTDR_IN_nid[a4] * la["PTDR_IN", i - 1, a4]) +
                       -(PTDR_IN_u[a4] * la["PTDR_IN", i - 1, a4]) +
                       -(PTDR_omega * la["PTDR_IN", i - 1, a4]) +
                        # Accumulate
                        (surv3 * la["PTDR_IN", i - 1, a3]) +
                        (PTDR_v[a3] * (1 - PTDR_f[a3]) * la["PTDR_L", i - 1, a3]) +
                        (PTDR_lambda[i - 1, a3] * PTDR_p[a3] * (1 - PTDR_f[a3]) * PTDR_x * la["PTDR_R", i - 1, a3]) +
                        ((1 - PTDR_f[a3]) * PTDR_r[a3] * la["PTDR_R", i - 1, a3]) +
                        (PTDR_lambda[i - 1, a3] * PTDR_p[a3] * (1 - PTDR_f[a3]) * PTDR_x * la["PTDS_R", i - 1, a3]) +
                        (PT_xi * PTDR_p[a3] * (1 - PTDR_f[a3]) * la["PTDS_T", i - 1, a3] * pt_nid_loss) +
                        (NT_xi * PTDR_p[a3] * (1 - PTDR_f[a3]) * la["NTDS_T", i - 1, a3] * nt_nid_loss) +
                        (PTDR_T_INphi_tau * la["PTDR_T_INphi", i - 1, a3]) +
                        (NTDR_T_INphi_tau * la["NTDR_T_INphi", i - 1, a3]) +
                        (la["NTDR_nid", i - 1, a3] * NTDR_nid_exit * nt_nid_loss) +
                        (la["PTDR_nid", i - 1, a3] * PTDR_nid_exit * pt_nid_loss) +
                       -(PTDR_IN_n[a3] *la["PTDR_IN", i - 1, a3]) +
                       -(PTDR_IN_kappa[a3] * la["PTDR_IN", i - 1, a3]) +
                       -(PTDR_IN_nid[a3] * la["PTDR_IN", i - 1, a3]) +
                       -(PTDR_IN_u[a3] * la["PTDR_IN", i - 1, a3]) +
                       -(PTDR_omega * la["PTDR_IN", i - 1, a3])

### Treatment Succeeding - Originating from Infectious DR-TB
la["PTDR_T_IIpsi", i, a4] <- (surv4 * la["PTDR_T_IIpsi", i - 1, a4]) +
                      (PTDR_II_kappa[a4] * PTDR_II_psi * la["PTDR_II", i - 1, a4]) +
                      (la["NTDR_mdt", i - 1, a4] * NTDR_mdt_exit * (1 - nt_mdt_loss) * PTDR_II_psi) +
                      (la["PTDR_mdt", i - 1, a4] * PTDR_mdt_exit * (1 - pt_mdt_loss) * PTDR_II_psi) +
                      (NT_xi * PTDR_p[a4] * PTDR_f[a4] * la["NTDS_T", i - 1, a4] * (1 - nt_mdt_loss) * PTDR_II_psi) +
                      (PT_xi * PTDR_p[a4] * PTDR_f[a4] * la["PTDS_T", i - 1, a4] * (1 - pt_mdt_loss) * PTDR_II_psi) +
                     -(PTDR_T_IIpsi_tau * la["PTDR_T_IIpsi", i - 1, a4]) +
                     -(PTDR_T_u[a4] * la["PTDR_T_IIpsi", i - 1, a4]) +
                      #Accumulate
                      (surv3 * la["PTDR_T_IIpsi", i - 1, a3]) +
                      (PTDR_II_kappa[a3] * PTDR_II_psi * la["PTDR_II", i - 1, a3]) +
                      (la["NTDR_mdt", i - 1, a3] * NTDR_mdt_exit * (1 - nt_mdt_loss) * PTDR_II_psi) +
                      (la["PTDR_mdt", i - 1, a3] * PTDR_mdt_exit * (1 - pt_mdt_loss) * PTDR_II_psi) +
                      (NT_xi * PTDR_p[a3] * PTDR_f[a3] * la["NTDS_T", i - 1, a3] * (1 - nt_mdt_loss) * PTDR_II_psi) +
                      (PT_xi * PTDR_p[a3] * PTDR_f[a3] * la["PTDS_T", i - 1, a3] * (1 - pt_mdt_loss) * PTDR_II_psi) +
                     -(PTDR_T_IIpsi_tau * la["PTDR_T_IIpsi", i - 1, a3]) +
                     -(PTDR_T_u[a3] * la["PTDR_T_IIpsi", i - 1, a3])

### Treatment Failing - Originating from Infectious DR-TB
la["PTDR_T_IIphi", i, a4] <- (surv4 * la["PTDR_T_IIphi", i - 1, a4]) +
                      (PTDR_II_kappa[a4] * PTDR_II_phi * la["PTDR_II", i - 1, a4]) +
                      (la["NTDR_mdt", i - 1, a4] * NTDR_mdt_exit * (1 - nt_mdt_loss) * PTDR_II_phi) +
                      (la["PTDR_mdt", i - 1, a4] * PTDR_mdt_exit * (1 - pt_mdt_loss) * PTDR_II_phi) +
                      (NT_xi * PTDR_p[a4] * PTDR_f[a4] * la["NTDS_T", i - 1, a4] * (1 - nt_mdt_loss) * PTDR_II_phi) +
                      (PT_xi * PTDR_p[a4] * PTDR_f[a4] * la["PTDS_T", i - 1, a4] * (1 - pt_mdt_loss) * PTDR_II_phi) +
                     -(PTDR_T_IIphi_tau * la["PTDR_T_IIphi", i - 1, a4]) +
                     -(PTDR_II_u[a4] * la["PTDR_T_IIphi", i - 1, a4]) +
                      #Accumulate
                      (surv3 * la["PTDR_T_IIphi", i - 1, a3]) +
                      (PTDR_II_kappa[a3] * PTDR_II_phi * la["PTDR_II", i - 1, a3]) +
                      (la["NTDR_mdt", i - 1, a3] * NTDR_mdt_exit * (1 - nt_mdt_loss) * PTDR_II_phi) +
                      (la["PTDR_mdt", i - 1, a3] * PTDR_mdt_exit * (1 - pt_mdt_loss) * PTDR_II_phi) +
                      (NT_xi * PTDR_p[a3] * PTDR_f[a3] * la["NTDS_T", i - 1, a3] * (1 - nt_mdt_loss) * PTDR_II_phi) +
                      (PT_xi * PTDR_p[a3] * PTDR_f[a3] * la["PTDS_T", i - 1, a3] * (1 - pt_mdt_loss) * PTDR_II_phi) +
                     -(PTDR_T_IIphi_tau * la["PTDR_T_IIphi", i - 1, a3]) +
                     -(PTDR_II_u[a3] * la["PTDR_T_IIphi", i - 1, a3])

### Treatment Succeeding - Originating from Non-Infectious DR-TB
la["PTDR_T_INpsi", i, a4] <- (surv4 * la["PTDR_T_INpsi", i - 1, a4]) +
                      (PTDR_IN_kappa[a4] * PTDR_IN_psi * la["PTDR_IN", i - 1, a4]) +
                      (la["NTDR_nid", i - 1, a4] * NTDR_nid_exit * (1 - nt_nid_loss) * PTDR_IN_psi) +
                      (la["PTDR_nid", i - 1, a4] * PTDR_nid_exit * (1 - pt_nid_loss) * PTDR_IN_psi) +
                      (NT_xi * PTDR_p[a4] * (1 - PTDR_f[a4]) * la["NTDS_T", i - 1, a4] * (1 - nt_nid_loss) * PTDR_IN_psi) +
                      (PT_xi * PTDR_p[a4] * (1 - PTDR_f[a4]) * la["PTDS_T", i - 1, a4] * (1 - pt_nid_loss) * PTDR_IN_psi) +
                     -(PTDR_T_INpsi_tau * la["PTDR_T_INpsi", i - 1, a4]) +
                     -(PTDR_T_u[a4] * la["PTDR_T_INpsi", i - 1, a4]) +
                      #Accumulate
                      (surv3 * la["PTDR_T_INpsi", i - 1, a3]) +
                      (PTDR_IN_kappa[a3] * PTDR_IN_psi * la["PTDR_IN", i - 1, a3]) +
                      (la["NTDR_nid", i - 1, a3] * NTDR_nid_exit * (1 - nt_nid_loss) * PTDR_IN_psi) +
                      (la["PTDR_nid", i - 1, a3] * PTDR_nid_exit * (1 - pt_nid_loss) * PTDR_IN_psi) +
                      (NT_xi * PTDR_p[a3] * (1 - PTDR_f[a3]) * la["NTDS_T", i - 1, a3] * (1 - nt_nid_loss) * PTDR_IN_psi) +
                      (PT_xi * PTDR_p[a3] * (1 - PTDR_f[a3]) * la["PTDS_T", i - 1, a3] * (1 - pt_nid_loss) * PTDR_IN_psi) +
                     -(PTDR_T_INpsi_tau * la["PTDR_T_INpsi", i - 1, a3]) +
                     -(PTDR_T_u[a3] * la["PTDR_T_INpsi", i - 1, a3])

### Treatment Failing - Originating from Non-Infectious DR-TB
la["PTDR_T_INphi", i, a4] <- (surv4 * la["PTDR_T_INphi", i - 1, a4]) +
                      (PTDR_IN_kappa[a4] * PTDR_IN_phi * la["PTDR_IN", i - 1, a4]) +
                      (la["NTDR_nid", i - 1, a4] * NTDR_nid_exit * (1 - nt_nid_loss) * PTDR_IN_phi) +
                      (la["PTDR_nid", i - 1, a4] * PTDR_nid_exit * (1 - pt_nid_loss) * PTDR_IN_phi) +
                      (NT_xi * PTDR_p[a4] * (1 - PTDR_f[a4]) * la["NTDS_T", i - 1, a4] * (1 - nt_nid_loss) * PTDR_IN_phi) +
                      (PT_xi * PTDR_p[a4] * (1 - PTDR_f[a4]) * la["PTDS_T", i - 1, a4] * (1 - pt_nid_loss) * PTDR_IN_phi) +
                     -(PTDR_T_INphi_tau * la["PTDR_T_INphi", i - 1, a4]) +
                     -(PTDR_II_u[a4] * la["PTDR_T_INphi", i - 1, a4]) +
                      #Accumulate
                      (surv3 * la["PTDR_T_INphi", i - 1, a3]) +
                      (PTDR_IN_kappa[a3] * PTDR_IN_phi * la["PTDR_IN", i - 1, a3]) +
                      (la["NTDR_nid", i - 1, a3] * NTDR_nid_exit * (1 - nt_nid_loss) * PTDR_IN_phi) +
                      (la["PTDR_nid", i - 1, a3] * PTDR_nid_exit * (1 - pt_nid_loss) * PTDR_IN_phi) +
                      (NT_xi * PTDR_p[a3] * (1 - PTDR_f[a3]) * la["NTDS_T", i - 1, a3] * (1 - nt_nid_loss) * PTDR_IN_phi) +
                      (PT_xi * PTDR_p[a3] * (1 - PTDR_f[a3]) * la["PTDS_T", i - 1, a3] * (1 - pt_nid_loss) * PTDR_IN_phi) +
                     -(PTDR_T_INphi_tau * la["PTDR_T_INphi", i - 1, a3]) +
                     -(PTDR_II_u[a3] * la["PTDR_T_INphi", i - 1, a3])

### Resolved - DR TB
la["PTDR_R", i, a4] <- (surv4 * la["PTDR_R", i - 1, a4]) +
                       (PTDR_II_n[a4] * la["PTDR_II", i - 1, a4]) +
                       (PTDR_IN_n[a4] * la["PTDR_IN", i - 1, a4]) +
                       (NTDR_T_IIpsi_tau * la["NTDR_T_IIpsi", i - 1, a4]) +
                       (NTDR_T_INpsi_tau * la["NTDR_T_INpsi", i - 1, a4]) +
                       (PTDR_T_INpsi_tau * la["PTDR_T_INpsi", i - 1, a4]) +
                       (PTDR_T_IIpsi_tau * la["PTDR_T_IIpsi", i - 1, a4]) +
                      -(PTDR_lambda[i - 1, a4] * PTDR_p[a4] * PTDR_x * la["PTDR_R", i - 1, a4]) +
                      -(PTDR_r[a4] * la["PTDR_R", i - 1, a4]) +
                      -(PTDS_lambda[i - 1, a4] * PTDS_p[a4] * PTDS_x * la["PTDR_R", i - 1, a4]) +
                       # Accumulate
                       (surv3 * la["PTDR_R", i - 1, a3]) +
                       (PTDR_II_n[a3] * la["PTDR_II", i - 1, a3]) +
                       (PTDR_IN_n[a3] * la["PTDR_IN", i - 1, a3]) +
                       (NTDR_T_IIpsi_tau * la["NTDR_T_IIpsi", i - 1, a3]) +
                       (NTDR_T_INpsi_tau * la["NTDR_T_INpsi", i - 1, a3]) +
                       (PTDR_T_INpsi_tau * la["PTDR_T_INpsi", i - 1, a3]) +
                       (PTDR_T_IIpsi_tau * la["PTDR_T_IIpsi", i - 1, a3]) +
                      -(PTDR_lambda[i - 1, a3] * PTDR_p[a3] * PTDR_x * la["PTDR_R", i - 1, a3]) +
                      -(PTDR_r[a3] * la["PTDR_R", i - 1, a3]) +
                      -(PTDS_lambda[i - 1, a3] * PTDS_p[a3] * PTDS_x * la["PTDR_R", i - 1, a3])

## Misdiagnosed and Treated - Infectious
la["PTDR_mdt", i, a4] <-  (surv4 * la["PTDR_mdt", i - 1, a4]) +
                          # Misdiagnosed and treated
                          (PTDR_II_mdt[a4] * la["PTDR_II", i - 1, a4]) +
                          # Exit outwards
                          -(la["PTDR_mdt", i - 1, a4] * PTDR_mdt_exit) +
                          -(PTDR_mdt_u[a4] * la["PTDR_mdt", i - 1, a4]) +
                          (surv3 * la["PTDR_mdt", i - 1, a3]) +
                          # Misdiagnosed and treated
                          (PTDR_II_mdt[a3] * la["PTDR_II", i - 1, a3]) +
                          # Exit outwards
                          -(la["PTDR_mdt", i - 1, a3] * PTDR_mdt_exit) +
                          -(PTDR_mdt_u[a3] * la["PTDR_mdt", i - 1, a3])

## Misdiagnosed and Treated - Non-Infectious
la["PTDR_nid", i, a4] <-  (surv4 * la["PTDR_nid", i - 1, a4]) +
                          # Misdiagnosed and treated
                          (PTDR_IN_nid[a4] * la["PTDR_IN", i - 1, a4]) +
                          # Exit outwards
                          -(la["PTDR_nid", i - 1, a4] * PTDR_nid_exit) +
                          -(PTDR_nid_u[a4] * la["PTDR_nid", i - 1, a4]) +
                          (surv3 * la["PTDR_mdt", i - 1, a3]) +
                          # Misdiagnosed and treated
                          (PTDR_IN_nid[a3] * la["PTDR_IN", i - 1, a3]) +
                          # Exit outwards
                          -(la["PTDR_nid", i - 1, a3] * PTDR_nid_exit) +
                          -(PTDR_nid_u[a3] * la["PTDR_nid", i - 1, a3])
