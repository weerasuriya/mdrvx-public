
# rm(surv)
surv <- 1 - u[a0]

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
