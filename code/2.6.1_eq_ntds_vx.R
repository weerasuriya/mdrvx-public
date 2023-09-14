
#rm(surv)
surv <- 1 - u[a0]

## Never Treated Drug Sensitive TB
## v_NTDS-Latent
la["v_NTDS_L", i, a1] <- (surv * la["v_NTDS_L", i - 1, a0]) +
                       ((1 - v_NTDS_p[a0]) * v_NTDS_lambda[i - 1, a0] * la["v_S", i - 1, a0]) +
                       ((1 - v_NTDS_p[a0]) * v_NTDS_lambda[i - 1, a0] * v_NTDS_x * la["v_NTDR_L", i - 1, a0]) +
                      -(v_NTDS_v[a0]  * la["v_NTDS_L", i - 1, a0]) +
                      -(v_NTDS_x * v_NTDS_lambda[i - 1, a0] * v_NTDS_p[a0] * la["v_NTDS_L", i - 1, a0]) +
                      -(v_NTDR_x * v_NTDR_lambda[i - 1, a0] * la["v_NTDS_L", i - 1, a0])

## Infectious v_NTDS
la["v_NTDS_II", i, a1] <- (surv * la["v_NTDS_II", i - 1, a0]) +
                        ## Conversion of Non-infectious to Infectious
                        (v_NTDS_omega  * la["v_NTDS_IN", i - 1, a0]) +
                        ## Deconstructed
                        (v_NTDS_p[a0] * v_NTDS_f[a0] * v_NTDS_lambda[i - 1, a0] * la["v_S", i - 1, a0]) +
                        (v_NTDS_p[a0] * v_NTDS_f[a0] * v_NTDS_lambda[i - 1, a0] * v_NTDS_x * la["v_NTDR_L", i - 1, a0]) +
                        (v_NTDS_p[a0] * v_NTDS_f[a0] * v_NTDS_lambda[i - 1, a0] * v_NTDS_x * la["v_NTDS_L", i - 1, a0]) +
                        ## End deconstruct
                        ## Reactivation from Latent
                        (v_NTDS_v[a0]  * v_NTDS_f[a0] * la["v_NTDS_L", i - 1, a0]) +
                        ## Reinfection of Resolved - DS
                        (v_NTDS_lambda[i - 1, a0] * v_NTDS_p[a0] * v_NTDS_f[a0] * v_NTDS_x * la["v_NTDS_R", i - 1, a0]) +
                        ## Reactivation from Resolved
                        (v_NTDS_f[a0] * v_NTDS_r[a0]  * la["v_NTDS_R", i - 1, a0]) +
                        ## Reinfection of Resolved - DR
                        (v_NTDS_lambda[i - 1, a0] * v_NTDS_p[a0] * v_NTDS_f[a0] * v_NTDS_x * la["v_NTDR_R", i - 1, a0]) +
                        ## Natural Cure
                       -(v_NTDS_II_n[a0]  * la["v_NTDS_II", i - 1, a0]) +
                        ## Case detection
                       -(v_NTDS_II_kappa[a0] * la["v_NTDS_II", i - 1, a0]) +
                        ## TB Death
                       -(v_NTDS_II_u[a0] * la["v_NTDS_II", i - 1, a0])

## Non-Infectious v_NTDS
la["v_NTDS_IN", i,a1] <- (surv * la["v_NTDS_IN", i-1,a0]) +
                       ## Deconstructed
                       (v_NTDS_p[a0] * (1 - v_NTDS_f[a0]) * v_NTDS_lambda[i - 1, a0] * la["v_S", i - 1, a0]) +
                       (v_NTDS_p[a0] * (1 - v_NTDS_f[a0]) * v_NTDS_lambda[i - 1, a0] * v_NTDS_x * la["v_NTDR_L", i - 1, a0]) +
                       (v_NTDS_p[a0] * (1 - v_NTDS_f[a0]) * v_NTDS_lambda[i - 1, a0] * v_NTDS_x * la["v_NTDS_L", i - 1, a0]) +
                       ## End deconstruct
                       ## Reactivation from Latent
                       (v_NTDS_v[a0]  * (1 - v_NTDS_f[a0]) * la["v_NTDS_L", i - 1, a0]) +
                       ## Reinfection of Resolved - DS
                       (v_NTDS_lambda[i - 1, a0] * v_NTDS_p[a0] * (1 - v_NTDS_f[a0]) * v_NTDS_x * la["v_NTDS_R", i - 1,a0]) +
                       ## Reactivation from Resolved
                       ((1 - v_NTDS_f[a0]) * v_NTDS_r[a0]  * la["v_NTDS_R", i-1,a0]) +
                       ## Reinfection of Resolved - DR
                       (v_NTDS_lambda[i - 1, a0] * v_NTDS_p[a0] * (1 - v_NTDS_f[a0]) * v_NTDS_x * la["v_NTDR_R", i - 1, a0]) +
                       ## Natural Cure
                      -(v_NTDS_IN_n[a0]  * la["v_NTDS_IN", i - 1, a0]) +
                       ## Case Detection
                      -(v_NTDS_IN_kappa[a0] * la["v_NTDS_IN", i - 1, a0]) +
                       ## TB Death
                      -(v_NTDS_IN_u[a0] * la["v_NTDS_IN", i - 1, a0]) +
                       ## Conversion from Non-infectious to Infectious
                      -(v_NTDS_omega  * la["v_NTDS_IN", i - 1, a0])

## v_NTDS on treatment
la["v_NTDS_T", i,a1]  <- (surv * la["v_NTDS_T", i-1,a0]) +
                       ## Treatment initiation from Infectious
                       (v_NTDS_II_kappa[a0]*la["v_NTDS_II", i-1,a0]) +
                       ## Treatment initiation from Non-infectious
                       (v_NTDS_IN_kappa[a0]*la["v_NTDS_IN", i-1,a0]) +
                       ## Resistance acquisition
                      -(v_NT_xi * la["v_NTDS_T", i-1,a0]) +
                       ## Treatment failures
                      -(v_NTDS_phi * la["v_NTDS_T", i-1,a0]) +
                       ## Treatment successes
                      -(v_NTDS_psi * la["v_NTDS_T", i-1,a0]) +
                       ## Death on treatment
                      -(v_NTDS_T_u[a0] * la["v_NTDS_T", i-1,a0])

# v_NTDS Resolved
la["v_NTDS_R", i, a1] <- (surv * la["v_NTDS_R", i - 1, a0]) +
                       ## Natural cure Infectious
                       (v_NTDS_II_n[a0]  * la["v_NTDS_II", i - 1, a0]) +
                       ## Natural cure Non-infectious
                       (v_NTDS_IN_n[a0]  * la["v_NTDS_IN", i - 1, a0]) +
                       ## Reinfection of Resolved - DS
                      -(v_NTDS_lambda[i - 1, a0] * v_NTDS_p[a0] * v_NTDS_x * la["v_NTDS_R", i - 1, a0]) +
                       ## Reactivation of Resolved
                      -(v_NTDS_r[a0]  * la["v_NTDS_R", i - 1, a0]) +
                       ## Reinfection of Resolved - DR
                      -(v_NTDR_lambda[i - 1, a0] * v_NTDR_p[a0] * v_NTDR_x * la["v_NTDS_R", i - 1, a0])
