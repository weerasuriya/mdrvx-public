
# rm(surv)
surv <- 1 - u[a0]

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
                        ## NTDS Treatment Failures
                        (NTDS_phi * la["NTDS_T", i - 1, a0] * (1 - ini_prop[current_step, a0])) +
                        ## PTDS Treatment Failures
                        (PTDS_phi * la["PTDS_T", i - 1, a0] * (1 - ini_prop[current_step, a0])) +
                        ## Natural cure
                       -(PTDS_II_n[a0] * la["PTDS_II", i - 1, a0]) +
                        ## Case detection
                       -(PTDS_II_kappa[a0] * la["PTDS_II", i - 1, a0]) +
                        ## TB death
                       -(PTDS_II_u[a0] * la["PTDS_II", i - 1, a0])

## Non-Infectious PTDS
la["PTDS_IN", i, a1] <- (surv * la["PTDS_IN", i - 1, a0]) +
                        ## NTDS Treatment Failures
                        (NTDS_phi * la["NTDS_T", i - 1, a0] * ini_prop[current_step, a0]) +
                        ## PTDS Treatment Failures
                        (PTDS_phi * la["PTDS_T", i - 1, a0] * ini_prop[current_step, a0]) +
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
