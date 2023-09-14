

#rm(surv)
surv <- 1 - u[a0]

## Previously Treated Drug Sensitive TB

## Infectious v_PTDS
la["v_PTDS_II", i, a1] <- (surv * la["v_PTDS_II", i - 1, a0]) +
                        ## Conversion of non-infectious to infectious
                        (v_PTDS_omega * la["v_PTDS_IN", i - 1, a0]) +
                        ## Reinfection of Resolved - DS
                        (v_PTDS_lambda[i - 1, a0] * v_PTDS_p[a0] * v_PTDS_f[a0] * v_PTDS_x * la["v_PTDS_R", i - 1, a0]) +
                        ## Reactivation from Resolved
                        (v_PTDS_f[a0] * v_PTDS_r[a0] * la["v_PTDS_R", i - 1, a0]) +
                        ## Reinfection of Resolved - DR
                        (v_PTDS_lambda[i - 1, a0] * v_PTDS_p[a0] * v_PTDS_f[a0] * v_PTDS_x * la["v_PTDR_R", i - 1, a0]) +
                        ## NTDS Treatment Failures
                        (v_NTDS_phi * la["v_NTDS_T", i - 1, a0] * (1 - ini_prop[current_step, a0])) +
                        ## PTDS Treatment Failures
                        (v_PTDS_phi * la["v_PTDS_T", i - 1, a0] * (1 - ini_prop[current_step, a0])) +
                        ## Natural cure
                       -(v_PTDS_II_n[a0] * la["v_PTDS_II", i - 1, a0]) +
                        ## Case detection
                       -(v_PTDS_II_kappa[a0] * la["v_PTDS_II", i - 1, a0]) +
                        ## TB death
                       -(v_PTDS_II_u[a0] * la["v_PTDS_II", i - 1, a0])

## Non-Infectious v_PTDS
la["v_PTDS_IN", i, a1] <- (surv * la["v_PTDS_IN", i - 1, a0]) +
                        ## NTDS Treatment Failures
                        (v_NTDS_phi * la["v_NTDS_T", i - 1, a0] * ini_prop[current_step, a0]) +
                        ## PTDS Treatment Failures
                        (v_PTDS_phi * la["v_PTDS_T", i - 1, a0] * ini_prop[current_step, a0]) +
                        ## Reinfection of Resolved - DS
                        (v_PTDS_lambda[i - 1, a0] * v_PTDS_p[a0] * (1 - v_PTDS_f[a0]) * v_PTDS_x * la["v_PTDS_R", i - 1, a0]) +
                        ## Reactivation from Resolved
                        ((1 - v_PTDS_f[a0]) * v_PTDS_r[a0] * la["v_PTDS_R", i - 1, a0]) +
                        ## Reinfection of Resolved - DR
                        (v_PTDS_lambda[i - 1, a0] * v_PTDS_p[a0] * (1 - v_PTDS_f[a0]) * v_PTDS_x * la["v_PTDR_R", i - 1, a0]) +
                        ## Natural cure
                       -(v_PTDS_IN_n[a0] * la["v_PTDS_IN", i - 1, a0]) +
                        ## Case detection
                       -(v_PTDS_IN_kappa[a0] * la["v_PTDS_IN", i - 1, a0]) +
                        ## TB death
                       -(v_PTDS_IN_u[a0] * la["v_PTDS_IN", i - 1, a0]) +
                        ## Conversion from non-infectious to infectious
                       -(v_PTDS_omega * la["v_PTDS_IN", i - 1, a0])

## v_PTDS on treatment
la["v_PTDS_T", i, a1] <-  (surv * la["v_PTDS_T", i - 1, a0]) +
                        ## Treatment initiation from Infectious
                        (v_PTDS_II_kappa[a0] * la["v_PTDS_II", i - 1, a0]) +
                        ## Treatment initiation from Non-infectious
                        (v_PTDS_IN_kappa[a0] * la["v_PTDS_IN", i - 1, a0]) +
                        ## Resistance acquisition
                       -(v_PT_xi * la["v_PTDS_T", i - 1, a0]) +
                        ## Treatment failures
                       -(v_PTDS_phi * la["v_PTDS_T", i - 1, a0]) +
                        ## Treatment successes
                       -(v_PTDS_psi * la["v_PTDS_T", i - 1, a0]) +
                        ## Death on treatment
                       -(v_PTDS_T_u[a0] * la["v_PTDS_T", i - 1, a0])

## v_PTDS Resolved
la["v_PTDS_R", i, a1] <- (surv * la["v_PTDS_R", i - 1, a0]) +
                       ## Natural cure - Infectious
                       (v_PTDS_II_n[a0] * la["v_PTDS_II", i - 1, a0]) +
                       ## Natural cure - Non-infectious
                       (v_PTDS_IN_n[a0] * la["v_PTDS_IN", i - 1, a0]) +
                       ## v_NTDS Treatment Success
                       (v_NTDS_psi * la["v_NTDS_T", i - 1, a0]) +
                       ## v_PTDS Treatment Success
                       (v_PTDS_psi * la["v_PTDS_T", i - 1, a0]) +
                       ## Reinfection of Resolved - to DS
                      -(v_PTDS_lambda[i - 1, a0] * v_PTDS_p[a0] * v_PTDS_x * la["v_PTDS_R", i - 1, a0]) +
                       ## Reactivation of Resolved
                      -(v_PTDS_r[a0] * la["v_PTDS_R", i - 1, a0]) +
                       ## Reinfection of Resolved - to DR
                      -(v_PTDR_lambda[i - 1, a0] * v_PTDR_p[a0] * v_PTDR_x * la["v_PTDS_R", i - 1, a0])
