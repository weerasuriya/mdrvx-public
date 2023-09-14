
#rm(surv)
surv <- 1 - u[a0]

# Susceptibles
la["v_S", i, a1] <-   (surv * la["v_S", i - 1, a0]) +
                   -(v_NTDS_lambda[i - 1, a0] * la["v_S", i - 1, a0]) +
                   -(v_NTDR_lambda[i - 1, a0] * la["v_S", i - 1, a0])
