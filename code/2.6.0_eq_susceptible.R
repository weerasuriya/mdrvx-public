
#rm(surv)
surv <- 1 - u[a0]

# Susceptibles
la["S", i, a1] <-   (surv * la["S", i - 1, a0]) +
                   -(NTDS_lambda[i - 1, a0] * la["S", i - 1, a0]) +
                   -(NTDR_lambda[i - 1, a0] * la["S", i - 1, a0])
