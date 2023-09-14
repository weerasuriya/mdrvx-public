# Unit cost sampler
library(here)
library(data.table)

setwd(here("code"))

# Read in cost-ranges
ca <- fread("../data/para_cost_ranges.csv")
ca[, range := max - min]
ca <- split(ca, by = 'param')

# Sample generator
sampler <- function(x) {
    rl <- list()
    rl$ds_dx_cost <- rnorm(n = x, mean = ca$ds_dx_cost$mlv, ca$ds_dx_cost$sd)
    rl$ds_tx_cost <- rnorm(n = x, mean = ca$ds_tx_cost$mlv, ca$ds_tx_cost$sd)
    rl$dr_dx_cost <- rl$ds_dx_cost * 1.2
    rl$dst_cost   <- rnorm(n = x, mean = ca$dst_cost$mlv, ca$dst_cost$sd)
    rl$vx_r       <- ca$vx_r$min + (ca$vx_r$range * rbeta(n = x, shape1 = ca$vx_r$alpha, shape2 = (ca$vx_r$beta_p * ca$vx_r$alpha)))
    rl$vx_c       <- ca$vx_c$min + (ca$vx_c$range * rbeta(n = x, shape1 = ca$vx_c$alpha, shape2 = (ca$vx_c$beta_p * ca$vx_c$alpha)))
    rl$dr_r1_cost <- ca$dr_r1_cost$min + (ca$dr_r1_cost$range * rbeta(n = x, shape1 = ca$dr_r1_cost$alpha, shape2 = (ca$dr_r1_cost$beta_p * ca$dr_r1_cost$alpha)))
    rl$dr_r2_cost <- ca$dr_r2_cost$min + (ca$dr_r2_cost$range * rbeta(n = x, shape1 = ca$dr_r2_cost$alpha, shape2 = (ca$dr_r2_cost$beta_p * ca$dr_r2_cost$alpha)))
    rl$dr_sc_cost <- ca$dr_sc_cost$min + (ca$dr_sc_cost$range * rbeta(n = x, shape1 = ca$dr_sc_cost$alpha, shape2 = (ca$dr_sc_cost$beta_p * ca$dr_sc_cost$alpha)))
    return(rl)
}

set.seed(123)
fwrite(as.data.table(sampler(1000)), file = "../data/para_cost_sampled.csv")


