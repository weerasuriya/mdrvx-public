# This script is run per fitted parameter set

# Set up ------------------------------------------------------------------
{
   library(here)
   # library(tidyverse)
   library(data.table)
   library(matrixStats)
   library(digest)
   library(R.utils)
   # library(foreach)
   # library(doParallel)

   setwd(here("code"))

   # Starting conditions
   startcon <- list(
      dt = 0.25,
      year1 = 1900,
      yearend = 2099,
      prop_col = 1,
      bgu = 1,
      fert = 1,
      init_dir = here("output")
   )

   # Import relevant data for model run
   source("1_import.R")

   # Import mainloop function
   source("build/2_main_rebuild.R")

   country = "CN"
}
## Epi Calculation Functions ---------------------------------------------------
epi_sum <- function(sy = 2027, eyl = list(2030, 2035, 2050), t = "bl", epidf, excol = c("para_hash", "scen", "Year")) {
   tb <- paste0(t, ".", as.character(match.call()$epidf[[3]]))
   op <- epidf[, .(tb = rowSums(.SD)), .SDcol = -excol][, Year := epidf$Year]
   setnames(op, old = "tb", new = tb)
   opl <- lapply(eyl, function(x) {
      setnames(op[Year>=2027 & Year<=x, .(sum(.SD)), .SDcol=tb], old = "V1", new = tb)[, c("Year") := .(..x)]
   })
   opdt <- rbindlist(opl)
   setattr(opdt, name = "var", value = tb)
   opdt
}

diff_epi_sum <- function(bl, vx) {
   vn <- paste0("diff_", gsub(x = attributes(bl)$var, pattern = "^.*\\.", replacement = ""))
   opdt <- data.table(
      Year = bl$Year,
      vn = bl[[attributes(bl)$var]] - vx[[attributes(vx)$var]]
   )
   setnames(opdt, old = "vn", new = vn)
   melt(opdt, id.vars = 'Year', variable.name = 'endpoint', value.name = 'diff_endpoint')
}

pe_adjust <- function(pe, t = 'bl') {
   ir <- melt(pe, id.vars = 'para_hash', variable.name = 'endpoint', value.name = paste0(t, '_result'))
   ir[, endpoint := as.character(endpoint)]
   ir[, Year := gsub(pattern = '.*(\\d{4})', x = ir$endpoint, replacement = '\\1')]
   ir[, endpoint := gsub(pattern = '(.*)_\\d{4}', x = ir$endpoint, replacement = '\\1')]
   setattr(ir, "var", value = t)
   ir[, !c("para_hash")]
}

diff_pe <- function(bl, vx) {
   ir <- merge(bl, vx)
   ir[, diff_endpoint := .((bl_result - vx_result)/bl_result)]
   ir$bl_result = NULL
   ir$vx_result = NULL
   setcolorder(ir,neworder = "Year")
   ir

}

diff_initrx <- function(sy = 2027, eyl = list(2030, 2035, 2050), bldf, vxdf, excol = c("para_hash", "scen", "Year")) {
   cbdf <- copy(bldf)
   cbdf[, diff_DS := DSTB_initRx - vxdf[, DSTB_initRx]]
   cbdf[, diff_DR := DRTB_initRx - vxdf[, DRTB_initRx]]
   cbdf[, c("DSTB_initRx","DRTB_initRx", "scen") := .(NULL, NULL, NULL) ]
   op <- lapply(eyl, function(x) {
      cbdf[Year >= sy & Year <= x, lapply(.SD, sum), .SD = c('diff_DS', 'diff_DR'), by = "para_hash"][, Year := ..x]
   })
   opl <- rbindlist(op)
   opl
}

## Cost Calculation Functions ------------------------------------------
vx_cost <- function(sy = 2027, eyl = list(2030, 2035, 2050), ivxdf, dr = 0.03) {
   op <- lapply(eyl, function(x) {
      rbind(ivxdf[Year >= sy & Year <= x, lapply(.SD, sum), .SD = !c("Year")][, `:=`(c("Year", "discount"), .(..x, FALSE))], ivxdf[Year >= sy & Year <= x, lapply(.SD/(1 + dr)^(.I - 1), sum), .SD = !c("Year")][, `:=`(c("Year", "discount"), .(..x, TRUE))])
   })
   opdf <- rbindlist(op)
   melt(opdf, id.var = c("Year", "discount"), variable.name = "vx_profile", value.name = "vaccine_cost")
}

cost_calc <- function(sy = 2027, eyl = list(2030, 2035, 2050), bl.cost.df, vx.cost.df, dr = 0.03, vn = 'tbrx') {
   opl <- lapply(eyl, function(x) {
      op = data.table(Year = bl.cost.df$Year, diff_tbrx = vx.cost.df[[vn]] - bl.cost.df[[vn]])
      rbind(op[Year >= sy & Year <= x, .(diff_tbrx = sum(diff_tbrx))][, `:=`(c("Year","discount"), .(..x, FALSE))],
      op[Year >= sy & Year <= x, .(diff_tbrx = sum(diff_tbrx/(1 + dr)^(.I - 1)))][, `:=`(c("Year", "discount"), .(..x, TRUE))])
   })
   opdt <- rbindlist(opl)
   setnames(opdt, old = 'diff_tbrx', new = paste0('diff_', vn))
   opdt
}

daly_calc <- function(sy = 2027, eyl = c(2030, 2035, 2050), excol = c("Year", "para_hash", "scen"), mort_df, yld_df, ledf_df = ledf, disc = 0.03, daly_wt = dalywt, t) {
   ey = max(eyl)
   stopifnot(ncol(mort_df[, -..excol]) == ncol(ledf_df[, !("Year")]))
   yll  <- mort_df[Year >= sy & Year <= ey, -..excol] * ledf_df[Year >= sy & Year <= ey, !("Year")] * 1e+03
   yld  <- yld_df[Year >= sy & Year <= ey, -..excol] * dalywt
   daly <- (yll + yld)[, `:=`(Year, mort_df[Year >= sy & Year <= ey, Year])][, .(udiscDALYs = rowSums(.SD)), by = Year][, `:=`(discDALYs, udiscDALYs/(1 + disc)^(.I - 1))]
   nm   <- paste0(t, "cumDALY")
   opl  <- lapply(eyl, function(x) {
      rbind(daly[Year >= sy & Year <= x, .(cumDALY = sum(udiscDALYs))][, `:=`(c("Year", "discount"), .(..x, FALSE))],
      daly[Year >= sy & Year <= x, .(cumDALY = sum(discDALYs))][, `:=`(c("Year", "discount"), .(..x, TRUE))])
   })
   opdf <- rbindlist(opl)
   setnames(opdf, old = "cumDALY", new = paste0(t, "_cumDALY"))
   opdf
}

bimp <- function(sy = 2027, ey = 2050, vxca, blca, vn = c("Year", "para_hash", "scen", "ds_dx", "dr_dx", "dst", "ds_tx", "dr_tx", "tbrx", "ccdc_tbrx")) {
   diffca <- rbind(vxca[, ..vn], blca[, ..vn])[, lapply(.SD, diff), by = c('Year', 'para_hash', 'scen')]
   diffca[, id := 'diff']
   mca <- rbind(vxca[, c(id = 'vx', .SD), .SDcol = vn], blca[, c(id = 'bl', .SD), .SDcol = vn], diffca)
   mca[, c("vax_hash", "cost_hash") := .(vax_hash, cost_hash)]
   mca
}


# Read External Files -----------------------------------------------------
# Source compiled function
cmp_mainloop   <- readRDS("../code/2.1_main_rebuild.R.cmp")

# Read in relevant external data required for baseline vs vaccine runs, cost calculations etc
# Life expectancy table for YLL calculations
ledf           <- fread("../data/unpd_le_cleaned_CN_v2.csv")

# Subsampled fitted parameter sets (x1000) and costs
parameter_sets <- fread("../data/fitted_parameters/split_failure_1K.csv")
cost_sets      <- fread("../data/para_cost_sampled.csv")


# Initial Value Set up ----------------------------------------------------
# DALY weights +/- range
dalywt    <- 0.333
dalywt_hi <- 0.454
dalywt_lo <- 0.224

# Discount Rate
disc_rate <- 0.03

# Set up possible vaccine profiles
{
   vx_start_year <- 2027
   coverageR     <- 0.8
   coverageM     <- 0.7
   effI          <- c(0)
   effD          <- c(0.3, 0.5, 0.7, 0.9)
   duration      <- c(5, 10, 100)
   mass_interval <- c(10)
   ageM          <- 10
   ageR          <- 9
   vxtype        <- c(1, 2, 3)
}

# Create grid of all possible vaccine profiles
para_vax_combn <- data.table(expand.grid(
   vx_start_year = vx_start_year,
   coverageR = coverageR,
   coverageM = coverageM,
   effI = effI,
   effD =  effD,
   duration = duration,
   mass_interval = mass_interval,
   ageM = ageM,
   ageR = ageR,
   vxtype = vxtype)
)

para_vax_combn[, vax_hash := apply(.SD, MARGIN = 1, FUN = function(x) strtrim(digest(x), 3))]


# System Check ------------------------------------------------------------
# Check current host and store SGE_TASK_ID if appropriate
if (system("uname -n", intern = T) == "archer") {
   op            <- file.path(here("output/trajectories/"))
   grid_task_int <- 1
   #runs <- nrow(para_vax_combn)
   runs <- 2
   # cl <- makeCluster(8)
   # registerDoParallel(cl)
} else {
   op            <- file.path(Sys.getenv("OPPATH"))
   grid_task_int <- as.numeric(Sys.getenv("SGE_TASK_ID"))
   # Set DT threads for cluster - disable internal parallelism in data.tables to prevent cluster multicore use
   setDTthreads(1)
   runs <- nrow(para_vax_combn)
   # Set up parallel processing cluster
   # cl <- makeCluster(8)
   # registerDoParallel(cl)
}

# Set up parameter set and pick row
para_variable <- as.list(parameter_sets[grid_task_int, ])
para_cost     <- as.list(cost_sets[grid_task_int, ])
print(para_cost$ds_tx_cost)
print(para_cost$dr_tx_cost)


# Generate parameter hash
para_hash <- para_variable$para_hash
cost_hash <- para_cost$cost_hash

rm(para_vax)

# Baseline Runs - No Vaccine ----------------------------------------------
print(para_cost$ds_tx_cost)
print(para_cost$dr_tx_cost)

sq <- cmp_mainloop(
   data = data,
   mode = 3,
   para_static = para_static,
   para_variable = para_variable,
   para_cost = para_cost,
   startcon = startcon,
   vaccine = 2,
   para_vax = NULL,
   transmission = 1,
   mdrac = 1970 ,
   rx = 1,
   cost = 1,
   scen = 1
)

print(para_cost$ds_tx_cost)
print(para_cost$dr_tx_cost)

pol <- cmp_mainloop(
   data = data,
   mode = 3,
   para_static = para_static,
   para_variable = para_variable,
   para_cost = para_cost,
   startcon = startcon,
   vaccine = 2,
   para_vax = NULL,
   transmission = 1,
   mdrac = 1970 ,
   rx = 1,
   cost = 1,
   scen = 2
)

ep_list <- list()

# Vaccine Profiles Loop ---------------------------------------------------

for(vax in seq_len(runs)) {

   {
      library(data.table)
      library(digest)
      setDTthreads(1)

      # Select vaccine profile
      para_vax <- as.list(para_vax_combn[vax, ])
      para_vax$vx_regimens <- c(10, 30)

      # vax hash
      vax_hash <- para_vax$vax_hash

      print(para_cost$ds_tx_cost)
      print(para_cost$dr_tx_cost)

      # Run vaccine model
      vx_pol <- cmp_mainloop(
         data = data,
         mode = 3,
         para_static = para_static,
         para_variable = para_variable,
         para_cost = para_cost,
         startcon = startcon,
         vaccine = 1,
         para_vax = para_vax,
         transmission = 1,
         mdrac = 1970 ,
         rx = 1,
         cost = 1,
         scen = 2
      )

      print(para_cost$ds_tx_cost)
      print(para_cost$dr_tx_cost)

      vx_sq <- cmp_mainloop(
         data = data,
         mode = 3,
         para_static = para_static,
         para_variable = para_variable,
         para_cost = para_cost,
         startcon = startcon,
         vaccine = 1,
         para_vax = para_vax,
         transmission = 1,
         mdrac = 1970 ,
         rx = 1,
         cost = 1,
         scen = 1
      )

      print(para_cost$ds_tx_cost)
      print(para_cost$dr_tx_cost)

      #endpoint years
      ep_yrs                      <- c(2030, 2035, 2050)

      # holding list
      epi_sq                      <- list()
      cost_sq                     <- list()
      cost_sq$para_vax            <- para_vax
      cost_sq$para_vax$vx_regimens <- NULL
      cost_sq$diff_tbrx           <- cost_calc(bl.cost.df = sq$difflist$ann_cost_array, vx.cost.df = vx_sq$difflist$ann_cost_array, vn = 'tbrx')
      cost_sq$diff_ccdc_tbrx      <- cost_calc(bl.cost.df = sq$difflist$ann_cost_array, vx.cost.df = vx_sq$difflist$ann_cost_array, vn = 'ccdc_tbrx')
      cost_sq$vx_cost             <- vx_cost(ivxdf = vx_sq$vxa[,.SD , .SDcols=patterns('Year|costT')])

      cost_sq$bl_dstb_daly        <- daly_calc(sy = 2027, eyl = ep_yrs, mort_df = sq$difflist$ann_dstb_mort_raw, yld_df = sq$difflist$ann_dstb_dalycmp, t = "dstb_bl")
      cost_sq$vx_dstb_daly        <- daly_calc(sy = 2027, eyl = ep_yrs, mort_df = vx_sq$difflist$ann_dstb_mort_raw, yld_df = vx_sq$difflist$ann_dstb_daly, t = "dstb_vx")
      cost_sq$diff_dstb_daly      <- merge(cost_sq$bl_dstb_daly, cost_sq$vx_dstb_daly)[, dstb_diff_daly := dstb_bl_cumDALY - dstb_vx_cumDALY]
      cost_sq$bl_drtb_daly        <- daly_calc(sy = 2027, eyl = ep_yrs, mort_df = sq$difflist$ann_drtb_mort_raw, yld_df = sq$difflist$ann_drtb_daly, t = "drtb_bl")
      cost_sq$vx_drtb_daly        <- daly_calc(sy = 2027, eyl = ep_yrs, mort_df = vx_sq$difflist$ann_drtb_mort_raw, yld_df = vx_sq$difflist$ann_drtb_daly, t = "drtb_vx")
      cost_sq$diff_drtb_daly      <- merge(cost_sq$bl_drtb_daly, cost_sq$vx_drtb_daly)[, drtb_diff_daly := drtb_bl_cumDALY - drtb_vx_cumDALY]
      cost_sq$bl_daly             <- daly_calc(sy = 2027, eyl = ep_yrs, mort_df = sq$difflist$ann_alltb_mort_raw, yld_df = sq$difflist$ann_dalycmp, t = "bl")
      cost_sq$vx_daly             <- daly_calc(sy = 2027, eyl = ep_yrs, mort_df = vx_sq$difflist$ann_alltb_mort_raw, yld_df = vx_sq$difflist$ann_dalycmp, t = "vx")
      cost_sq$diff_daly           <- merge(cost_sq$bl_daly, cost_sq$vx_daly)[, diff_daly := bl_cumDALY - vx_cumDALY]

      cost_sq$master_cost         <- Reduce(f = merge, x = list(cost_sq$diff_tbrx, cost_sq$diff_ccdc_tbrx, cost_sq$vx_cost, cost_sq$diff_daly, cost_sq$diff_dstb_daly, cost_sq$diff_drtb_daly))
      cost_sq$master_cost[, c("total_cost", "ccdc_total_cost") := .(diff_tbrx + vaccine_cost, diff_ccdc_tbrx + vaccine_cost)]
      cost_sq$master_cost[, c("ICER" , "CCDC_ICER") := .(total_cost/diff_daly, ccdc_total_cost/diff_daly)]
      cost_sq$master_cost[, c("costT", "regimen") := tstrsplit(vx_profile, "_")]

      cost_sq$master_cost[ , c("country", "scen", "para_hash", "cost_hash") := .(country, "sq", para_hash, cost_hash)]
      cost_sq$master_cost         <- cbind(cost_sq$master_cost, as.data.table(cost_sq$para_vax))
      cost_sq$bimp                <- bimp(vxca = vx_sq$difflist$ann_cost_array, blca = sq$difflist$ann_cost_array)
      cost_sq$vbimp               <- vx_sq$vxa[, c("vax_hash", "cost_hash") := .(vax_hash, cost_hash)]

      epi_sq$diff_pe              <- diff_pe(bl = pe_adjust(pe = sq$point_estimates, t = "bl"), vx = pe_adjust(pe = vx_sq$point_estimates, t = "vx"))
      epi_sq$diff_alltb_inc       <- diff_epi_sum(bl = epi_sum(epidf = sq$difflist$ann_alltb_inc_raw, t = "bl"), vx = epi_sum(epidf = vx_sq$difflist$ann_alltb_inc_raw, t = "vx"))
      epi_sq$diff_alltb_mort      <- diff_epi_sum(bl = epi_sum(epidf = sq$difflist$ann_alltb_mort_raw, t = "bl"), vx = epi_sum(epidf = vx_sq$difflist$ann_alltb_mort_raw, t = "vx"))
      epi_sq$diff_drtb_inc        <- diff_epi_sum(bl = epi_sum(epidf = sq$difflist$ann_drtb_inc_raw, t = "bl"), vx = epi_sum(epidf = vx_sq$difflist$ann_drtb_inc_raw, t = "vx"))
      epi_sq$diff_drtb_mort       <- diff_epi_sum(bl = epi_sum(epidf = sq$difflist$ann_drtb_mort_raw, t = "bl"), vx = epi_sum(epidf = vx_sq$difflist$ann_drtb_mort_raw, t = "vx"))
      epi_sq$diff_initRx          <- diff_epi_sum(bl = epi_sum(epidf = sq$difflist$ann_drtb_mort_raw, t = "bl"), vx = epi_sum(epidf = vx_sq$difflist$ann_drtb_mort_raw, t = "vx"))
      epi_sq$master_epi           <- rbindlist(epi_sq)
      epi_sq$para_vax             <- para_vax
      epi_sq$para_vax$vx_regimens   <- NULL
      epi_sq$master_epi[ , c("country", "scen", "para_hash", "cost_hash") := .(country, "sq", para_hash, cost_hash)]
      epi_sq$diff_initrx          <- diff_initrx(bldf = sq$difflist$ann_initRx, vxdf = vx_sq$difflist$ann_initRx)[, vxp:= para_vax$vax_hash]
      epi_sq$master_epi           <- cbind(epi_sq$master_epi, as.data.table(epi_sq$para_vax))

      # holding list
      epi_pol                     <- list()
      cost_pol                    <- list()
      cost_pol$para_vax           <- para_vax
      cost_pol$para_vax$vx_regimens <- NULL
      cost_pol$diff_tbrx          <- cost_calc(bl.cost.df = pol$difflist$ann_cost_array, vx.cost.df = vx_pol$difflist$ann_cost_array, vn = 'tbrx')
      cost_pol$diff_ccdc_tbrx     <- cost_calc(bl.cost.df = pol$difflist$ann_cost_array, vx.cost.df = vx_pol$difflist$ann_cost_array, vn = 'ccdc_tbrx')
      cost_pol$vx_cost            <- vx_cost(ivxdf = vx_pol$vxa[,.SD , .SDcols=patterns('Year|costT')])

      cost_pol$bl_dstb_daly       <- daly_calc(sy = 2027, eyl = ep_yrs, mort_df = pol$difflist$ann_dstb_mort_raw, yld_df = pol$difflist$ann_dstb_dalycmp, t = "dstb_bl")
      cost_pol$vx_dstb_daly       <- daly_calc(sy = 2027, eyl = ep_yrs, mort_df = vx_pol$difflist$ann_dstb_mort_raw, yld_df = vx_pol$difflist$ann_dstb_daly, t = "dstb_vx")
      cost_pol$diff_dstb_daly     <- merge(cost_pol$bl_dstb_daly, cost_pol$vx_dstb_daly)[, dstb_diff_daly := dstb_bl_cumDALY - dstb_vx_cumDALY]
      cost_pol$bl_drtb_daly       <- daly_calc(sy = 2027, eyl = ep_yrs, mort_df = pol$difflist$ann_drtb_mort_raw, yld_df = pol$difflist$ann_drtb_daly, t = "drtb_bl")
      cost_pol$vx_drtb_daly       <- daly_calc(sy = 2027, eyl = ep_yrs, mort_df = vx_pol$difflist$ann_drtb_mort_raw, yld_df = vx_pol$difflist$ann_drtb_daly, t = "drtb_vx")
      cost_pol$diff_drtb_daly     <- merge(cost_pol$bl_drtb_daly, cost_pol$vx_drtb_daly)[, drtb_diff_daly := drtb_bl_cumDALY - drtb_vx_cumDALY]
      cost_pol$bl_daly            <- daly_calc(sy = 2027, eyl = ep_yrs, mort_df = pol$difflist$ann_alltb_mort_raw, yld_df = pol$difflist$ann_dalycmp, t = "bl")
      cost_pol$vx_daly            <- daly_calc(sy = 2027, eyl = ep_yrs, mort_df = vx_pol$difflist$ann_alltb_mort_raw, yld_df = vx_pol$difflist$ann_dalycmp, t = "vx")
      cost_pol$diff_daly          <- merge(cost_pol$bl_daly, cost_pol$vx_daly)[, diff_daly := bl_cumDALY - vx_cumDALY]

      cost_pol$master_cost        <- Reduce(f = merge, x = list(cost_pol$diff_tbrx, cost_pol$diff_ccdc_tbrx, cost_pol$vx_cost, cost_pol$diff_daly, cost_pol$diff_dstb_daly, cost_pol$diff_drtb_daly))
      cost_pol$master_cost[, c("total_cost", "ccdc_total_cost") := .(diff_tbrx + vaccine_cost, diff_ccdc_tbrx + vaccine_cost)]
      cost_pol$master_cost[, c("ICER" , "CCDC_ICER") := .(total_cost/diff_daly, ccdc_total_cost/diff_daly)]
      cost_pol$master_cost[, c("costT", "regimen") := tstrsplit(vx_profile, "_")]
      cost_pol$bimp               <- bimp(vxca = vx_pol$difflist$ann_cost_array, blca = pol$difflist$ann_cost_array)

      cost_pol$master_cost[ , c("country", "scen", "para_hash", "cost_hash") := .(country, "pol", para_hash, cost_hash)]
      cost_pol$master_cost        <- cbind(cost_pol$master_cost, as.data.table(cost_pol$para_vax))
      cost_pol$vbimp              <- vx_pol$vxa[, c("vax_hash", "cost_hash") := .(vax_hash, cost_hash)]

      epi_pol$diff_pe             <- diff_pe(bl = pe_adjust(pe = pol$point_estimates, t = "bl"), vx = pe_adjust(pe = vx_pol$point_estimates, t = "vx"))
      epi_pol$diff_alltb_inc      <- diff_epi_sum(bl = epi_sum(epidf = pol$difflist$ann_alltb_inc_raw, t = "bl"), vx = epi_sum(epidf = vx_pol$difflist$ann_alltb_inc_raw, t = "vx"))
      epi_pol$diff_alltb_mort     <- diff_epi_sum(bl = epi_sum(epidf = pol$difflist$ann_alltb_mort_raw, t = "bl"), vx = epi_sum(epidf = vx_pol$difflist$ann_alltb_mort_raw, t = "vx"))
      epi_pol$diff_drtb_inc       <- diff_epi_sum(bl = epi_sum(epidf = pol$difflist$ann_drtb_inc_raw, t = "bl"), vx = epi_sum(epidf = vx_pol$difflist$ann_drtb_inc_raw, t = "vx"))
      epi_pol$diff_drtb_mort      <- diff_epi_sum(bl = epi_sum(epidf = pol$difflist$ann_drtb_mort_raw, t = "bl"), vx = epi_sum(epidf = vx_pol$difflist$ann_drtb_mort_raw, t = "vx"))
      epi_pol$master_epi          <- rbindlist(epi_pol)
      epi_pol$para_vax            <- para_vax
      epi_pol$para_vax$vx_doses   <- NULL
      epi_pol$para_vax$vx_price   <- NULL
      epi_pol$diff_initrx         <- diff_initrx(bldf = pol$difflist$ann_initRx, vxdf = vx_pol$difflist$ann_initRx)[, vxp:= para_vax$vax_hash]
      epi_pol$master_epi[ , c("country", "scen", "para_hash", "cost_hash") := .(country, "pol", para_hash, cost_hash)]
      epi_pol$master_epi          <- cbind(epi_pol$master_epi, as.data.table(epi_sq$para_vax))

      op_epi   <- rbind(epi_sq$master_epi, epi_pol$master_epi)
      op_cost  <- rbind(cost_sq$master_cost, cost_pol$master_cost)
      op_bimp  <- rbind(cost_sq$bimp, cost_pol$bimp)
      op_vbimp <- rbind(cost_sq$vbimp, cost_pol$vbimp)
      op_initrx <- rbind(epi_sq$diff_initrx, epi_pol$diff_initrx)

      rm(epi_pol, epi_sq, cost_pol, cost_sq)
      gc(full = T)
      ep_list[[vax]] <- list(op_epi, op_cost, op_bimp, op_vbimp, op_initrx)
      # list(op_bimp, op_vbimp)
   }
}

# stopCluster(cl)

{
   epi   <- rbindlist(lapply(ep_list, function(x) x[[1]]))
   cost  <- rbindlist(lapply(ep_list, function(x) x[[2]]))
   bimp  <- rbindlist(lapply(ep_list, function(x) x[[3]]))
   vbimp <- rbindlist(lapply(ep_list, function(x) x[[4]]))
   initrx <- rbindlist(lapply(ep_list, function(x) x[[5]]))
}

# bimp  <- rbindlist(lapply(ep_list, function(x) x[[1]]))
# vbimp <- rbindlist(lapply(ep_list, function(x) x[[2]]))

{
   cost$vx_profile = NULL
   cost$costT = NULL
   cost$vx_start_year = NULL
   cost$ageM = NULL
   cost$ageR = NULL
   cost$effI = NULL
   cost$mass_interval = NULL
   cost$coverageM = NULL
   cost$coverageR = NULL
   epi$vx_start_year = NULL
   epi$ageM = NULL
   epi$ageR = NULL
   epi$effI = NULL
   epi$mass_interval = NULL
   epi$coverageM = NULL
   epi$coverageR = NULL
}

{
   fwrite(epi, file = file.path(op, paste0("epi_vxe_", para_hash, ".csv")))
   fwrite(cost, file = file.path(op, paste0("cost_vxe_", para_hash, ".csv")))
   fwrite(bimp, file = file.path(op, paste0("bimp_vxe_", para_hash, ".csv")))
   fwrite(vbimp, file = file.path(op, paste0("vbimp_vxe_", para_hash, ".csv")))
   fwrite(initrx, file = file.path(op, paste0("initrx_vxe_", para_hash, ".csv")))

}
