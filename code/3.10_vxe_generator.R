#!/usr/bin/Rscript --vanilla
# Second generation vaccine endpoint generator
start <- proc.time()
library(here)
library(digest)
library(data.table)
library(logger)
library(argparse, quietly = TRUE)

if (system("uname -n", intern = T) == "archer") {
  op <- here("output", "vxe", "CN", "local")
  try(if (!dir.exists(op)) {
    dir.create(path = op, recursive = TRUE)
  })
  grid_task_int <- 1
  shell_dt <- format(Sys.time(), format = "%Y-%m-%d-%H%M")
  ext_log_level <- as.name("INFO")
} else {
  shell_dt <- system("date -u +%Y-%m-%d-%H%M -r $HOME/MDR-Aeras/mdrvx/code/3.6_vx_runs_caller.sh", intern = T)
  op <- here("output", "vxe", "CN", shell_dt)
  try(if (!dir.exists(op)) {
    dir.create(path = op, recursive = TRUE)
  })
  grid_task_int <- as.numeric(Sys.getenv("SGE_TASK_ID"))
  ext_log_level <- as.name(Sys.getenv("MCMC_LOG_LEVEL"))
  
  # Set DT threads for cluster - disable internal parallelism in data.tables to prevent cluster multicore use
  setDTthreads(1)
}

log_formatter(formatter_sprintf)
logger <- layout_glue_generator(format = '[{time}]\t{grid_task_int}\t{level}\t{msg}')
log_layout(logger)
log_appender(appender_tee(file = here("output", "logs", "vxe", sprintf("vxe_%s.log", grid_task_int))))
log_info("START VXE_GENERATOR SCRIPT")
log_info("OUTPUT PATH: %s", op)

tryCatch(
  expr = {
    parser <- ArgumentParser()
    parser$add_argument("-u", type = "logical", default=FALSE)
    parser$add_argument("-f", type = "character")
    parser$add_argument("-t", type = "integer", default=10)
    args <- parser$parse_args()
    log_threshold(ext_log_level)
    log_info('COMMAND LINE OPTS PARSED')
  }, error = function(e) {
    log_error("COMMAND LINE OPTS ERROR: %s", e$message)
    stop(e$call)
  }
)

tryCatch(
  expr = {

    init_dir <- setwd(here("code"))
    log_debug("WD: %s", init_dir)

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
    log_debug("STARTCON: %s", startcon)
  }
)

tryCatch(
  expr = {
    # Import relevant data for model run
    source("1_import.R")
    rm(para_cost)
    rm(para_variable)

    country = "CN"
    # Source compiled function
    cmp_mainloop   <- readRDS("../code/2.1_main_rebuild.R.cmp")

    # Read in relevant external data required for baseline vs vaccine runs, cost calculations etc
    # Life expectancy table for YLL calculations
    ledf           <- melt(fread("../data/unpd_le_cleaned_CN_v2.csv"), id.vars = "Year", variable.factor = F, variable.name = "MGrp", value.name = "LE")
    setkey(ledf, Year, MGrp)

    # Subsampled fitted parameter sets (x1000) and costs
    parameter_sets <- fread("../data/fitted_parameters/april_2020_1K_china.csv")
    cost_sets      <- fread("../data/para_cost_sampled.csv")
  }, error = function(e) {
    log_error("Error importing external data: %s", e$message)
    stop(e$call)
  }
)

# Settable arguments

# DALY weights +/- range
dalywt    <- 0.333
log_trace("DALY WEIGHT: %s", dalywt)
# Discount Rate
disc_rate <- 0.03
log_trace("DISCOUNT RATE: %s", disc_rate)

# Set up possible vaccine profiles
vx_start_year <- 2027
coverageR     <- 0.8
coverageM     <- 0.7
effI          <- c(0)
effD          <- c(0.3, 0.5, 0.7, 0.9)
duration      <- c(5, 10)
mass_interval <- c(10)
ageM          <- 10
ageR          <- 9
vxtype        <- c(1, 2, 3)

# Create grid of all possible vaccine profiles
para_vax_combn <- data.table(expand.grid(
# para_vax_combn_master <- data.table(expand.grid(
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

# Create a subgrid for main manuscript vax endpoints
# para_vax_combn <- para_vax_combn_master[effD == 0.5 & duration == 10 &  vxtype != 2]

tryCatch(
  expr = {
    para_vax_combn[, vax_hash := apply(.SD, MARGIN = 1, FUN = function(x) strtrim(digest(x), 3))]
    npvx <- nrow(para_vax_combn)
    log_trace("VACCINE COMBINATIONS OK")
  }, error = function(e) {
    log_error("VACCINE COMBINATIONS: %s", e$message)
    stop(e$call)
  }
)

tryCatch(
  expr = {
    para_variable <- as.list(parameter_sets[grid_task_int, ])
    log_trace("PARA_VARIABLE: %s", para_variable)
    para_cost     <- as.list(cost_sets[grid_task_int, ])
    log_trace("PARA_COST: %s", para_cost)
    para_hash <- para_variable$para_hash
    cost_hash <- para_cost$cost_hash
    log_info("PARA_HASH: %s", para_hash)
    log_info("COST_HASH: %s", cost_hash)
  }, error = function(e) {
    log_error("ERROR SETTING PARA_VAR/PARA_COST: %s", e$message)
    stop(e$call)
  }
)

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

log_info("BASELINE - SQ - DONE")

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

log_info("BASELINE - POL - DONE")

# Settable analysis params
ep_years        <- c('2030', '2035', '2050')
cumul_output    <- c("ann_alltb_inc_raw", "ann_alltb_mort_raw", "ann_drtb_inc_raw", "ann_drtb_mort_raw")
retained_tables <- list(
  icer = list(),
  averted_burden = list(),
  averted_treatment = list(),
  pe_reduction = list()
)

for (vax in seq_len(nrow(para_vax_combn))) {
# for (vax in 1:2) {
  log_info("START VX SEQUENCES")
  # Set up vaccine
  para_vax <- as.list(para_vax_combn[vax, ])
  para_vax$vx_regimens <- c(10, 30)
  
  # vax hash
  vax_hash <- para_vax$vax_hash
  vax_char <- list(effD = para_vax$effD, vxtype = para_vax$vxtype, duration = para_vax$duration)
  log_debug("VAX CHAR: %s: %s", names(vax_char), vax_char)
  log_info("VX-SQ >> %s/%s >> %s >> START", vax, npvx, vax_hash)

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

  log_info("VX-SQ >> %s/%s >> %s >> CHECKPOINT: MODEL", vax, npvx, vax_hash)

  op_sq <- mapply(FUN = merge, x = sq$mdifflist, y = vx_sq$mdifflist, SIMPLIFY = FALSE, MoreArgs = list(suffixes = c(".bl", ".vx")))
  log_debug("VX-SQ >> %s/%s >> %s >> CHECKPOINT: MERGE OUTPUT OK", vax, npvx, vax_hash)
  
  lapply(op_sq, function(x) {
    invisible(x[, MVal.diff := MVal.vx - MVal.bl])
  })
  
  log_debug("VX-SQ >> %s/%s >> %s >> CHECKPOINT: DIFF OUTPUT OK", vax, npvx, vax_hash)
  
  # Interim calculations lists
  ic_sq <- list()

  # DALY calculations
  tryCatch(
    expr = {
      ic_sq$YLL       <- ledf[Year>=vx_start_year][op_sq$ann_alltb_mort_raw[Year>=vx_start_year], on = c("Year", "MGrp")]
      ic_sq$YLL       <- ic_sq$YLL[, YLL := MVal.diff * 1e3 * LE]
      ic_sq$YLL       <- ic_sq$YLL[, .(YLL = sum(YLL)), by = .(Year, scen, para_hash)]
      setkey(ic_sq$YLL, Year, scen, para_hash)
      ic_sq$YLD       <- op_sq$ann_dalycmp[Year>=vx_start_year][, YLD := MVal.diff * dalywt][, .(YLD = sum(YLD)), by = .(Year, scen, para_hash)]
      ic_sq$DALY      <- ic_sq$YLD[ic_sq$YLL][, DALY_nd := YLD + YLL][Year >= vx_start_year][, DALY_d := DALY_nd / (1 + disc_rate)^(.I-1)]

      # Cost calculations
      ic_sq$vx_cost   <- melt(vx_sq$vxa, id.vars = c("Year", "scen", "para_hash"), variable.name = "MGrp", value.name = "MVal")
      setkey(ic_sq$vx_cost, Year, scen, para_hash)
      ic_sq$merged_cost      <- merge(op_sq$ann_cost_array[Year>=vx_start_year], ic_sq$vx_cost, by = c("Year", "scen", "para_hash"), suffixes = c(".bl", ".vx"), allow.cartesian = T)
      ic_sq$cost      <- ic_sq$merged_cost[MGrp.bl == "ccdc_tbrx"][grep(pattern = "costT.*", x = MGrp.vx)]
      ic_sq$cost[, cost_nd := MVal + MVal.diff][, .(Year, scen, para_hash, MGrp.vx, cost_nd)]
      ic_sq$cost      <- ic_sq$cost[MGrp.vx == 'costT_USD10', cost_d := cost_nd/(1 + disc_rate)^(.I-1)]
      ic_sq$cost      <- ic_sq$cost[MGrp.vx == 'costT_USD30', cost_d := cost_nd/(1 + disc_rate)^(.I-1)]

      # ICER calculations
      ic_sq$icer_calc <- merge(x = ic_sq$cost, y = ic_sq$DALY[, !c('YLD', 'YLL')])
      ic_sq$icer      <- rbindlist(lapply(ep_years, function(x) {
        ic_sq$icer_calc[Year %between% c(vx_start_year, x), ][, .(cost_nd = sum(cost_nd), cost_d = sum(cost_d), DALY_nd = sum(DALY_nd), DALY_d = sum(DALY_d)), by = .(scen, para_hash, MGrp.vx)][, TH:= ..x]
      }))
      ic_sq$icer[, c("icer_nd", "icer_d") := .(-cost_nd/DALY_nd, -cost_d/DALY_d)]
      ic_sq$icer <- melt(ic_sq$icer, id.vars = c("scen", "para_hash", "MGrp.vx", "TH"), measure.vars = c("icer_d", "icer_nd"), variable.name = "MGrp", value.name = "MVal")
      log_debug("VX-SQ >> %s/%s >> %s >> CHECKPOINT: HE SEQUENCE OK", vax, npvx, vax_hash)
    }, error = function(e) {
      log_error("VX-SQ >> %s/%s >> %s >> HE SEQUENCE: %s", vax, npvx, vax_hash, e$message)
    }
  )

  # Point estimate calculations
  tryCatch(
    expr = {
      ic_sq$pe_reduction <- merge(x = sq$point_estimates, y = vx_sq$point_estimates, suffixes = c(".bl", ".vx"))[, MVal.diff := (MVal.bl - MVal.vx)/MVal.bl]
      log_debug("VX-SQ >> %s/%s >> %s >> CHECKPOINT: POINT EST SEQ OK", vax, npvx, vax_hash)
    }, error = function(e) {
      log_error("VX-SQ >> %s/%s >> %s >> POINT EST SEQ: %s", vax, npvx, vax_hash, e$message)
    }
  )

  # Summed averted cases and deaths
  tryCatch(
    expr = {
      ic_sq$averted_burden <- rbindlist(
        lapply(ep_years, function(y) {
          rbindlist(
            lapply(cumul_output, function(x) {
              op_sq[[x]][Year %between% c(vx_start_year, y)][, .(MVal = sum(MVal.diff)), by = .(scen, para_hash)][, c("endpoint", "endpoint_year") := .(..x, ..y)]
            }))
        }))
      log_debug("VX-SQ >> %s/%s >> %s >> CHECKPOINT: AVERT BURD SEQ OK", vax, npvx, vax_hash)
    }, error = function(e) {
      log_debug("VX-SQ >> %s/%s >> %s >> AVERT BURD SEQ: %s", vax, npvx, vax_hash, e$message)
    }
  )

  # Averted treatment
  tryCatch(
    expr = {
      ic_sq$averted_treatment <- rbindlist(lapply(ep_years, function(x) {
        op_sq$ann_initRx[Year %between% c(vx_start_year, x), .(MVal = sum(MVal.diff)), by = .(scen, para_hash, MGrp)][, endpoint_year := ..x]
      }))
      log_debug("VX-SQ >> %s/%s >> %s >> CHECKPOINT: AVERT TREAT SEQ OK", vax, npvx, vax_hash)
    }, error = function(e) {
      log_debug("VX-SQ >> %s/%s >> %s >> AVERT TREAT SEQ: %s", vax, npvx, vax_hash, e$message)
    }
  )

  log_info("VX-POL >> %s/%s >> %s >> START", vax, npvx, vax_hash)
  
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
  
  log_info("VX-POL >> %s/%s >> %s >> CHECKPOINT: MODEL", vax, npvx, vax_hash)
  
  op_pol <- mapply(FUN = merge, x = pol$mdifflist, y = vx_pol$mdifflist, SIMPLIFY = FALSE, MoreArgs = list(suffixes = c(".bl", ".vx")))
  log_debug("VX-POL >> %s/%s >> %s >> CHECKPOINT: MERGE OUTPUT OK", vax, npvx, vax_hash)
  
  lapply(op_pol, function(x) {
    invisible(x[, MVal.diff := MVal.vx - MVal.bl])
  })

  log_debug("VX-POL >> %s/%s >> %s >> CHECKPOINT: DIFF OUTPUT OK", vax, npvx, vax_hash)
  
  # Interim calculations lists
  ic_pol <- list()

  # DALY calculations
  tryCatch(
    expr = {
      ic_pol$YLL       <- ledf[Year>=vx_start_year][op_pol$ann_alltb_mort_raw[Year>=vx_start_year], on = c("Year", "MGrp")]
      ic_pol$YLL       <- ic_pol$YLL[, YLL := MVal.diff * 1e3 * LE]
      ic_pol$YLL       <- ic_pol$YLL[, .(YLL = sum(YLL)), by = .(Year, scen, para_hash)]
      setkey(ic_pol$YLL, Year, scen, para_hash)
      ic_pol$YLD       <- op_pol$ann_dalycmp[Year>=vx_start_year][, YLD := MVal.diff * dalywt][, .(YLD = sum(YLD)), by = .(Year, scen, para_hash)]
      ic_pol$DALY      <- ic_pol$YLD[ic_pol$YLL][, DALY_nd := YLD + YLL][Year >= vx_start_year][, DALY_d := DALY_nd / (1 + disc_rate)^(.I-1)]

      # Cost calculations
      ic_pol$vx_cost   <- melt(vx_pol$vxa, id.vars = c("Year", "scen", "para_hash"), variable.name = "MGrp", value.name = "MVal")
      setkey(ic_pol$vx_cost, Year, scen, para_hash)
      ic_pol$merged_cost      <- merge(op_pol$ann_cost_array[Year>=vx_start_year], ic_pol$vx_cost, by = c("Year", "scen", "para_hash"), suffixes = c(".bl", ".vx"), allow.cartesian = T)
      ic_pol$cost      <- ic_pol$merged_cost[MGrp.bl == "ccdc_tbrx"][grep(pattern = "costT.*", x = MGrp.vx)]
      ic_pol$cost[, cost_nd := MVal + MVal.diff][, .(Year, scen, para_hash, MGrp.vx, cost_nd)]
      ic_pol$cost      <- ic_pol$cost[MGrp.vx == 'costT_USD10', cost_d := cost_nd/(1 + disc_rate)^(.I-1)]
      ic_pol$cost      <- ic_pol$cost[MGrp.vx == 'costT_USD30', cost_d := cost_nd/(1 + disc_rate)^(.I-1)]


      # ICER calculations
      ic_pol$icer_calc <- merge(x = ic_pol$cost, y = ic_pol$DALY[, !c('YLD', 'YLL')])
      ic_pol$icer      <- rbindlist(lapply(ep_years, function(x) {
        ic_pol$icer_calc[Year %between% c(vx_start_year, x), ][, .(cost_nd = sum(cost_nd), cost_d = sum(cost_d), DALY_nd = sum(DALY_nd), DALY_d = sum(DALY_d)), by = .(scen, para_hash, MGrp.vx)][, TH:= ..x]
      }))
      ic_pol$icer[, c("icer_nd", "icer_d") := .(-cost_nd/DALY_nd, -cost_d/DALY_d)]
      ic_pol$icer <- melt(ic_pol$icer, id.vars = c("scen", "para_hash", "MGrp.vx", "TH"), measure.vars = c("icer_d", "icer_nd"), variable.name = "MGrp", value.name = "MVal")
      log_debug("VX-POL >> %s/%s >> %s >> CHECKPOINT: HE SEQUENCE OK", vax, npvx, vax_hash)
    }, error = function(e) {
      log_error("VX-POL >> %s/%s >> %s >> HE SEQUENCE: %s", vax, npvx, vax_hash, e$message)
    }
  )

  # Point estimate calculations
  tryCatch(
    expr = {
      ic_pol$pe_reduction <- merge(x = pol$point_estimates, y = vx_pol$point_estimates, suffixes = c(".bl", ".vx"))[, MVal.diff := (MVal.bl - MVal.vx)/MVal.bl]
      log_debug("VX-POL >> %s/%s >> %s >> CHECKPOINT: POINT EST SEQ OK", vax, npvx, vax_hash)
    }, error = function(e) {
      log_error("VX-POL >> %s/%s >> %s >> POINT EST SEQ: %s", vax, npvx, vax_hash, e$message)
    }
  )

  # Summed averted cases and deaths
  tryCatch(
    expr = {
      ic_pol$averted_burden <- rbindlist(
        lapply(ep_years, function(y) {
          rbindlist(
            lapply(cumul_output, function(x) {
              op_pol[[x]][Year %between% c(vx_start_year, y)][, .(MVal = sum(MVal.diff)), by = .(scen, para_hash)][, c("endpoint", "endpoint_year") := .(..x, ..y)]
            }))
        }))
      log_debug("VX-POL >> %s/%s >> %s >> CHECKPOINT: AVERT BURD SEQ OK", vax, npvx, vax_hash)
    }, error = function(e) {
      log_debug("VX-POL >> %s/%s >> %s >> AVERT BURD SEQ: %s", vax, npvx, vax_hash, e$message)
    }
  )

  # Averted treatment
  tryCatch(
    expr = {
      ic_pol$averted_treatment <- rbindlist(lapply(ep_years, function(x) {
        op_pol$ann_initRx[Year %between% c(vx_start_year, x), .(MVal = sum(MVal.diff)), by = .(scen, para_hash, MGrp)][, endpoint_year := ..x]
      }))
      log_debug("VX-POL >> %s/%s >> %s >> CHECKPOINT: AVERT TREAT SEQ OK", vax, npvx, vax_hash)
    }, error = function(e) {
      log_debug("VX-POL >> %s/%s >> %s >> AVERT TREAT SEQ: %s", vax, npvx, vax_hash, e$message)
    }
  )
  
  # Bind SQ and POL results together and assign out of loop
  lapply(names(retained_tables), function(x) {
    retained_tables[[x]][[vax]] <<- rbind(ic_sq[[x]], ic_pol[[x]])[, c(names(vax_char)) := .(vax_char$effD, vax_char$vxtype, vax_char$duration)]
  })

}

log_info("VACCINE SIMULATION COMPLETE")

# Finalise retained tables
tryCatch(
  expr = {
    lapply(names(retained_tables), function(x) {
      f_retained_table <- rbindlist(retained_tables[[x]])
      invisible(fwrite(f_retained_table, file = file.path(op, sprintf("%s_%s.csv", x, para_hash))))
      log_debug("WRITE OUT: %s", x)
    })
  }, error = function(e) {
    log_error("WRITE OUT ERROR: %s", e$message)
  }
)

end <- proc.time()
el <- end - start
log_info("ELAPSED: %s seconds", el["elapsed"])
log_info("END OF SCRIPT")
