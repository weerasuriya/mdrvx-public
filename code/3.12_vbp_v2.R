#!/usr/bin/Rscript --vanilla
# Second generation vaccine endpoint generator
start <- proc.time()
library(here)
library(digest)
library(data.table)
library(logger)
library(argparse, quietly = TRUE)

if (system("uname -n", intern = T) == "archer") {
  op <- here("output", "vbp", "CN", "local")
  try(if (!dir.exists(op)) {
    dir.create(path = op, recursive = TRUE)
  })
  grid_task_int <- 10
  shell_dt <- format(Sys.time(), format = "%Y-%m-%d-%H%M")
  ext_log_level <- as.name("INFO")
} else {
  shell_dt <- system("date -u +%Y-%m-%d-%H%M -r $HOME/MDR-Aeras/mdrvx/code/3.14_vbp_caller.sh", intern = T)
  op <- here("output", "vbp", "CN", shell_dt)
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
log_appender(appender_tee(file = here("output", "logs", "vbp", sprintf("vbp_%s.log", grid_task_int))))
log_info("START VBP_GENERATOR SCRIPT")
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
coverageR     <- seq(0, 1, 0.05)
coverageM     <- seq(0, 1, 0.05)
effI          <- c(0)
effD          <- c(0.5)
duration      <- c(10)
mass_interval <- 10
ageM          <- 60
ageR          <- 9
vxtype        <- c(3)
widthM        <- c(10)

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
  vxtype = vxtype,
  widthM = widthM)
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
    try(expr = log_info("PARA_HASH: %s", para_hash))
    try(expr = log_info("COST_HASH: %s", cost_hash))
    log_error("ERROR SETTING PARA_VAR/PARA_COST: %s", e$message)
    stop(e$call)
  }
)

tryCatch(
  expr = {
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
  },
  error = function(e) {
    log_error("BASELINE - SQ: %s", e$message)
  }
)

tryCatch(
  expr = {
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
  },
  error = function(e) {
    log_error("BASELINE - POL: %s", e$message)
  }
)

# Settable analysis params
ep_years        <- c('2050')
cumul_output    <- c("ann_alltb_inc_raw", "ann_alltb_mort_raw", "ann_drtb_inc_raw", "ann_drtb_mort_raw")
retained_tables <- list(
  vbp = list(),
  averted_burden = list(),
  averted_treatment = list(),
  pe_reduction = list()
)

# Comparison year
comp_yr <- 2018

for (vax in seq_len(nrow(para_vax_combn))) {
  # for (vax in 1:2) {
  log_info("START VX SEQUENCES")
  # Set up vaccine
  para_vax <- as.list(para_vax_combn[vax, ])
  para_vax$vx_regimens <- c(1,2)
  
  # vax hash
  vax_hash <- para_vax$vax_hash
  
  vax_char <- list(effD = para_vax$effD,
                   vxtype = para_vax$vxtype,
                   duration = para_vax$duration,
                   ageR = para_vax$ageR,
                   ageM = para_vax$ageM,
                   widthM = para_vax$widthM,
                   covM = para_vax$coverageM,
                   covR = para_vax$coverageR,
                   country = "China")
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
    # MVALDIFF is VACCINE - BASELINE
    invisible(x[, MVal.diff := MVal.vx - MVal.bl])
  })
  
  log_debug("VX-SQ >> %s/%s >> %s >> CHECKPOINT: DIFF OUTPUT OK", vax, npvx, vax_hash)
  
  # Interim calculations lists
  ic_sq <- list()
  tmp_sq <- list()
  
  # VBP calculations
  tryCatch(
    expr = {
      # DALY calculations
      ic_sq$YLL       <- ledf[Year>=comp_yr][op_sq$ann_alltb_mort_raw[Year>=comp_yr], on = c("Year", "MGrp")]
      ic_sq$YLL       <- ic_sq$YLL[, YLL := MVal.diff * 1e3 * LE]
      ic_sq$YLL       <- ic_sq$YLL[, .(YLL = sum(YLL)), by = .(Year, scen, para_hash)]
      setkey(ic_sq$YLL, Year, scen, para_hash)
      ic_sq$YLD       <- op_sq$ann_dalycmp[Year>=comp_yr][, YLD := MVal.diff * dalywt][, .(YLD = sum(YLD)), by = .(Year, scen, para_hash)]
      ic_sq$DALY      <- ic_sq$YLD[ic_sq$YLL][, DALY_nd := YLD + YLL][Year >= comp_yr][, DALY_d := DALY_nd / (1 + disc_rate)^(.I-1)]
      ic_sq$aDALY     <- ic_sq$DALY[, .(Year, scen, para_hash, aDALY_d = -DALY_d, aDALY_nd = -DALY_nd)]
      
      # Total costs
      # Should be: Baseline TB cost - Vaccine TB cost - Vaccine Delivery Cost
      # Vaccine TB cost - Baseline TB cost calculated as MVal.Diff - needs NEGATION
      
      # Delivery costs should be calculated in model
      tmp_sq$vx_cost       <- melt(vx_sq$vxa, id.vars = c("Year", "scen", "para_hash"), variable.name = "MGrp", value.name = "MVal")
      ic_sq$vx_cost        <- tmp_sq$vx_cost[grep(pattern = "deliv_costT_USD1", x = MGrp)][Year>=comp_yr]
      setkey(ic_sq$vx_cost, Year, scen, para_hash)
      
      ic_sq$prog_cost      <- op_sq$ann_cost_array[Year>=comp_yr & MGrp == 'ccdc_tbrx']
      setkey(ic_sq$prog_cost, Year, scen, para_hash)
      
      tmp_sq$merged_cost   <- merge(ic_sq$prog_cost, ic_sq$vx_cost, suffixes = c(".prog", ".vx"), by = c("Year", "scen", "para_hash"))
      # Calculate total cost
      tmp_sq$cost          <- tmp_sq$merged_cost[, .(Year = Year, scen = scen, para_hash = para_hash, totalcost_nd = (-MVal.diff)-MVal)]
      tmp_sq$cost[, totalcost_d := totalcost_nd/(1 + disc_rate)^(.I-1)]
      tmp_sq$costbenefit   <- merge(tmp_sq$cost, ic_sq$aDALY, by = c("Year", "scen", "para_hash"))
      
      # Number of vaccinations
      tmp_sq$dd_cost        <- vx_sq$vxa[Year>=comp_yr, .(Year, scen, para_hash, deliv_costR_USD1, deliv_costM_USD1)]
      tmp_sq$dd_cost[, c("delivR_d", "delivM_d") := .(deliv_costR_USD1/(1 + disc_rate)^(.I-1), deliv_costM_USD1/(1 + disc_rate)^(.I-1))]
      tmp_sq$dd_cost <- tmp_sq$dd_cost[, .(delivR_nd = sum(deliv_costR_USD1),
                                           delivR_d = sum(delivR_d),
                                           delivM_nd = sum(deliv_costM_USD1),
                                           delivM_d = sum(delivM_d)), by = .(scen, para_hash)]
      
      # Number of vaccinations
      setkey(vx_sq$ann_imms, Year, scen, para_hash)
      tmp_sq$vaccinations <- vx_sq$ann_imms[Year>=comp_yr, .(Year, scen, para_hash, vaccinations = Total, vaccinationsR = Routine, vaccinationsM = Mass)]
      tmp_sq$vbp_long <- merge(tmp_sq$costbenefit, tmp_sq$vaccinations, by = c('Year', 'scen', 'para_hash'))
      tmp_sq$vbp <- tmp_sq$vbp_long[, .(totalcost_d = sum(totalcost_d),
                                        totalcost_nd = sum(totalcost_nd),
                                        aDALY_d = sum(aDALY_d),
                                        aDALY_nd = sum(aDALY_nd),
                                        vaccinations = sum(vaccinations),
                                        vaccinationsR = sum(vaccinationsR),
                                        vaccinationsM = sum(vaccinationsM)), by = .(para_hash, scen)]
      
      ic_sq$vbp <- merge(tmp_sq$vbp, tmp_sq$dd_cost)
      
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
  tmp_pol <- list()
  
  # VBP calculations
  tryCatch(
    expr = {
      # DALY calculations
      ic_pol$YLL       <- ledf[Year>=comp_yr][op_pol$ann_alltb_mort_raw[Year>=comp_yr], on = c("Year", "MGrp")]
      ic_pol$YLL       <- ic_pol$YLL[, YLL := MVal.diff * 1e3 * LE]
      ic_pol$YLL       <- ic_pol$YLL[, .(YLL = sum(YLL)), by = .(Year, scen, para_hash)]
      setkey(ic_pol$YLL, Year, scen, para_hash)
      ic_pol$YLD       <- op_pol$ann_dalycmp[Year>=comp_yr][, YLD := MVal.diff * dalywt][, .(YLD = sum(YLD)), by = .(Year, scen, para_hash)]
      ic_pol$DALY      <- ic_pol$YLD[ic_pol$YLL][, DALY_nd := YLD + YLL][Year >= comp_yr][, DALY_d := DALY_nd / (1 + disc_rate)^(.I-1)]
      ic_pol$aDALY     <- ic_pol$DALY[, .(Year, scen, para_hash, aDALY_d = -DALY_d, aDALY_nd = -DALY_nd)]
      
      # Total costs
      # Should be: Baseline TB cost - Vaccine TB cost - Vaccine Delivery Cost
      # Vaccine TB cost - Baseline TB cost calculated as MVal.Diff - needs NEGATION
      
      # Delivery costs should be calculated in model
      tmp_pol$vx_cost       <- melt(vx_pol$vxa, id.vars = c("Year", "scen", "para_hash"), variable.name = "MGrp", value.name = "MVal")
      ic_pol$vx_cost        <- tmp_pol$vx_cost[grep(pattern = "deliv_costT_USD1", x = MGrp)][Year>=comp_yr]
      setkey(ic_pol$vx_cost, Year, scen, para_hash)
      
      ic_pol$prog_cost      <- op_pol$ann_cost_array[Year>=comp_yr & MGrp == 'ccdc_tbrx']
      setkey(ic_pol$prog_cost, Year, scen, para_hash)
      
      tmp_pol$merged_cost   <- merge(ic_pol$prog_cost, ic_pol$vx_cost, suffixes = c(".prog", ".vx"), by = c("Year", "scen", "para_hash"))
      # Calculate total cost
      tmp_pol$cost          <- tmp_pol$merged_cost[, .(Year = Year, scen = scen, para_hash = para_hash, totalcost_nd = (-MVal.diff)-MVal)]
      tmp_pol$cost[, totalcost_d := totalcost_nd/(1 + disc_rate)^(.I-1)]
      tmp_pol$costbenefit   <- merge(tmp_pol$cost, ic_pol$aDALY, by = c("Year", "scen", "para_hash"))
      
      # Number of vaccinations
      tmp_pol$dd_cost        <- vx_pol$vxa[Year>=comp_yr, .(Year, scen, para_hash, deliv_costR_USD1, deliv_costM_USD1)]
      tmp_pol$dd_cost[, c("delivR_d", "delivM_d") := .(deliv_costR_USD1/(1 + disc_rate)^(.I-1), deliv_costM_USD1/(1 + disc_rate)^(.I-1))]
      tmp_pol$dd_cost <- tmp_pol$dd_cost[, .(delivR_nd = sum(deliv_costR_USD1),
                                             delivR_d = sum(delivR_d),
                                             delivM_nd = sum(deliv_costM_USD1),
                                             delivM_d = sum(delivM_d)), by = .(scen, para_hash)]
      
      # Number of vaccinations
      setkey(vx_pol$ann_imms, Year, scen, para_hash)
      tmp_pol$vaccinations <- vx_pol$ann_imms[Year>=comp_yr, .(Year, scen, para_hash, vaccinations = Total, vaccinationsR = Routine, vaccinationsM = Mass)]
      tmp_pol$vbp_long <- merge(tmp_pol$costbenefit, tmp_pol$vaccinations, by = c('Year', 'scen', 'para_hash'))
      tmp_pol$vbp <- tmp_pol$vbp_long[, .(totalcost_d = sum(totalcost_d),
                                          totalcost_nd = sum(totalcost_nd),
                                          aDALY_d = sum(aDALY_d),
                                          aDALY_nd = sum(aDALY_nd),
                                          vaccinations = sum(vaccinations),
                                          vaccinationsR = sum(vaccinationsR),
                                          vaccinationsM = sum(vaccinationsM)), by = .(para_hash, scen)]
      
      ic_pol$vbp <- merge(tmp_pol$vbp, tmp_pol$dd_cost)
      
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
    retained_tables[[x]][[vax]] <<- rbind(ic_sq[[x]], ic_pol[[x]])[, c(names(vax_char)) := vax_char]
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
