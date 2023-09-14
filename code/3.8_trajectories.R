library(here)
library(digest)
library(data.table)
library(logger)

# Set working directory to code directory
init_dir <- setwd(here("code"))

# Check current host and store SGE_TASK_ID if appropriate
if (system("uname -n", intern = T) == "archer") {
  op <- here("output", "trajectories", "CN", "singletons")
  grid_task_int <- 1
  ext_log_level <- as.name("INFO")
  shell_dt <- Sys.time()
} else {
  # Check datestamp of cluster shell controller file
  shell_dt <- system("date -u +%Y-%m-%d-%H%M -r $HOME/MDR-Aeras/mdrvx/code/3.9_traj_cluster_job.sh", intern = T)
  
  op <- here("output", "trajectories", "CN", shell_dt)
  try(if (!dir.exists(op)) {
    dir.create(path = op, recursive = TRUE)
  })
  grid_task_int <- as.numeric(Sys.getenv("SGE_TASK_ID"))
  ext_log_level <- as.name(Sys.getenv("MCMC_LOG_LEVEL"))
  # Set DT threads for cluster - disable internal parallelism in data.tables to prevent cluster multicore use
  setDTthreads(1)
}

# ---------------------
# Set up logging
log_formatter(formatter_sprintf)
logger <- layout_glue_generator(format = '[{time}]\t{level}\t{grid_task_int}\t{msg}')
log_layout(logger)
log_appender(appender_tee(file = here("output", "logs", "trajectories", sprintf("trajectories_CN_%s.log", grid_task_int))))
log_threshold(ext_log_level)

# Initial information
log_info("START TRAJECTORIES JOB")
log_info("SESSION INFO: %s", sessionInfo()[["R.version"]]$version.string)
log_info("OUTPUT PATH: %s", op)
log_info("SHELL TIMESTAMP: %s", shell_dt)
log_info("R-LIBPATH: %s", .libPaths())

# External data import
tryCatch(
  expr = {
    # Import relevant data for model run
    source("1_import.R")
    rm(para_variable)
    
    # Parameter sets
    parameter_sets <- fread("../data/fitted_parameters/april_2020_1K_china.csv")
    cost_sets      <- fread("../data/para_cost_sampled.csv")
    
    # Source compiled function
    cmp_mainloop <- readRDS("../code/2.1_main_rebuild.R.cmp")
    
  }, error = function(e) {
    log_error("Import: %s", e$message)
    stop(e$call)
  }
)

tryCatch(
  expr = {
    if (grid_task_int == 1 | nrow(parameter_sets) == 1) {
      para_variable <- as.list(parameter_sets[1, ])
      para_cost     <- as.list(cost_sets[1, ])
    } else if (grid_task_int <= nrow(parameter_sets)){
      para_variable <- as.list(parameter_sets[grid_task_int, ])
      para_cost     <- as.list(cost_sets[grid_task_int, ])
    } else if (grid_task_int > nrow(parameter_sets)) {
      para_variable <- as.list(parameter_sets[nrow(parameter_sets) %% grid_task_int, ])
      para_cost     <- as.list(cost_sets[nrow(parameter_sets) %% grid_task_int, ])
    }
    log_info("PARA_HASH: %s", para_variable$para_hash)
    log_info("PARA_COST: %s", para_cost$cost_hash)
  }, error = function(e) {
    try(expr = log_info("PARA_HASH: %s", para_variable$para_hash))
    try(expr = log_info("PARA_COST: %s", para_cost$cost_hash))
    log_error("Parameter set up: %s", e$message)
    stop(e$call)
  }
)

# Set common mainloop arguments
margs <- list()

# Starting conditions
margs$startcon <- startcon <- list(
  dt = 0.25,
  year1 = 1900,
  yearend = 2099,
  prop_col = 1,
  bgu = 1,
  fert = 1,
  opp = op
)

margs$mode <- 4
margs$cost <- 1

# Vaccine trajectories
log_debug("COMMON_ARGS: %s", margs)

ms <- list()
ms$aglo <- ms$aglp <- list(
  All = c(0:99),
  A09 = c(0:9),
  A1019 = c(10:19),
  A2029 = c(20:29),
  A3039 = c(30:39),
  A4049 = c(40:49),
  A5059 = c(50:59),
  A6069 = c(60:69),
  A7079 = c(70:79),
  A8089 = c(80:89),
  A9099 = c(90:99)
)

# Baseline trajectories
# Baseline trajectories
# Status Quo
tryCatch(
  expr = {
    results_sq <- cmp_mainloop(
      data = data,
      mode = margs$mode,
      para_static = para_static,
      para_variable = para_variable,
      para_cost = para_cost,
      startcon = startcon,
      vaccine = 0,
      para_vax = NULL,
      transmission = 1,
      mdrac = 1970,
      rx = 1,
      cost = margs$cost,
      scen = 1,
      modespec = ms
    )
    log_info("BASELINE - SQ - DONE")
  }, error = function(e) {
    log_error("BASELINE - SQ: %s", as.character(e$message))
  }
)

# Policy
tryCatch(
  expr = {
    results_pol <- cmp_mainloop(
      data = data,
      mode = margs$mode,
      para_static = para_static,
      para_variable = para_variable,
      para_cost = para_cost,
      startcon = startcon,
      vaccine = 0,
      para_vax = NULL,
      transmission = 1,
      mdrac = 1970,
      rx = 1,
      cost = margs$cost,
      scen = 2,
      modespec = ms
    )
    log_info("BASELINE - POLICY - DONE")
  }, error = function(e) {
    log_error("BASELINE - POLICY: %s", as.character(e$message))
  }
)

rm(results_sq)
rm(results_pol)
gc(full = T)

# Vaccine trajectories
# for (vxt in c(1,2,3)) {
#   
#   log_info("START VXTYPE %s", vxt)
#   para_vax$effD <- 0.5
#   para_vax$duration <- 10
#   para_vax$vxtype <- vxt
#   para_vax$vx_regimens <- 10
#   
#   tryCatch(
#     expr = {
#       vx_results_sq <- cmp_mainloop(
#         data = data,
#         mode = margs$mode,
#         para_static = para_static,
#         para_variable = para_variable,
#         para_cost = para_cost,
#         startcon = startcon,
#         vaccine = 1,
#         para_vax = para_vax,
#         transmission = 1,
#         mdrac = 1970,
#         rx = 1,
#         cost = margs$cost,
#         scen = 1
#       )
#       log_info("VXTYPE%s SQ RUN DONE", para_vax$vxtype)
#     }, error = function(e) {
#       log_error("VXTYPE %s: %s", para_vax$vxtype, e$message)
#       stop(e$call)
#     }
#   )
# 
#   tryCatch(
#     expr = {
#       vx_results_pol <- cmp_mainloop(
#         data = data,
#         mode = margs$mode,
#         para_static = para_static,
#         para_variable = para_variable,
#         para_cost = para_cost,
#         startcon = startcon,
#         vaccine = 1,
#         para_vax = para_vax,
#         transmission = 1,
#         mdrac = 1970,
#         rx = 1,
#         cost = margs$cost,
#         scen = 2
#       )
#       log_info("VXTYPE%s POL RUN DONE", para_vax$vxtype)
#     }, error = function(e) {
#       log_error("VXTYPE %s: %s", para_vax$vxtype, e$message)
#       stop(e$call)
#     }
#   )
#   
#   rm(vx_results_sq)
#   rm(vx_results_pol)
#   gc(full = T)
#   
#   log_info("START VXTYPE %s", vxt)
#   para_vax$effD <- 0.5
#   para_vax$duration <- 10
#   para_vax$vxtype <- vxt
#   para_vax$vx_regimens <- 10
#   para_vax$widthM <- 10
#
#   tryCatch(
#     expr = {
#       vx_results_sq <- cmp_mainloop(
#         data = data,
#         mode = margs$mode,
#         para_static = para_static,
#         para_variable = para_variable,
#         para_cost = para_cost,
#         startcon = startcon,
#         vaccine = 1,
#         para_vax = para_vax,
#         transmission = 1,
#         mdrac = 1970,
#         rx = 1,
#         cost = margs$cost,
#         scen = 1,
#         modespec = ms
#       )
#       log_info("VXTYPE%s SQ RUN DONE", para_vax$vxtype)
#     }, error = function(e) {
#       log_error("VXTYPE %s - SQ: %s", para_vax$vxtype, e$message)
#     }
#   )
#
#   tryCatch(
#     expr = {
#       vx_results_pol <- cmp_mainloop(
#         data = data,
#         mode = margs$mode,
#         para_static = para_static,
#         para_variable = para_variable,
#         para_cost = para_cost,
#         startcon = startcon,
#         vaccine = 1,
#         para_vax = para_vax,
#         transmission = 1,
#         mdrac = 1970,
#         rx = 1,
#         cost = margs$cost,
#         scen = 2,
#         modespec = ms
#       )
#       log_info("VXTYPE %s POL RUN DONE", para_vax$vxtype)
#     }, error = function(e) {
#       log_error("VXTYPE %s - POL: %s", para_vax$vxtype, e$message)
#     }
#   )
#
#   rm(vx_results_sq) 
#   rm(vx_results_pol)
#   gc(full = T)
#
# }
log_info("END OF SCRIPT")

