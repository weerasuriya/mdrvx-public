# **Manual** Rejection-sampler Controller script for cluster
# Each run of this script is an SGE task

# Clear all
rm(list = ls())

# Libraries
library(here)
library(logger)
library(data.table)
library(jsonlite)

# Set working directory to code directory
setwd(here("code"))

# Check datestamp of cluster shell controller file
shell_dt <- system("date -u +%Y-%m-%d-%H%M -r 3.4_MCMC_cluster_job.sh", intern = T)

# Check current host and store SGE_TASK_ID if appropriate
if (system("uname -n", intern = T) == "archer") {
  op <- here("output", "mcmc_sampler", "CN")
  meta_folder <- here("output", "mcmc_sampler", "CN")
  grid_task_int <- 1
  ext_log_level <- as.name("INFO")
  easyabc_libloc <- NULL
} else {
  op <- here("output", "mcmc_sampler", "CN", shell_dt)
  try(if (!dir.exists(op)) {
    dir.create(path = op, recursive = TRUE)
  })
  meta_folder <- here("output", "mcmc_sampler", "CN", shell_dt)
  grid_task_int <- as.numeric(Sys.getenv("SGE_TASK_ID"))
  ext_log_level <- as.name(Sys.getenv("MCMC_LOG_LEVEL"))
  easyabc_libloc <- "/home/lsh1604836/R/x86_64-library/3.5"
  # Set DT threads for cluster - disable internal parallelism in data.tables to prevent cluster multicore use
  setDTthreads(1)
}

# Set up logging
log_formatter(formatter_sprintf)
logger <- layout_glue_generator(format = '[{time}]\t{level}\t{grid_task_int}\t{msg}')
log_layout(logger)
log_appender(appender_tee(file = here("output", "logs", "mcmc", sprintf("mcmc_CN_%s.log", grid_task_int))))
log_threshold(ext_log_level)

# Initial information
log_info("SESSION INFO: %s", sessionInfo()[["R.version"]]$version.string)
log_info("OUTPUT PATH: %s", op)
log_info("SHELL TIMESTAMP: %s", shell_dt)
log_info("R-LIBPATH: %s", .libPaths())

# Secondary Libraries
tryCatch(
  expr = {
    library(digest)
    library(EasyABC, lib.loc = easyabc_libloc)
    log_info("EASYABC LOCATION: %s", easyabc_libloc)
    log_info("EASYABC VERSION: %s", packageVersion("EasyABC"))
    library(coda)
  }, error = function(e) {
    log_error("SECONDARY LIBRARY LOAD FAILURE: %s", e)
    stop()
  }
)

tryCatch(
  expr = {
    # Import data
    source ("1_import.R", local = TRUE)
    # Read in compiled mainloop function, assign to .GlobalEnv
    cmp_mainloop <- readRDS("2.1_main_rebuild.R.cmp")
    log_debug("DATA AND FUNCTION IMPORT")
  }, error = function(e) {
    log_error("DATA AND FUNCTION IMPORT: %s", e)
    stop()
  }
)

# Starting conditions
startcon = list(
  dt = 0.25,
  year1 = 1900,
  yearend = 2099,
  prop_col = 1,
  bgu = 1,
  fert = 1
)

# Other switches
mode         <- 1
vaccine      <- 0
transmission <- 1
rx           <- 1
mdrac        <- 1970

tryCatch(
  expr = {
    # Remove para_variable
    pvnames <- names(para_variable)
    rm(para_variable)
    log_debug("PVNAMES: %s", pvnames)
    # Import fit targets
    target_ranges <- as.list(split(fread("../data/fit_targets/fit_targets.csv", header = TRUE), by = "Parameter"))
    log_debug("TARGET RANGE: %s", target_ranges)
    target_names <- names(target_ranges)
    log_debug("TARGET NAME: %s" , target_names)
    # Prior range import
    para_ranges <- fread("../data/para_ranges.csv")
    log_debug("PARA RANGE: %s", para_ranges)
    #Prior ranges
    prior_ranges <- lapply(split(para_ranges, f = para_ranges$Name), function(x) c("unif", x$low, x$high))
    prop_frac <- 1/50
    prop_frac_small <- 1/100
    prop_frac_ultra <- 1/200
    log_info("PROP_FRAC: %s", prop_frac)
    log_info("PROP_FRAC_SMALL: %s", prop_frac_small)
    log_info("PROP_FRAC_ULTRA: %s", prop_frac_ultra)
    prop_ranges <- vapply(split(para_ranges, f = para_ranges$Name), function(x) x$high-x$low, FUN.VALUE = numeric(1)) * prop_frac
    prop_ranges['CDRscale'] <- para_ranges[Name == "CDRscale", high - low] * prop_frac_small
    prop_ranges['chr'] <- para_ranges[Name == "chr", high - low] * prop_frac_small
    prop_ranges['uiscaleC'] <- para_ranges[Name == "uiscaleC", high - low] * prop_frac_ultra
    prop_ranges['uiscaleE'] <- para_ranges[Name == "uiscaleE", high - low] * prop_frac_ultra
    prop_ranges['uiscaleA'] <- para_ranges[Name == "uiscaleA", high - low] * prop_frac_ultra
    prop_ranges['velderly'] <- para_ranges[Name == "velderly", high - low] * prop_frac_small
    log_debug("PROP RANGE: %s", prop_ranges)
    log_debug("PRIOR RANGE: %s", prior_ranges)
  }, error = function(e) {
    log_error("TARGET / PARAM / PRIOR RANGE: %s", e)
    stop()
  }
)

# Set variables
# Chain length (abc_n), write out interval (wo_interval) and target
abc_n       <- 5e5
log_info("CHAIN LENGTH: %s", abc_n)
wo_interval <- if (system("uname -n", intern = T) == "archer") {10} else {100}
log_info("WO_INT: %s", wo_interval)
target      <- 26
log_info("TARGET: %s", target)
seed_search = TRUE
log_info("SEED SEARCH SETTING: %s", seed_search)

# Read seed values and set initial parameter
tryCatch(
  expr = {
    seed_file <- "seed_search_26.csv"
    seed_path <- here("data", "mcmc_seeds", "China", seed_file)
    param_seeds <- fread(seed_path, header = TRUE)
    log_info("SEED IMPORT: %s", seed_file)
    if (grid_task_int == 1 | nrow(param_seeds) == 1) {
      init_param <- param_seeds[1, ]
    } else if (grid_task_int <= nrow(param_seeds)){
      init_param <- param_seeds[grid_task_int, ]
    } else if (grid_task_int > nrow(param_seeds)) {
      init_param <- param_seeds[nrow(param_seeds) %% grid_task_int, ]
    }
    log_debug("INIT PARAM: %s", init_param)
    log_info("INIT PARAM: %s", init_param$unique_id)
    init_param <- unlist(init_param[, !"unique_id"])
  }, error = function(e) {
    log_error("SEED IMPORT: %s", e$message)
    stop(e$call)
  }
)


# Target summary statistic
tss             <- c(1,0)

# Global counter
gc              <- 1

# Try counter
tc              <- 1

# Initialse matrix to store hits
hitarray        <- list()

# Initialise matrix to store **breaker** parameter
brk_param_array <- list()

# Accepted hitarray
ahitarray       <- list()

# Initialise max hits counter
max_hits        <- 0

# Call ABC_mcmccore function
# Metadata writer
run_meta_file <- file.path(meta_folder, sprintf("%s_run_metadata.json", shell_dt))
run_meta <- list()
run_meta$prop_frac <- prop_frac
run_meta$target <- target
run_meta$seed_file <- seed_file
run_meta$seed_search <- seed_search
run_meta$timestamp <- shell_dt
run_meta$notes <-
  'Genesis seed search. 3 independent seeds. Target 17.
  PF - 1/50; PFS = 1/100
'

if (!file.exists(run_meta_file)) {
  write_json(x = run_meta, pretty = T, path = run_meta_file)
}

# Call ABC_mcmccore function
source("3.3_ABC_mcmc_core.R")

# Define restart buffer
buffer = function() {
  log_info("RESTART")
  withRestarts(
    expr = tryCatch(
      expr = {
        stopifnot(exists("buffer_param"))
        abc_mcmc_result <- invisible(ABC_mcmc(method = "Marjoram_original",
                                              model = ABC_mcmccore,
                                              prior = prior_ranges,
                                              summary_stat_target = tss,
                                              proposal_range = prop_ranges,
                                              n_rec = abc_n,
                                              verbose = TRUE,
                                              n_between_sampling = 1,
                                              acceptance = TRUE,
                                              rejection = TRUE,
                                              init_param = buffer_param))
      },
      error = function(e) {
        log_error("FAILED-RESTART: %s", e)
        invokeRestart("buffer")
      }
    ),
    buffer = buffer)
}

# ABC mcmc function call
try(
  {
    withRestarts(
      expr = tryCatch(
        expr = {
          log_info("MAIN-START")
          abc_mcmc_result <- invisible(ABC_mcmc(
            method = "Marjoram_original",
            model = ABC_mcmccore,
            prior = prior_ranges,
            summary_stat_target = tss,
            proposal_range = prop_ranges,
            n_rec = abc_n,
            verbose = TRUE,
            n_between_sampling = 1,
            acceptance = TRUE,
            rejection = TRUE,
            init_param = init_param))
        },
        error = function(e) {
          log_error("FAILED-MAIN: %s", e)
          invokeRestart("buffer")
        }
      ),
      buffer = buffer
    )
    log_info("FINISH ABC_MCMC FUNCTION")
  }
)

log_info("END OF SCRIPT >> GC: %s >> TC: %s >> G/T: %3.2f >> MAX HITS: %s", gc, tc, gc/tc, max_hits)

# Save results
# write.csv(ret_param_array, paste0("../output/", grid_task_int, "_run_ret_parameter_table.csv"), col.names = FALSE)
# write.csv(hitarray, paste0("../output/", grid_task_int, "_run_hit_table.csv"), col.names = FALSE)
# write.csv(brk_param_array, paste0("../output/", grid_task_int, "_run_brk_parameter_table.csv"), col.names = FALSE)
# write.csv(aret_param_array, paste0("../output/", grid_task_int, "_run_a_ret_parameter_table.csv"), col.names = FALSE)
# write.csv(ahitarray, paste0("../output/", grid_task_int, "_run_a_hit_table.csv"), col.names = FALSE)


# # Core function rebuilder
# library(compiler)
# source("3.3_ABC_mcmc_core.R")
# cmp_ABC_mcmccore <- cmpfun(ABC_mcmccore)
# saveRDS(cmp_ABC_mcmccore, "3.3_ABC_mcmc_core.R.cmp")
