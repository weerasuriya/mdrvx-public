# **Manual** Rejection-sampler Controller script for cluster
# Each run of this script is an SGE task

# Clear all
rm(list = ls())

# Libraries
library(here)
library(digest)
library(data.table)
library(jsonlite)

# Set working directory to code directory
setwd(here("code"))

# Read in compiled mainloop function, assign to .GlobalEnv
cmp_mainloop <- readRDS("2.1_main_rebuild.R.cmp")

# Starting conditions
startcon <- read_json("../data/startcon.json")

# Import data
source ("1_import.R", local = TRUE)

# Other switches
mode = 1
vaccine = 0
transmission = 1
rx = 1
mdrac = 1970

# Import fit targets
target_ranges <- fread("../data/fit_targets/fit_targets.csv", header = TRUE)

# Target range pairs
range_pairs   <- split(target_ranges,by = "Parameter")

# Target names
target_names  <- target_ranges$Parameter

# Range check
range_check   <- function(x, val) {
   return({val >= range_pairs[[x]]$LL & val <= range_pairs[[x]]$UL})
}

# Remove para_variable
#rm(para_variable)

# Prior range import
para_ranges            <- fread("../data/para_ranges.csv", header = TRUE)

# number of ABC runs
abc_n                  <- 1e+6
min_accept             <- 10000
accepted_param_counter <- 0

# Initialse matrix to store hits
hitarray               <- list()

# Initialise matrix to store **breaker** parameter
brk_param_array        <- list()

# Store current SGE_TASK_ID
# Check current host and store SGE_TASK_ID if appropriate
if (system("uname -n", intern = T) == "archer") {
  op <- file.path("..", "output", "random_sampler")
  grid_task_int <- 1
} else {
  grid_task_int <- as.numeric(Sys.getenv("SGE_TASK_ID"))
  folder_num <- sprintf(fmt = "%05d", grid_task_int %/% 1000)
  op <- file.path(Sys.getenv("OPPATH"), folder_num)
  try(if (!dir.exists(op)) {
    dir.create(path = op, recursive = TRUE)
  })
  # Set DT threads for cluster - disable internal parallelism in data.tables to
  # prevent cluster multicore use
  setDTthreads(1)
}

# Generate list of parameter names and ranges
prl <- split(para_ranges, f = para_ranges$Name)

# Uniform prior generator with constraints
unif_prior_gen <- function(prl) {
   pv               <- lapply(prl, function(x) runif(1, min = x$low, max = x$high))
   # Constraints
   pv$nelderly      <- runif(1, min = prl$nelderly$low, max = pv$n)
   pv$pelderly      <- runif(1, min = pv$padult, max = prl$pelderly$high)
   pv$velderly      <- runif(1, min = pv$vadult, max = prl$pelderly$high)
   pv$felderly      <- runif(1, min = prl$felderly$low, max = pv$fadult)
   pv$relderly      <- runif(1, min = pv$radult, max = prl$relderly$high)
   pv$uiscaleC      <- runif(1, min = pv$uiscaleA, max = prl$uiscaleC$high)
   pv$uiscaleE      <- runif(1, min = pv$uiscaleA, max = prl$uiscaleE$high)
   pv$CDRscaleE     <- runif(1, min = prl$CDRscaleE$low, max = pv$CDRscale)
   return(pv)
}

# Ongoing maximum hit counter
max_hit <- 0

# Set up loop of random sampled runs
for (i in 1:abc_n) {

   # Generate random parameter set
   para_variable         <- unif_prior_gen(prl)

   # Unique ID
   unique_id <- digest(para_variable)

   tryCatch(
      {
         results <<- cmp_mainloop(data,
            mode = 1,
            para_static = para_static,
            para_variable = para_variable,
            para_cost = para_cost,
            startcon = startcon,
            vaccine = 0,
            transmission = 1,
            mdrac = 1970,
            rx = 1)

            hits <- sapply(target_names, function(x) range_check(x, results[[x]]), USE.NAMES = T)
            print(paste("Hits: ", sum(hits)))
            browser
            hitarray[[i]] <- c(c(total = sum(hits), task = grid_task_int, unique_id = unique_id), as.list(hits), para_variable)
            if (sum(hits) > max_hit) {max_hit <<- sum(unlist(hits))}
            accepted_param_counter <- accepted_param_counter + 1

         },
         error = function(e) {
            brk_param_array[[i]] <<- c(list(e$call$rtype, task = grid_task_int, unique_id = unique_id), para_variable)
         }
      )

      if (exists("para_variable")) rm(para_variable)
      if (exists("results")) rm(results)

      if (i%%1000 == 0 | i == abc_n | accepted_param_counter == min_accept + 1) {
         hitarray        <- rbindlist(hitarray)
         brk_param_array <- rbindlist(brk_param_array)
         fwrite(hitarray, file.path(op, paste0(grid_task_int, "_run_hit_table.csv")), col.names = FALSE, sep = ",", append = TRUE, row.names = FALSE, logical01 = T)
         fwrite(brk_param_array, file.path(op, paste0(grid_task_int, "_run_brk_parameter_table.csv")), col.names = FALSE, sep = ",", append = TRUE, row.names = FALSE)
         # Reset matrices
         hitarray        <- list()
         brk_param_array <- list()
      }

      if (accepted_param_counter == (min_accept + 1)) break

   }

   # if (
   #    any(is.na(unlist(results))) |
   #    any(is.nan(unlist(results))) |
   #    results == "Reject-ntu" |
   #    results == "Reject-gen" |
   #    results == "Reject"
   # ) {
   #    # If core model rejects or has absurd results, then reject ABC run as well.
   #    if (any(is.na(results)) | any(is.nan(results))) {
   #       rtype = 3
   #    } else if (results == "Reject-gen") {
   #       rtype = 2
   #    } else if (results == "Reject-ntu") {
   #       rtype = 1
   #    } else {
   #       rtype = 4
   #    }
   #
   # } else {
   #    hits <- mapply(function(res, min, max) {return(res <= max && res >= min)}, res = results, min = target_ranges$LL, max = target_ranges$UL, SIMPLIFY = F)
   #    hitarray[[i]] <- c(list(total = sum(unlist(hits)), task = grid_task_int, unique_id = unique_id), hits, para_variable)
   # }
   #
   # rm(para_variable, results)
   # if (exists("rtype")) {
   #    rm(rtype)
   # }
   #
   # if (sum(unlist(hits)) > max_hit) {max_hit <<- sum(unlist(hits))}
   #
   # print(paste("Max Hits:", max_hit))
   #
   # print(paste("i = ", i))

hitarray
