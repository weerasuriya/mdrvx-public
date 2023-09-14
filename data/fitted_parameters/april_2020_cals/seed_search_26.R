library(data.table)
library(here)
library(digest)
library(jsonlite)
setwd(here())



session <- list()
session$raw_master <- "data/fitted_parameters/april_2020_cals/2020-03-10-2258-results.csv"
session$prng_seed <- 123
session$md5_raw_master <- digest(object = raw_master, algo = 'md5', file = T)
raw_26 <- fread("data/fitted_parameters/april_2020_cals/2020-03-10-2258-results.csv")[hit_sum==26]
set.seed(session$prng_seed)
seeds <- raw_26[sample(x = 1:nrow(raw_26), size = 1000, replace = F)]
session$outfile <- "data/fitted_parameters/april_2020_cals/seed_search_26.csv"
fwrite(file = session$outfile, x = seeds)
session$md5_outfile <- digest(object = session$outfile, file = T, algo = "md5")
session$note <- "Using output of ascending MCMC, generate 1000 seeds for full parameter exploration in China"
write_json(x = session, path = "data/fitted_parameters/april_2020_cals/seed_search.json", pretty = T)
