library(data.table)
library(here)
library(digest)
library(jsonlite)
setwd(here())

session <- list()
session$raw_file <- "data/fitted_parameters/april_2020_cals/2020-04-17-1840-results.csv"
session$seed <- 123
session$outfile <- "data/fitted_parameters/april_2020_cals/april_2020_1K_china.csv"
session$md5_raw_file <- digest(file = T, object = session$raw_file, algo = "md5")

raw <- fread(session$raw_file, key = "unique_id")
raw <- unique(raw, by = "unique_id")

set.seed(session$seed)

subsampled <- raw[sample(x = nrow(raw), size = 1000, replace = F), ]

fwrite(x = subsampled, file = session$outfile)

session$run_time <- format(Sys.time(), "%Y-%m-%d-%H%M")
session$md5_outfile <- digest(file = T, object = session$outfile, algo = "md5")

write_json(x = session, path = "data/fitted_parameters/april_2020_cals/april_2020_1K_metadata.json", pretty = T)
