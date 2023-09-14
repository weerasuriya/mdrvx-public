#!/usr/bin/Rscript --vanilla
# Second generation trajectory aggregator

suppressPackageStartupMessages(library(here))
library(digest)
library(data.table)
library(logger)
library(argparse, quietly = TRUE)
suppressPackageStartupMessages(library(zip))

if (!system("uname -n", intern = T) == "archer") {
  setDTthreads(1)
}

log_formatter(formatter_sprintf)
logger <- layout_glue_generator(format = '[{time}]\t{level}\t{msg}')
log_layout(logger)
log_appender(appender_tee(file = here("output", "logs", "traj_aggrs.log")))
log_info("START TRAJ_AGGREGATOR SCRIPT")

tryCatch(
  expr = {
    parser <- ArgumentParser()
    parser$add_argument("-f", type = "character", default = "singletons")
    parser$add_argument("-d", type = "character", default="INFO")
    parser$add_argument("-s", type = "character", default=TRUE)
    parser$add_argument("-m", type = "character", default=TRUE)
    args <- parser$parse_args()
    log_threshold(as.symbol(args$d))
    log_info('COMMAND LINE OPTS PARSED')
  }, error = function(e) {
    log_error("COMMAND LINE OPTS ERROR: %s", e$message)
    stop(e$call)
  }
)

tryCatch(
  expr = {
    # Read files and generate stems
    setwd(here("output", "trajectories", "CN", args[["f"]]))
    log_info("WD: %s", getwd())
    log_info("READ TRAJECTORY FILENAMES - START")
    files <- list.files(pattern = "^.*\\.csv$")
    stems <- unique(gsub(pattern = "(vx.__.*)__\\d_\\w{32}.csv", x = files, replacement = "\\1"))
    #stems_r <- unique(gsub(pattern = "(vx.)-(.*)-\\d_\\w{32}.csv", x = files, replacement = "\\1_\\2"))
    log_info("READ TRAJECTORY FILENAMES - OK")
  }, error = function(e) {
    log_error("ERROR STEMS: %s", e$message)
    stop(e$call)
  }
)

tryCatch(
  expr = {
    log_info("READ-RBIND - START")
    rbound <- lapply(stems, function(x) {
      rdt <- rbindlist(
        lapply(
          dir(pattern = x, include.dirs = F), fread, data.table = TRUE
        )
      )
      setkeyv(rdt, c("Year", "para_hash", "scen", "MGrp", "vxtype"))
      #gc(full = T)
    })
    rbound <- setNames(rbound, stems)
    log_info("READ-RBIND - OK")
  }, error = function(e) {
    log_error("ERROR READ-RBIND")
  }
)

if (args[['s']]) {
  tryCatch(
    expr = {
      log_info("SUBTRAJ - START")
      if (!dir.exists("subtraj")) dir.create("subtraj", recursive = T)
      subtraj <- list()
      for (nm in stems) {
        log_debug("START SUBTRAJ CALC: %s", nm)
        colnms <- colnames(rbound[[nm]])
        colnms <- colnms[colnms != "para_hash"]
        colnms <- colnms[colnms != "MVal"]
        subtraj[[nm]] <- rbound[[nm]][, .(min = min(MVal), med = median(MVal), max = max(MVal)), by = colnms]
        fwrite(x = subtraj[[nm]], file = sprintf("subtraj/%s.csv", nm))
        log_debug("WRITE OUT SUBTRAJ: %s", nm)
      }
      log_info("SUBTRAJ WRITE OUT - OK")
    }, error = function(e) {
      log_error("SUBTRAJ ERROR: %s", e$message)
    }
  )
}

tryCatch(
  expr = {
    log_info("STACKED OUTPUT - START")
    bl_stems <- grep(pattern = "^vx0.*", x = stems, value = T)
    uniq_stems <- gsub(x = bl_stems, pattern = "^vx.__(.*)$", replacement = "\\1")
    
    stacked_output <- lapply(uniq_stems, function(x) {
      # Construct dt names
      vset <- as.list(paste0("vx", c(0, 1,2,3), "__", x))
      vset_bound <- rbindlist(lapply(vset, function(x) rbound[[x]]))
      vset_bound
    })
    log_debug("STACKED OUTPUT - RBIND OK")
    stacked_output <- setNames(stacked_output, uniq_stems)
    
    stacked_output <- lapply(names(stacked_output), function(x) {
      colnms <- colnames(stacked_output[[x]])
      colnms <- colnms[colnms != "para_hash"]
      colnms <- colnms[colnms != "MVal"]
      stacked_output[[x]][, .(min = min(MVal), med = median(MVal), max = max(MVal)), by = colnms]
    })
    
    log_debug("STACKED OUTPUT - MIN/MED/MAX OK")
    stacked_output <- setNames(stacked_output, uniq_stems)
    
    if (!dir.exists("stacked")) dir.create("stacked", recursive = T)
    for (nm in names(stacked_output)) {
      fwrite(x = stacked_output[[nm]], file = sprintf("stacked/%s.csv", nm))
      log_debug("WRITE OUT STACKED: %s", nm)
    }
    log_info("STACKED OUTPUT - WRITE OUT OK")
  }, error = function(e) {
    log_error("STACKED OUTPUT: %s", e$message)
  }
)

tryCatch(
  expr = {
    log_info("MERGE SEQUENCE - START")
    bl_stems <- grep(pattern = "^vx0.*", x = stems, value = T)
    vx_stems <- sapply(c(1,2,3), function(x) gsub(pattern = "(vx)\\d(.*)", x = bl_stems, replacement = sprintf("\\1%s\\2", x)), simplify = T)
    log_debug("GENERATE MERGE STEMS - OK")
    
    merged_trajectories <- lapply(vx_stems, function(vxn) {
      blname <- gsub(pattern = "(vx)\\d(.*)", replacement = "\\10\\2", x = vxn)
      log_debug("MERGE %s- START", blname)
      op <- merge(x = rbound[[blname]][Year>=2027], y = rbound[[vxn]][Year>=2027], allow.cartesian = T, by = c("Year", "para_hash", "scen", "MGrp"), suffixes = c(".bl", ".vx"))
      op[, MVal.diff := MVal.vx - MVal.bl]
      log_debug("MERGE %s - OK", blname)
      op
      
    })
    merged_trajectories <- setNames(merged_trajectories, vx_stems)
    
    log_info("BL-VX TRAJECTORIES MERGE MMM - START")
    
    merged_trajectories <- lapply(names(merged_trajectories), function(x) {
      merged_trajectories[[x]][, .(min = min(MVal.diff), med = median(MVal.diff), max = max(MVal.diff)), by = c("Year", "scen", "MGrp", "vxtype.vx")]
    })
    
    merged_trajectories <- setNames(merged_trajectories, vx_stems)
    log_info("BL-VX TRAJECTORIES MERGE MMM - OK")
    
    log_info('BL-VX TRAJECTORIES MERGE - OK')
  }, error = function(e) {
    log_error("BL-VX TRAJECTORY MERGE: %s", e$message)
  }
)

if (args[['m']]) {
  tryCatch(
    expr = {
      log_info("MERGE WRITE - START")
      if (!dir.exists("merged")) dir.create("merged", recursive = T)
      for (nm in names(merged_trajectories)) {
        fwrite(x = merged_trajectories[[nm]], file = sprintf("merged/%s.csv", nm))
        log_debug("WRITE OUT MERGED: %s", nm)
      }
      log_info("MERGED WRITE - OK")
    }, error = function(e) {
      log_error("MERGED WRITE ERROR: %s", e$message)
    }
  )
}

tryCatch(
  expr = {
    try({
      log_info("ZIP MERGED - START")
      zipr(zipfile = sprintf("%s_merged.zip", basename(args[["f"]])), files = dir("merged", full.names = T))
    })
    try({
      log_info("ZIP STACKED - START")
      zipr(zipfile = sprintf("%s_stacked.zip", basename(args[["f"]])), files = dir("stacked", full.names = T))
    })
    try({
      log_info("ZIP SUBTRAJ - START")
      zipr(zipfile = sprintf("%s_subtraj.zip", basename(args[["f"]])), files = dir("subtraj", full.names = T))
    })
  }, warning = function(w) {
    log_error("ZIP: %s", w$message)
  }
)

log_info("END OF SCRIPT")
