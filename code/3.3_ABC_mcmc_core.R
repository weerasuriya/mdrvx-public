# Core ABC rejection sampler

# Input is a **vector** containing the para_variable values

ABC_mcmccore <- function(para_variable_vector) {

  # TRY COUNTER
  tc <<- tc + 1
  log_debug("TRY COUNTER: %s", tc)

   tryCatch(expr = {
      para_variable <- as.list(para_variable_vector)
      names(para_variable) <- pvnames
      unique_id <- digest(para_variable)
      log_trace("PARA_VARIABLE: %s", para_variable)
   }, error = function(e) {
      log_error("PARA VARIABLE: %s", e)
      stop()
   })


  # Constraint test
  if (para_variable$uiscaleC < para_variable$uiscaleA |
    para_variable$uiscaleE < para_variable$uiscaleA |
    para_variable$felderly > para_variable$fadult |
    para_variable$nelderly > para_variable$n |
    para_variable$pelderly < para_variable$padult |
    para_variable$velderly < para_variable$vadult |
    para_variable$relderly < para_variable$radult |
    para_variable$CDRscale < para_variable$CDRscaleE) {
    # Reject based on constraints
    log_debug("FAILED PARAMETER CONSTRAINTS")
    return(invisible(c(0, 0)))
  }

  # Run model
  results <- tryCatch({
    cmp_mainloop(
      data,
      mode = 1,
      para_static,
      para_variable,
      para_cost,
      startcon,
      vaccine = 0,
      para_vax,
      transmission = 1,
      mdrac = 1970,
      rx = 1,
      cost = 0,
      scen = 1)
  }, error = function(e) {
    return(invisible(e))
  })

  log_trace("RUN RESULT: %s", results)

  # Analyse the result - determine where to store and return value
  if (any(class(results)=="simpleError")) {
    brk_param_array[[gc]] <- c(list(rtype = results$call$type, task = grid_task_int, unique_id = unique_id), para_variable)
    log_debug("RUN FAILED: %s", results$message)
    return(invisible(c(0, 0)))
  } else {
    hits <- vapply(target_names, function(x) {return(target_ranges[[x]]$LL <= results[[x]] & results[[x]] <= target_ranges[[x]]$UL)}, FUN.VALUE = logical(length = 1))
    hitarray[[gc]] <<- c(list(total = sum(hits), task = grid_task_int, unique_id = unique_id), hits, para_variable)

    log_debug("HITSUM: %s", hitarray[[gc]]$total)

    if (hitarray[[gc]]$total >= target) {
      ahitarray[[gc]] <<- c(list(total = sum(hits), task = grid_task_int, unique_id = unique_id), hits, para_variable)
      ret_int <- 1
    } else {
      ret_int <- 0
    }
    log_debug("GC: %s", gc)

    if (sum(hits) > max_hits) {
      log_info("MAX HIT INCREASE: %s -> %s", max_hits, sum(hits))
      max_hits <<- sum(hits)
    } else if (sum(hits) == target & seed_search == TRUE) {
      log_info("SEED SEARCH: %s", sum(hits))
    }

    if (gc%%wo_interval == 0 | gc == abc_n) {

      hitarray <<- rbindlist(hitarray[!sapply(hitarray, is.null)])
      fwrite(hitarray, file.path(op, paste0(grid_task_int, "_run_hit_table.csv")), col.names = FALSE, sep = ",", append = TRUE, row.names = FALSE, logical01 = T)

      if (length(ahitarray) >= 1) {
        ahitarray <<- rbindlist(ahitarray[!sapply(ahitarray, is.null)])
        fwrite(ahitarray, file.path(op, paste0(grid_task_int, "_run_ahit_table.csv")), col.names = FALSE, sep = ",", append = TRUE, row.names = FALSE, logical01 = T)
      }
      if (length(brk_param_array) >= 1) {
        brk_param_array <<- rbindlist(brk_param_array[!sapply(brk_param_array, is.null)])
        fwrite(brk_param_array, file.path(op, paste0(grid_task_int, "_run_brk_parameter_table.csv")), col.names = FALSE, sep = ",", append = TRUE, row.names = FALSE)
      }

      # Reset matrices
      hitarray <<- list()
      brk_param_array <<- list()
      ahitarray <<- list()
      log_info("CHAIN LENGTH / WRITE OUT: %s >> TRY LENGTH: %s >> G/T: %3.2f", gc, tc, gc/tc)
    }

    gc <<- gc + 1
    #assign("buffer_param", para_variable, pos = .GlobalEnv)
    buffer_param <<- para_variable

    if (ret_int == 1) {
      rm(ret_int)
      return(invisible(c(1, 0)))
    } else if (ret_int == 0) {
      return(invisible(c(0, 0)))
    }

  }
}
