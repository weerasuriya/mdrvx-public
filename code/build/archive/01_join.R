mainloop <- function(data, mode = 0, para_static, para_variable, startcon, vaccine, transmission = 1, mdrac = 1, rx = 1) {
  time.start <- Sys.time()
  require(here)
  setwd(here("code"))

  list2env(data, environment(), hash = TRUE)
  list2env(para_static, environment(), hash = TRUE)
  list2env(para_variable, environment(), hash = TRUE)
  list2env(startcon, environment(), hash = TRUE)

  # Set up age classes - children, adults, younger adults, elderly
  chiyrs    <- 15
  aduyrs    <- 50
  yaduyrs   <- 40
  eldyrs    <- (max_age - chiyrs - aduyrs)

  # Global MDR switch. If mdrac==0, then disable MDR acquisition in this run.
  if (mdrac == 0) {
    if (exists("xi") == TRUE) {
      rm(xi)
    }
    assign("NT_xi", 0, envir = environment())
    assign("PT_xi", 0, envir = environment())
    mdrac_yr <- year1
  } else (
    mdrac_yr <- mdrac
  )

