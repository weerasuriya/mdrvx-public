if (vaccine == 1) {
  vx_startstep                                 <- ((vx_start_year - year1) * (1/dt)) + 1

  if (ageM == -1) widthM <- -1
  if (ageR != -1 & ageM != -1) {ageGap <- ageM - ageR}
  if (widthM == -1) {
    widthM <- (age_cls - 1) - ageM
  }

  # Set up vaccination and waning matrices vaccinated
  vaccinated <- vaccinatedR <- vaccinatedM <- notvaccinated <- wanes <- matrix(0, steps, age_cls)

  # Set up timing If the duration of protection is greater than mass campaign
  # interval, then set duration = mass campaign interval
  if (duration > mass_interval) {
    mass_interval                              <- duration
  }
  # Vector of **step numbers** where routine vaccination occurs
  # Vector of **step numbers** where mass vaccination occurs
  if (ageR != -1) {

    rtvxsteps                                    <- seq(vx_startstep, steps, (1/dt))
  ## Vaccine matrices - insert ## Routine vaccination
  vaccinated[rtvxsteps, (ageR + 1)]            <- coverageR
  vaccinatedR[rtvxsteps, (ageR + 1)]           <- coverageR
  }
  # Mass vaccination campaigns
  if (ageM != -1) {
    # Vector of **step numbers** where mass vaccination occurs
    massvxsteps                                  <- seq(vx_startstep, steps, (mass_interval/dt))

    # Mass vaccination campaigns
    vaccinated[massvxsteps, (ageM + 1):(ageM + 1 + widthM)]  <- coverageM
    vaccinatedM[massvxsteps, (ageM + 1):(ageM + 1 + widthM)] <- coverageM
  }


  if (duration < (yearend - vx_start_year)) {

    if (ageR != -1) {
    ## Waning matrices - insert ## Waning of routinly vaccinated
      wanes[seq((vx_startstep + (duration/dt)), steps, 1/dt), ageR + 1 + duration] <- 1
    }

    if (ageM != -1) {
    # Waning of mass vaccinated
    massexitsteps <- massvxsteps + (duration/dt)
      wanes[massexitsteps[massexitsteps <= steps], (ageM + 1 + duration):min((ageM + 1 + widthM + duration), age_cls)] <- 1
    # Drop routine vaccinees who no longer would wane out as they are revaccinated in
      # If both routine and mass campaigns activated
      if (ageR != -1) {
    # mass campaigns
    dropsteps <- c()

    for (st in 1:length(massvxsteps)) {
          # dropsteps <- c(dropsteps, (massvxsteps[st] + 1):min(steps, (massvxsteps[st]  - (ageM-ageR)/dt + ((duration)/dt))))
          dropsteps <- c(dropsteps, (massvxsteps[st] + 1):min(steps, (massvxsteps[st]  + ((duration - ageGap)/dt))))
    }
        if (ageR + duration >= ageM) {
    wanes[dropsteps, ageR + 1 + duration] <- 0
        }
      }

    }
  } else {
    # Lifelong vaccination - no waning
    wanes                                        <- matrix(0, steps, age_cls)
  }

  # Convert vaccine matrices to compartment specific matrices depending on vaccine
  # type
  if (vxtype == 1) {
    vaxmatS <- vaxmatL <- vaxmatR <- vaccinated

  } else if (vxtype == 2) {
    vaxmatS <- vaccinated
    vaxmatR <- vaxmatL <- notvaccinated

  } else if (vxtype == 3) {
    vaxmatS <- notvaccinated
    vaxmatR <- vaxmatL <- vaccinated

  }

  # Extra matrix to represent _vaccinations_ rather than transitions
  vaxmatI <- vaccinated
}
