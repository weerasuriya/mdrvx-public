# Calculate cumulative discounted costs over a given time frame
disc_cumul_cost <- function(dtab, int_year = 2027, end_year = 2099, discount = 0.03, ccol = 2) {

  costVector <- dtab[Year>=int_year & Year<=end_year, 2]
  discounted_costs <- costVector
  for (i in seq_len(length(nrow(costVector)))) {
    discounted_costs[i, ] <- costVector[i, ]/((1 + discount)^(i - 1))
  }
  return(sum(discounted_costs))
}

# Years lost due to disability calculator [YLD]
yld_calculator <- function(table, daly_wt, startyr, endyr, discount = FALSE, discount_rate = 0.035) {

  # Convert to data.table
  require(data.table)
  table <- data.table(table)

  # Extract annual data
  yld_annual <- rowSums(table[, -1]) * daly_wt
  colnames(table)[1] <- "Year"

  # Isolate start and end indices
  yld_start_index <- table[Year == startyr, which = TRUE]
  yld_end_index <- table[Year == endyr, which = TRUE]

  # Isolate relevant index range
  yld_range <- c(yld_start_index:yld_end_index)

  # Subset annual table to relevant range
  yld_annual <- yld_annual[yld_range]

  # Check if discounting enabled
  if (discount == FALSE) {
    return(sum(yld_annual))
  } else if (discount == TRUE) {
    discounted_yll <- yld_annual
    for (i in seq_len(length(yld_annual))) {
      discounted_yll[i] <- yld_annual[i]/(1 + discount_rate)^(i - 1)
    }
    return(sum(discounted_yll))
  }
}

# Years of Life Lost due to Mortality [YLL] yld_calculator
yll_calculator <- function(mort, le, startyr, endyr, discount = FALSE, discount_rate = NULL) {

  # Convert to data.table
  require(data.table)
  mort <- data.table(mort)
  le <- data.table(le)

  # Extract relevant mortality years and ages
  mort_start_index <- mort[Year == startyr, which = TRUE]
  mort_end_index <- mort[Year == endyr, which = TRUE]
  mort_range <- c(mort_start_index:mort_end_index)
  le_start_index <- le[Year == startyr, which = TRUE]
  le_end_index <- le[Year == endyr, which = TRUE]
  le_range <- c(le_start_index:le_end_index)

  # calculate YLL
  yll <- mort[mort_range, -1] * le[le_range, -1] * 1000

  # Discount if appropriate and return
  if (discount == FALSE) {
    return(sum(yll))
  } else if (discount == TRUE) {
    discounted_yll <- yll
    for (i in seq_len(nrow(yll))) {
      discounted_yll[i, ] <- yll[i, ]/((1 + discount_rate)^(i - 1))
    }
    return(sum(discounted_yll))
  }
}

# Define bodger function for dealing with multiple model trend data
bodger <- function(table) {
  require(tidyverse)
  require(magrittr)
  output_table <- NULL
  desc <- c("median", "min", "max")
  agegroups <- colnames(table)[-c(1, 2)]
  grp <- 1
  for (grp in 1:length(agegroups)) {
    med <- table[, c(1:2, grp + 2)] %>% spread(key = Year, value = agegroups[grp]) %>%
                                select(-(Hash)) %>%
                                summarise_all(median) %>%
                                add_column(.before = 1, agegroups[grp])

    minima <- table[, c(1:2, grp + 2)] %>% spread(key = Year, value = agegroups[grp]) %>%
                                select(-(Hash)) %>%
                                summarise_all(min) %>%
                                add_column(.before = 1, agegroups[grp])

    maxima <- table[, c(1:2, grp + 2)] %>% spread(key = Year, value = agegroups[grp]) %>%
                                select(-(Hash)) %>%
                                summarise_all(max) %>%
                                add_column(.before = 1, agegroups[grp])
    output_table <- bind_rows(output_table,med,minima, maxima)
  }
  output_table <- add_column(.data = output_table, .after = 1, rep(desc, length(agegroups)))
  colnames(output_table)[1] <- "AgeGrp"
  colnames(output_table)[2] <- "Stat"
  return(output_table)
}
