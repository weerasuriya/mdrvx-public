# Controller file for MDR-TB/Vaccine model.
# This file clears the environment, imports and executes the relevant other files.

# Clear working environment
rm(list=ls())
options(max.print=999999)
# Set user: Chathika=0, Rebecca=2, Finn=3 [needs to add in directory entry]

#   ____________________________________________________________________________
#   Set Path                                                                ####

C = 0

if (C == 0) {
  home <- "/home/chathika/Documents/PhD/MDR-Aeras/mdrvx"
}
if (C == 1) {
  home <- "/home/lsh355020/India_Gates/"
}

if (C == 2) {
  home <- '/Users/lsh355020/R/Chathika/MDR_China'
}

if (C==3) {
  home <- "C:/Users/cfmcquaid/Documents/R/AerasMDR/mdrvx"
}

setwd(home)

#   __________________ #< 8e1c7738c8da40bbb8f30c8bf11a888b ># __________________
#   Starting Conditions                                                     ####
startcon <- list(
  'dt'      = 0.25,
  'year1'   = 1900,
  'yearend' = 2099,
  "prop_col" = 1,
  'bgu' = 1,
  'fert' =1)

#   __________________ #< 40092d9ad0de565b0959b620f1d714cc ># __________________
#   Import External Data                                                    ####

source("1_import.R")

setwd(home)

#   __________________ #< 0924481b301e70d9ac6dfb624d3bc490 ># __________________
#   Initialise Mainloop                                                     ####

source("2_main.R")

#   __________________ #< 61d35efabacf18196a403f463e41908d ># __________________
#   Execute                                                                 ####
results <- mainloop(para_static, para_variable, startcon, vaccine=0,transmission=1,mdrac=0)
