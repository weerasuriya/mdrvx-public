## ******* READ IN DATA Read in data from csv files (Stored in Data folder)
library(data.table)
library(jsonlite)
# List containing data objects
data = list()
country_code = "CN"
## contact matrix - this is the data from Read 2014 (made reciprocal)
data$myneta           <- data.matrix(read.csv("../data/myneta2.csv", header = FALSE)[, 2:101])

## Births - births from 1900-2099, absolute numbers, per annum
data$births           <- data.matrix(fread("../data/demographics/demographics_births.csv", select = c("Year", country_code)))

## Deaths/mortality
#data$death_rate      <- data.matrix(read.csv("../data/modified_death_rate_WPP2019.csv", header = TRUE))
data$death_rate       <- data.matrix(fread("../data/demographics/demographics_mortality.csv")[country==country_code][, country := NULL])

## Population by age and year, from 1900-2099. Values for 1900-1949 are copied
## from 1950-1999.
#data$total_population <- data.matrix(read.csv("../data/data_total_population_1900-2099_WPP2019.csv", header = T)[, 2:101])
data$total_population <- data.matrix(fread("../data/demographics/demographics_total_population.csv")[country==country_code][, c("country", "Year") := NULL])

# Import population classes which needed to be summed and retained during model
# run.
data$popclass         <- read.csv("../data/calc_popclass.csv", header = TRUE)

# Initial proportions of population with various TB stages [i.e. to be
# distributed among the SIRetc epidemiologic compartments]. Assumed the same
# across all ages.
data$init_tb_prop     <- read.csv("../data/data_prop.csv", header = TRUE, row.names = 1)

# Import Parameters para_static - input variables which do not change / not sampled / pre-scaling
para_static           <- read_json("../data/para_static.json")

# para_variable - input variables which are sampled / scaling / fitting factors etc
df_para_variable      <- read.csv("../data/para_variable.csv", header = TRUE)
para_variable         <- setNames(as.list(df_para_variable[, 2]), df_para_variable[, 1])
rm(df_para_variable)

# para_cost - cost model values
para_cost             <- read_json("../data/para_cost.json")

# para_vax - vaccine related parameters
para_vax              <- read_json("../data/para_vax.json", simplifyVector = F)

# para_range - ranges for para_var members
para_ranges           <- data.frame(read.csv("../data/para_ranges.csv", header = TRUE))

## Treatment Success Rate (=> leads to psi terms)
data$psi_pp           <- fread("../data/data_psi_pp.csv", header = T)

# Weight for China MDR-TB Rx - R1/2
data$mdrtb_rx         <- fread("../data/data_mdrtb_rx.csv")

## AGES
data$widthage         <- 1 # Take individual age classes
data$max_age          <- 99 # Maximum age
data$age_cls          <- ceiling((data$max_age + 1)/data$widthage) # Number of age classes for matrix
age_nms               <- c(as.character(0:data$max_age)) # Names for age classes, assuming a 0-year old start

# DST probability scale up
dstp_scaleup	      <- fread("../data/data_dstp_scaleup.csv")
