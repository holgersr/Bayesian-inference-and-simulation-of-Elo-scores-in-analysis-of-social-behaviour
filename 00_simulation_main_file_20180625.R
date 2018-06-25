## #######################
## Simulation main file ##
## #######################
library("rstan")
library("viridis")
library("loo")

source("00b_functions_ml_approach.R")
source("00b_simulation_functions.R")

today <- format(Sys.Date(), "%Y%m%d")

## choose (umcomment) one of the following:
#source("01_simulation_1.R")
#source("02_simulation_2.R")
#source("03_simulation_3.R")
#source("04_simulation_4.R")
#source("05_simulation_5.R")