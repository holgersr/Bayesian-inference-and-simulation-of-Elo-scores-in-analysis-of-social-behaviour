## #######################
## Simulation main file ##
## #######################
library("rstan")
library("viridis")
library("loo")

source("00b_functions_ml_approach_20180625.R")
source("00b_simulation_functions_20180625.R")

today <- format(Sys.Date(), "%Y%m%d")

## choose (umcomment) one of the following:
#source("01_simulation_1_20180625.R")
#source("02_simulation_2_20180625.R")
#source("03_simulation_3_20180625.R")
#source("04_simulation_4_20180625.R")
#source("05_simulation_5_20180625.R")
