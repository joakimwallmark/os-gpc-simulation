options(scipen = 999)
source("functions/get-abs-distance.R")
source("functions/tg-density.R")
library(mirt)
library(sn)
library(foreach)
library(doParallel)
library(patchwork)
library(ggplot2)
library(TestGardener)
source("functions/make-dataList.R") # this needs to be last, after library(TestGardener) for added functionality

# Load models -------------------------------------------------------------
data_string <- "natmat18"
data_string <- "swesat14"
load(paste("simulation-data/sim-", data_string, "-mirt-mod.RData", sep = ""))
load(paste("simulation-data/sim-", data_string, "-tg-mod.RData", sep = ""))
load(paste("simulation-data/sim-", data_string, "-item-data.RData", sep = ""))
load(paste("simulation-data/sim-", data_string, "-scenarios.RData", sep = ""))

os_cycles <- 10
max_surp <- 20
iterations <- 1000

# For testing
# iterations <- 3
# scenarios <- scenarios[c(15,13)]
# scenarios <- scenarios[c(9,11)]
# scenarios <- scenarios[1:2]
