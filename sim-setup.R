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
load(paste("simulation-data/sim-swesat14-mirt-mod.RData", sep = ""))
load(paste("simulation-data/sim-swesat14-tg-mod.RData", sep = ""))
load(paste("simulation-data/sim-swesat14-item-data.RData", sep = ""))
load(paste("simulation-data/sim-swesat14-scenarios.RData", sep = ""))

os_cycles <- 10
max_surp <- 20
iterations <- 1000
