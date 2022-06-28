# Loading libraries and functions
library(estmeansd)
library(doParallel)
library(foreach)
library(metafor)
library(fdrtool)

load('../data/simulation_settings.RData')

run_sims_malevel(batches = all_batches_ma, 
                 batchsize = batchsize_ma, 
                 seeds = seeds_ma, 
                 nboot = 1000)
