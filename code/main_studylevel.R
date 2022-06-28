# Loading libraries and functions
library(estmeansd)
library(doParallel)
library(foreach)
library(metafor)
library(fdrtool)
library(metaBLUE)

load('../data/simulation_settings.RData')

# QE and BC with and without bootstrap
res_est <- run_sims_studylevel(batches = all_batches_studylevel1,
                               batchsize = batchsize_studylevel1,
                               seeds = seeds_studylevel1,
                               methods_bootstrap = c('qe', 'bc', 'mln'),
                               nboot = 1e3)

# Getting true values by Monte Carlo method
res_truth <- run_sims_studylevel(batches = all_batches_studylevel2,
                                 batchsize = batchsize_studylevel2,
                                 seeds = seeds_studylevel2,
                                 methods_bootstrap = c(),
                                 nboot = 0)
