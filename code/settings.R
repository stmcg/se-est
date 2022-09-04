rm(list = ls())
source('helper_functions_simulation.R')
source('mln.R')

################################################################################
## Study-level simulations
################################################################################

all_dist_studylevel <- c('normal', 'log-normal', 'log-normal', 'log-normal', 'halfnormal')
all_par1_studylevel <- c(5, 5, 5, 5, 0.1)
all_par2_studylevel <- c(1, 0.25, 0.5, 1, NA)
all_n_studylevel <- c(50, 250, 1000)
all_scenarios_studylevel <- c('S1', 'S2', 'S3')
sim_scenarios_studylevel <- as.data.frame(matrix(NA, ncol = 5))
colnames(sim_scenarios_studylevel) <- c('scenario', 'dist', 'n', 'par1', 'par2')
i <- 1
for (scenario in all_scenarios_studylevel){
  for (j in 1:length(all_dist_studylevel)){
    for (n in all_n_studylevel){
      sim_scenarios_studylevel[i, 'scenario'] <- scenario 
      sim_scenarios_studylevel[i, 'dist'] <- all_dist_studylevel[j] 
      sim_scenarios_studylevel[i, 'n'] <- n
      sim_scenarios_studylevel[i, 'par1'] <- all_par1_studylevel[j]
      sim_scenarios_studylevel[i, 'par2'] <- all_par2_studylevel[j]
      i <- i + 1
    }
  }
}

all_batches_studylevel1 <- 1:20
batchsize_studylevel1 <- 50

all_batches_studylevel2 <- 1:20
batchsize_studylevel2 <- 5000


################################################################################
## Meta-analytic level simulations
################################################################################

all_dist_ma <- c('log-normal')
all_par1_ma <- 5
all_par2_ma <- 0.25
all_n_studies_ma <- c(9, 30)
all_scenarios_ma <- c('S1', 'S2', 'S3')
all_prop_medians <- c(0, 1/3, 2/3, 1)
all_tausq_ma <- c(6, 0)
all_batches_ma <- 1:40
batchsize_ma <- 25


sim_scenarios_ma <- as.data.frame(matrix(NA, ncol = 7))
colnames(sim_scenarios_ma) <- c('dist', 'par1', 'par2', 'n_studies', 'scenario', 'tausq', 'prop_medians')
i <- 1
for (tausq in all_tausq_ma){
  for (scenario in all_scenarios_ma){
    for (prop_medians in all_prop_medians){
      for (n_studies in all_n_studies_ma){
        if (!((prop_medians == 0) & (scenario %in% c('S2', 'S3')))){
          sim_scenarios_ma[i, 'dist'] <- all_dist_ma
          sim_scenarios_ma[i, 'par1'] <- all_par1_ma
          sim_scenarios_ma[i, 'par2'] <- all_par2_ma
          sim_scenarios_ma[i, 'tausq'] <- tausq
          
          sim_scenarios_ma[i, 'n_studies'] <- n_studies
          sim_scenarios_ma[i, 'scenario'] <- scenario 
          sim_scenarios_ma[i, 'prop_medians'] <- prop_medians 
        }
        i <- i + 1
      }
    }
  }
}
sim_scenarios_ma <- sim_scenarios_ma[complete.cases(sim_scenarios_ma),]


################################################################################
## Getting seeds
################################################################################
set.seed(1234)
seed_init <- sample.int(2^30, size = 3)
seeds_studylevel1 <- get_ma_seeds(nscenarios = nrow(sim_scenarios_studylevel), n_all_batches = length(all_batches_studylevel1), seed_init = seed_init[1])
seeds_studylevel2 <- get_ma_seeds(nscenarios = nrow(sim_scenarios_studylevel), n_all_batches = length(all_batches_studylevel2), seed_init = seed_init[2])
seeds_ma <- get_ma_seeds(nscenarios = nrow(sim_scenarios_ma), n_all_batches = length(all_batches_ma), seed_init = seed_init[3])

save.image(file = '../data/simulation_settings.RData')
