library(estmeansd)
library(doParallel)
library(foreach)
library(metafor)
library(fdrtool)
library(data.table)

load('../data/simulation_settings.RData')

M <- nrow(sim_scenarios_ma)
nMonteCarlo <- 1e6
registerDoParallel(cores=10)

set.seed(1)
seeds <- sample.int(2^30, size = 10)
res <- foreach(m = 1:(M / 2)) %dopar% {
  set.seed(seeds[m])
  rowind <- m * 2 - 1
  dist <- sim_scenarios_ma[rowind, 'dist']
  par1 <- sim_scenarios_ma[rowind, 'par1']
  par2 <- sim_scenarios_ma[rowind, 'par2']
  scenario <- sim_scenarios_ma[rowind, 'scenario']
  tausq <- sim_scenarios_ma[rowind, 'tausq']
  prop_medians <- sim_scenarios_ma[rowind, 'prop_medians']
  
  sample_sizes <- round(runif(n = nMonteCarlo, min = 100, max = 500))
  medians <- rep(FALSE, times = nMonteCarlo)
  if (prop_medians > 0){
    medians[1:round(prop_medians * nMonteCarlo)] <- TRUE
  }
  re <- rnorm(nMonteCarlo, mean = 0, sd = sqrt(tausq))

  res_qe <- res_bc <- res_mln <- rep(NA, times = nMonteCarlo)
  for (i in 1:nMonteCarlo){
    temp <- sim_data(dist = dist, par1 = par1, par2 = par2, n = sample_sizes[i], medians = medians[i], re = re[i])
    if (medians[i]) {
      res_qe[i] <- apply_method(method = 'qe', quants = temp, scenario = scenario, n = sample_sizes[i], nboot = 0)$mean_est
      res_bc[i] <- apply_method(method = 'bc', quants = temp, scenario = scenario, n = sample_sizes[i], nboot = 0)$mean_est
      res_mln[i] <- apply_method(method = 'mln', quants = temp, scenario = scenario, n = sample_sizes[i], nboot = 0)$mean_est
    } else {
      res_qe[i] <- res_bc[i] <- res_mln[i] <- temp[1]
    }
  }
  return(list(qe = tausq / var(res_qe), bc = tausq / var(res_bc), mln = tausq / var(res_mln)))
}

res <- rbindlist(res)
res_table <- res[rep(1:10, each = 2),]

save(res_table, file = '../results/true_I2.RData')
