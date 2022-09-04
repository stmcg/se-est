run_sims_studylevel <- function(batches, batchsize, seeds, methods_bootstrap, nboot){
  ncores <- n_batches <- length(batches)
  if (ncores > 1){
    registerDoParallel(cores=ncores)
  }
  
  for (m in 1:nrow(sim_scenarios_studylevel)){
    seeds_batch <-seeds[[m]]
    dist <- sim_scenarios_studylevel[m, 'dist']
    par1 <- sim_scenarios_studylevel[m, 'par1']
    par2 <- sim_scenarios_studylevel[m, 'par2']
    n <- sim_scenarios_studylevel[m, 'n']
    scenario <- sim_scenarios_studylevel[m, 'scenario']
    
    temp <- foreach (t = 1:n_batches) %dopar%{
      set.seed(seeds_batch[batches[t]])
      res_qe <- res_bc <- res_mln <- res_hozo <- res_wan <- res_luo <- 
        res_shi <- res_yang <- 
        data.frame(selected_dist = rep(NA, batchsize),
                   mean_est = rep(NA, batchsize),
                   se_est_naive = rep(NA, batchsize), 
                   se_est_boot = rep(NA, batchsize))
      
      for (i in 1:batchsize){
        dat <- sim_data(dist = dist, par1 = par1, par2 = par2, n = n, medians = TRUE, re = 0)
        
        if ('qe' %in% methods_bootstrap){
          fit.qe <- apply_method(method = 'qe', quants = dat, scenario = scenario, n = n, nboot = nboot)
        } else {
          fit.qe <- apply_method(method = 'qe', quants = dat, scenario = scenario, n = n, nboot = 0)
        }
        
        if ('bc' %in% methods_bootstrap){
          fit.bc <- apply_method(method = 'bc', quants = dat, scenario = scenario, n = n, nboot = nboot)
        } else {
          fit.bc <- apply_method(method = 'bc', quants = dat, scenario = scenario, n = n, nboot = 0)
        }
        if ('mln' %in% methods_bootstrap){
          fit.mln <- apply_method(method = 'mln', quants = dat, scenario = scenario, n = n, nboot = nboot)
        } else {
          fit.mln <- apply_method(method = 'mln', quants = dat, scenario = scenario, n = n, nboot = 0)
        }
        
        if (scenario == 'S1' & dist == 'log-normal'){
          res_hozo[i, 'mean_est'] <- dat[3]
          res_hozo[i, 'se_est_naive'] <- sqrt(((dat[1] - 2 * dat[3] + dat[5])^2 / 4 + (dat[5] - dat[1])^2) / 12) / sqrt(n)
          
          res_wan[i, 'mean_est'] <- dat[3]
          res_wan[i, 'se_est_naive'] <- Wan.std(X = dat[c(1, 3, 5)], n = n, type = scenario)$sigmahat / sqrt(n)
          
          res_luo[i, 'mean_est'] <- Luo.mean(X = dat[c(1, 3, 5)], n = n, type = scenario)$muhat
          res_luo[i, 'se_est_naive'] <- Wan.std(X = dat[c(1, 3, 5)], n = n, type = scenario)$sigmahat / sqrt(n)
          
          C <- stats::qnorm((n - 0.375)/(n + 0.25))
          ls_mean <- Luo.mean(X = log(dat[c(1, 3, 5)]), n = n, type = scenario)$muhat
          ls_var <- (log(dat[5]) - log(dat[1]))^2 / (2 * C)^2 / (1.01 + 0.25 / log(n)^2)
          sigma4 <- ((log(dat[5]) - log(dat[1]))/ (2 * C))^4 / (1 + 2.23 / (log(n)^2))
          res_shi[i, 'mean_est'] <- exp(ls_mean + ls_var / 2) / (1 + 0.565 / n * ls_var + 0.37 / n * sigma4)
          res_shi[i, 'se_est_naive'] <- 
            sqrt(exp(2 * ls_mean + 2 * ls_var) / (1 + 2.26 / n * ls_var + 5.92 / n * sigma4) - 
                   exp(2 * ls_mean + ls_var) / (1 + 2.26 / n * ls_var + 1.48 / n * sigma4)) / sqrt(n)
          
          res_yang[i, 'mean_est'] <- BLUE_s(X = dat[c(1, 3, 5)], n = n, type = scenario)$muhat
          res_yang[i, 'se_est_naive'] <- BLUE_s(X = dat[c(1, 3, 5)], n = n, type = scenario)$sigmahat / sqrt(n)
        }
        
        
        res_qe[i, 'selected_dist'] <- fit.qe$selected_dist
        
        res_qe[i, 'mean_est'] <- fit.qe$mean_est
        res_bc[i, 'mean_est'] <- fit.bc$mean_est
        res_mln[i, 'mean_est'] <- fit.mln$mean_est
        
        res_qe[i, 'se_est_naive'] <- fit.qe$se_est_naive
        res_bc[i, 'se_est_naive'] <- fit.bc$se_est_naive
        res_mln[i, 'se_est_naive'] <- fit.mln$se_est_naive
        
        res_qe[i, 'se_est_boot'] <- fit.qe$se_est_boot
        res_bc[i, 'se_est_boot'] <- fit.bc$se_est_boot
        res_mln[i, 'se_est_boot'] <- fit.mln$se_est_boot
        
      }
      truth_mc <- nboot == 0
      if (truth_mc){
        # We seek to get true estimates by MC
        batchres <- list(res_qe = res_qe$mean_est, res_bc = res_bc$mean_est, 
                         res_mln = res_mln$mean_est, res_hozo = res_hozo$mean_est, 
                         res_wan = res_wan$mean_est, res_luo = res_luo$mean_est, 
                         res_shi = res_shi$mean_est, res_yang = res_yang$mean_est)
      } else {
        batchres <- list(res_qe = res_qe, res_bc = res_bc, res_mln = res_mln, 
                         res_hozo = res_hozo, res_wan = res_wan, res_luo = res_luo, 
                         res_shi = res_shi, res_yang = res_yang)
      }
      save(batchres, file = paste0('../results/StudyLevel_SimScenario=', m, 'Batch=', batches[t], 'Truth=', truth_mc, '.RData'))
    }
  }
}


run_sims_malevel <- function(batches, batchsize, seeds, nboot){
  ncores <- n_batches <- length(batches)
  if (ncores > 1){
    registerDoParallel(cores=ncores)
  }
  
  for (m in 1:nrow(sim_scenarios_ma)){
    seeds_batch <-seeds[[m]]
    dist <- sim_scenarios_ma[m, 'dist']
    par1 <- sim_scenarios_ma[m, 'par1']
    par2 <- sim_scenarios_ma[m, 'par2']
    n_studies <- sim_scenarios_ma[m, 'n_studies']
    scenario <- sim_scenarios_ma[m, 'scenario']
    prop_medians <- sim_scenarios_ma[m, 'prop_medians']
    tausq <- sim_scenarios_ma[m, 'tausq']
    random_effect <- tausq > 0
    
    temp <- foreach (t = 1:n_batches) %dopar%{
      set.seed(seeds_batch[batches[t]])
      res_qe_naive <- res_qe_boot <- 
        res_bc_naive <- res_bc_boot <- 
        res_mln_naive <- res_mln_boot <- 
        data.frame(pooled_mean_est = rep(NA, batchsize),
                   pooled_mean_ci_lb = rep(NA, batchsize),
                   pooled_mean_ci_ub = rep(NA, batchsize),
                   tausq_est = rep(NA, batchsize), 
                   tausq_ci_lb = rep(NA, batchsize), 
                   tausq_ci_ub = rep(NA, batchsize), 
                   I2 = rep(NA, batchsize))
      
      for (i in 1:batchsize){
        dat_ma <- sim_data_ma(dist = dist, par1 = par1, par2 = par2, n_studies = n_studies, 
                              prop_medians = prop_medians, tausq = tausq)
        yi_qe <- yi_bc <- yi_mln <- 
          vi_qe_naive <- vi_qe_boot <- 
          vi_bc_naive <- vi_bc_boot <- 
          vi_mln_naive <- vi_mln_boot <- rep(NA, n_studies)
        
        for (study in 1:n_studies){
          temp <- dat_ma[[study]]
          if (temp$median == TRUE) {
            fit_qe_study <- apply_method(method = 'qe', quants = temp$dat, scenario = scenario, n = temp$n, nboot = nboot)
            fit_bc_study <- apply_method(method = 'bc', quants = temp$dat, scenario = scenario, n = temp$n, nboot = nboot)
            fit_mln_study <- apply_method(method = 'mln', quants = temp$dat, scenario = scenario, n = temp$n, nboot = nboot)
            
            yi_qe[study] <- fit_qe_study$mean_est
            yi_bc[study] <- fit_bc_study$mean_est
            yi_mln[study] <- fit_mln_study$mean_est
            
            vi_qe_naive[study] <- fit_qe_study$se_est_naive^2
            vi_qe_boot[study] <- fit_qe_study$se_est_boot^2
            vi_bc_naive[study] <- fit_bc_study$se_est_naive^2
            vi_bc_boot[study] <- fit_bc_study$se_est_boot^2
            vi_mln_naive[study] <- fit_mln_study$se_est_naive^2
            vi_mln_boot[study] <- fit_mln_study$se_est_boot^2
          } else {
            yi_qe[study] <- yi_bc[study] <- yi_mln[study] <- temp$dat[1]
            vi_qe_naive[study] <- vi_qe_boot[study] <- 
              vi_bc_naive[study] <- vi_bc_boot[study] <- 
              vi_mln_naive[study] <- vi_mln_boot[study] <- temp$dat[2]^2 / temp$n
          } 
        }
        
        if (random_effect){
          method <- 'REML'
        } else {
          method <- 'FE'
        }
        fit_qe_naive <- rma.uni_robust(yi = yi_qe, vi = vi_qe_naive, method = method)
        fit_qe_boot <- rma.uni_robust(yi = yi_qe, vi = vi_qe_boot, method = method)
        fit_bc_naive <- rma.uni_robust(yi = yi_bc, vi = vi_bc_naive, method = method)
        fit_bc_boot <- rma.uni_robust(yi = yi_bc, vi = vi_bc_boot, method = method)
        fit_mln_naive <- rma.uni_robust(yi = yi_mln, vi = vi_mln_naive, method = method)
        fit_mln_boot <- rma.uni_robust(yi = yi_mln, vi = vi_mln_boot, method = method)
        
        mycol <- c('pooled_mean_est', 'pooled_mean_ci_lb', 'pooled_mean_ci_ub')
        res_qe_naive[i, mycol] <- extract_pooled_mean(fit_qe_naive)
        res_qe_boot[i, mycol] <- extract_pooled_mean(fit_qe_boot)
        res_bc_naive[i, mycol] <- extract_pooled_mean(fit_bc_naive)
        res_bc_boot[i, mycol] <- extract_pooled_mean(fit_bc_boot)
        res_mln_naive[i, mycol] <- extract_pooled_mean(fit_mln_naive)
        res_mln_boot[i, mycol] <- extract_pooled_mean(fit_mln_boot)
        
        mycol <- c('tausq_est', 'tausq_ci_lb', 'tausq_ci_ub')
        if (random_effect){
          res_qe_naive[i, mycol] <- extract_tausq(fit_qe_naive)
          res_qe_boot[i, mycol] <- extract_tausq(fit_qe_boot)
          res_bc_naive[i, mycol] <- extract_tausq(fit_bc_naive)
          res_bc_boot[i, mycol] <- extract_tausq(fit_bc_boot)
          res_mln_naive[i, mycol] <- extract_tausq(fit_mln_naive)
          res_mln_boot[i, mycol] <- extract_tausq(fit_mln_boot)
          
          res_qe_naive[i, 'I2'] <- fit_qe_naive$I2
          res_qe_boot[i, 'I2'] <- fit_qe_boot$I2
          res_bc_naive[i, 'I2'] <- fit_bc_naive$I2
          res_bc_boot[i, 'I2'] <- fit_bc_boot$I2
          res_mln_naive[i, 'I2'] <- fit_mln_naive$I2
          res_mln_boot[i, 'I2'] <- fit_mln_boot$I2
        } else {
          res_qe_naive[i, mycol] <- res_qe_boot[i, mycol] <- res_bc_naive[i, mycol] <- 
            res_bc_boot[i, mycol] <- res_mln_naive[i, mycol] <- res_mln_boot[i, mycol] <- 0
          
          res_qe_naive[i, 'I2'] <- res_qe_boot[i, 'I2'] <- res_bc_naive[i, 'I2'] <- 
            res_bc_boot[i, 'I2'] <- res_mln_naive[i, 'I2'] <- res_mln_boot[i, 'I2'] <- 0
        }
      }
      batchres <- list(res_qe_naive = res_qe_naive, res_qe_boot = res_qe_boot, 
                       res_bc_naive = res_bc_naive, res_bc_boot = res_bc_boot, 
                       res_mln_naive = res_mln_naive, res_mln_boot = res_mln_boot)
      save(batchres, file = paste0('../results/SimScenario=', m, 'Batch=', batches[t], '.RData'))
    }
  }
}

# Simulate parametric bootstrap data
sim_boot_data <- function(method, fit, n){
  if (method == 'qe'){
    selected.dist <- names(which.min(fit$values))
    if (selected.dist == 'normal'){
      dat <- rnorm(n = n, mean = fit$norm.par[1], sd = fit$norm.par[2])
    } else if (selected.dist == 'log-normal'){
      dat <- rlnorm(n = n, meanlog = fit$lnorm.par[1], sdlog = fit$lnorm.par[2])
    } else if (selected.dist == 'gamma'){
      dat <- rgamma(n = n, shape = fit$gamma.par[1], rate = fit$gamma.par[2])
    } else if (selected.dist == 'weibull'){
      dat <- rweibull(n = n, shape = fit$weibull.par[1], scale = fit$weibull.par[2])
    } else if (selected.dist == 'beta'){
      dat <- rbeta(n = n, shape1 = fit$beta.par[1], shape2 = fit$beta.par[2])
    }
  } else if (method == 'bc' | method == 'mln'){
    dat_trans <- rnorm(n = n, mean = fit$location, sd = fit$scale)
    dat <- estmeansd:::inv.smooth.bc.transform(transformed.vals = dat_trans, lambda = fit$shape)
  }
  all_quants <- unname(quantile(dat, probs = c(0, 0.25, 0.5, 0.75, 1)))
  return(all_quants)
}

# Apply qe or bc method without bootstrap
apply_method_helper <- function(method, quants, scenario, n){
  min.val <- quants[1]
  q1.val <- quants[2]
  med.val <- quants[3]
  q3.val <- quants[4]
  max.val <- quants[5]
  
  if (scenario == 'S1'){
    q1.val <- NA; q3.val <- NA
  } else if (scenario == 'S2'){
    min.val <- NA; max.val <- NA
  }
  
  if (method == 'qe'){
    fit <- suppressWarnings(
      qe.mean.sd(min.val = min.val, q1.val = q1.val, med.val = med.val, 
                 q3.val = q3.val, max.val = max.val, n = n))
  } else if (method == 'qe fit'){
    fit <- suppressWarnings(
      qe.fit(min.val = min.val, q1.val = q1.val, med.val = med.val, 
             q3.val = q3.val, max.val = max.val, n = n)
    )
  } else if (method == 'bc'){
    fit <- suppressWarnings(
      bc.mean.sd(min.val = min.val, q1.val = q1.val, med.val = med.val, 
                 q3.val = q3.val, max.val = max.val, n = n)
    )
  } else if (method == 'mln'){
    fit <- suppressWarnings(
      mln.mean.sd(min.val = min.val, q1.val = q1.val, med.val = med.val, 
                 q3.val = q3.val, max.val = max.val, n = n)
    )
  }
  return(fit)
}

# Apply qe or bc method with possible bootstrap
apply_method <- function(method, quants, scenario, n, nboot){
  if (scenario == 'S1' | scenario == 'S3'){
    shift.param1 <- ifelse(quants[1] <= 0, 0.5 - quants[1], 0)
  } else if (scenario == 'S2'){
    shift.param1 <- ifelse(quants[2] <= 0, 0.5 - quants[2], 0)
  }
  quants <- quants + shift.param1
  
  fit.temp <- apply_method_helper(method = method, quants = quants, scenario = scenario, n = n)
  mean_est <- fit.temp$est.mean - shift.param1
  se_est_naive <- fit.temp$est.sd / sqrt(n)
  if (method == 'qe'){
    selected_dist <- fit.temp$selected.dist
  } else {
    selected_dist <- NA
  }
  
  if (nboot > 0){
    if (method == 'qe'){
      fit <- apply_method_helper(method = 'qe fit', quants = quants, scenario = scenario, n = n)
    } else if (method == 'bc' | method == 'mln'){
      fit <- fit.temp
    }
    res_boot <- rep(NA, nboot)
    for (b in 1:nboot){
      dat_boot <- sim_boot_data(method = method, fit = fit, n = n)
      if (scenario == 'S1' | scenario == 'S3'){
        shift.param2 <- ifelse(dat_boot[1] <= 0, 0.5 - dat_boot[1], 0)
      } else if (scenario == 'S2'){
        shift.param2 <- ifelse(dat_boot[2] <= 0, 0.5 - dat_boot[2], 0)
      }
      dat_boot <- dat_boot + shift.param2
      fit_boot <- apply_method_helper(method = method, quants = dat_boot, scenario = scenario, n = n)
      res_boot[b] <- fit_boot$est.mean - shift.param2
    }
    se_est_boot <- sd(res_boot)
  } else {
    se_est_boot <- NA
  }
  return(list(mean_est = mean_est, se_est_naive = se_est_naive, se_est_boot = se_est_boot, 
              selected_dist = selected_dist))
}

# Simulate data under given distribution. Returns the 5 number summary or mean/sd
sim_data <- function(dist, par1, par2, n, medians, re = 0){
  if (dist == 'normal'){
    dat <- rnorm(n = n, mean = par1, sd = par2) + re
  } else if (dist == 'log-normal'){
    dat <- rlnorm(n = n, meanlog = par1, sdlog = par2) + re
  } else if (dist == 'halfnormal'){
    dat <- rhalfnorm(n = n, theta = par1) + re
  }
  if (medians){
    res <- unname(quantile(dat, probs = c(0, 0.25, 0.5, 0.75, 1)))
  } else {
    res <- c(mean(dat), sd(dat))
  }
  return(res)
}

sim_data_ma <- function(dist, par1, par2, n_studies, prop_medians, tausq){
  res <- vector(mode = "list", length = n_studies)
  sample_sizes <- round(runif(n = n_studies, min = 100, max = 500))
  medians <- rep(FALSE, times = n_studies)
  if (prop_medians > 0){
    medians[1:round(prop_medians * n_studies)] <- TRUE
  }
  if (tausq > 0){
    re <- rnorm(n_studies, mean = 0, sd = sqrt(tausq))
  } else {
    re <- rep(0, n_studies)
  }
  for (study in 1:n_studies){
    temp <- sim_data(dist = dist, par1 = par1, par2 = par2, n = sample_sizes[study], medians = medians[study], re = re[study])
    res[[study]] <- list(dat = temp, n = sample_sizes[study], medians = medians[study], re = re[study])
  }
  return(res)
}


get_table_studylevel <- function(res_est, res_truth, measure){
  ind <- 1:nrow(sim_scenarios_studylevel)
  
  
  res <- data.frame(dist = rep(NA, length(ind)), 
                    par1 = rep(NA, length(ind)), 
                    par2 = rep(NA, length(ind)), 
                    n = rep(NA, length(ind)), 
                    qe_naive = rep(NA, length(ind)), 
                    qe_boot = rep(NA, length(ind)), 
                    bc_naive = rep(NA, length(ind)), 
                    bc_boot = rep(NA, length(ind)), 
                    mln_naive = rep(NA, length(ind)), 
                    mln_boot = rep(NA, length(ind))) 

  row_ind <- 1
  for (i in ind){
    res_est_temp <- res_est[[i]]; truth <- res_truth[i,]
    res[row_ind, 'dist'] <- sim_scenarios_studylevel[i, 'dist']
    res[row_ind, 'par1'] <- sim_scenarios_studylevel[i, 'par1']
    res[row_ind, 'par2'] <- sim_scenarios_studylevel[i, 'par2']
    res[row_ind, 'n'] <- sim_scenarios_studylevel[i, 'n']
    
    if (measure == 'mre'){
      res[row_ind, 'qe_naive'] <- get_mre(res_est_temp$res_qe$se_est_naive, truth['qe'])
      res[row_ind, 'qe_boot'] <- get_mre(res_est_temp$res_qe$se_est_boot, truth['qe'])
      res[row_ind, 'bc_naive'] <- get_mre(res_est_temp$res_bc$se_est_naive, truth['bc'])
      res[row_ind, 'bc_boot'] <- get_mre(res_est_temp$res_bc$se_est_boot, truth['bc'])
      res[row_ind, 'mln_naive'] <- get_mre(res_est_temp$res_mln$se_est_naive, truth['mln'])
      res[row_ind, 'mln_boot'] <- get_mre(res_est_temp$res_mln$se_est_boot, truth['mln'])
      row_ind <- row_ind + 1
    } else if (measure == 'are'){
      res[row_ind, 'qe_naive'] <- get_are(res_est_temp$res_qe$se_est_naive, truth['qe'])
      res[row_ind, 'qe_boot'] <- get_are(res_est_temp$res_qe$se_est_boot, truth['qe'])
      res[row_ind, 'bc_naive'] <- get_are(res_est_temp$res_bc$se_est_naive, truth['bc'])
      res[row_ind, 'bc_boot'] <- get_are(res_est_temp$res_bc$se_est_boot, truth['bc'])
      res[row_ind, 'mln_naive'] <- get_are(res_est_temp$res_mln$se_est_naive, truth['mln'])
      res[row_ind, 'mln_boot'] <- get_are(res_est_temp$res_mln$se_est_boot, truth['mln'])
      row_ind <- row_ind + 1
    } else if (measure == 'rmse'){
      res[row_ind, 'qe_naive'] <- get_rmse(res_est_temp$res_qe$se_est_naive, truth['qe'])
      res[row_ind, 'qe_boot'] <- get_rmse(res_est_temp$res_qe$se_est_boot, truth['qe'])
      res[row_ind, 'bc_naive'] <- get_rmse(res_est_temp$res_bc$se_est_naive, truth['bc'])
      res[row_ind, 'bc_boot'] <- get_rmse(res_est_temp$res_bc$se_est_boot, truth['bc'])
      res[row_ind, 'mln_naive'] <- get_rmse(res_est_temp$res_mln$se_est_naive, truth['mln'])
      res[row_ind, 'mln_boot'] <- get_rmse(res_est_temp$res_mln$se_est_boot, truth['mln'])
      row_ind <- row_ind + 1
    }
  }
  return(res)
}



get_table_ma <- function(res_ma, dist, par1, par2, scenario, tausq, truth_table_I2){
  ind <- which(sim_scenarios_ma$dist == dist & 
                 sim_scenarios_ma$par1 == par1 & 
                 sim_scenarios_ma$par2 == par2 & 
                 sim_scenarios_ma$scenario == scenario & 
                 sim_scenarios_ma$tausq == tausq)
  truth_pooled_mean <- get_true_pooled_mean(dist = dist, par1 = par1, par2 = par2)
  truth_tausq <- tausq
  
  pooled_mean_bias <- pooled_mean_var <- tausq_bias <- tausq_var <- pooled_mean_cov <- tausq_cov <- I2_bias <- I2_median_bias <- matrix(NA, nrow = length(ind), ncol = 8)
  colnames(pooled_mean_bias) <- colnames(pooled_mean_var) <- colnames(tausq_bias) <- colnames(tausq_var) <- colnames(pooled_mean_cov) <- colnames(tausq_cov) <- colnames(I2_bias) <- colnames(I2_median_bias) <- 
    c('K', 'p', 'QE', 'QE Boot', 'BC', 'BC Boot', 'MLN', 'MLN Boot')
  
  row_ind <- 1
  for (i in ind){
    t <- res_ma[[i]]
    truth_I2 <- truth_table_I2[i, ]

    pooled_mean_bias[row_ind, 'K'] <- pooled_mean_var[row_ind, 'K'] <- pooled_mean_cov[row_ind, 'K'] <- 
      tausq_bias[row_ind, 'K'] <- tausq_var[row_ind, 'K'] <- tausq_cov[row_ind, 'K'] <- 
      I2_bias[row_ind, 'K'] <- I2_median_bias[row_ind, 'K'] <- sim_scenarios_ma[i, ]$n_studies
    
    pooled_mean_bias[row_ind, 'p'] <- pooled_mean_var[row_ind, 'p'] <- pooled_mean_cov[row_ind, 'p'] <- 
      tausq_bias[row_ind, 'p'] <- tausq_var[row_ind, 'p'] <- tausq_cov[row_ind, 'p'] <- 
      I2_bias[row_ind, 'p'] <- I2_median_bias[row_ind, 'p'] <- sim_scenarios_ma[i, ]$prop_medians

    pooled_mean_bias[row_ind, 'QE'] <- get_bias(t$res_qe_naive$pooled_mean_est, truth_pooled_mean)
    pooled_mean_bias[row_ind, 'QE Boot'] <- get_bias(t$res_qe_boot$pooled_mean_est, truth_pooled_mean)
    pooled_mean_bias[row_ind, 'BC'] <- get_bias(t$res_bc_naive$pooled_mean_est, truth_pooled_mean)
    pooled_mean_bias[row_ind, 'BC Boot'] <- get_bias(t$res_bc_boot$pooled_mean_est, truth_pooled_mean)
    pooled_mean_bias[row_ind, 'MLN'] <- get_bias(t$res_mln_naive$pooled_mean_est, truth_pooled_mean)
    pooled_mean_bias[row_ind, 'MLN Boot'] <- get_bias(t$res_mln_boot$pooled_mean_est, truth_pooled_mean)

    pooled_mean_var[row_ind, 'QE'] <- var(t$res_qe_naive$pooled_mean_est)
    pooled_mean_var[row_ind, 'QE Boot'] <- var(t$res_qe_boot$pooled_mean_est)
    pooled_mean_var[row_ind, 'BC'] <- var(t$res_bc_naive$pooled_mean_est)
    pooled_mean_var[row_ind, 'BC Boot'] <- var(t$res_bc_boot$pooled_mean_est)
    pooled_mean_var[row_ind, 'MLN'] <- var(t$res_mln_naive$pooled_mean_est)
    pooled_mean_var[row_ind, 'MLN Boot'] <- var(t$res_mln_boot$pooled_mean_est)

    tausq_bias[row_ind, 'QE'] <- get_bias(t$res_qe_naive$tausq_est, truth_tausq)
    tausq_bias[row_ind, 'QE Boot'] <- get_bias(t$res_qe_boot$tausq_est, truth_tausq)
    tausq_bias[row_ind, 'BC'] <- get_bias(t$res_bc_naive$tausq_est, truth_tausq)
    tausq_bias[row_ind, 'BC Boot'] <- get_bias(t$res_bc_boot$tausq_est, truth_tausq)
    tausq_bias[row_ind, 'MLN'] <- get_bias(t$res_mln_naive$tausq_est, truth_tausq)
    tausq_bias[row_ind, 'MLN Boot'] <- get_bias(t$res_mln_boot$tausq_est, truth_tausq)

    tausq_var[row_ind, 'QE'] <- var(t$res_qe_naive$tausq_est)
    tausq_var[row_ind, 'QE Boot'] <- var(t$res_qe_boot$tausq_est)
    tausq_var[row_ind, 'BC'] <- var(t$res_bc_naive$tausq_est)
    tausq_var[row_ind, 'BC Boot'] <- var(t$res_bc_boot$tausq_est)
    tausq_var[row_ind, 'MLN'] <- var(t$res_mln_naive$tausq_est)
    tausq_var[row_ind, 'MLN Boot'] <- var(t$res_mln_boot$tausq_est)

    pooled_mean_cov[row_ind, 'QE'] <- get_cov(t$res_qe_naive$pooled_mean_ci_lb, t$res_qe_naive$pooled_mean_ci_ub, truth_pooled_mean)
    pooled_mean_cov[row_ind, 'QE Boot'] <- get_cov(t$res_qe_boot$pooled_mean_ci_lb, t$res_qe_boot$pooled_mean_ci_ub, truth_pooled_mean)
    pooled_mean_cov[row_ind, 'BC'] <- get_cov(t$res_bc_naive$pooled_mean_ci_lb, t$res_bc_naive$pooled_mean_ci_ub, truth_pooled_mean)
    pooled_mean_cov[row_ind, 'BC Boot'] <- get_cov(t$res_bc_boot$pooled_mean_ci_lb, t$res_bc_boot$pooled_mean_ci_ub, truth_pooled_mean)
    pooled_mean_cov[row_ind, 'MLN'] <- get_cov(t$res_mln_naive$pooled_mean_ci_lb, t$res_mln_naive$pooled_mean_ci_ub, truth_pooled_mean)
    pooled_mean_cov[row_ind, 'MLN Boot'] <- get_cov(t$res_mln_boot$pooled_mean_ci_lb, t$res_mln_boot$pooled_mean_ci_ub, truth_pooled_mean)

    tausq_cov[row_ind, 'QE'] <- get_cov(t$res_qe_naive$tausq_ci_lb, t$res_qe_naive$tausq_ci_ub, truth_tausq)
    tausq_cov[row_ind, 'QE Boot'] <- get_cov(t$res_qe_boot$tausq_ci_lb, t$res_qe_boot$tausq_ci_ub, truth_tausq)
    tausq_cov[row_ind, 'BC'] <- get_cov(t$res_bc_naive$tausq_ci_lb, t$res_bc_naive$tausq_ci_ub, truth_tausq)
    tausq_cov[row_ind, 'BC Boot'] <- get_cov(t$res_bc_boot$tausq_ci_lb, t$res_bc_boot$tausq_ci_ub, truth_tausq)
    tausq_cov[row_ind, 'MLN'] <- get_cov(t$res_mln_naive$tausq_ci_lb, t$res_mln_naive$tausq_ci_ub, truth_tausq)
    tausq_cov[row_ind, 'MLN Boot'] <- get_cov(t$res_mln_boot$tausq_ci_lb, t$res_mln_boot$tausq_ci_ub, truth_tausq)

    I2_bias[row_ind, 'QE'] <- get_bias(t$res_qe_naive$I2, truth_I2['qe'])
    I2_bias[row_ind, 'QE Boot'] <- get_bias(t$res_qe_boot$I2, truth_I2['qe'])
    I2_bias[row_ind, 'BC'] <- get_bias(t$res_bc_naive$I2, truth_I2['bc'])
    I2_bias[row_ind, 'BC Boot'] <- get_bias(t$res_bc_boot$I2, truth_I2['bc'])
    I2_bias[row_ind, 'MLN'] <- get_bias(t$res_mln_naive$I2, truth_I2['mln'])
    I2_bias[row_ind, 'MLN Boot'] <- get_bias(t$res_mln_boot$I2, truth_I2['mln'])
    
    I2_median_bias[row_ind, 'QE'] <- get_median_bias(t$res_qe_naive$I2, truth_I2['qe'])
    I2_median_bias[row_ind, 'QE Boot'] <- get_median_bias(t$res_qe_boot$I2, truth_I2['qe'])
    I2_median_bias[row_ind, 'BC'] <- get_median_bias(t$res_bc_naive$I2, truth_I2['bc'])
    I2_median_bias[row_ind, 'BC Boot'] <- get_median_bias(t$res_bc_boot$I2, truth_I2['bc'])
    I2_median_bias[row_ind, 'MLN'] <- get_median_bias(t$res_mln_naive$I2, truth_I2['mln'])
    I2_median_bias[row_ind, 'MLN Boot'] <- get_median_bias(t$res_mln_boot$I2, truth_I2['mln'])
    
    row_ind <- row_ind + 1
  }
  
  new_order <- order(pooled_mean_bias[, 'K'], pooled_mean_bias[, 'p'])
  pooled_mean_bias <- pooled_mean_bias[new_order,]
  pooled_mean_var <- pooled_mean_var[new_order,]
  tausq_bias <- tausq_bias[new_order,]
  tausq_var <- tausq_var[new_order,]
  pooled_mean_cov <- pooled_mean_cov[new_order,]
  tausq_cov <- tausq_cov[new_order,]
  I2_bias <- I2_bias[new_order,]
  I2_median_bias <- I2_median_bias[new_order,]
  
  return(list(pooled_mean_bias = pooled_mean_bias, pooled_mean_var = pooled_mean_var, 
              tausq_bias = tausq_bias, tausq_var = tausq_var, 
              pooled_mean_cov = pooled_mean_cov, tausq_cov = tausq_cov, 
              I2_bias = I2_bias, I2_median_bias = I2_median_bias))
}

get_true_pooled_mean <- function(dist = dist, par1 = par1, par2 = par2){
  if (dist == 'log-normal'){
    res <- exp(par1 + par2^2 / 2)
  }
  return(res)
}

rma.uni_robust <- function(yi, vi, method){
  if (any(is.na(vi))){
    res <- NULL
  } else {
    res <- try(metafor::rma.uni(yi = yi, vi = vi, method = method), silent = TRUE)
    if ('try-error' %in% class(res)){
      res <- metafor::rma.uni(yi = yi, vi = vi, method = 'DL')
    }
  }
  return(res)
}

get_ma_seeds <- function(nscenarios, n_all_batches, seed_init){
  set.seed(seed_init)
  myseeds <- vector("list", nscenarios)
  for (s in 1:nscenarios){
    myseeds[[s]] <- sample.int(2^30, size = n_all_batches)
  }
  return(myseeds)
}

unpack_res <- function(studylevel, truth_mc = FALSE){
  if (studylevel){
    if (truth_mc){
      batches <- all_batches_studylevel2
    } else {
      batches <- all_batches_studylevel1
    }
    M <- nrow(sim_scenarios_studylevel)
    prefix <- 'StudyLevel_'; suffix <- paste0('Truth=', truth_mc)
  } else {
    batches <- all_batches_ma
    M <- nrow(sim_scenarios_ma)
    prefix <- ''; suffix <- ''
  }
  res <- vector("list", M)
  for (m in 1:M){
    for (t in 1:length(batches)){
      load(paste0('../results/', prefix, 'SimScenario=', m, 'Batch=', batches[t], suffix, '.RData'))
      if (t == 1){
        temp <- batchres 
      } else {
        for (listitem in 1:length(temp)){
          temp[[listitem]] <- rbind(temp[[listitem]], batchres[[listitem]])
        }
      }
    }
    res[[m]] <- temp
  }
  return(res)
}

get_truth_table <- function(res){
  M <- nrow(sim_scenarios_studylevel)
  res_table <- matrix(NA, nrow = M, ncol = 8)
  colnames(res_table) <- c('qe', 'bc', 'mln', 'hozo', 'wan', 'luo', 'shi', 'yang')
  for (m in 1:M){
    res_table[m, 'qe'] <- sd(res[[m]]$res_qe)
    res_table[m, 'bc'] <- sd(res[[m]]$res_bc)
    res_table[m, 'mln'] <- sd(res[[m]]$res_mln)
    res_table[m, 'hozo'] <- sd(res[[m]]$res_hozo)
    res_table[m, 'wan'] <- sd(res[[m]]$res_wan)
    res_table[m, 'luo'] <- sd(res[[m]]$res_luo)
    res_table[m, 'shi'] <- sd(res[[m]]$res_shi)
    res_table[m, 'yang'] <- sd(res[[m]]$res_yang)
  }
  return(res_table)
}

extract_pooled_mean <- function(fit){
  if (!is.null(fit)){
    res <- c(fit$beta, fit$ci.lb, fit$ci.ub)
  } else {
    res <- c(NA, NA, NA)
  }
  return(res)
}

extract_tausq <- function(fit){
  if (!is.null(fit)){
    res <- confint(fit)[[1]]['tau^2', ]
  } else {
    res <- c(NA, NA, NA)
  }
  return(res)
}

get_cov <- function(lb, ub, truth){
  return(mean(lb <= truth & ub >= truth))
}

get_are <- function(ests, truth){
  return(round(100 * mean((ests - truth) / truth), 0))
}

get_bias <- function(ests, truth){
  return(round(mean(ests - truth), 4))
}

get_median_bias <- function(ests, truth){
  return(round(median(ests - truth), 4))
}

get_mre <- function(ests, truth){
  return(round(100 * median((ests - truth) / truth), 0))
}

get_rmse <- function(ests, truth){
  return(round(sqrt(mean(ests - truth)^2), 4))
}

