rm(list = ls())

source('mln.R')
source('helper_functions_simulation.R')
source('helper_functions_application.R')

library('estmeansd')
library('metafor')
library('doParallel')
library('foreach')

vars_mortality <- c(
  "Age", "O2SatNoOx", "RespRate", "Hgb", "Leukocyte", "Lymphocyte", "Neutrophil", 
  "Platelets", "APTT", "DDimer", "Fibrinogen", "INR", "Prothrombin", "ALAT", 
  "Albumin", "ASAT", "LDH", "BUN", "Creatinine", "CRP", "IL6", "Procalcitonin", 
  "BNP", "CK", "CKMB", "TroponinI"
)
nboot <- 1000
registerDoParallel(cores=9)

set.seed(1234)
seeds <- sample.int(2^30, size = length(vars_mortality))

options(warn = 1)

temp <- foreach (i = 1:length(vars_mortality)) %dopar%{
  set.seed(seeds[i])
  varname <- vars_mortality[i]
  load(paste0('../data/dat.', varname, '_Mortality.Rda'))
  df <- get(paste0('dat.', varname, '_Mortality'))
  
  # Data processing
  df <- remove_ties(df = df)
  df$bowley.g1 <- (df$q1.g1 + df$q3.g1 - 2 * df$med.g1) / (df$q3.g1 - df$q1.g1)
  df$bowley.g2 <- (df$q1.g2 + df$q3.g2 - 2 * df$med.g2) / (df$q3.g2 - df$q1.g2)
  df <- df[!is.na(df$n.g1) & !is.na(df$n.g2), ]
  ind_smallsize <- (!is.na(df$bowley.g1) & df$n.g1 < 10) | 
    (!is.na(df$bowley.g2) & df$n.g2 < 10)
  df <- df[(is.na(df$bowley.g1) | df$n.g1 >= 10) & (is.na(df$bowley.g2) | df$n.g2 >= 10),]
  ind_skewness <- (is.na(df$bowley.g1) | abs(df$bowley.g1) < 0.75) & 
    (is.na(df$bowley.g2) | abs(df$bowley.g2) < 0.75)
  df <- df[ind_skewness,]
  print(paste0('Removed ', sum(ind_smallsize, na.rm = TRUE), ' for ', varname, ' due to small sample size'))
  print(paste0('Removed ', sum(!ind_skewness, na.rm = TRUE), ' for ', varname, ' due to skewness'))
  print(paste0('Revised ', sum(df$any_ties, na.rm = TRUE), ' for ', varname, ' due to ties'))
  if (varname == 'DDimer'){
    df <- df[!df$author %in% c('Yan, Xisheng (JMV)', 'Xu, Bo (J. Inf.)'),]
  } else if (varname == 'CRP'){
    df <- df[!df$author %in% c('Yan, Xisheng (JMV)'),]
  } else if (varname == 'TroponinI'){
    df <- df[!df$author %in% c('Yan, Xisheng (JMV)'),]
  }
  
  res_tab <- data.frame(outcome = units[units$Outcome == varname,]$Outcome.Name, 
                        n = NA, 
                        n.medians = NA, 
                        est.qe.naive = NA, ci.lb.qe.naive = NA, ci.ub.qe.naive = NA, 
                        est.qe.boot = NA, ci.lb.qe.boot = NA, ci.ub.qe.boot = NA, 
                        est.bc.naive = NA, ci.lb.bc.naive = NA, ci.ub.bc.naive = NA, 
                        est.bc.boot = NA, ci.lb.bc.boot = NA, ci.ub.bc.boot = NA, 
                        est.mln.naive = NA, ci.lb.mln.naive = NA, ci.ub.mln.naive = NA, 
                        est.mln.boot = NA, ci.lb.mln.boot = NA, ci.ub.mln.boot = NA, 
                        tau2.qe.naive = NA, 
                        tau2.qe.boot = NA, 
                        tau2.bc.naive = NA, 
                        tau2.bc.boot = NA, 
                        tau2.mln.naive = NA, 
                        tau2.mln.boot = NA, 
                        I2.qe.naive = NA, 
                        I2.qe.boot = NA, 
                        I2.bc.naive = NA, 
                        I2.bc.boot = NA, 
                        I2.mln.naive = NA, 
                        I2.mln.boot = NA)
  if (nrow(df) >= 6){
    res <- mean_analysis(df, nboot = nboot)
    
    res_tab[1, 'n'] <- res[[1]]$k
    res_tab[1, 'n.medians'] <- paste0(round(res$prop_medians * res[[1]]$k), 
                                      ' (' ,round(100 * res$prop_medians), '%)')
    
    res_tab[1, 'est.qe.naive'] <- res$qe_naive$beta
    res_tab[1, 'ci.lb.qe.naive'] <- res$qe_naive$ci.lb
    res_tab[1, 'ci.ub.qe.naive'] <- res$qe_naive$ci.ub
    
    res_tab[1, 'est.qe.boot'] <- res$qe_boot$beta
    res_tab[1, 'ci.lb.qe.boot'] <- res$qe_boot$ci.lb
    res_tab[1, 'ci.ub.qe.boot'] <- res$qe_boot$ci.ub
    
    res_tab[1, 'est.bc.naive'] <- res$bc_naive$beta
    res_tab[1, 'ci.lb.bc.naive'] <- res$bc_naive$ci.lb
    res_tab[1, 'ci.ub.bc.naive'] <- res$bc_naive$ci.ub
    
    res_tab[1, 'est.bc.boot'] <- res$bc_boot$beta
    res_tab[1, 'ci.lb.bc.boot'] <- res$bc_boot$ci.lb
    res_tab[1, 'ci.ub.bc.boot'] <- res$bc_boot$ci.ub
    
    res_tab[1, 'est.mln.naive'] <- res$mln_naive$beta
    res_tab[1, 'ci.lb.mln.naive'] <- res$mln_naive$ci.lb
    res_tab[1, 'ci.ub.mln.naive'] <- res$mln_naive$ci.ub
    
    res_tab[1, 'est.mln.boot'] <- res$mln_boot$beta
    res_tab[1, 'ci.lb.mln.boot'] <- res$mln_boot$ci.lb
    res_tab[1, 'ci.ub.mln.boot'] <- res$mln_boot$ci.ub
    
    res_tab[1, 'tau2.qe.naive'] <- res$qe_naive$tau2
    res_tab[1, 'tau2.qe.boot'] <- res$qe_boot$tau2
    res_tab[1, 'tau2.bc.naive'] <- res$bc_naive$tau2
    res_tab[1, 'tau2.bc.boot'] <- res$bc_boot$tau2
    res_tab[1, 'tau2.mln.naive'] <- res$mln_naive$tau2
    res_tab[1, 'tau2.mln.boot'] <- res$mln_boot$tau2
    
    res_tab[1, 'I2.qe.naive'] <- res$qe_naive$I2
    res_tab[1, 'I2.qe.boot'] <- res$qe_boot$I2
    res_tab[1, 'I2.bc.naive'] <- res$bc_naive$I2
    res_tab[1, 'I2.bc.boot'] <- res$bc_boot$I2
    res_tab[1, 'I2.mln.naive'] <- res$mln_naive$I2
    res_tab[1, 'I2.mln.boot'] <- res$mln_boot$I2
    
    df_with_results <- res$df
    save(df_with_results, file = paste0('../results/', varname, '.RData'))
  }
  return(res_tab)
}

res_tab <- do.call(rbind.data.frame, temp)

save.image('../results/application.RData')
