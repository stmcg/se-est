## Loading unit/category link
units <- read.csv('../data/Units_Categories.csv')

## Functions

remove_ties <- function(df, inflation_factor = 0.025){
  for (i in 1:nrow(df)){
    any_ties <- FALSE
    if (!is.na(df[i, 'q1.g1'])){
      if (df[i, 'q1.g1'] == df[i, 'med.g1']){
        df[i, 'med.g1'] <- df[i, 'med.g1'] * (1 + inflation_factor)
        df[i, 'q3.g1'] <- df[i, 'q3.g1'] * (1 + inflation_factor)
        any_ties <- TRUE
      }
      if (df[i, 'med.g1'] == df[i, 'q3.g1']){
        df[i, 'q3.g1'] <- df[i, 'q3.g1'] * (1 + inflation_factor)
        any_ties <- TRUE
      }
    }
    if (!is.na(df[i, 'q1.g2'])){
      if (df[i, 'q1.g2'] == df[i, 'med.g2']){
        df[i, 'med.g2'] <- df[i, 'med.g2'] * (1 + inflation_factor)
        df[i, 'q3.g2'] <- df[i, 'q3.g2'] * (1 + inflation_factor)
        any_ties <- TRUE
      }
      if (df[i, 'med.g2'] == df[i, 'q3.g2']){
        df[i, 'q3.g2'] <- df[i, 'q3.g2'] * (1 + inflation_factor)
        any_ties <- TRUE
      }
    }
    df[i, 'any_ties'] <- any_ties
  }
  return(df)
}

## Applies transformation-based methods
mean_analysis <- function(df, outcome_name = 'NA', nboot){
  method <- 'REML'
  
  for (i in 1:nrow(df)){
    scenario_g1 <- metamedian:::get.scenario(q1.val = df[i,'q1.g1'], 
                                             med.val = df[i,'med.g1'], 
                                             q3.val = df[i,'q3.g1'], 
                                             mean.val = df[i,'mean.g1'], 
                                             sd.val = df[i,'sd.g1'])
    if (scenario_g1 == 'S2'){
      quants <- unname(unlist(c(df[i, c('q1.g1')], df[i, c('q1.g1', 'med.g1', 'q3.g1')], NA)))
      qe <- apply_method(method = 'qe', quants = quants, scenario = 'S2', n = df[i, 'n.g1'], nboot = nboot)
      bc <- apply_method(method = 'bc', quants = quants, scenario = 'S2', n = df[i, 'n.g1'], nboot = nboot)
      mln <- apply_method(method = 'mln', quants = quants, scenario = 'S2', n = df[i, 'n.g1'], nboot = nboot)
      df[i, 'mean.g1_qe'] <- qe$mean_est
      df[i, 'se.g1_qe_naive'] <- qe$se_est_naive
      df[i, 'se.g1_qe_boot'] <- qe$se_est_boot
      df[i, 'mean.g1_bc'] <- bc$mean_est
      df[i, 'se.g1_bc_naive'] <- bc$se_est_naive
      df[i, 'se.g1_bc_boot'] <- bc$se_est_boot
      df[i, 'mean.g1_mln'] <- mln$mean_est
      df[i, 'se.g1_mln_naive'] <- mln$se_est_naive
      df[i, 'se.g1_mln_boot'] <- mln$se_est_boot
    } else {
      df[i, 'mean.g1_qe'] <- df[i, 'mean.g1_bc'] <- df[i, 'mean.g1_mln'] <- df[i, 'mean.g1']
      df[i, 'se.g1_qe_naive'] <- df[i, 'se.g1_qe_boot'] <- 
        df[i, 'se.g1_bc_naive'] <- df[i, 'se.g1_bc_boot'] <- 
        df[i, 'se.g1_mln_naive'] <- df[i, 'se.g1_mln_boot'] <- 
        df[i, 'sd.g1'] / sqrt(df[i, 'n.g1'])
    }
    scenario_g2 <- metamedian:::get.scenario(q1.val = df[i,'q1.g2'], 
                                             med.val = df[i,'med.g2'], 
                                             q3.val = df[i,'q3.g2'], 
                                             mean.val = df[i,'mean.g2'], 
                                             sd.val = df[i,'sd.g2'])
    if (scenario_g2 == 'S2'){
      quants <- unname(unlist(c(df[i, c('q1.g2')], df[i, c('q1.g2', 'med.g2', 'q3.g2')], NA)))
      qe <- apply_method(method = 'qe', quants = quants, scenario = 'S2', n = df[i, 'n.g2'], nboot = nboot)
      bc <- apply_method(method = 'bc', quants = quants, scenario = 'S2', n = df[i, 'n.g2'], nboot = nboot)
      mln <- apply_method(method = 'mln', quants = quants, scenario = 'S2', n = df[i, 'n.g2'], nboot = nboot)
      df[i, 'mean.g2_qe'] <- qe$mean_est
      df[i, 'se.g2_qe_naive'] <- qe$se_est_naive
      df[i, 'se.g2_qe_boot'] <- qe$se_est_boot
      df[i, 'mean.g2_bc'] <- bc$mean_est
      df[i, 'se.g2_bc_naive'] <- bc$se_est_naive
      df[i, 'se.g2_bc_boot'] <- bc$se_est_boot
      df[i, 'mean.g2_mln'] <- mln$mean_est
      df[i, 'se.g2_mln_naive'] <- mln$se_est_naive
      df[i, 'se.g2_mln_boot'] <- mln$se_est_boot
    } else {
      df[i, 'mean.g2_qe'] <- df[i, 'mean.g2_bc'] <- df[i, 'mean.g2_mln'] <- df[i, 'mean.g2']
      df[i, 'se.g2_qe_naive'] <- df[i, 'se.g2_qe_boot'] <- 
        df[i, 'se.g2_bc_naive'] <- df[i, 'se.g2_bc_boot'] <- 
        df[i, 'se.g2_mln_naive'] <- df[i, 'se.g2_mln_boot'] <- 
        df[i, 'sd.g2'] / sqrt(df[i, 'n.g2'])
    }
    
    df[i, 'median_indicator'] <- scenario_g1 == 'S2' & scenario_g2 == 'S2'
  }
  
  qe_naive <- rma.uni(yi = df$mean.g1_qe - df$mean.g2_qe, 
                      vi = df$se.g1_qe_naive^2 + df$se.g2_qe_naive^2, control = list(maxiter=200))
  qe_boot <- rma.uni(yi = df$mean.g1_qe - df$mean.g2_qe, 
                     vi = df$se.g1_qe_boot^2 + df$se.g2_qe_boot^2, control = list(maxiter=200))
  bc_naive <- rma.uni(yi = df$mean.g1_bc - df$mean.g2_bc, 
                      vi = df$se.g1_bc_naive^2 + df$se.g2_bc_naive^2, control = list(maxiter=200))
  bc_boot <- rma.uni(yi = df$mean.g1_bc - df$mean.g2_bc, 
                     vi = df$se.g1_bc_boot^2 + df$se.g2_bc_boot^2, control = list(maxiter=200))
  mln_naive <- rma.uni(yi = df$mean.g1_mln - df$mean.g2_mln, 
                       vi = df$se.g1_mln_naive^2 + df$se.g2_mln_naive^2, control = list(maxiter=200))
  mln_boot <- rma.uni(yi = df$mean.g1_mln - df$mean.g2_mln, 
                      vi = df$se.g1_mln_boot^2 + df$se.g2_mln_boot^2, control = list(maxiter=200))
  res <- list(qe_naive = qe_naive, qe_boot = qe_boot, 
              bc_naive = bc_naive, bc_boot = bc_boot, 
              mln_naive = mln_naive, mln_boot = mln_boot, 
              prop_medians = mean(df$median_indicator), 
              df = df)
  
  return(res)
}


