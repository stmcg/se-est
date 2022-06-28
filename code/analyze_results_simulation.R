rm(list = ls())

library('xtable')
library('fdrtool')
load('../data/simulation_settings.RData')

res_studylevel_est <- unpack_res(studylevel = TRUE, truth_mc = FALSE)
res_studylevel_truth <- unpack_res(studylevel = TRUE, truth_mc = TRUE)
res_studylevel_truth_small <- get_truth_table(res = res_studylevel_truth)
res_ma <- unpack_res(studylevel = FALSE)
load('../results/true_I2.RData')

## Simple Illustration
ind <- which(sim_scenarios_studylevel$dist == 'log-normal' & 
               sim_scenarios_studylevel$par1 == 5 & 
               sim_scenarios_studylevel$par2 == 0.25 & 
               sim_scenarios_studylevel$n == 1000 & 
               sim_scenarios_studylevel$scenario == 'S1')
res_est_temp <- res_studylevel_est[[ind]]; res_truth_temp <- res_studylevel_truth_small[ind,]
naive_se_bestcase <- sqrt( (exp(0.25^2) - 1) * exp(2 * 5 + 0.25^2) ) / sqrt(1000)

pdf('../results/fig/Illustration_Naive.pdf', height = 3, width = 6)

par(mfrow = c(1, 3), mar=c(3.1, 4.1, 3.1, 2.1))
ylim <- c(1, 3.9)
boxplot(res_est_temp$res_qe$se_est_naive, ylab = 'Estimated SE', 
        main = 'QE Approach', pch = 20, ylim = ylim)
abline(h = res_truth_temp['qe'], lty = 2, col = 'red')

boxplot(res_est_temp$res_bc$se_est_naive, ylab = 'Estimated SE', 
        main = 'BC Approach', pch = 20, ylim = ylim)
abline(h = res_truth_temp['bc'], lty = 2, col = 'red')

boxplot(res_est_temp$res_mln$se_est_naive, ylab = 'Estimated SE', 
        main = 'MLN Approach', pch = 20, ylim = ylim)
abline(h = res_truth_temp['mln'], lty = 2, col = 'red')
dev.off()



pdf('../results/fig/Illustration_Boot.pdf', height = 3, width = 6)

par(mfrow = c(1, 3), mar=c(3.1, 4.1, 3.1, 2.1))
ylim <- c(1, 5.5); names <- c('Naïve', 'Bootstrap')
boxplot(res_est_temp$res_qe$se_est_naive, 
        res_est_temp$res_qe$se_est_boot, ylab = 'Estimated SE', 
        main = 'QE Approach', pch = 20, ylim = ylim, 
        names = names)
abline(h = res_truth_temp['qe'], lty = 2, col = 'red')

boxplot(res_est_temp$res_bc$se_est_naive, 
        res_est_temp$res_bc$se_est_boot, ylab = 'Estimated SE', 
        main = 'BC Approach', pch = 20, ylim = ylim, 
        names = names)
abline(h = res_truth_temp['bc'], lty = 2, col = 'red')

boxplot(res_est_temp$res_mln$se_est_naive, 
        res_est_temp$res_mln$se_est_boot, ylab = 'Estimated SE', 
        main = 'MLN Approach', pch = 20, ylim = ylim, 
        names = names)
abline(h = res_truth_temp['mln'], lty = 2, col = 'red')
dev.off()



pdf('../results/fig/Illustration_OtherMethods.pdf', height = 6.5, width = 6)

par(mfrow = c(3, 2), mar=c(3.1, 4.1, 3.1, 2.1))
ylim <- c(1, 4)
boxplot(res_est_temp$res_hozo$se_est_naive, ylab = 'Estimated SE', 
        main = 'Hozo et al.', pch = 20, ylim = ylim)
abline(h = res_truth_temp['hozo'], lty = 2, col = 'red')

ylim <- c(0.95, 2.5)
boxplot(res_est_temp$res_wan$se_est_naive, ylab = 'Estimated SE', 
        main = 'Wan et al.', pch = 20, ylim = ylim)
abline(h = res_truth_temp['wan'], lty = 2, col = 'red')

boxplot(res_est_temp$res_luo$se_est_naive, ylab = 'Estimated SE', 
        main = 'Luo et al.', pch = 20, ylim = ylim)
abline(h = res_truth_temp['luo'], lty = 2, col = 'red')

boxplot(res_est_temp$res_shi$se_est_naive, ylab = 'Estimated SE', 
        main = 'Shi et al.', pch = 20, ylim = ylim)
abline(h = res_truth_temp['shi'], lty = 2, col = 'red')

boxplot(res_est_temp$res_yang$se_est_naive, ylab = 'Estimated SE', 
        main = 'Yang et al.', pch = 20, ylim = ylim)
abline(h = res_truth_temp['yang'], lty = 2, col = 'red')

dev.off()

sd(rlnorm(5e6, 5, 0.25)) / sqrt(1000)




pdf('../results/fig/Distributions.pdf', height = 3.5, width = 8)
par(mfrow = c(1, 3), mar = c(5.1, 4.1, 4.1, 2.1))

x_lnorm <- seq(from = 0, to = 500, length.out = 1e4)
y_lnorm_1 <- dlnorm(x_lnorm, 5, 0.25)  
y_lnorm_2 <- dlnorm(x_lnorm, 5, 0.5)  
y_lnorm_3 <- dlnorm(x_lnorm, 5, 1)  
plot(x_lnorm, y_lnorm_1, type = 'l', main = 'Log-Normal', xlab = 'x', ylab = 'Density', 
     col = '#E41A1C', lwd = 2)
lines(x_lnorm, y_lnorm_2, col = '#377EB8', lwd = 2)
lines(x_lnorm, y_lnorm_3, col = '#4DAF4A', lwd = 2)

x_norm <- seq(from = 0, to = 10, length.out = 1e4)
y_norm <- dnorm(x_norm, 5, 1)
plot(x_norm, y_norm, type = 'l', main = 'Normal', xlab = 'x', ylab = 'Density', 
     lwd = 2)

x_hnorm <- seq(from = 0, to = 50, length.out = 1e4)
y_hnorm <- dhalfnorm(x_hnorm, 0.1)
plot(x_hnorm, y_hnorm, type = 'l', main = 'Half-Normal', xlab = 'x', ylab = 'Density', 
     lwd = 2)
dev.off()



## Study-level tables

res_studylevel_mre <- get_table_studylevel(res_est = res_studylevel_est, res_truth = res_studylevel_truth_small, measure = 'mre')
res_studylevel_are <- get_table_studylevel(res_est = res_studylevel_est, res_truth = res_studylevel_truth_small, measure = 'are')
res_studylevel_rmse <- get_table_studylevel(res_est = res_studylevel_est, res_truth = res_studylevel_truth_small, measure = 'rmse')

xtable(res_studylevel_mre, digits = 0)
xtable(res_studylevel_are, digits = 0)
xtable(res_studylevel_rmse)


# Meta-analysis tables
res_ma_tables_S1 <- get_table_ma(res_ma = res_ma, dist = 'log-normal', par1 = 5, par2 = 0.25, scenario = 'S1', tausq = 6, truth_table_I2 = 100 * as.matrix(res_table))
res_ma_tables_S2 <- get_table_ma(res_ma = res_ma, dist = 'log-normal', par1 = 5, par2 = 0.25, scenario = 'S2', tausq = 6, truth_table_I2 = 100 * as.matrix(res_table))
res_ma_tables_S3 <- get_table_ma(res_ma = res_ma, dist = 'log-normal', par1 = 5, par2 = 0.25, scenario = 'S3', tausq = 6, truth_table_I2 = 100 * as.matrix(res_table))

xtable(res_ma_tables_S1$tausq_bias, digits = 2)
xtable(res_ma_tables_S2$tausq_bias, digits = 2)
xtable(res_ma_tables_S3$tausq_bias, digits = 2)

xtable(res_ma_tables_S1$tausq_var, digits = 2)
xtable(res_ma_tables_S2$tausq_var, digits = 2)
xtable(res_ma_tables_S3$tausq_var, digits = 2)

xtable(res_ma_tables_S1$tausq_cov, digits = 2)
xtable(res_ma_tables_S2$tausq_cov, digits = 2)
xtable(res_ma_tables_S3$tausq_cov, digits = 2)

xtable(res_ma_tables_S1$pooled_mean_bias, digits = 2)
xtable(res_ma_tables_S2$pooled_mean_bias, digits = 2)
xtable(res_ma_tables_S3$pooled_mean_bias, digits = 2)

xtable(res_ma_tables_S1$pooled_mean_var, digits = 2)
xtable(res_ma_tables_S2$pooled_mean_var, digits = 2)
xtable(res_ma_tables_S3$pooled_mean_var, digits = 2)

xtable(res_ma_tables_S1$pooled_mean_cov, digits = 2)
xtable(res_ma_tables_S2$pooled_mean_cov, digits = 2)
xtable(res_ma_tables_S3$pooled_mean_cov, digits = 2)

xtable(res_ma_tables_S1$I2_bias, digits = 2)
xtable(res_ma_tables_S2$I2_bias, digits = 2)
xtable(res_ma_tables_S3$I2_bias, digits = 2)

xtable(res_ma_tables_S1$I2_median_bias, digits = 2)
xtable(res_ma_tables_S2$I2_median_bias, digits = 2)
xtable(res_ma_tables_S3$I2_median_bias, digits = 2)

