rm(list = ls())

load('../results/application.RData')
library('stringr')
library('xtable')

res_tab <- res_tab[!is.na(res_tab$n),]
res_tab$est.ci.qe.naive <- paste0(round(res_tab$est.qe.naive, 2), ' [', round(res_tab$ci.lb.qe.naive, 2), ', ', round(res_tab$ci.ub.qe.naive, 2), ']')
res_tab$est.ci.qe.boot <- paste0(round(res_tab$est.qe.boot, 2), ' [', round(res_tab$ci.lb.qe.boot, 2), ', ', round(res_tab$ci.ub.qe.boot, 2), ']')
res_tab$est.ci.bc.naive <- paste0(round(res_tab$est.bc.naive, 2), ' [', round(res_tab$ci.lb.bc.naive, 2), ', ', round(res_tab$ci.ub.bc.naive, 2), ']')
res_tab$est.ci.bc.boot <- paste0(round(res_tab$est.bc.boot, 2), ' [', round(res_tab$ci.lb.bc.boot, 2), ', ', round(res_tab$ci.ub.bc.boot, 2), ']')
res_tab$est.ci.mln.naive <- paste0(round(res_tab$est.mln.naive, 2), ' [', round(res_tab$ci.lb.mln.naive, 2), ', ', round(res_tab$ci.ub.mln.naive, 2), ']')
res_tab$est.ci.mln.boot <- paste0(round(res_tab$est.mln.boot, 2), ' [', round(res_tab$ci.lb.mln.boot, 2), ', ', round(res_tab$ci.ub.mln.boot, 2), ']')
res_tab[, c('I2.qe.naive.r', 'I2.qe.boot.r', 'I2.bc.naive.r', 'I2.bc.boot.r', 'I2.mln.naive.r', 'I2.mln.boot.r')] <- 
  round(res_tab[, c('I2.qe.naive', 'I2.qe.boot', 'I2.bc.naive', 'I2.bc.boot', 'I2.mln.naive', 'I2.mln.boot')])

xtable(res_tab[c('outcome', 'n', 'n.medians', 'est.ci.mln.naive', 'est.ci.mln.boot')]) # Main text
print(xtable(res_tab[c('outcome', 'est.ci.qe.naive', 'est.ci.qe.boot')]), include.rownames = F) # Appendix
print(xtable(res_tab[c('outcome', 'est.ci.bc.naive', 'est.ci.bc.boot')]), include.rownames = F) # Appendix
xtable(res_tab[c('outcome', 'tau2.qe.naive', 'tau2.qe.boot', 'tau2.bc.naive', 'tau2.bc.boot', 'tau2.mln.naive', 'tau2.mln.boot')], digits = 2) # Appendix
xtable(res_tab[c('outcome', 'I2.qe.naive.r', 'I2.qe.boot.r', 'I2.bc.naive.r', 'I2.bc.boot.r', 'I2.mln.naive.r', 'I2.mln.boot.r')], digits = 0) # Appendix

# Change in tau^2
summary((res_tab$tau2.qe.boot - res_tab$tau2.qe.naive) / res_tab$tau2.qe.naive, na.rm = T)
summary((res_tab$tau2.bc.boot - res_tab$tau2.bc.naive) / res_tab$tau2.bc.naive, na.rm = T)
summary((res_tab$tau2.mln.boot - res_tab$tau2.mln.naive) / res_tab$tau2.mln.naive, na.rm = T)

# Change in I^2
summary(res_tab$I2.qe.naive - res_tab$I2.qe.boot, na.rm = T)
summary(res_tab$I2.bc.naive - res_tab$I2.bc.boot, na.rm = T)
summary(res_tab$I2.mln.naive - res_tab$I2.mln.boot, na.rm = T)


# Before/after plot of I^2

pdf('../results/fig/I2.pdf', width = 12, height = 12)

plot_res <- res_tab
labs <- gsub("(.*)\\(.*", "\\1", plot_res$outcome)
labs <- gsub(' ', '', labs)
labs[labs == 'RespiratoryRate'] <- 'Resp. Rate'
labs[labs == 'SpO2-withoutO2'] <- 'SpO2'
pos1 <- 0.9; pos2 <- 2.1
shift <- 0.15

# QE Panel
par(mfrow = c(1, 3), mar=c(5.1, 5.6, 4.1, 2.1))
cex.main <- 2.5
cex.axis <- cex.lab <- 2; cex.names <- 1.25
n <- nrow(plot_res)
plot(rep(pos1, n), plot_res[1:n,]$I2.qe.naive, xlim = c(0, 3), ylim = c(0, 100), 
     main = 'QE Approach', ylab = expression(hat(I)^2), xaxt = 'n', xlab = '', 
     cex.main = cex.main, cex.lab = cex.lab, cex.axis	= cex.axis, pch = 19)
axis(side = 1, at = c(pos1, pos2), labels = c('Naïve', 'Bootstrap'), cex.axis = cex.axis, line = 1, tick = F)
axis(side = 1, at = c(pos1, pos2), labels = c('', ''))
points(rep(pos2, n), plot_res[1:n,]$I2.qe.boot, pch = 19)
d1 <- c(1.9, # Age
        2.75, # SpO2 - without O2
        -0.4, # Respiratory Rate
        0, #Hemoglobin
        -0.95, #Leukocyte
        -0.4, # Lymphocyte
        2, #Neutrophil
        0.2, #Platelets
        -3.1, # APTT
        0.8, # D-Dimer
        0, #Fibrinogen
        -0.7, #INR
        0, #Prothrombin
        0, #ALAT
        1.3, # Albumin
        -4.1, # ASAT
        -1.85, #LDH
        -0.3, #BUN
        0.6, #Creatinine
        -0.5, #CRP
        0, #IL-6
        1.75, #PCT
        -0.7, #CK
        3 #CK-MB
)
d2 <- c(0.15, # Age
        0.5, # SpO2 - without O2
        0.6, # Respiratory Rate
        0, #Hemoglobin
        -0.4, #Leukocyte
        0.4, # Lymphocyte
        1.25, #Neutrophil
        -0.6, #Platelets
        0.4, # APTT
        0, # D-Dimer
        0, #Fibrinogen
        0, #INR
        0, #Prothrombin
        0, #ALAT
        -0.25, # Albumin
        -0.95, # ASAT
        -0.3, #LDH
        0.25, #BUN
        0, #Creatinine
        0.45, #CRP
        0, #IL-6
        0, #PCT
        -0.4, #CK
        0.5 #CK-MB
)
y1 <- plot_res[1:n,]$I2.qe.naive + d1
y2 <- plot_res[1:n,]$I2.qe.boot + d2
for (i in 1:n){
  lines(c(pos1, pos2), c(plot_res[i, 'I2.qe.naive'], plot_res[i, 'I2.qe.boot']),
        type = 'l', lty = 1)
  lines(c(pos1 - shift - 0.02, pos1), c(y1[i], plot_res[i, 'I2.qe.naive']),
        type = 'l', lty = 2, col = 'blue')
  lines(c(pos2, pos2 +  shift + 0.02), c(plot_res[i, 'I2.qe.boot'], y2[i]),
        type = 'l', lty = 2, col = 'blue')
}
text(rep(pos1 - shift, n), y1, labels = labs, pos = 2, xpd=NA, cex = cex.names)
text(rep(pos2 + shift, n), y2, labels = labs, pos = 4, xpd=NA, cex = cex.names)

# BC Panel
plot(rep(pos1, n), plot_res[1:n,]$I2.bc.naive, xlim = c(0, 3), ylim = c(0, 100), 
     main = 'BC Approach', ylab = expression(hat(I)^2), xaxt = 'n', xlab = '', 
     cex.main = cex.main, cex.lab = cex.lab, cex.axis	= cex.axis, pch = 19)
axis(side = 1, at = c(pos1, pos2), labels = c('Naïve', 'Bootstrap'), cex.axis = cex.axis, line = 1, tick = F)
axis(side = 1, at = c(pos1, pos2), labels = c('', ''))
points(rep(pos2, n), plot_res[1:n,]$I2.bc.boot, pch = 19)
d3 <- c(3, # Age
        3.2, # SpO2 - without O2
        0, # Respiratory Rate
        -0.3, #Hemoglobin
        -0.65, #Leukocyte
        -1.6, # Lymphocyte
        1.25, #Neutrophil
        0, #Platelets
        -2, # APTT
        -0.1, # D-Dimer
        0.3, #Fibrinogen
        -3.2, #INR
        0, #Prothrombin
        0, #ALAT
        2, # Albumin
        -2.3, # ASAT
        -4.2, #LDH
        0.6, #BUN
        1.8, #Creatinine
        0, #CRP
        0, #IL-6
        1.3, #PCT
        1, #CK
        2.25 #CK-MB
)
d4 <- c(0, # Age
        -0.3, # SpO2 - without O2
        0, # Respiratory Rate
        -0.6, #Hemoglobin
        0, #Leukocyte
        -0.5, # Lymphocyte
        0.3, #Neutrophil
        0, #Platelets
        0.15, # APTT
        0, # D-Dimer
        0.6, #Fibrinogen
        0, #INR
        0, #Prothrombin
        0, #ALAT
        0.8, # Albumin
        0, # ASAT
        0.5, #LDH
        -0.5, #BUN
        0.75, #Creatinine
        0, #CRP
        0, #IL-6
        -0.5, #PCT
        0, #CK
        0 #CK-MB
)
y3 <- plot_res[1:n,]$I2.bc.naive + d3
y4 <- plot_res[1:n,]$I2.bc.boot + d4
for (i in 1:n){
  lines(c(pos1, pos2), c(plot_res[i, 'I2.bc.naive'], plot_res[i, 'I2.bc.boot']),
        type = 'l', lty = 1)
  lines(c(pos1 - shift - 0.02, pos1), c(y3[i], plot_res[i, 'I2.bc.naive']),
        type = 'l', lty = 2, col = 'blue')
  lines(c(pos2, pos2 +  shift + 0.02), c(plot_res[i, 'I2.bc.boot'], y4[i]),
        type = 'l', lty = 2, col = 'blue')
}
text(rep(pos1 - shift, n), y3, labels = labs, pos = 2, xpd=NA, cex = cex.names)
text(rep(pos2 + shift, n), y4, labels = labs, pos = 4, xpd=NA, cex = cex.names)

# MLN Panel
plot(rep(pos1, n), plot_res[1:n,]$I2.mln.naive, xlim = c(0, 3), ylim = c(0, 100), 
     main = 'MLN Approach', ylab = expression(hat(I)^2), xaxt = 'n', xlab = '', 
     cex.main = cex.main, cex.lab = cex.lab, cex.axis	= cex.axis, pch = 19)
axis(side = 1, at = c(pos1, pos2), labels = c('Naïve', 'Bootstrap'), cex.axis = cex.axis, line = 1, tick = F)
axis(side = 1, at = c(pos1, pos2), labels = c('', ''))
points(rep(pos2, n), plot_res[1:n,]$I2.mln.boot, pch = 19)
d5 <- c(2, # Age
        0, # SpO2 - without O2
        0, # Respiratory Rate
        0, #Hemoglobin
        0.75, #Leukocyte
        2.95, # Lymphocyte
        1.85, #Neutrophil
        0, #Platelets
        -3.4, # APTT
        -2.2, # D-Dimer
        0.5, #Fibrinogen
        -0.95, #INR
        0, #Prothrombin
        -0.5, #ALAT
        1.1, # Albumin
        -2.5, # ASAT
        -0.25, #LDH
        0, #BUN
        -0.8, #Creatinine
        0, #CRP
        0, #IL-6
        0, #PCT
        -1.15, #CK
        0 #CK-MB
)
d6 <- c(1.85, # Age
        0.7, # SpO2 - without O2
        -0.25, # Respiratory Rate
        0.6, #Hemoglobin
        -0.3, #Leukocyte
        0, # Lymphocyte
        1.75, #Neutrophil
        0, #Platelets
        0.75, # APTT
        -1, # D-Dimer
        0, #Fibrinogen
        0.15, #INR
        0, #Prothrombin
        0, #ALAT
        -0.35, # Albumin
        0, # ASAT
        0.25, #LDH
        0.5, #BUN
        -0.45, #Creatinine
        0, #CRP
        0, #IL-6
        0, #PCT
        -1, #CK
        0 #CK-MB
)
y5 <- plot_res[1:n,]$I2.mln.naive + d5
y6 <- plot_res[1:n,]$I2.mln.boot + d6
for (i in 1:n){
  lines(c(pos1, pos2), c(plot_res[i, 'I2.mln.naive'], plot_res[i, 'I2.mln.boot']),
        type = 'l', lty = 1)
  lines(c(pos1 - shift - 0.02, pos1), c(y5[i], plot_res[i, 'I2.mln.naive']),
        type = 'l', lty = 2, col = 'blue')
  lines(c(pos2, pos2 +  shift + 0.02), c(plot_res[i, 'I2.mln.boot'], y6[i]),
        type = 'l', lty = 2, col = 'blue')
}
text(rep(pos1 - shift, n), y5, labels = labs, pos = 2, xpd=NA, cex = cex.names)
text(rep(pos2 + shift, n), y6, labels = labs, pos = 4, xpd=NA, cex = cex.names)
dev.off()


# Additional IL6 results
load("../results/IL6.RData")

df_with_results$mean.qe.naive <- df_with_results$mean.g1_qe - df_with_results$mean.g2_qe
df_with_results$se.qe.naive <- sqrt(df_with_results$se.g1_qe_naive^2 + df_with_results$se.g2_qe_naive^2)
df_with_results$se.qe.boot <- sqrt(df_with_results$se.g1_qe_boot^2 + df_with_results$se.g2_qe_boot^2)

df_with_results$mean.bc.naive <- df_with_results$mean.g1_bc - df_with_results$mean.g2_bc
df_with_results$se.bc.naive <- sqrt(df_with_results$se.g1_bc_naive^2 + df_with_results$se.g2_bc_naive^2)
df_with_results$se.bc.boot <- sqrt(df_with_results$se.g1_bc_boot^2 + df_with_results$se.g2_bc_boot^2)

df_with_results$mean.mln.naive <- df_with_results$mean.g1_mln - df_with_results$mean.g2_mln
df_with_results$se.mln.naive <- sqrt(df_with_results$se.g1_mln_naive^2 + df_with_results$se.g2_mln_naive^2)
df_with_results$se.mln.boot <- sqrt(df_with_results$se.g1_mln_boot^2 + df_with_results$se.g2_mln_boot^2)

res.naive <- metafor::rma.uni(yi = mean.mln.naive, sei = se.mln.naive, data = df_with_results, method = 'REML')
res.boot <- metafor::rma.uni(yi = mean.mln.naive, sei = se.mln.boot, data = df_with_results, method = 'REML')

df_with_results[, c('author', 'mean.qe.naive', 'se.qe.naive', 'se.qe.boot')]
df_with_results[, c('author', 'mean.bc.naive', 'se.bc.naive', 'se.bc.boot')]

print(xtable(
  cbind(df_with_results[, c('author', 'n.g1', 'n.g2', 'mean.mln.naive', 'se.mln.naive', 'se.mln.boot')], weights(res.naive), weights(res.boot)))
      , include.rownames = F)
