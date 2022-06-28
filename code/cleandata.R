## Loading helper functions and data
rm(list = ls())
load('../data/Age.RData')
load('../data/Blood.RData')
load('../data/Clinical.RData')
load('../data/Coagulation.RData')
load('../data/Heart.RData')
load('../data/Infections.RData')
load('../data/Kidney.RData')
load('../data/Liver.RData')

# Getting the relevant outcome names
all_objects <- ls()
vars_mortality <- c(
  "Age", "O2SatNoOx", "RespRate", "Hgb", "Leukocyte", "Lymphocyte", "Neutrophil", 
  "Platelets", "APTT", "DDimer", "Fibrinogen", "INR", "Prothrombin", "ALAT", 
  "Albumin", "ASAT", "LDH", "BUN", "Creatinine", "CRP", "IL6", "Procalcitonin", 
  "BNP", "CK", "CKMB", "TroponinI"
  )
vars_icu <- c(
  "Age", "RespRate", "Hgb", "Leukocyte", "Lymphocyte", "Neutrophil", 
  "Platelets", "APTT", "DDimer", "Prothrombin", "ALAT", 
  "Albumin", "ASAT", "LDH", "BUN", "Creatinine", "CRP", "Procalcitonin", 
  "CK", "CKMB", "TroponinI"
)

# Step 1: Clean and save data sets
setwd('../data')
for (analysis_num in c(1, 3)){
  if (analysis_num == 1){
    varnames <- vars_mortality
    analysis_type <- 'Mortality'
  } else if (analysis_num == 3){
    varnames <- vars_icu
    analysis_type <- 'ICU'
  }
  for (varname in varnames){
    temp <- get(paste0(varname, analysis_num, 'clean'))
    if (varname == 'DDimer' & analysis_num == 1){
      temp <- temp[!(temp$ID %in% c('7CAQ5CRL')),]
    } else if (varname == 'CRP' & analysis_num == 1){
      temp <- temp[!(temp$ID %in% c('VN52PHAB-VAL', 'VN52PHAB')),]
    }
    temp <- temp[, c('author', 
                     'n.g1', 'q1.g1', 'med.g1', 'q3.g1', 'mean.g1', 'sd.g1', 
                     'n.g2', 'q1.g2', 'med.g2', 'q3.g2', 'mean.g2', 'sd.g2')]
    mystr <- paste0('dat.', varname, '_', analysis_type)
    assign(mystr, temp)
    save(list = mystr, file = paste0(mystr, '.Rda'))
  }
}
