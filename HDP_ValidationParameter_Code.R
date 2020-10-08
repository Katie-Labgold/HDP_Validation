##---------------------------------------------------------------------------##
## HDP Validation Parameter Code
## K. Labgold
## Last Updated: 09/28/2020
## 
## Code to calculate the Sensitivity, Specificity, PPV, and NPV of
## hypertensive disorders during pregnancy (HDP) identified in hospital 
## discharge record ICD-10 codes, compared with medical record (gold-standard).
##---------------------------------------------------------------------------##


####---- 0. Load Data & Packages ----####

library(dplyr)
library(binom)


hdp_analysis <- readRDS("hdp_analysis") # Data not publicly available



####---- 1. Write Function to Calculate Validation Parameters ----####

# Write a function - compare hospital discharge (HDD) to medical record (MR)
validation.param.function <- function(hdd_var, mr_var) {
                                a_v <- sum(hdp_analysis[,hdd_var] == 1 & hdp_analysis[,mr_var] == 1) # TP
                                b_v <- sum(hdp_analysis[,hdd_var] == 1 & hdp_analysis[,mr_var] == 0) # FP
                                c_v <- sum(hdp_analysis[,hdd_var] == 0 & hdp_analysis[,mr_var] == 1) # FN
                                d_v <- sum(hdp_analysis[,hdd_var] == 0 & hdp_analysis[,mr_var] == 0) # TN
  
  
                              # Sensitivity ---
                              se <- paste0(round((a_v / (a_v + c_v))*100, 1), "%") # TP / (TP + FN)
                              
                              se_low <- round((binom.confint(a_v, (a_v + c_v), 
                                                             conf.level = 0.95, methods = "exact")$lower)*100, 1)
                              
                              se_up <- round((binom.confint(a_v, (a_v + c_v), 
                                                            conf.level = 0.95, methods = "exact")$upper)*100, 1)
                              
                              # Specificity ---
                              sp <- paste0(round((d_v / (b_v + d_v))*100, 1), "%") # TN / (TN + FP)
                              
                              sp_low <- round((binom.confint(d_v, (b_v + d_v), 
                                                             conf.level = 0.95, methods = "exact")$lower)*100, 1)
                              
                              sp_up <- round((binom.confint(d_v, (b_v + d_v), 
                                                            conf.level = 0.95, methods = "exact")$upper)*100, 1) 
                              
                              # PPV ---              
                              ppv <- paste0(round((a_v / (a_v + b_v))*100, 1),"%") # TP / (TP + FP)
                              
                              ppv_low <- round((binom.confint(a_v, (a_v + b_v), 
                                                              conf.level = 0.95, methods = "exact")$lower)*100, 1)
                              
                              ppv_up <- round((binom.confint(a_v, (a_v + b_v), 
                                                             conf.level = 0.95, methods = "exact")$upper)*100, 1)
                              
                              # NPV ---
                              npv <- paste0(round((d_v / (c_v + d_v))*100, 1), "%") # TN / (TN + FN)
                              
                              npv_low <- round((binom.confint(d_v, (c_v + d_v), 
                                                              conf.level = 0.95, methods = "exact")$lower)*100, 1)
                              
                              npv_up <- round((binom.confint(d_v, (c_v + d_v), 
                                                             conf.level = 0.95, methods = "exact")$upper)*100, 1)
                              
                              # Return all output --
                              return(list(se, se_low, se_up,
                                          sp, sp_low, sp_up,
                                          ppv, ppv_low, ppv_up, 
                                          npv, npv_low, npv_up))
}




####---- 2. Calculate Validation Parameters using function from #1 ----####

t.out1 <- as.data.frame(matrix(nrow = 9, ncol = 14)) # Empty table, 9 rows, 14 columns

names(t.out1) <- c("hdd_var", "mr_var", # Name table columns
                   "Sensitivity", "Se_Lower", "Se_Upper",
                   "Specificity", "Sp_Lower", "Sp_Upper",
                   "PPV", "PPV_Lower", "PPV_Upper",
                   "NPV", "NPV_Lower", "NPV_Upper")

t.out1$hdd_var <- c("DX_AllHTN", "DX_AllHTN_noUnsp", "DX_HELLP", # Vector with HDD var names
                    "DX_Eclampsia", "DX_SeverePRE", "DX_MildPRE", 
                    "DX_SIPRE", "DX_ChronicHTN", "DX_GestHTN")

t.out1$mr_var <- c("IDPREG_AllHTN_NoUnspec", "IDPREG_AllHTN_NoUnspec", # Vector with MR var names - same order as HDD
                   "IDPREG_HELLP", "IDPREG_Eclampsia", "IDPREG_SeverePRE", 
                   "IDPREG_MildPRE", "IDPREG_SIPRE","IDPREG_ChronicHTN", "IDPREG_GestHTN")

for (i in (1:9)) { # Fill the table 
      t.out1[i, c(3:14)] <- validation.param.function(t.out1$hdd_var[i], t.out1$mr_var[i])
}

#View(t.out1)




####---- 3. Function to Calculate Prevalence and Bias ----####

bias.function <- function(hdd_var, mr_var) {
                      a_v <- sum(hdp_analysis[,hdd_var] == 1 & hdp_analysis[,mr_var] == 1)
                      b_v <- sum(hdp_analysis[,hdd_var] == 1 & hdp_analysis[,mr_var] == 0)
                      c_v <- sum(hdp_analysis[,hdd_var] == 0 & hdp_analysis[,mr_var] == 1)
                      d_v <- sum(hdp_analysis[,hdd_var] == 0 & hdp_analysis[,mr_var] == 0)
                      
                      total_HDD <- a_v + b_v
                      total_MR <- a_v + c_v
                      n <- a_v + b_v + c_v + d_v
                      
                      # Prevalence HDD--- 
                      prev_HDD <- (total_HDD /(n))*100 #total HDD / all
                      
                      prev_HDD_round <- paste0(round(prev_HDD, 1),"%") # Rounding to one decimal places
                      
                      
                      prev_HDD_low <- round((binom.confint(total_HDD, n, 
                                                      conf.level = 0.95, methods = "exact")$lower)*100, 1)
                      
                      prev_HDD_up <- round((binom.confint(total_HDD, n, 
                                                     conf.level = 0.95, methods = "exact")$upper)*100, 1)
                      
                      # Prevalence MR---
                      prev_MR <- (total_MR /(n))*100 #total MR / all
                      
                      prev_MR_round <- paste0(round(prev_MR, 1),"%") # Rounding to one decimal places
                      
                      prev_MR_low <- round((binom.confint(total_MR, n, 
                                                           conf.level = 0.95, methods = "exact")$lower)*100, 1)
                      
                      prev_MR_up <- round((binom.confint(total_MR, n, 
                                                          conf.level = 0.95, methods = "exact")$upper)*100, 1)
                      
                      # Expected Bias---
                      
                      PR_Bias <- round(prev_MR/prev_HDD, 1)
                      
                      PD_Bias <- paste0(round(prev_MR - prev_HDD, 1), "%")
                      
                      case_diff <- total_MR - total_HDD
                      
                      # Return output---
                      return(list(total_HDD, total_MR,
                                  prev_HDD_round, prev_HDD_low, prev_HDD_up,
                                  prev_MR_round, prev_MR_low, prev_MR_up,
                                  PR_Bias, PD_Bias, case_diff))
}



####---- 4. Calculate prevalence and expected bias using function from #3 ----####

t.out2 <- as.data.frame(matrix(nrow = 9, ncol = 13)) # Empty table, 9 rows, 13 columns

names(t.out2) <- c("hdd_var", "mr_var", # Name table columns
                   "total_HDD", "total_MR",
                   "prev_HDD", "prev_HDD_low", "prev_HDD_up",
                   "prev_MR", "prev_MR_low", "prev_MR_up",
                   "PR_Bias", "PD_Bias","case_diff")

t.out2$hdd_var <- c("DX_AllHTN", "DX_AllHTN_noUnsp", "DX_HELLP", # Vector with HDD var names
                    "DX_Eclampsia", "DX_SeverePRE", "DX_MildPRE", 
                    "DX_SIPRE", "DX_ChronicHTN", "DX_GestHTN")

t.out2$mr_var <- c("IDPREG_AllHTN_NoUnspec", "IDPREG_AllHTN_NoUnspec", # Vector with MR var names - same order as HDD
                   "IDPREG_HELLP", "IDPREG_Eclampsia", "IDPREG_SeverePRE", 
                   "IDPREG_MildPRE", "IDPREG_SIPRE","IDPREG_ChronicHTN", 
                   "IDPREG_GestHTN")

for (i in (1:13)) { # Fill the table
  t.out2[i, c(3:13)] <- bias.function(t.out2$hdd_var[i], t.out2$mr_var[i])
}

#View(t.out2)
