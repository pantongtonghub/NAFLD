#### Forward code --------------
library(tidyverse)
library(data.table)
library(TwoSampleMR)
Gut_microbiota_IDs <- paste0('ebi-a-GCST9001', 6908:7118)  
exp_data <- extract_instruments(outcomes = Gut_microbiota_IDs, p1 = 1e-05)
out_data <- extract_outcome_data(snps = exp_data$SNP,outcomes = 'ebi-a-GCST90091033')
dat <- harmonise_data(exposure_dat = exp_data,outcome_dat = out_data)
dat$F_statistic <- round(dat$beta.exposure ^ 2 / dat$se.exposure ^ 2,2)
res <- mr(dat)
res_hete <- TwoSampleMR::mr_heterogeneity(dat)
res_plei <- TwoSampleMR::mr_pleiotropy_test(dat)
res_leaveone <- TwoSampleMR::mr_leaveoneout(dat)
res <- generate_odds_ratios(res)
res$estimate <- paste0(
  format(round(res$or,2),nsmall = 2),'(',
  format(round(res$or_lci95,2),nsmall = 2),'-',
  format(round(res$or_uci95,2),nsmall = 2),')')
rio::export(res,'All_results for forward MR.xlsx')
res_ivw <- res %>%
  filter(., method == 'Inverse variance weighted') %>%
  filter(., pval < 0.05) %>%
  select(., exposure, nsnp, or, estimate, pval)
rio::export(res_ivw, 'IVW for forward MR.xlsx')


#### Reverse code --------------
library(tidyverse)
library(data.table)
library(TwoSampleMR)
Gut_microbiota_IDs <- paste0('ebi-a-GCST9001', 6908:7118)  
exp_data <- extract_instruments(outcomes = 'ebi-a-GCST90091033',p1 = 5e-06)
out_data <- extract_outcome_data(snps = exp_data$SNP,outcomes = Gut_microbiota_IDs)
dat <- harmonise_data(exposure_dat = exp_data,outcome_dat = out_data)
dat$r.exposure <- get_r_from_bsen(b = dat$beta.exposure,
                                  dat$se.exposure,
                                  dat$samplesize.exposure)
dat$r.outcome <- get_r_from_bsen(b = dat$beta.outcome,
                                 dat$se.outcome,
                                 dat$samplesize.outcome)
dat <- steiger_filtering(dat)
dat <- subset(dat, dat$steiger_dir == 'TRUE')
dat$F_statistic <- round(dat$beta.exposure ^ 2 / dat$se.exposure ^ 2,2)
res <- mr(dat)
res_hete <- TwoSampleMR::mr_heterogeneity(dat)
res_plei <- TwoSampleMR::mr_pleiotropy_test(dat)
res_leaveone <- TwoSampleMR::mr_leaveoneout(dat)
res <- generate_odds_ratios(res)
res$estimate <- paste0(
  format(round(res$or,2),nsmall = 2),'(',
  format(round(res$or_lci95,2),nsmall = 2),'-',
  format(round(res$or_uci95,2),nsmall = 2),')')
rio::export(res,'All_results for reverse MR.xlsx')
res_ivw <- res %>%
  filter(., method == 'Inverse variance weighted') %>%
  filter(., pval < 0.05) %>%
  select(., outcome, nsnp, or, estimate, pval)
rio::export(res_ivw, 'IVW for reverse MR.xlsx')

