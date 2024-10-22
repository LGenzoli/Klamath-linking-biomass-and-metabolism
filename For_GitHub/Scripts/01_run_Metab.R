########################################################
#
# Here, I calculate stream metabolism for 11 sites for ~3.5 mo on the Klamath River.
# Script created on 8Feb2021 by L. Genzoli
# This is final metabolism run compiled and QA-ed DO data from all 11 sites
#
# SV and OR data was updated after the compilation, so new files are brought in 
#   to run those sites, but no tribal DO data is posted publicly (SV, OR, WE, KAT),
#   and must be accessed by contacting the Karuk and Yurok Tribes.
#
# For each site, I saved 4 files, but in the end, I ended up using the RDA file in
#   in the next script(02_Extract_Metab_Data.R) to extract parameters and CIs. This 
#   script will run for the non-tribal data sites, but would need to have additional 
#   folders in the project to save the output. So, my metabolism output data really 
#   starts in the next script. Below are files originally saved in this script: 
#   1) RDA file
#   2) csv of the metab params
#   3) csv of the modeled DO (to compare to measured)
#   4) stock plot of the measured vs. modeled data
#
# Script relies on: oxy_20210212_with_EvenTS_03.csv
# Script will take a very long time to run!
#
#
#
#
########################################################



## Packages:
library(tidyverse)
library(lubridate)
library(rstan)
library(streamMetabolizer)

setwd("/Users/laurelgenzoli/Dropbox/2018-2022_PHD/2019_Klamath/For_GitHub/Data")
oxy_data <- read.csv("oxy_20210212_with_EvenTS_03.csv")%>%
  select(-1)%>%
  filter(!is.na(site))

unique(oxy_data$site)
str(oxy_data)
head(oxy_data)

##-----Make Metabolism Data Frame and check with a plot:
########################################################

oxy <- oxy_data %>% 
  mutate(solar.time = as.POSIXct(solar.time, tz = "UTC"))%>%
  filter(QA != "fouled")%>%
  select(-QA, -bp)

ggplot(oxy, aes(x = solar.time, y = DO.obs))+
  geom_point(size = 0.2, alpha = 0.1) +
  facet_wrap(~site, ncol = 2)+
  theme_classic()

## quick check that DOsat looks good in the fixed data
oxy%>% filter(site == "WE")%>%
  slice(4500:5000)%>%
  ggplot(aes(x = solar.time, y = DO.sat), alpha = 0.5)+
  geom_point(size = 0.5)+
  geom_line(aes(x = solar.time, y = DO.obs), alpha = 0.3)


#######################################
## LOAD OR DATA that got updated 

OR_data <- read_csv("OR_2019_UPDATE.csv")




## names and specs will be the same for each site as follows:
## pool k "normal" (bobs specs from email)
##############################################
## for all sites except KAT:

nn_name_poolNormal <- mm_name(type = "bayes", 
                              pool_K600 = "normal", 
                              err_obs_iid = T, err_proc_iid =T)
nn_specs_poolNormal <- specs(nn_name_poolNormal, 
                             K600_daily_meanlog_meanlog=2.25, 
                             K600_daily_meanlog_sdlog=0.75, 
                             K600_daily_sdlog_sigma=0.02, 
                             burnin_steps=1000, saved_steps=1000) 


## for KAT use 0.56, 0.75 (mean and sd)
kat_nn_name_poolNormal <- mm_name(type = "bayes", 
                                  pool_K600 = "normal", 
                                  err_obs_iid = T, err_proc_iid =T)
kat_nn_specs_poolNormal <- specs(kat_nn_name_poolNormal, 
                                 K600_daily_meanlog_meanlog=0.56, 
                                 K600_daily_meanlog_sdlog=0.75, 
                                 K600_daily_sdlog_sigma=0.02, 
                                 burnin_steps=1500, saved_steps=1500) 




## Run I5 with  pool_K600 = "normal"
###############################################
I5 <- oxy %>% filter(site == "I5") %>% select(-site, -discharge)

mm_I5 <- metab(nn_specs_poolNormal, data=I5)
plot_metab_preds(mm_I5)
I5_params <- get_params(mm_I5)
I5_specs <- get_specs(mm_I5)
I5_mod_DO <- predict_DO(mm_I5)

ggplot(I5_params, aes(y = K600.daily, x = ER.daily))+geom_point()

saveRDS(mm_I5, file = "full_model_output/I5_pool_K600_normal_evenTS.RDS") 
write.csv(I5_params, "params/I5_pool_K600_normal_evenTS.csv")
write.csv(I5_mod_DO, "mod_DO/I5_pool_K600_normal_evenTS.csv")

pdf("mod_DO_plot/I5_pool_K600_normal_evenTS.pdf", width = 160, height = 14)
plot_DO_preds(mm_I5)
dev.off()


## Run TH with  pool_K600 = "normal"
###############################################

## TH
###############################################
TH <- oxy %>% filter(site == "TH") %>% select(-site, -discharge)

mm_TH <- metab(nn_specs_poolNormal, data=TH)
plot_metab_preds(mm_TH)
TH_params <- get_params(mm_TH)
TH_specs <- get_specs(mm_TH)
TH_mod_DO <- predict_DO(mm_TH)

ggplot(TH_params, aes(y = K600.daily, x = ER.daily))+geom_point()

saveRDS(mm_TH, file = "full_model_output/TH_pool_K600_normal_evenTS.RDS") 

write.csv(TH_params, "params/TH_pool_K600_normal_evenTS.csv")
write.csv(TH_mod_DO, "mod_DO/TH_pool_K600_normal_evenTS.csv")

pdf("mod_DO_plot/TH_pool_K600_normal_evenTS.pdf", width = 160, height = 14)
plot_DO_preds(mm_TH)
dev.off()
###############################################



## Run ABC with  pool_K600 = "normal"
###############################################

ABC <- oxy %>% filter(site == "ABC") %>% select(-site, -discharge)

mm_ABC <- metab(nn_specs_poolNormal, data=ABC)
plot_metab_preds(mm_ABC)
ABC_params <- get_params(mm_ABC)
ABC_specs <- get_specs(mm_ABC)
ABC_mod_DO <- predict_DO(mm_ABC)

ggplot(ABC_params, aes(y = K600.daily, x = ER.daily))+geom_point()

saveRDS(mm_ABC, file = "full_model_output/ABC_pool_K600_normal_evenTS.RDS") 

write.csv(ABC_params, "params/ABC_pool_K600_normal_evenTS.csv")
write.csv(ABC_mod_DO, "mod_DO/ABC_pool_K600_normal_evenTS.csv")

pdf("mod_DO_plot/ABC_pool_K600_normal_evenTS.pdf", width = 160, height = 14)
plot_DO_preds(mm_ABC)
dev.off()
###############################################


## Run BB with  pool_K600 = "normal"
###############################################

BB <- oxy %>% filter(site == "BB") %>% select(-site, -discharge)

mm_BB <- metab(nn_specs_poolNormal, data=BB)
plot_metab_preds(mm_BB)
BB_params <- get_params(mm_BB)
BB_specs <- get_specs(mm_BB)
BB_mod_DO <- predict_DO(mm_BB)

ggplot(BB_params, aes(y = K600.daily, x = ER.daily))+geom_point()

saveRDS(mm_BB, file = "full_model_output/BB_pool_K600_normal_evenTS.RDS") 

write.csv(BB_params, "params/BB_pool_K600_normal_evenTS.csv")
write.csv(BB_mod_DO, "mod_DO/BB_pool_K600_normal_evenTS.csv")

pdf("mod_DO_plot/BB_pool_K600_normal_evenTS.pdf", width = 160, height = 14)
plot_DO_preds(mm_BB)
dev.off()
###############################################



## Run RP with  pool_K600 = "normal"
###############################################

RP <- oxy %>% filter(site == "RP") %>% select(-site, -discharge)

mm_RP <- metab(nn_specs_poolNormal, data=RP)
plot_metab_preds(mm_RP)
RP_params <- get_params(mm_RP)
RP_specs <- get_specs(mm_RP)
RP_mod_DO <- predict_DO(mm_RP)

ggplot(RP_params, aes(y = K600.daily, x = ER.daily))+geom_point()

saveRDS(mm_RP, file = "full_model_output/RP_pool_K600_normal_evenTS.RDS") 

write.csv(RP_params, "params/RP_pool_K600_normal_evenTS.csv")
write.csv(RP_mod_DO, "mod_DO/RP_pool_K600_normal_evenTS.csv")

pdf("mod_DO_plot/RP_pool_K600_normal_evenTS.pdf", width = 160, height = 14)
plot_DO_preds(mm_RP)
dev.off()
###############################################



## Run SV with  pool_K600 = "normal"
###############################################

## Here, I'm grabbing data from the big metabolism analysis
## There was missing data in the OG DO data

setwd('/Users/laurelgenzoli/Dropbox/2018-2022_PHD')
sv <- read_csv("2023_Klamath_Metab/Klamath-Metabolism/DATA/O2_data_for_metab/SV_2005_2022_30min_QA.csv")%>%
  mutate(year = year(date))%>%
  mutate(jday = yday(date))%>%
  filter(year == 2019)%>%
  filter(jday %in% c(167:275))%>%
  select(-date, -discharge.daily, - year, -jday)
  
  sv%>%
  ggplot(aes(x = solar.time, y = DO.obs), alpha = 0.5)+
  geom_point(size = 0.5)

 SV <- sv 
 setwd("/Users/laurelgenzoli/Dropbox/2018-2022_PHD/2019_Klamath/Analysis/data/Metab_Output_2021/pool_Norm_EvenTS_03")
 

#SV <- oxy %>% filter(site == "SV") %>% select(-site, -discharge)

mm_SV <- metab(nn_specs_poolNormal, data=SV)
plot_metab_preds(mm_SV)
SV_params <- get_params(mm_SV)
SV_specs <- get_specs(mm_SV)
SV_mod_DO <- predict_DO(mm_SV)

ggplot(SV_params, aes(y = GPP.daily, x = date))+
  geom_point()

saveRDS(mm_SV, file = "full_model_output/SV_pool_K600_normal_evenTS_UPDATE.RDS") 

write.csv(SV_params, "params/SV_pool_K600_normal_evenTS_UPDATE.csv")
write.csv(SV_mod_DO, "mod_DO/SV_pool_K600_normal_evenTS_UPDATE.csv")

pdf("mod_DO_plot/SV_pool_K600_normal_evenTS_UPDATE.pdf", width = 160, height = 14)
plot_DO_preds(mm_SV)
dev.off()
###############################################


## Run HC with  pool_K600 = "normal"
###############################################

HC <- oxy %>% filter(site == "HC") %>% select(-site, -discharge)

mm_HC <- metab(nn_specs_poolNormal, data=HC)
plot_metab_preds(mm_HC)
HC_params <- get_params(mm_HC)
HC_specs <- get_specs(mm_HC)
HC_mod_DO <- predict_DO(mm_HC)

ggplot(HC_params, aes(y = K600.daily, x = ER.daily))+geom_point()

saveRDS(mm_HC, file = "full_model_output/HC_pool_K600_normal_evenTS.RDS") 

write.csv(HC_params, "params/HC_pool_K600_normal_evenTS.csv")
write.csv(HC_mod_DO, "mod_DO/HC_pool_K600_normal_evenTS.csv")

pdf("mod_DO_plot/HC_pool_K600_normal_evenTS.pdf", width = 160, height = 14)
plot_DO_preds(mm_HC)
dev.off()
###############################################


## Run OMR with  pool_K600 = "normal"
###############################################

OMR <- oxy %>% filter(site == "OMR") %>% select(-site, -discharge)

mm_OMR <- metab(nn_specs_poolNormal, data=OMR)
plot_metab_preds(mm_OMR)
OMR_params <- get_params(mm_OMR)
OMR_specs <- get_specs(mm_OMR)
OMR_mod_DO <- predict_DO(mm_OMR)

ggplot(OMR_params, aes(y = K600.daily, x = ER.daily))+geom_point()

saveRDS(mm_OMR, file = "full_model_output/OMR_pool_K600_normal_evenTS.RDS") 

write.csv(OMR_params, "params/OMR_pool_K600_normal_evenTS.csv")
write.csv(OMR_mod_DO, "mod_DO/OMR_pool_K600_normal_evenTS.csv")

pdf("mod_DO_plot/OMR_pool_K600_normal_evenTS.pdf", width = 160, height = 14)
plot_DO_preds(mm_OMR)
dev.off()
###############################################



## Run OR with  pool_K600 = "normal"
## I updated the data being used (OR_data) and added _UPDATE to the file name
###############################################

OR <- OR_data 

mm_OR <- metab(nn_specs_poolNormal, data=OR)
plot_metab_preds(mm_OR)
OR_params <- get_params(mm_OR)
OR_specs <- get_specs(mm_OR)
OR_mod_DO <- predict_DO(mm_OR)

ggplot(OR_params, aes(y = K600.daily, x = ER.daily))+geom_point()

saveRDS(mm_OR, file = "full_model_output/OR_pool_K600_normal_evenTS_UPDATE.RDS") 

write.csv(OR_params, "params/OR_pool_K600_normal_evenTS_UPDATE.csv")
write.csv(OR_mod_DO, "mod_DO/OR_pool_K600_normal_evenTS_UPDATE.csv")

pdf("mod_DO_plot/OR_pool_K600_normal_evenTS_UPDATE.pdf", width = 160, height = 14)
plot_DO_preds(mm_OR)
dev.off()
###############################################



## Run WE with  pool_K600 = "normal"
###############################################

WE <- oxy %>% filter(site == "WE") %>% select(-site, -discharge)

mm_WE <- metab(nn_specs_poolNormal, data=WE)
plot_metab_preds(mm_WE)
WE_params <- get_params(mm_WE)
WE_specs <- get_specs(mm_WE)
WE_mod_DO <- predict_DO(mm_WE)

ggplot(WE_params, aes(y = K600.daily, x = ER.daily))+geom_point()

saveRDS(mm_WE, file = "full_model_output/WE_pool_K600_normal_evenTS.RDS") 

write.csv(WE_params, "params/WE_pool_K600_normal_evenTS.csv")
write.csv(WE_mod_DO, "mod_DO/WE_pool_K600_normal_evenTS.csv")

pdf("mod_DO_plot/WE_pool_K600_normal_evenTS.pdf", width = 160, height = 14)
plot_DO_preds(mm_WE)
dev.off()
###############################################



## Run KAT with  pool_K600 = "normal"
###############################################

KAT <- oxy %>% filter(site == "KAT") %>% select(-site, -discharge)

mm_KAT <- metab(kat_nn_specs_poolNormal, data=KAT)
plot_metab_preds(mm_KAT)
KAT_params <- get_params(mm_KAT)
KAT_specs <- get_specs(mm_KAT)
KAT_mod_DO <- predict_DO(mm_KAT)

ggplot(KAT_params, aes(y = K600.daily, x = ER.daily))+geom_point()

saveRDS(mm_KAT, file = "full_model_output/KAT_pool_K600_normal_evenTS.RDS") 
write.csv(KAT_params, "params/KAT_pool_K600_normal_evenTS.csv")
write.csv(KAT_mod_DO, "mod_DO/KAT_pool_K600_normal_evenTS.csv")

pdf("mod_DO_plot/KAT_pool_K600_normal_evenTS.pdf", width = 160, height = 14)
plot_DO_preds(mm_KAT)
dev.off()

