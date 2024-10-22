##
## Extract .csv and .rda files from the .RDS streamMetabolizer output files
##  that are huge. RDA files are not included in data release due to size, so 
##  this script doesn't run unless you were to re-run the metabolism script,
##  saving the RDA files, then using this script to extract them. Instead, I just
##  present the extracted data in a .csv and .rda file.
##
##
##
##  These metabolism output files are the basis for the analysis (along with 
##    biomass).
##
#####################################################

library(tidyverse)
library(streamMetabolizer)

setwd('/Users/laurelgenzoli/Dropbox/2018-2022_PHD/For_GitHub/Data')


## Read in .RDS File from metabolism run:

## ABC:

ABC <- readRDS("full_model_output/ABC_pool_K600_normal_evenTS.RDS")
abc <- get_params(ABC, uncertainty = "ci")
abc.rda <- get_fit(ABC)[-2]
write_csv(abc, "final_metab_params/abc.csv")
save(abc.rda, file = "final_metab_rda/abc.rda")


## BB:

BB <- readRDS("Analysis/data/Metab_Output_2021/pool_Norm_EvenTS_03/full_model_output/BB_pool_K600_normal_evenTS.RDS")
bb <- get_params(BB, uncertainty = "ci")
bb.rda <- get_fit(BB)[-2]
write_csv(bb, "final_metab_params/bb.csv")
save(bb.rda, file = "final_metab_rda/bb.rda")


## HC:

HC <- readRDS("full_model_output/HC_pool_K600_normal_evenTS.RDS")
hc <- get_params(HC, uncertainty = "ci")
hc.rda <- get_fit(HC)[-2]
write_csv(hc, "final_metab_params/hc.csv")
save(hc.rda, file = "Data/final_metab_rda/hc.rda")


## I5:

I5 <- readRDS("full_model_output/I5_pool_K600_normal_evenTS.RDS")
i5 <- get_params(I5, uncertainty = "ci")
i5.rda <- get_fit(I5)[-2]
write_csv(i5, "final_metab_params/i5.csv")
save(i5.rda, file = "final_metab_rda/i5.rda")


## KAT:

KAT <- readRDS("full_model_output/KAT_pool_K600_normal_evenTS.RDS")
kat <- get_params(KAT, uncertainty = "ci")
kat.rda <- get_fit(KAT)[-2]
write_csv(kat, "final_metab_params/kat.csv")
save(kat.rda, file = "final_metab_rda/kat.rda")


## OMR:

OMR <- readRDS("full_model_output/OMR_pool_K600_normal_evenTS.RDS")
omr <- get_params(OMR, uncertainty = "ci")
omr.rda <- get_fit(OMR)[-2]
write_csv(omr, "final_metab_params/omr.csv")
save(omr.rda, file = "final_metab_rda/omr.rda")

## OR:

OR <- readRDS("full_model_output/OR_pool_K600_normal_evenTS.RDS")
or <- get_params(OR, uncertainty = "ci")
or.rda <- get_fit(OR)[-2]
write_csv(or, "final_metab_params/or.csv")
save(or.rda, file = "final_metab_rda/or.rda")


## RP:

RP <- readRDS("full_model_output/RP_pool_K600_normal_evenTS.RDS")
rp <- get_params(RP, uncertainty = "ci")
rp.rda <- get_fit(RP)[-2]
write_csv(rp, "final_metab_params/rp.csv")
save(rp.rda, file = "final_metab_rda/rp.rda")


## SV:

SV <- readRDS("full_model_output/SV_pool_K600_normal_evenTS.RDS")
sv <- get_params(SV, uncertainty = "ci")
sv.rda <- get_fit(SV)[-2]
write_csv(sv, "final_metab_params/sv.csv")
save(sv.rda, file = "final_metab_rda/sv.rda")




## TH:

TH <- readRDS("full_model_output/TH_pool_K600_normal_evenTS.RDS")
th <- get_params(TH, uncertainty = "ci")
th.rda <- get_fit(TH)[-2]
write_csv(th, "final_metab_params/th.csv")
save(th.rda, file = "final_metab_rda/th.rda")

## WE:

WE <- readRDS("full_model_output/WE_pool_K600_normal_evenTS.RDS")
we <- get_params(WE, uncertainty = "ci")
we.rda <- get_fit(WE)[-2]
write_csv(we, "final_metab_params/we.csv")
save(we.rda, file = "final_metab_rda/we.rda")
