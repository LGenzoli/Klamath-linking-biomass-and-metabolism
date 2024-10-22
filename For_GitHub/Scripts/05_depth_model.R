##model deoth
#use2 models, variable intercept and fixed slope
#variable slope and variable intercept


library(brms)
library(loo)
library(tidybayes)
library(tidyverse)
library(ggplot2)


setwd('/Users/laurelgenzoli/Dropbox/2018-2022_PHD/2019_Klamath/2023_PUB')
klam_hg<- read.csv("Data/LongProfiles/r_curves_out.csv")

head(klam_hg)

## Variable intercept, fixed slope
hg_fixed_out<-brm(logZ ~  logQ + (1|site), data=klam_hg)

## Don't need to run this one; bob did the model 
## Variable Slope, Variable Intercept
hg_variable_out<-brm(logZ~ logQ + (1 + logQ|site) , data=klam_hg)

## Model Comparisons:
loo_fixed<-loo(hg_fixed_out)
loo_variable<- loo(hg_variable_out)
loo_compare(loo_fixed, loo_variable)
###looks like "best" model is fixed intercept


##look at random efects
summary(hg_fixed_out)

z.site.intercept <- ranef(hg_fixed_out) %>%
  data.frame()%>%
  select(1)%>%
  rownames_to_column("site")
write_csv(z.site.intercept, "Data/LongProfiles/z_site_intercept_out.csv")

fsteps<- hg_fixed_out %>% 
  spread_draws(b_Intercept, r_site[site,])


head(fsteps)

fsteps$intercept<- fsteps$b_Intercept + fsteps$r_site ## addding random effect to group level effect

ggplot(fsteps, aes(x=as.factor(site), y=intercept  )) + 
  geom_violin()+
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  xlab("site")+
  ylab("log intercept")+
  theme_classic()



###same for variable slope here, but we ignore the intercepts because these will covary
ranef(hg_variable_out)
vsteps<- hg_variable_out %>% 
  spread_draws(b_logQ, r_site[site,])


ggplot(vsteps, aes(x=as.factor(site), y=b_logQ  )) + 
  geom_violin()+
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  xlab("site")+
  ylab("slope random effect")+
  theme_classic()


depth<- exp(-1.04+ ranef[site=="I5"] + logQ*0.36) /3.28

