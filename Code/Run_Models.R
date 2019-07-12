###############################################################################
#                                                                     March '19
#          Discrete Choice Stan Model - RBT Prey Selection
#                             Turbidity
#
#  Notes:
#  * Changed the order of the prey to; c(2,4,5,7,19,21,1)
#  * Dropped the first size bin for all prey 
#    - changes to formatter, and subsequent code
#  * See 'Model_Turb_V1_all_no_worms.r' - Previous version
#  * See 'Model_Turb_V1' - Previous version
#    - U:\STAN_Mod_Turb_2
#  * Latest version has no worms
#
###############################################################################
rm(list = ls(all = TRUE))
library(dplyr)
library(rstan)
library(here)
# here()

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# Sys.setenv(LOCAL_CPPFLAGS = '-march=native')  # can't use this with the larger model...Crash!

#-----------------------------------------------------------------------------#

# source("U:/Desktop/Fish_Git/DiscreteChoice/Code/Functions.R"
source(paste0(getwd(),"/Code/Functions_Turb.R"))

#-----------------------------------------------------------------------------#
data.in = model.set.up.no.worms.turb(model.name = "Length")

list2env(data.in, globalenv())

params <- c("mu_sp_all", "mu_sz",
            "beta_sp_all", "beta_sp_turb_all", "beta_sz_turb", "beta_sz",
            "ts_sp_eff", "ts_sz_eff",
            "beta_f_sz_p_sp", "beta_f_sz_p_sz", "beta_f_sz_p_sp_all",
            "beta_f_sz_sp", "beta_f_sz_sz", "beta_f_sz_sp_all")

## MCMC settings
ni <- 500   # 500iter ~ 7.7h, with new tuning ~5.4h
nt <- 1
nb <- 250
nc <- 3

# m <- stan_model(paste0(getwd(), "/Stan_Code/",'Model_Turb_V1.stan'))

setwd("U:/Desktop/Fish_Git/DCTurbidity/Stan_Code")

fit <- stan("Model_Turb_v1.stan", data = data.in,
# fit <- stan(m, data = data.in,
            pars = params,
            # init = inits,
            # control = list(adapt_delta = .9),
            control = list(max_treedepth = 12),
            # control = list(max_treedepth = 15, adapt_delta = .925),
            chains = nc, thin = nt, iter = ni, warmup = nb)

save.image("U:/Desktop/Fish_Git/DCTurbidity/working_runs/Model_Turb_v1_Length_more_parms_500iter.RData")


#-----------------------------------------------------------------------------#


