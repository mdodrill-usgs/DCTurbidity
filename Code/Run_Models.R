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
source(paste0(getwd(),"/Code/Functions.R"))

#-----------------------------------------------------------------------------#
data.in = model.set.up.no.worms.turb(model.name = "Length")

list2env(data.in, globalenv())

params <- c("beta_sz", "mu_sz")# "ts_sz_eff",# "t_sz_eff", "s_sz_eff",
# "beta_sp_all", "mu_sp_all", "ts_sp_eff", #"t_sp_eff", "s_sp_eff",
# "p", "p_a",
# "gamma_st",
# "beta_sp_turb_all", "beta_sz_turb")


## MCMC settings
ni <- 100   
nt <- 1
nb <- 50
nc <- 3

# m <- stan_model(paste0(getwd(), "/Stan_Code/",'Model_Turb_V1.stan'))

setwd("U:/Desktop/Fish_Git/DCTurbidity/Stan_Code")

fit <- stan("Model_Turb_v1_test.stan", data = data.in,
            # fit <- stan(m, data = data.in, 
            pars = params,
            # init = inits,
            # control = list(adapt_delta = .9),
            # control = list(max_treedepth = 14),
            control = list(max_treedepth = 15, adapt_delta = .925),
            chains = nc, thin = nt, iter = ni, warmup = nb)











