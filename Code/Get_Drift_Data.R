###############################################################################
#                                                                        May 19
#     Discrete Choice Stan Model - Get Drift Data to add as covariate
#
#  Notes:
#  * Stole this from:
#    - U:\Desktop\FB_Git\Diet\NO_Fullness\Data_Prep.R
#  
#  To Do:
#  * 
#
###############################################################################
# library(devtools)
# install_github(repo = 'jmuehlbauer-usgs/R-packages', subdir = 'foodbase')

library(foodbase)
library(dplyr)

#-----------------------------------------------------------------------------#
NO_trips = c('GC20120419','GC20120705','GC20120913','GC20130110',
             'GC20130404','GC20130625','GC20130912','GC20140109',
             'GC20140403','GC20140626','GC20140911','GC20150108')

# drift.tmp = readDB(gear = "Drift", type = "Sample", updater = TRUE)
drift.tmp = readDB(gear = "Drift", type = "Sample", updater = FALSE)

drift.tmp.2 = drift.tmp[which(drift.tmp$TripID %in% NO_trips),]
# all gear type 4 for NO trips

drift.dat = sampspec(samp = drift.tmp.2, species = "Big9")

#--------------------------------------
# cut down the sample table
keep.col.2 = c("BarcodeID",
               "TripID",
               # "Region",
               # "Reach",
               "Date",
               "SampleNumber",
               "RiverMile",
               "DepthTotal",
               "DepthSample",
               "DepthIntegrated",
               # "GearID",
               "TimeDay", 
               "TimeBegin", 
               "TimeElapsed",
               "Volume")
# "Notes")

ltl.samp = drift.dat$Samples[,which(names(drift.dat$Samples) %in% keep.col.2)] 

ltl.samp$no.site = ifelse(ltl.samp$RiverMile <= 0, 'I',
                          ifelse(ltl.samp$RiverMile >= 17.22 & ltl.samp$RiverMile <= 20.58, 'II', 
                                 ifelse(ltl.samp$RiverMile >= 37.57 & ltl.samp$RiverMile <= 42.11, 'III',
                                        ifelse(ltl.samp$RiverMile > 59 & ltl.samp$RiverMile < 62.5, 'IVa',
                                               ifelse(ltl.samp$RiverMile > 63 & ltl.samp$RiverMile < 65.6, 'IVb',
                                                      ifelse(ltl.samp$RiverMile > 65.6, 'downstream', 'other'))))))

# Only data for the NO sites
ltl.samp.2 = ltl.samp[which(ltl.samp$no.site %in% c("I", "II", "III", "IVa", "IVb")),]

#--------------------------------------
# only the total mass

# counts...
# drift.dat$Specimens$CountTotal = rowSums(drift.dat$Specimens[,3:23])
# mass...
drift.dat$Biomass$MassTotal = rowSums(drift.dat$Specimens[,3:23])
ltl.specs = drift.dat$Biomass[,c(1,2,24)]

#--------------------------------------
# merge the sample table with the specimen counts

no.drift.dat = left_join(ltl.samp.2, ltl.specs, by = "BarcodeID")

no.drift.dat$SpeciesID = as.character(no.drift.dat$SpeciesID)

#--------------------------------------
# should it be, sum across taxa for a sample, then average? (this gives the same avg mass as below)

# get the avg. mass per trip, site, & taxa
t.dat = dplyr::select(no.drift.dat, no.trip = TripID, no.site,
                      taxa = SpeciesID, mass = MassTotal) %>%
  filter(no.site %in% c("I", "II", "III", "IVa", "IVb")) %>%
  group_by(no.trip, no.site, taxa) %>%
  summarise(avg.mass = mean(mass))

# sum across taxa
dat3 = group_by(t.dat, no.trip, no.site) %>%
  summarise(drift.mass = sum(avg.mass))


#-----------------------------------------------------------------------------#
# write.table(dat3, file = "Drift_Data_2019_05_28_mass.txt", row.names = F, sep = "\t")
