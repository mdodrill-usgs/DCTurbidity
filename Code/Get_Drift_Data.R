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
# need to cut out the worms here...
# counts...
drift.dat.count = drift.dat$Specimens

drift.dat.count.2 = drift.dat.count[which(drift.dat.count$SpeciesID != "OLIG"),]

drift.dat.count.2$CountTotal = rowSums(drift.dat.count.2[,3:23])


# need to cut out the worms here...
# mass...
drift.dat.mass = drift.dat$Biomass

drift.dat.mass.2 = drift.dat.mass[which(drift.dat.mass$SpeciesID != "OLIG"),]

drift.dat.mass.2$MassTotal = rowSums(drift.dat.mass.2[,3:23])


ltl.specs = drift.dat.mass.2[,c(1,2,24)]

#--------------------------------------
# merge the sample table with the specimen counts

no.drift.dat = left_join(ltl.samp.2, ltl.specs, by = "BarcodeID")

no.drift.dat$SpeciesID = as.character(no.drift.dat$SpeciesID)

no.drift.dat$SpeciesID.2 = ifelse(no.drift.dat$SpeciesID == "CHIP", "CHIA", no.drift.dat$SpeciesID)
no.drift.dat$SpeciesID.2 = ifelse(no.drift.dat$SpeciesID == "SIMP", "SIMA", no.drift.dat$SpeciesID.2)

#--------------------------------------
# sum up the life-stage mass estimates
t.dat = dplyr::select(no.drift.dat, ID = BarcodeID, no.trip = TripID, no.site, taxa = SpeciesID.2, mass = MassTotal) %>% 
  filter(no.site %in% c("I", "II", "III", "IVa", "IVb")) %>%
  group_by(ID, no.trip, no.site, taxa) %>%
  summarise(mass.2 = sum(mass))

# get the avg. mass per trip, site, & taxa
t.dat.2 = group_by(t.dat, no.trip, no.site, taxa) %>% 
          summarise(avg.mass = mean(mass.2))


#-----------------------------------------------------------------------------#
write.table(t.dat.2, file = "Drift_Data_2019_05_31_mass_by_taxa.txt", row.names = F, sep = "\t")
