###############################################################################
#                                                                      March 19
#          Discrete Choice Stan Model - RBT Prey Selection
#
#  Notes:
#  * Stole some of this from the DiscreteChoice functions:
#    - 'Functions.R'
#  * Stole some from: 'Model_Turb_V1_all_no_worms.R'
#  
#  To Do:
#  * 
#
###############################################################################

#-----------------------------------------------------------------------------#
model.set.up.no.worms.turb = function(model.name = NULL){
  #-----------------------------------------------------------------------------#
  # read in the drift data
  
  if(is.null(model.name)){
    message("Need a model name...")
  } else {
    tmp.A = read.table(paste0(getwd(), "/Data/Drift_A.txt"), sep = "\t", header = TRUE)
    colnames(tmp.A)[3:80] = as.character(1:78) 
    
    # cut out worms
    tmp.A = tmp.A[,-c((ncol(tmp.A) - 18):ncol(tmp.A))]
    
    # only by taxa, not size
    tmp.A.w = read.table(paste0(getwd(), "/Data/Drift_A_w.txt"), sep = "\t", header = TRUE)
    
    # cut out worms
    tmp.A.w = tmp.A.w[,-9]
    
    # these should match
    if(identical(tmp.A[,1:2], tmp.A.w[,1:2]) == FALSE){
      message("Stop - Problem with the data")
    }
    #-----------------------------------------------------------------------------#
    # read in the diet data
    
    y.tmp2 = read.table(paste0(getwd(), "/Data/Diet.txt"), sep = "\t", header = TRUE)
    colnames(y.tmp2)[6:83] = as.character(1:78) 
    
    # cut out the worms
    y.tmp2 = y.tmp2[,-c(60:78)]
    
    # only by taxa, not size
    w.tmp2 = read.table(paste0(getwd(), "/Data/Diet_w.txt"), sep = "\t", header = TRUE)
    
    # cut out the worms
    w.tmp2 = w.tmp2[,-10]
    
    # these should match
    if(identical(y.tmp2[,1:3], w.tmp2[,1:3]) == FALSE){
      message("Stop - Problem with the data")
    } 
    
    #-----------------------------------------------------------------------------#
    # format/match/check data
    
    drift.ts = paste(tmp.A[,1], tmp.A[,2])
    diet.ts = sort(unique(paste(y.tmp2[,1], y.tmp2[,2])))
    
    # subset only the drift that matches the diet
    A = as.matrix(tmp.A[match(diet.ts, drift.ts),3:ncol(tmp.A)])
    w.a = as.matrix(tmp.A.w[match(diet.ts, drift.ts),3:ncol(tmp.A.w)])
    Nst = as.numeric(nrow(A))  # number of site & trips
    
    # diet
    y.tmp = y.tmp2[order(paste(y.tmp2[,1], y.tmp2[,2])),]
    w.tmp = w.tmp2[order(paste(w.tmp2[,1], w.tmp2[,2])),]
    
    y.tmp <<- y.tmp  # writes this to the global
    
    if(identical(y.tmp[,1:3], w.tmp[,1:3]) == FALSE){
      message("Stop - Problem with the data")
    }
    y.in = as.matrix(y.tmp[,6:ncol(y.tmp)])
    w.in = as.matrix(w.tmp[,4:ncol(w.tmp)])
    
    #-----------------------------------------------------------------------------#
    # format some variables for the model
    Nind = nrow(y.in)   # Number of individuals
    Nsp = 6             # number of prey taxa 
    upper = c(7, 15, 10, 13, 8, 12) - 1
    Nspsz = sum(upper)  # number of taxa and sizes
    spsz = c(1:Nspsz)
    idx_ts_ind = as.numeric(as.factor(paste(y.tmp$trip, y.tmp$site)))
    
    spp = sort(as.numeric(rep(c(1:Nsp), upper)))
    
    # new varible for the fixed paramaters in the drift portion of the model 
    spp2 = sort(as.numeric(rep(c(1:Nsp), upper-1)))
    
    # better way to do this? 
    idx = c(rep(1, upper[1]),
            rep(upper[1] + 1, upper[2]),
            rep(upper[1] + upper[2] + 1, upper[3]),
            rep(upper[1] + upper[2] + upper[3] + 1, upper[4]),
            rep(upper[1] + upper[2] + upper[3] + upper[4] + 1, upper[5]),
            rep(upper[1] + upper[2] + upper[3] + upper[4] + upper[5] + 1, upper[6]))
    
    idx2 = c(rep(upper[1], upper[1]),
             rep(upper[1] + upper[2], upper[2]),
             rep(upper[1] + upper[2] + upper[3], upper[3]),
             rep(upper[1] + upper[2] + upper[3] + upper[4], upper[4]),
             rep(upper[1] + upper[2] + upper[3] + upper[4] + upper[5], upper[5]),
             rep(upper[1] + upper[2] + upper[3] + upper[4] + upper[5] + upper[6], upper[6]))
    
    out = NULL
    for(i in 1:Nsp){
      out[[i]] = seq(2, upper[i] + 1)  
    }
    
    sz = do.call(c,out)
    
    len.idx = sz
    len.idx <<- len.idx
    
    u.idx = unique(idx)
    u.idx2 = unique(idx2)
    
    # need these (not for the current model)
    # Ns = 5
    # Nt = length(unique(y.tmp[,1]))
    # trip = as.numeric(y.tmp[,1])
    # site = as.numeric(y.tmp[,2])
    
    #-----------------------------------------------------------------------------#
    # fish size covariate (one sample has no size, set to the mean)
    t.sz = y.tmp[,4] 
    
    t.sz[which(t.sz == 2287)] = NA
    
    library(arm)
    fish_sz = scale(t.sz, center = T, scale = T)
    
    # fish_sz = t.sz - mean(t.sz, na.rm = T) 
    
    fish_sz[is.na(fish_sz)] = 0
    
    fish_sz = fish_sz[,1]
    
    #-----------------------------------------------------------------------------#
    # need emp_a (drift prop for each taxa across all the data)
    tmp.a = data.frame(spp = spp,
                       sums = colSums(A))
    
    tmp.a.2 = group_by(tmp.a, spp) %>%
      summarize(my.sum = sum(sums))
    
    a.tmp.w = colSums(w.a)
    
    tmp.a.all = tmp.a.2$my.sum + a.tmp.w
    
    emp_a_2 = tmp.a.all / sum(tmp.a.all)
    
    emp_a_2 <<- emp_a_2
    
    # emp.dat2 = data.frame(taxa = name.key$name,
    #                       emp = emp_a_2)
    #-----------------------------------------------------------------------------#
    if(model.name == "Length"){
      # calc the average size of prey in the drift (probably a better way :/)
      temp = data.frame(sz = sz,
                        sums = colSums(A))
      
      temp.2 = group_by(temp, sz) %>%
        summarize(my.sum = sum(sums))
      
      tmp.vec = list()
      
      for(i in 1:dim(temp.2)[1]){
        tmp.vec[[i]] = rep(as.character(temp.2[i,1]), times = as.numeric(temp.2[i,2]))
      }
      
      all.len = do.call('c', tmp.vec)
      
      avg.log.measure = mean(log(as.numeric(all.len)))
    }
    
    #-----------------------------------------------------------------------------#
    if(model.name == "Width"){
      # width calc (see C:\Users\mdodrill\Desktop\FB_DOWN\Analysis\IMAGE\Image_V2.R)
      # wrote using 'dput' on 'out.w'
      
      fit.dat.w = structure(list(V1 = c("SIMA", "SIML", "CHIA", "CHIL", "GAM", "NZMS", "LUMB"),
                                 V2 = c(-0.712264417633861, -1.585825466801, -1.33393917278193,
                                        -1.80644981811059, -1.34684364955096, -0.295195222715158, 
                                        -1.49260182871226),
                                 V3 = c(0.589397273790521, 0.800389362070142, 0.614460012015383,
                                        0.563246021013028, 0.890778931760663, 0.424033122560842, 
                                        0.137416988547813)),
                            row.names = c(NA, -7L), class = "data.frame")
      # cut out worms
      fit.dat.w = fit.dat.w[-7,]
      
      name.key = data.frame(num = c(1:7),
                            abr.name = c("NZMS", "GAM", "SIMA", "SIML",
                                         "CHIA", "CHIL", "LUMB"),
                            name = c("Snails", "Gammarus",
                                     "Black Fly Adult", "Black Fly Larva",
                                     "Midge Adult", "Midge Larva","Worms"))
      name.key = name.key[-7,]
      
      tmp.ar2 = list()
      
      fit.c = vector(length = Nsp)
      
      del = vector(length = Nsp)
      
      area.df = list()
      
      for(i in 1:Nsp){
        sub = fit.dat.w[which(fit.dat.w$V1 == name.key[i,2]),]
        
        fit.c[i] = sub[,3]
        
        tmp.sz = sz[which(spp == i)]
        
        tmp.ar = sub[,2] + sub[,3] * log(tmp.sz)
        
        # del[i] = sub[,2] + sub[,3] * avg.log.measure  #log area for average size
        
        tmp.ar2[[i]] = exp(tmp.ar)
        
        area.df[[i]] = data.frame(Taxa = name.key[i,3],
                                  Size = tmp.sz,
                                  width = exp(tmp.ar))
      }
      
      dat.width = do.call('rbind', area.df)
      
      avg.log.measure = mean(log(dat.width$width))
      measure = dat.width$width
    }
    #-----------------------------------------------------------------------------#
    if(model.name == "Area"){
      # area calc (see C:\Users\mdodrill\Desktop\FB_DOWN\Analysis\IMAGE\Image_V2.R)
      # wrote using 'dput'
      
      fit.dat = structure(list(V1 = c("SIMA", "SIML", "CHIA", "CHIL", "GAM", "NZMS", "LUMB"),
                               V2 = c(0.224460916479209, -1.78662907514783, -1.34743919150295,
                                      -2.15070453700794, -1.45525542402561, -0.74863545618833,
                                      -1.85570798696119),
                               V3 = c(0.94958327948254, 1.79612553310148, 1.54243803018606,
                                      1.67798706868551, 1.89261925347815, 1.65868421232674,
                                      1.30725930407769)),
                          row.names = c(NA, -7L), class = "data.frame")
      
      fit.dat = fit.dat[-7,]
      
      name.key = data.frame(num = c(1:7),
                            abr.name = c("NZMS", "GAM", "SIMA", "SIML",
                                         "CHIA", "CHIL", "LUMB"),
                            name = c("Snails", "Gammarus",
                                     "Black Fly Adult", "Black Fly Larva",
                                     "Midge Adult", "Midge Larva","Worms"))
      name.key = name.key[-7,]
      
      tmp.ar2 = list()
      
      fit.c = vector(length = Nsp)
      
      del = vector(length = Nsp)
      
      area.df = list()
      
      for(i in 1:Nsp){
        sub = fit.dat[which(fit.dat$V1 == name.key[i,2]),]
        
        fit.c[i] = sub[,3]
        
        tmp.sz = sz[which(spp == i)]
        
        tmp.ar = sub[,2] + sub[,3] * log(tmp.sz)
        
        # del[i] = sub[,2] + sub[,3] * avg.log.measure  #log area for average size
        
        tmp.ar2[[i]] = exp(tmp.ar)
        
        area.df[[i]] = data.frame(Taxa = name.key[i,3],
                                  Size = tmp.sz,
                                  area = exp(tmp.ar))
      }
      
      # area = do.call('c', tmp.ar2)
      area = do.call('rbind', area.df)
      avg.log.measure = mean(log(area$area))
      measure = area$area
      
    }    
    #-----------------------------------------------------------------------------#
    if(model.name == "Mass"){
      # estimate mass (see foodbase package species list for regression parameters)
      # species list to estimate biomass
      sp.tmp = read.table(paste0(getwd(), "/Data/SpeciesList_2019_02_07.txt"), sep = "\t", header = TRUE)
      
      sp.tmp.2 = sp.tmp[which(sp.tmp$SpeciesID %in% c("CHIL", "CHIA", "SIML", "SIMA",
                                                      "GAMM", "NZMS")),]
      sp.tmp.3 = sp.tmp.2[,which(names(sp.tmp.2) %in% c("SpeciesID", "RegressionA", "RegressionB"))]
      
      sp.key = data.frame(num = c(2,4,5,7,19,21),
                          name = c("NZMS", "GAMM", "SIMA", "SIML", "CHIA", "CHIL"))
      sp.key$A = sp.tmp.3[match(sp.key$name, sp.tmp.3$SpeciesID),]$RegressionA
      sp.key$B = sp.tmp.3[match(sp.key$name, sp.tmp.3$SpeciesID),]$RegressionB
      
      out.list = list() 
      
      for(i in 1:Nsp){
        tmp.sz = sz[which(spp == i)]
        
        out.list[[i]] = sp.key[i,"A"] * tmp.sz^sp.key[i,"B"]
        
      }
      
      # mass
      measure = do.call(c,out.list)
      avg.log.measure = mean(log(measure))
      
      # check / look at regressions
      # test = as.data.frame(cbind(sz, as.factor(spp), mass))
      # windows()
      # p = ggplot(test, aes(x = sz, y = mass, group = spp)) +
      #     geom_point(aes(color = spp)) + 
      #     geom_line(aes(color = spp)) 
      # p
      
      
      #-----------------------------------------------------------------------------#  
    }
    
    # dirichlet params 
    alpha = rep(.2, Nsp)
    
    # dummy matrix
    X <- as.matrix(model.matrix(~ as.factor(spp) - 1))   
    
    idx_first = unique(idx)  # position of the first bin for each taxa
    tmpper = seq(1:Nspsz)
    not_first = tmpper[which(!tmpper %in% idx_first)]
    #-----------------------------------------------------------------------------#
    # add this in above 
    if(model.name == "Length"){
      measure = sz
    }
    #-----------------------------------------------------------------------------#
    # variables specific to the turbidity model
    d.turb = read.table(paste0(getwd(), "/Data/NO_Turb_ts_summary.txt"), sep = "\t", header = TRUE)
    
    # make sure the order is the same 
    turb.ts = paste(d.turb[,1], d.turb[,2])
    
    # cut out trips to match the diet
    d.turb2 = d.turb[match(diet.ts, turb.ts),]
    
    tmp.turb = scale(log(d.turb2$ts.mean), center = T, scale = T)
      
    tmp.turb.2 = tmp.turb[,1] 
    
    # expand to match the ts order of diet (at the individual level, this is so
    # I can use the .* operator in stan)
    turb = tmp.turb.2[match(paste(y.tmp[,1], y.tmp[,2]), paste(d.turb2[,1], d.turb2[,2]))]
    
    
    # Rosenfeld and Taylor 2003
    tmp.turb.pmr = 100 - (44.8 * log10(d.turb2$ts.mean + 1))
    
    tmp.turb.pmr = ifelse(tmp.turb.pmr < 0, 0, tmp.turb.pmr)
    
    tmp.turb.pmr.2 = scale(tmp.turb.pmr, center = T, scale = T)
    
    turb.pmr = tmp.turb.pmr.2[,1]
    
    #-----------------------------------------------------------------------------#
    # # Drift Data  (the mass is in mg/m3)
    # d.drift = read.table(paste0(getwd(), "/Data/Drift_Data_2019_05_28_mass.txt"),
    #                      sep = "\t", header = TRUE)
    # drift.ts = paste(d.drift[,1], d.drift[,2])
    # d.drift.2 = d.drift[match(diet.ts, drift.ts),]
    # 
    # tmp.drift = scale(d.drift.2$drift.mass, center = T, scale = T)
    # 
    # tmp.drift.2 = tmp.drift[,1] # convert from a matrix to vector
    # 
    # # expand to match the ts order of diet (at the individual level, this is so
    # # I can use the .* operator in stan)
    # mass = tmp.drift.2[match(paste(y.tmp[,1], y.tmp[,2]), paste(d.drift.2[,1], d.drift.2[,2]))]
    # 
    #-----------------------------------------------------------------------------#
    d.drift = read.table(paste0(getwd(), "/Data/Drift_Data_2019_05_31_mass_by_taxa.txt"),
                         sep = "\t", header = TRUE)
    
    fkn.key = data.frame(num = c(1:6),
                         abrev = c("NZMS", "GAMM", "SIMA", "SIML", "CHIA", "CHIL"))
    
    dat.out = list()
  
    for(i in 1:6){
      sub = d.drift[d.drift$taxa == fkn.key[i,2],]
      sub.ts = paste(sub[,1], sub[,2], sep = " ") 
      
      sub.tmp = sub[match(diet.ts, sub.ts),4]
      sub.tmp[is.na(sub.tmp)] = 0
      
      sub.tmp.tmp = scale(sub.tmp, center = T, scale = T)
      
      dat.out[[i]] = sub.tmp.tmp[match(paste(y.tmp[,1], y.tmp[,2]), sub.ts)]
    }
    
    d.dat.7 = do.call(cbind, dat.out) 
    
    # there are some non matches that result in NA, I think these are from when
    # a taxa doesn't occur in a given site, trip (I hope ;)
    d.dat.7[is.na(d.dat.7)] = 0
    
    mass = d.dat.7
    
    
    
    #-----------------------------------------------------------------------------#
    # RBT density 
    d.den = read.table(paste0(getwd(), "/Data/RBT_Density.txt"),
                         sep = "\t", header = TRUE)
    
    d.den.ts = paste(d.den[,2], d.den[,1])
    
    d.den.2 = d.den[match(diet.ts, d.den.ts),]
    
    tmp.den = scale(d.den.2$density, center = T, scale = T)
    
    tmp.den.2 = tmp.den[,1] # convert from a matrix to vector
    
    trout = tmp.den.2[match(paste(y.tmp[,1], y.tmp[,2]), paste(d.den.2[,2], d.den.2[,1]))]
    #-----------------------------------------------------------------------------#
    # # add in a covariate - prey concentration, overall not by taxa
    # # version b has 0's and the both levels of depth integrated
    # d.dat = read.csv(paste0(getwd(), "/Data/Drift_data_2019_03_19b.csv"), stringsAsFactors = FALSE)
    # 
    # # for the turbidity/this version, cut out the worms
    # 
    # d.dat.tmp = d.dat[which(d.dat$taxa != "Worms"),]
    # 
    # d.keep = c("DriftSampleID", "new.name", "new.Rtot", "TripID", "site")
    # d.dat.2 = d.dat.tmp[,which(names(d.dat.tmp) %in% d.keep)]
    # 
    # # d.dat.3 = group_by(d.dat.2, TripID, site, new.name) %>% 
    # #           summarise(mean = mean(new.Rtot)) %>%
    # #           ungroup() %>%
    # #           as.data.frame()
    # 
    # # for each sample, sum to get total concentration acrosss taxa
    # d.tmp.2 = group_by(d.dat.2, DriftSampleID, TripID, site) %>% 
    #   summarise(sum.Rtot = sum(new.Rtot))
    # # get the mean conc. for each trip & site
    # d.dat.3 = group_by(d.tmp.2, TripID, site) %>% 
    #   summarise(mean = mean(sum.Rtot))
    # 
    # 
    # dat.ts = paste(d.dat.3$TripID, d.dat.3$site)
    # 
    # out = d.dat.3[match(diet.ts, dat.ts),]
    # 
    # tmp.conc = scale(out$mean, center = T, scale = T)
    # 
    # conc = tmp.conc[match(paste(y.tmp[,1], y.tmp[,2]), paste(out[,1], out[,2]))]
    #-----------------------------------------------------------------------------#
    # add in a covariate - prey concentration, for each taxa...
    # version b has 0's and the both levels of depth integrated
    d.dat = read.csv(paste0(getwd(), "/Data/Drift_data_2019_03_19b.csv"), stringsAsFactors = FALSE)
    
    # take out the worms for this version
    d.dat.tmp = d.dat[which(d.dat$taxa != "Worms"),]
    
    d.keep = c("new.name", "new.Rtot", "TripID", "site")
    d.dat.2 = d.dat.tmp[,which(names(d.dat.tmp) %in% d.keep)]
    
    d.dat.3 = group_by(d.dat.2, TripID, site, new.name) %>% 
      summarise(mean = mean(new.Rtot)) %>%
      ungroup() %>%
      as.data.frame()
    
    
    # dat.ts = paste(d.dat.3[,1], d.dat.3[,2], sep = " ") 
    # as.matrix(tmp.A[match(diet.ts, drift.ts),3:ncol(tmp.A)])
    
    # name.key.7 = c(2,4,5,7,19,21,1)
    name.key.7 = c(2,4,5,7,19,21) # no worms
    
    dat.out = list()
    for(i in 1:length(name.key.7)){
      sub = d.dat.3[d.dat.3$new.name == name.key.7[i],]
      sub.ts = paste(sub[,1], sub[,2], sep = " ") 
      
      # dat.out[[i]] = sub[match(diet.ts, sub.ts),4]
      sub.tmp = sub[match(diet.ts, sub.ts),4]
      sub.tmp[is.na(sub.tmp)] = 0
      
      # dat.out[[i]] = sub.tmp - mean(sub.tmp)  # center on the mean
      sub.tmp.tmp = scale(sub.tmp, center = T, scale = T)
      
      dat.out[[i]] = sub.tmp.tmp[match(paste(y.tmp[,1], y.tmp[,2]), sub.ts)]
    }
    
    d.dat.4 = do.call(cbind, dat.out) 
    
    # there are some non matches that result in NA, I think these are from when
    # a taxa doesn't occur in a given site, trip (I hope ;)
    d.dat.4[is.na(d.dat.4)] = 0
    
    conc = d.dat.4
    
    
    #-----------------------------------------------------------------------------#
    # Competition index
    NO_trips = c('GC20120419','GC20120705','GC20120913','GC20130110',
                 'GC20130404','GC20130625','GC20130912','GC20140109',
                 'GC20140403','GC20140626','GC20140911','GC20150108')
    
    c.dat = read.csv(paste0(getwd(), "/Data/Comp_Idx.csv"), stringsAsFactors = FALSE)
    
    c.dat$no.trip = rep(NO_trips, 5)
    
    c.dat.2 = c.dat[match(diet.ts, paste(c.dat[,4], c.dat[,1])),]
    
    comp = as.vector(scale(c.dat.2$Comp, center = T, scale = T))
    #-----------------------------------------------------------------------------#
    
    data.in = list(Nspsz = Nspsz, Nst = Nst, Nsp = Nsp, Nind = Nind,
                   sp = spp, sp_idx2 = spp2,# spsz = spsz,
                   a = A, w_a = w.a, y = y.in, w = w.in,
                   idx = idx, idx2 = idx2, idx_first = idx_first, not_first = not_first,
                   idx_ts_ind = idx_ts_ind,
                   alpha = alpha, a_Nsz = upper,   
                   X = X, 
                   #  emp_a = emp_a_2,
                   sz = measure,  u_idx = u.idx, u_idx2 = u.idx2,
                   avg_log_len = avg.log.measure,
                   turb = turb, mass = mass, trout = trout, conc = conc,
                   fish_sz = fish_sz,
                   turb_pmr = turb.pmr, comp = comp)  
    
    return(data.in)  
  }
}
#-----------------------------------------------------------------------------#






