// This model is based on "Model_Turb_v1.stan" (folder STAN_Mod_Turb_2) which includes random
// trip and random site effects for both the intercept and slope
// Fixed the first coef for each taxa in the drift portion of the model 
// Fixed the intercept for worms in the diet portion of the model
// Added derived parms: mu_sp_all; beta_sp_all
// includes log(sz)

data {
  int Nspsz;           // Number of taxa & size bins      
  int Nst;             // Number of sites & trips
  int Nsp;             // Number of taxa
  int Nind;            // Number of individuals
  
  int idx_ts_ind[Nind];// Index matching individuals to trip, sites

  int u_idx[Nsp];      // Diet
  int u_idx2[Nsp];     // Diet
  int idx[Nspsz];      // Drift
  int idx2[Nspsz];     // Drift
  
  int idx_first[Nsp];  //Drift index for the first position of each taxa in the size bin taxa matrix (i.e., lprod)  
  int sp_idx2[Nspsz - Nsp];  //Drift new index for each taxa and size bin - 1, for each taxa
  int not_first[Nspsz - Nsp]; //Drift oposite of idx_first
  
  vector[Nspsz] sz;        // sizes - diet
  vector[Nind] fish_sz;    // centered & scaled fork length of fish
  real avg_log_len;        // average log length (size measure)
  
  int y[Nind, Nspsz];      // diet matrix 
  int w[Nind, Nsp];        // diet matrix - extra
  matrix[Nspsz, Nsp] X;    // dummy
  
  int a[Nst, Nspsz];       // drift matrix
  int w_a[Nst, Nsp];       // drift matrix - extra
  int sp[Nspsz];
  
  vector<lower = 0> [Nsp]alpha;  // dirichlet params 
  
  vector[Nind] turb;      // log turbidity, centered & scaled, expanded to match individual diet dims
} 

parameters {
  // Drift
  real<lower = 0> sig_mu;
  vector<lower = 0>[Nsp] sig_lprod;
  vector[Nspsz - Nsp] mu_lprod;
  matrix [Nst, Nspsz - Nsp]lprod;
  simplex[Nsp] tmp_ps_a[Nst];
  
  // Diet
  vector[Nsp-1] mu_sp;
  vector<lower = 0>[Nsp-1] sig_sp_trip;
  vector[Nsp] mu_sz;
  vector<lower = 0>[Nsp] sig_sz_trip;
  matrix[Nst, Nsp-1]ts_sp_eff;
  matrix[Nst, Nsp]ts_sz_eff;
  vector[Nsp-1] beta_sp_turb;
  vector[Nsp] beta_sz_turb;
  vector[Nsp] beta_f_sz_p_sz;
  vector[Nsp - 1] beta_f_sz_p_sp;
}

transformed parameters {
  // Drift
  matrix[Nst,Nspsz] eprod_a;
  matrix[Nst,Nspsz] n_eprod_a;
  matrix[Nst,Nspsz] p_a;
  matrix[Nst,Nsp] ps_a;
  matrix[Nst,Nspsz] sum_eprod_a;
  matrix[Nst,Nspsz] fix_lprod;
  
  // Diet
  matrix[Nind, Nspsz] p;
  matrix[Nind, Nsp] ps;
  matrix[Nind, Nspsz] beta1;
  matrix[Nind, Nspsz] eprod_u;
  vector[Nind] sum_eprod_u;
  matrix[Nind, Nsp] beta_sz;
  matrix[Nind, Nsp] fix_beta_sp; 

  // drift ////////////////////////////////////////////
  for(i in 1:Nst){                 
     ps_a[i,] = tmp_ps_a[i]'; 
  } 
  
  for(i in 1:Nst){
    for(j in 1:Nsp){
      fix_lprod[i,idx_first[j]] = 0.0;
    }
  }
  
  for(i in 1:Nst){
    for(j in 1:Nspsz - Nsp){
      fix_lprod[i,not_first[j]] = lprod[i,j];
    }
  }
  
    for(j in 1:Nspsz){
      eprod_a[,j] = exp(fix_lprod[,j]);
    }
    
    for(i in 1:Nst){       // error if I remove the lines below from the i loop
      for(j in 1:Nspsz){
      sum_eprod_a[i,j] = sum(eprod_a[i,idx[j]:idx2[j]]);
    }
      
      for(j in 1:Nspsz){
        n_eprod_a[i,j] = eprod_a[i,j] / sum_eprod_a[i,j];   
        p_a[i,j] = n_eprod_a[i,j] * ps_a[i,sp[j]];   
      }
      
    }
    
  // diet ////////////////////////////////////////////
  for(i in 1:Nind){
    fix_beta_sp[i,Nsp] = 0.0;
  }
  
  for(k in 1:Nsp - 1){
    fix_beta_sp[,k] = mu_sp[k] + ts_sp_eff[idx_ts_ind[],k] + beta_sp_turb[k] * turb + 
                      beta_f_sz_p_sp[k] * turb .* fish_sz;
  }
  
  for(k in 1:Nsp){
      beta_sz[,k] = mu_sz[k] + ts_sz_eff[idx_ts_ind[],k] + beta_sz_turb[k] * turb +
                    beta_f_sz_p_sz[k] * turb .* fish_sz;
   }
  
  for(i in 1:Nind){
    for(k in 1:Nspsz){
      beta1[i,k] = dot_product(fix_beta_sp[i,], X[k,]) + 
                   dot_product(beta_sz[i,], X[k,]) * (log(sz[k]) - avg_log_len);
    }  
  }
  
  for(k in 1:Nspsz){
    eprod_u[,k] = p_a[idx_ts_ind[],k] .* exp(beta1[,k]);  // where the avail comes in 
  }
     
  for(i in 1:Nind){ 
    sum_eprod_u[i] = sum(eprod_u[i,]);
  }
     
  for(k in 1:Nspsz){
    p[,k] = eprod_u[,k] ./ sum_eprod_u;
  }  
  
  for(i in 1:Nind){    
     for(k in 1:Nsp){
       ps[i,k] = sum(p[i, u_idx[k]:u_idx2[k]]);
     }
  }
}

model {
   // drift ////////////////////////////////////////////
  for(i in 1:Nst){
    tmp_ps_a[i] ~ dirichlet(alpha[]);
  }
  
  sig_mu ~ normal(0, 10);  
  sig_lprod ~ normal(0, 10);
  mu_lprod ~ normal(0, sig_mu); 
  
  for(j in 1:Nspsz - Nsp){
    lprod[,j] ~ normal(mu_lprod[j], sig_lprod[sp_idx2[j]]);
  }
    
  for(i in 1:Nst){
    w_a[i,1:Nsp] ~ multinomial(ps_a[i,]');
    a[i,1:Nspsz] ~ multinomial(p_a[i,]');
  }
  
   // diet ////////////////////////////////////////////
   beta_f_sz_p_sz ~ normal(0,10);
   beta_f_sz_p_sp ~ normal(0,10);
   mu_sp ~ normal(0, 10);
   
   beta_sp_turb ~ normal(0, 10);
   beta_sz_turb ~ normal(0, 10);
   
   sig_sp_trip ~ normal(0, 10);
   sig_sz_trip ~ normal(0, 10);
   
   for(k in 1:Nsp-1){
     ts_sp_eff[,k] ~ normal(0, sig_sp_trip[k]);  
   }
   
   for(k in 1:Nsp){
     ts_sz_eff[,k] ~ normal(0, sig_sz_trip[k]);  
   }
   
   mu_sz ~ normal(0, 10);
    
  for(i in 1:Nind){         
    y[i,1:Nspsz] ~ multinomial(p[i,]');
    w[i,1:Nsp] ~ multinomial(ps[i,]');
  }
}

generated quantities {
  vector[Nsp] mu_sp_all;
  vector[Nsp] mu_sp_tmp;
  
  vector[Nsp] beta_t_tmp;
  vector[Nsp] beta_sp_turb_all;
  matrix[Nst, Nsp] beta_sp_all;
  
  vector[Nsp] beta_f_sz_p_sp_all;
  vector[Nsp] beta_f_sz_p_sp_tmp;

  // matrix[Nst, Nspsz] tmp_st;
  // matrix[Nst, Nsp] tmp_st_sum;
  // matrix[Nst, Nsp] gamma_st;

  ///////////////////////////////////
  mu_sp_tmp[Nsp] = 0.0;

  for(i in 1:Nsp-1){
    mu_sp_tmp[i] = mu_sp[i];
  }

  for(i in 1:Nsp){
    mu_sp_all[i] =  mu_sp_tmp[i] - mean(mu_sp_tmp[]);
  }

  /////////////////////////////////
  beta_t_tmp[Nsp] = 0.0;

  for(i in 1:Nsp-1){
    beta_t_tmp[i] = beta_sp_turb[i];
  }

  for(i in 1:Nsp){
    beta_sp_turb_all[i] =  beta_t_tmp[i] - mean(beta_t_tmp[]);
  }

  for(i in 1:Nst){
    for(j in 1:Nsp){
      beta_sp_all[i,j] = fix_beta_sp[i,j] - mean(fix_beta_sp[i,]);
    }
  }
  /////////////////////////////////
  beta_f_sz_p_sp_tmp[Nsp] = 0.0;
  
  for(i in 1:Nsp-1){
    beta_f_sz_p_sp_tmp[i] = beta_f_sz_p_sp[i];
  }
  
  for(i in 1:Nsp){
    beta_f_sz_p_sp_all[i] = beta_f_sz_p_sp_tmp[i] - mean(beta_f_sz_p_sp_tmp[]);
  }

  // ///////////////////////////////////
  // // Calc. gamma_st (predicted diet proportions) for each site and trip
  // for(j in 1:Nspsz){
  //   for(i in 1:Nst){
  //      tmp_st[i,j] = p_a[i,j] * exp(fix_beta_sp[i,sp[j]] + beta_sz[i, sp[j]] * log(sz[j]));
  //   }
  // }
  // 
  // for(i in 1:Nst){
  //   for(j in 1:Nsp){
  //     tmp_st_sum[i,j] = sum(tmp_st[i,u_idx[j]:u_idx2[j]]);
  //   }
  // }
  // 
  // for(i in 1:Nst){
  //   for(j in 1:Nsp){
  //     gamma_st[i,j] = tmp_st_sum[i,j] / sum(tmp_st_sum[i,]);
  //   }
  // }
}

