/* COSP model */

functions{
  matrix Block_COV(matrix coords, int[] ind_ls, real phi, int nb, vector hA,
    matrix Dist_M, real gamma){
    // coords: coords of the ALS variables
    // int_ls: records the first index for each plot
    // nb: number of observed response
    // phi: parameter
    // hA: weights (h_{li}) for each cell in the fine grid
    // gamma: hyperparameter in the tapering kernel
    
    matrix[nb, nb] C_B;
    real d_temp;
    
    for(i in 1:nb){
      for(j in i:nb){
        C_B[i, j] = 0.0;
        if(i == j){
          for(k in ind_ls[i]:(ind_ls[i+1]-1)){
            for(l in k:(ind_ls[j+1]-1)){
              d_temp = Dist_M[k, l];
              if(l == k){
                C_B[i, j] += exp(-phi * d_temp) * (hA[k]^2);
              }else{
                if(d_temp < gamma){
                    C_B[i, j] += 2 * exp(-phi * d_temp) * 
                        pow((1 - d_temp / gamma), 4) * (1 + 4 * d_temp / gamma) *
                          (hA[k] * hA[l]);
                }
              }
            }
          }
        }else{
          if(min(Dist_M[ind_ls[i]:(ind_ls[i+1]-1), ind_ls[j]:(ind_ls[j+1]-1)]) 
              < gamma){
             for(k in ind_ls[i]:(ind_ls[i+1]-1)){
                for(l in ind_ls[j]:(ind_ls[j+1]-1)){
                    d_temp = Dist_M[k, l];
                    if(d_temp < gamma){
                        C_B[i, j] += exp(-phi * d_temp) *
                           pow((1 - d_temp / gamma), 4) * (1 + 4 * d_temp / gamma) * 
                           (hA[k] * hA[l]);
                    }
                }
            }
          }
          C_B[j, i] = C_B[i, j];
        }
      }
    }
    return C_B;
  }
}

data {
  int<lower=1> na;
  int<lower=1> nb;
  int<lower=1> p;
  vector[nb] y;
  matrix[nb, p] HX;
  vector[nb] Dh;
  matrix[na, 2] gridA;
  vector[na] hA;
  int plotid_ind[nb + 1];
  vector[p] mu_beta;
  matrix[p, p] V_beta;
  real ap;
  real bp;
  real ss;
  real st;
  matrix[na, na] Dist_M;
  real gamma;
}

transformed data {
  matrix[p, p] invV_beta;
  vector[p] invVmu_beta;
  invV_beta = inverse_spd(V_beta);
  invVmu_beta = invV_beta * mu_beta;
}

parameters{
  real<lower = 0> sigma;
  real<lower = 0> tau;
  real<lower = 0> phi;
}

transformed parameters {
  real sigmasq = square(sigma);
  real tausq = square(tau);
}

model{
  matrix[nb, nb] C_B;
  
  matrix[nb, nb] COV;
  vector[nb] Mu;
  
  phi ~ gamma(ap, bp);
  sigma ~ normal(0, ss);
  tau ~ normal(0, st);
  C_B = Block_COV(gridA, plotid_ind, phi, nb, hA, Dist_M, gamma);
  
  COV = add_diag(sigmasq * C_B + HX * V_beta * HX', tausq * Dh);
  Mu = HX * mu_beta;
  
  y ~ multi_normal(Mu, COV);  
}

    
    
    
    
    
    
    
    
    
    
    
    
    
