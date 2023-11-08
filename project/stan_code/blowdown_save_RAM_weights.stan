/* COSP model */

functions{
  matrix Block_COV(matrix coords, int[] ind_ls, real phi, int nb, vector hA){
    // coords: coords of the ALS variables
    // int_ls: records the first index for each plot
    // nb: number of observed response
    // phi: parameter
    // hA: weights (h_{li}) for each cell in the fine grid
    
    matrix[nb, nb] C_B;
    
    for(i in 1:nb){
      for(j in i:nb){
        C_B[i, j] = 0.0;
        
        if(i == j){
          for(k in ind_ls[i]:(ind_ls[i+1]-1)){
            for(l in k:(ind_ls[j+1]-1)){
              if(l == k){
                C_B[i, j] += exp(-phi * distance(coords[k,:], coords[l,:])) * 
                          (hA[k]^2);
              }else{
                C_B[i, j] += 2 * exp(-phi * distance(coords[k,:], coords[l,:])) *
                          (hA[k] * hA[l]);
              }
            }
          }
        }else{
          for(k in ind_ls[i]:(ind_ls[i+1]-1)){
            for(l in ind_ls[j]:(ind_ls[j+1]-1)){
              C_B[i, j] += exp(-phi * distance(coords[k,:], coords[l,:])) *
                           (hA[k] * hA[l]);
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
  C_B = Block_COV(gridA, plotid_ind, phi, nb, hA);
  
  COV = add_diag(sigmasq * C_B + HX * V_beta * HX', tausq * Dh);
  Mu = HX * mu_beta;
  
  y ~ multi_normal(Mu, COV);  
}

    
    
    
    
    
    
    
    
    
    
    
    
    
