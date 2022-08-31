/* COSP model */

functions{
  matrix Block_COV(matrix coords, int[] ind_ls, real phi, int nb, vector hA){
    // coords: coords of the ALS variables
    // int_ls: records the first index for each plot
    // nb: number of observed response
    // phi: parameter
    // hA: weights (hij) for each coords
    
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
  
  real COSP_lpdf(vector y, matrix HX, matrix C_B, real sigmasq, real tausq,
  vector mu_beta, matrix V_beta, vector Dh, int nb){
    // y: vector of observed response
    // HX: averaged ALS variables
    // C_B: output of function Block_COV
    // sigmasq, tausq: parameters
    // mu_beta, V_beta: mean and covariance matrix of the prior of regression coefficient beta
    // Dh: scales of noise of the responses
    
    matrix[nb, nb] COV;
    vector[nb] Mu;
    
    COV = sigmasq * C_B + HX * V_beta * HX' + tausq * diag_matrix(Dh);
    Mu = HX * mu_beta;
    
    return multi_normal_lpdf(y | Mu, COV);
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
  //y ~ COSP(HX, C_B, sigmasq, tausq, mu_beta, V_beta, Dh, nb);
}

    
    
    
    
    
    
    
    
    
    
    
    
    
