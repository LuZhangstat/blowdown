/* Response NNGP model */

functions{
  matrix Block_COV(matrix coords, int[] ind_ls, real phi, int nb){
    // coords: coords of the ALS variables
    // int_ls: records the first index for each plot
    // nb: number of observed response
    // phi: parameter
    
    matrix[nb, nb] C_B;
    
    for(i in 1:nb){
      for(j in i:nb){
        C_B[i, j] = 0.0;
        
        if(i == j){
          for(k in ind_ls[i]:(ind_ls[i+1]-1)){
            for(l in k:(ind_ls[j+1]-1)){
              if(l == k){
                C_B[i, j] += exp(-phi * distance(coords[k,:], coords[l,:]));
              }else{
                C_B[i, j] += 2 * exp(-phi * distance(coords[k,:], coords[l,:]));
              }
            }
          }
          C_B[i, j] = C_B[i, j] / (ind_ls[i+1]-ind_ls[i])^2;
        }else{
          for(k in ind_ls[i]:(ind_ls[i+1]-1)){
            for(l in ind_ls[j]:(ind_ls[j+1]-1)){
              C_B[i, j] += exp(-phi * distance(coords[k,:], coords[l,:]));
            }
          }
          C_B[i, j]  = C_B[i, j] / ((ind_ls[i+1]-ind_ls[i]) * 
          (ind_ls[j+1]-ind_ls[j]));
          C_B[j, i] = C_B[i, j];
        }
      }
    }
    return C_B;
  }
  
  real COSP_lpdf(vector y, matrix HX, matrix C_B, real sigmasq, real tausq,
                 vector Dh, int nb, int p){
    // y: vector of observed response
    // HX: averaged ALS variables
    // C_B: output of function Block_COV
    // sigmasq, tausq: parameters
    // mu_beta, V_beta: mean and covariance matrix of the prior of regression coefficient beta
    // Dh: scales of noise of the responses
    
    matrix[nb, nb] V_star;
    matrix[nb, nb] V_star_CholL;
    matrix[p, p] inv_V_beta_star;
    vector[nb] inv_V_star_CholL_y;
    matrix[nb, p] inv_V_star_CholL_HX;
    matrix[p, p] V_beta_ast;
    vector[p] mu_beta_ast;
    
   
    V_star = add_diag(sigmasq * C_B, tausq * Dh);
    V_star_CholL = cholesky_decompose(V_star);
    inv_V_star_CholL_y = mdivide_left_tri_low(V_star_CholL, y);
    inv_V_star_CholL_HX = mdivide_left_tri_low(V_star_CholL, HX);
    V_beta_ast = inverse_spd(crossprod(inv_V_star_CholL_HX));
    mu_beta_ast = V_beta_ast * (inv_V_star_CholL_HX' * inv_V_star_CholL_y);
    
    return (multi_normal_cholesky_lpdf(y | rep_vector(0.0, nb), V_star_CholL) - 
        multi_normal_lpdf(mu_beta_ast | rep_vector(0.0, p), V_beta_ast));
  }
  
  // vector beta_omega_rng(vector y, matrix HX, matrix C_B, real sigmasq,
  //                       real tausq, vector invVmu_beta, matrix invV_beta, 
  //                       vector Dh, int nb, int p){
    //                         
    //   // recover posterior samples of beta and omega
    //   // p: number of predictors
    //   
    //   vector[nb+p] m;
    //   matrix[nb+p, nb+p] invM;
    //   matrix[nb+p, nb+p] M;
    //   matrix[nb, nb] invC_B;
    //   
    //   invC_B = inverse_spd(C_B);
    //   
    //   invM = crossprod((diag_matrix(sqrt(1 ./ (Dh * tausq))) * 
    //                     append_col(HX, identity_matrix(nb))));
    //   invM[1:p, 1:p] = invM[1:p, 1:p] + invV_beta;
    //   invM[(p+1):(nb+p), (p+1):(nb+p)] = invM[(p+1):(nb+p), (p+1):(nb+p)] +
    //         invC_B / sigmasq;
    //   
    //   M = inverse_spd(invM);
    //   m[1:p] = HX'* (diag_matrix(1 ./ (Dh * tausq)) * y) + invVmu_beta;
    //   m[(p+1):(p+nb)] = diag_matrix(1 ./ (Dh * tausq)) * y;
    //   
    //   return multi_normal_rng(M*m, M);
    // }
    
}

data {
  int<lower=1> na;
  int<lower=1> nb;
  int<lower=1> p;
  vector[nb] y;
  matrix[nb, p] HX;
  vector[nb] Dh;
  matrix[na, 2] gridA;
  int plotid_ind[nb + 1];
  real ap;
  real bp;
  real ss;
  real st;
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
  matrix[nb, nb] V_star;
  vector[nb] Mu;
  
  phi ~ gamma(ap, bp);
  sigma ~ normal(0, ss);
  tau ~ normal(0, st);
  C_B = Block_COV(gridA, plotid_ind, phi, nb);
  
  //V_star = add_diag(sigmasq * C_B, tausq * Dh);
  //y ~ multi_normal(Mu, COV);  
  //C_B =  identity_matrix(nb);
  y ~ COSP(HX, C_B, sigmasq, tausq, Dh, nb, p);
}

// generated quantities {
  //   vector[p + nb] beta_omega_B;
  //   {
    //     matrix[nb, nb] C_B;
    //     C_B = Block_COV(gridA, plotid_ind, phi, nb); 
    //   }
    //   beta_omega_B = beta_omega_rng(y, HX, C_B, sigmasq, tausq, invVmu_beta, 
    //                  invV_beta, Dh, nb, p);
    // }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
