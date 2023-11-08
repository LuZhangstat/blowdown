/* benchmark model */

functions{
  matrix COV_exp(matrix coords, real phi, int nb){
    // coords: coords of the grid
    // nb: number of the locations
    // phi: parameter
    
    matrix[nb, nb] C_B;
    
    for(i in 1:nb){
      for(j in i:nb){
        C_B[i, j] = exp(-phi * distance(coords[i,:], coords[j,:]));
        if(i != j){
            C_B[j, i] = C_B[i, j];
        }
      }
    }
    return C_B;
  }
  
}

data {
  int<lower=1> n;  // No. of obs
  int<lower=1> p;
  vector[n] y;
  matrix[n, p] X;
  matrix[n, 2] gridA;
  vector[p] mu_beta;
  matrix[p, p] V_beta;
  real ap;
  real bp;
  real ss;
  real st;
}

transformed data {
  vector[n] onesv;
  onesv = rep_vector(1.0, n);
}


parameters{
  real<lower = 0> sigma;
  real<lower = 0> tau;
  real<lower = 0> phi;
  vector[p] beta; // regression coef
}

transformed parameters {
  real sigmasq = square(sigma);
  real tausq = square(tau);
}

model{
  matrix[n, n] C_B;
  matrix[n, n] COV;
  //vector[nb] Muw;
  
  beta ~ multi_normal(mu_beta, V_beta);
  phi ~ gamma(ap, bp);
  sigma ~ normal(0, ss);
  tau ~ normal(0, st);
  
  C_B = COV_exp(gridA, phi, n);
  COV = add_diag(sigmasq * C_B, tausq * onesv);
  
  y ~ multi_normal(X * beta, COV);
}

    
    
    
    
    
    
    
    
    
    
    
    
    
