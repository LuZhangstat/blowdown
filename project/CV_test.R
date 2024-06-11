### CV test for real data analysis ###
rm(list = ls())
library(Matrix)
library(fields)
library(ggplot2)
source("./project/utils.R")
library(dplyr)
library(caret)

# Compute RMSPE
rmspe <- function(observed, predicted) {
  residuals <- observed - predicted
  squared_residuals <- residuals^2
  mean_squared_residuals <- mean(squared_residuals)
  sqrt(mean_squared_residuals)
}


plot.y <- read.csv("./data/outcome_y_plot.csv")
plot.y %>% glimpse

plot.grid.x <- read.csv("./data/predictor_x_grid.csv")
plot.grid.x %>% glimpse

## generate data on a coarser grid ##
# create the shared location for every 9 cells
plot.grid.x$coord.x2 <- (plot.grid.x$coord.x + 1 * (plot.grid.x$coord.x %% 3 == 1) - 
                           1 * (plot.grid.x$coord.x %% 3 == 0))
plot.grid.x$coord.y2 <- (plot.grid.x$coord.y + 1 * (plot.grid.x$coord.y %% 3 == 1) - 
                           1 * (plot.grid.x$coord.y %% 3 == 0))

plot.grid.x_lit <- plot.grid.x %>% 
  group_by(obs.plot.id, coord.x2, coord.y2) %>% 
  summarise(height.m = mean(height.m), n = n()) %>% 
  inner_join(plot.y, by = "obs.plot.id") %>% arrange(obs.plot.id) %>%
  select(sub.region, obs.plot.id, coord.x2, coord.y2, height.m, n) 
plot.grid.x_lit %>% glimpse()


## Compute the average height.m over each plot and summary with observation in one dataset
x_obs_bar <- plot.grid.x_lit %>% group_by(obs.plot.id) %>% 
  summarize(height.mean = sum(height.m * n) / sum(n)) %>% select(obs.plot.id, 
                                                                 height.mean) 
obs_xy <- plot.y %>% select(obs.plot.id, sub.region, volume.m3.ha) %>% 
  inner_join(x_obs_bar, by = "obs.plot.id") %>% 
  mutate(intercept = 1) %>%
  # select(volume.m3.ha, intercept, height.mean)
  mutate(x.Frohn = height.mean * (sub.region == "Frohn")) %>%
  mutate(x.Laas = height.mean * (sub.region == "Laas")) %>%
  mutate(x.Mauthen = height.mean * (sub.region == "Mauthen")) %>%
  mutate(x.Liesing = height.mean * (sub.region == "Liesing")) %>%
  mutate(x.Ploecken = height.mean * (sub.region == "Ploecken")) %>%
  select(volume.m3.ha, intercept, x.Frohn, x.Laas, x.Mauthen, x.Liesing, 
         x.Ploecken)

###################
## Model fitting ##
###################

## precomputation for the model fitting ##
HX = as.matrix(obs_xy[, -1]); y = obs_xy[, 1] # HX: the H_{BA} X in the paper
plot.grid.x_lit <- plot.grid.x_lit %>% group_by(obs.plot.id) %>% 
  mutate(counts = sum(n)) %>% mutate(weight = n / counts)
plot.grid.x_lit %>% glimpse()

Dh = plot.grid.x_lit %>% group_by(obs.plot.id) %>% # the sum_{i = 1}^{n_a} h^2_{li} 
  summarise(Dh = sum(weight^2)) %>% pull

## the coords of ALS variables 
grid.A <- plot.grid.x_lit %>% # select predictors in region Frohn
  select(coord.x2, coord.y2) # coordinates of predictors over regions with observation
plotid_ind <- c(order(plot.grid.x_lit$obs.plot.id)[      # record the first index of each plot
  !duplicated(plot.grid.x_lit$obs.plot.id)], nrow(grid.A) + 1)       

library(cmdstanr)
library(bayesplot)
file <- file.path(getwd(), "project/stan_code/blowdown_save_RAM_weights_tapering_phi_unif.stan")
mod <- cmdstan_model(file)


# cross-validation #
set.seed(123) 
k = 5
folds <- createFolds(1:nrow(obs_xy), k = k, list = TRUE, 
                     returnTrain = FALSE)
fit_flag = FALSE
if(fit_flag){
  for (i in 1:k) {
    cat(i, "\n")
    # Indices for the prediction set
    prediction_indices <- folds[[i]]
    # Indices for the training set
    training_indices <- setdiff(1:nrow(obs_xy), prediction_indices)
    
    # Create training and prediction sets
    HX_T <-  HX[training_indices, ]; y_T = y[training_indices]
    Dh_T <- Dh[training_indices]
    
    HXU_T <- HX[prediction_indices, ]
    DhU_T <- Dh[prediction_indices]
    
    ## the coords of ALS variables 
    plot.grid.x_lit_T = plot.grid.x_lit %>% 
      filter(obs.plot.id %in% training_indices)
    
    grid.A_T <- plot.grid.x_lit_T %>% select(coord.x2, coord.y2) # coordinates of predictors over regions with observation
    plotid_ind_T <- c(order(plot.grid.x_lit_T$obs.plot.id)[      # record the first index of each plot
      !duplicated(plot.grid.x_lit_T$obs.plot.id)], nrow(grid.A_T) + 1)       
    
    hA_T = plot.grid.x_lit_T$weight
    
    pred.grid.x_lit_T <- plot.grid.x_lit %>%   # store withheld data in pred.grid.x_lit_T
      filter(obs.plot.id %in% prediction_indices)
    
    grid.Au_T <- pred.grid.x_lit_T %>% select(coord.x2, coord.y2) # coordinates of predictors over regions for prediction
    predid_ind_T <- c(order(pred.grid.x_lit_T$obs.plot.id)[     # record the first index of each polygon
      !duplicated(pred.grid.x_lit_T$obs.plot.id)], nrow(grid.Au_T) + 1)
    
    na = nrow(grid.A_T); nb = length(y_T); p = ncol(HX_T);
    
    # compute and store the distance matrix #
    # cloest distance between plots in Mauthen and Ploecken is 1.876 km
    # furtherest distance amony responses in Frohn is 2.61 km
    gamma = 1.8
    Dist_M_T <- rdist(grid.A_T[, -1] / 1000, grid.A_T[, -1] / 1000) 
    
    
    #-------------------------- Set parameters of priors --------------------------#
    mu_beta = rep(0, p)     # mean vector in the Gaussian prior of beta
    V_beta = diag(p) * 1000    # covariance matrix in the Gaussian prior of beta
    ss = 20000       # scale parameter in the inverse gamma prior of sigma 
    st = 10000       # scale parameter in the inverse gamma prior of tau  
    ap = 3/500; bp = 3/0.1       # lower and upper bound of uniform prior of phi 
    
    
    data <- list(na = na, nb = nb, p = p, y = y_T, HX = HX_T, Dh = Dh_T, 
                 gridA = grid.A_T[, -1] / 1000, hA = hA_T,
                 plotid_ind = plotid_ind_T,
                 mu_beta = mu_beta, V_beta = V_beta,
                 ap = ap, bp = bp, ss = ss, st = st, Dist_M = Dist_M_T, 
                 gamma = gamma)
    
    fit <- mod$sample(
      data = data,
      seed = 1234,
      chains = 4,
      parallel_chains = 4,
      refresh = 50,
      save_warmup = TRUE,
      iter_warmup = 500,
      iter_sampling = 500,
      sig_figs = 18
    )
    #7439.5s
    filename = paste0("./results/fit_phi_unif_CV", i, "_v2.RDS")
    fit$save_object(file = filename)
    fit$summary()
    fit_draws <- fit$draws(inc_warmup = FALSE)
    
    pick_sample_id <- seq(5, 2000, by = 10) # pick 200 samples
    
    phi_ls <- c(fit_draws[, , "phi"])[pick_sample_id]
    sigmasq_ls <- c(fit_draws[, , "sigmasq"])[pick_sample_id]
    tausq_ls <- c(fit_draws[, , "tausq"])[pick_sample_id]
    coords_A = grid.A_T[, -1] /1000
    coords_AU = grid.Au_T[, -1] / 1000
    ind_ls_B = plotid_ind_T
    ind_ls_BU = predid_ind_T
    hA = plot.grid.x_lit_T$weight
    hAU = pred.grid.x_lit_T$weight
    
    ## recover posterior samples of beta and omega for both observed and unobserved plots ##
    t <- proc.time()
    beta_omega_sam <- sample_beta_omega_h_tapering(
      phi_ls, sigmasq_ls, tausq_ls, coords_A, coords_AU, hA, hAU, ind_ls_B, 
      ind_ls_BU, HX_T, Dh_T, y_T, mu_beta, V_beta, gamma, flat_prior = FALSE)
    proc.time() -t
    
    filename2 = paste0("./results/RDA_recov_sam_CV", i, "_v2.RDS")
    save(beta_omega_sam, file = filename2)
  }
}
# prediction returned by COS
y_pre_sam <- matrix(NA, nrow = nrow(obs_xy), ncol = 200)
for(i in 1:k){
  cat(i, "\n")
  # Indices for the prediction set
  prediction_indices <- folds[[i]]
  # Indices for the training set
  training_indices <- setdiff(1:nrow(obs_xy), prediction_indices)
  
  # Create training and prediction sets
  HX_T <-  HX[training_indices, ]; y_T = y[training_indices]
  Dh_T <- Dh[training_indices]
  
  HXU_T <- HX[prediction_indices, ]
  DhU_T <- Dh[prediction_indices]
  
  ## the coords of ALS variables 
  plot.grid.x_lit_T = plot.grid.x_lit %>% 
    filter(obs.plot.id %in% training_indices)
  
  grid.A_T <- plot.grid.x_lit_T %>% select(coord.x2, coord.y2) # coordinates of predictors over regions with observation
  plotid_ind_T <- c(order(plot.grid.x_lit_T$obs.plot.id)[      # record the first index of each plot
    !duplicated(plot.grid.x_lit_T$obs.plot.id)], nrow(grid.A_T) + 1)       
  
  hA_T = plot.grid.x_lit_T$weight
  
  pred.grid.x_lit_T <- plot.grid.x_lit %>%   # store withheld data in pred.grid.x_lit_T
    filter(obs.plot.id %in% prediction_indices)
  
  grid.Au_T <- pred.grid.x_lit_T %>% select(coord.x2, coord.y2) # coordinates of predictors over regions for prediction
  predid_ind_T <- c(order(pred.grid.x_lit_T$obs.plot.id)[     # record the first index of each polygon
    !duplicated(pred.grid.x_lit_T$obs.plot.id)], nrow(grid.Au_T) + 1)
  
  na = nrow(grid.A_T); nb = length(y_T); p = ncol(HX_T);
  
  # compute and store the distance matrix #
  # cloest distance between plots in Mauthen and Ploecken is 1.876 km
  # furtherest distance amony responses in Frohn is 2.61 km
  gamma = 1.8
  Dist_M_T <- rdist(grid.A_T[, -1] / 1000, grid.A_T[, -1] / 1000) 
  
  filename = paste0("./results/fit_phi_unif_CV", i, "_v2.RDS")
  CVRDS1 <- readRDS(filename)
  fit_draws <- CVRDS1$draws(inc_warmup = FALSE)
  
  pick_sample_id <- seq(5, 2000, by = 10) # pick 200 samples
  
  phi_ls <- c(fit_draws[, , "phi"])[pick_sample_id]
  sigmasq_ls <- c(fit_draws[, , "sigmasq"])[pick_sample_id]
  tausq_ls <- c(fit_draws[, , "tausq"])[pick_sample_id]
  coords_A = grid.A_T[, -1] /1000
  coords_AU = grid.Au_T[, -1] / 1000
  ind_ls_B = plotid_ind_T
  ind_ls_BU = predid_ind_T
  hA = plot.grid.x_lit_T$weight
  hAU = pred.grid.x_lit_T$weight
  filename2 = paste0("./results/RDA_recov_sam_CV", i, "_v2.RDS")
  load(filename2)
  
  yU_ls <- pred_sample_y(beta_omega_sam$beta_ls, beta_omega_sam$omega_BU_ls, 
                         tausq_ls, HXU_T, DhU_T)
  
  yU_pre <-  data.frame(pred.poly.id = unique(pred.grid.x_lit_T$obs.plot.id), 
                        yU_pm = rowMeans(yU_ls),
                        omega_BU_pm = colMeans(t(beta_omega_sam$omega_BU_ls) + 
                                                 beta_omega_sam$beta_ls[1,]))
  y_pre_sam[prediction_indices, ] <- yU_ls
}
y.hat <- apply(y_pre_sam, 1, mean)
y.var <- apply(y_pre_sam, 1, var)
y.sd <- apply(y_pre_sam, 1, sd)
y.95.bounds <- t(apply(y_pre_sam, 1, quantile, c(0.025, 0.975)))
cos <- c(sum((y.95.bounds[, 1] < y) & (y.95.bounds[, 2] > y)) / length(y),
         rmspe(y, y.hat), mean(abs(y.hat - y)),
         mean(crps_norm(y, y.hat, y.var)))

plot(y.hat, y)
abline(a = 0, b = 1)

# 
library(spBayes) ##Using spBayes_0.4-4 (prior versions have a bug with spSVC LOO)
library(scoringRules)
library(knitr)

set.seed(123)
## SVI.
starting.svi <- list("phi"=3/1, "sigma.sq"=20000, "tau.sq"=10000)
tuning.svi <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)
priors.svi <- list("phi.Unif"=list(3/500, 3/0.1),
                   "sigma.sq.IG"=list(2, 20000),
                   "tau.sq.IG"=c(2, 10000))

## Save LOO posterior predictive samples.
n.samples <- 25000
nonsp.CV.pred.samples <- matrix(NA, length(y), 500)
svi.CV.pred.samples <- matrix(NA, length(y), 0.25*n.samples/2 +1)

for(i in 1:k){
  
  y.loo <- y[-folds[[i]]]
  x.loo <- as.matrix(obs_xy[-folds[[i]], 3:7])
  coords.loo <- plot.y[-folds[[i]], 4:5]/1000
  
  
  y.ho <- y[folds[[i]]]
  x.ho <- as.matrix(obs_xy[folds[[i]], 3:7])
  coords.ho <- plot.y[folds[[i]], 4:5]/1000
  
  p <- 1+5
  
  ################################################
  ## Non-spaial model
  ################################################
  m.0 <- bayesLMRef(lm(y.loo ~ x.loo), n.samples = 500)
  
  ################################################
  ## Spatial varying intercept 
  ################################################
  cov.model <- "exponential"
  m.1 <- spSVC(y.loo~x.loo, coords=coords.loo, starting=starting.svi, svc.cols=1,
               tuning=tuning.svi, priors=priors.svi, cov.model=cov.model,
               n.samples=n.samples, verbose = FALSE)
  
  m.1 <- spRecover(m.1, start=floor(0.75*n.samples), thin=2,  verbose = FALSE)
  
  ################################################
  ## Prediction
  ################################################
  ##Non-spatial
  beta.samples <- m.0$p.beta.tauSq.samples[,1:p]
  tau.sq.samples <- m.0$p.beta.tauSq.samples[,p+1]
  
  nonsp.CV.pred.samples[folds[[i]], ] <- cbind(1,x.ho) %*% t(beta.samples) + 
    matrix(rnorm(length(tau.sq.samples) * length(folds[[i]]), 
                 sqrt(tau.sq.samples)), 
           nrow = length(folds[[i]]), ncol = length(tau.sq.samples), 
           byrow = TRUE)
  
  
  ##SVI
  svi.CV.pred.samples[folds[[i]], ] <- 
          spPredict(m.1, pred.covars=cbind(1,x.ho), 
                    pred.coords=coords.ho, 
                    verbose = FALSE)$p.y.predictive.samples
  
  print(i)
}

##save.image(file="sp_gp_x_val.RData")
##load("sp_gp_x_val.RData")


## Summarize LOO results.
y.hat <- apply(nonsp.CV.pred.samples, 1, mean)
y.var <- apply(nonsp.CV.pred.samples, 1, var)
y.sd <- apply(nonsp.CV.pred.samples, 1, sd)
y.95.bounds <- t(apply(nonsp.CV.pred.samples, 1, quantile, c(0.025, 0.975)))
non.sp <- c(sum((y.95.bounds[, 1] < y) & (y.95.bounds[, 2] > y)) / length(y),
            rmspe(y, y.hat), mean(abs(y.hat - y)),
            mean(crps_norm(y, y.hat, y.var)))

y.hat <- apply(svi.CV.pred.samples, 1, mean)
y.var <- apply(svi.CV.pred.samples, 1, var)
y.sd <- apply(svi.CV.pred.samples, 1, sd)
y.95.bounds <- t(apply(svi.CV.pred.samples, 1, quantile, c(0.025, 0.975)))
svi <- c(sum((y.95.bounds[, 1] < y) & (y.95.bounds[, 2] > y)) / length(y),
         rmspe(y, y.hat), mean(abs(y.hat - y)),
         mean(crps_norm(y, y.hat, y.var)))

dt <- t(as.data.frame(round(rbind(non.sp, svi, cos), 2)))
rownames(dt) <- c("CI cover", "RMSPE", "MAE", "CRPS")
kable(dt, "latex", col.names = c("Non-spatial", "SVI", "COS"), 
      booktabs = TRUE)

# \begin{tabular}{lrrr}
# \toprule
# & Non-spatial & SVI & COS\\
# \midrule
# CI cover & 0.35 & 0.97 & 0.90\\
# RMSPE & 225.91 & 163.47 & 165.80\\
# MAE & 187.14 & 133.48 & 130.56\\
# CRPS & 768.55 & 5566.66 & 5480.68\\
# \bottomrule
# \end{tabular}

