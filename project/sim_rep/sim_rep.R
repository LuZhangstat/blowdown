rm(list = ls())
library(dplyr)
library(Matrix)
library(fields)
library(ggplot2)
source("./project/utils.R")

## read in arguments from command line ##
arg <- commandArgs(trailingOnly = FALSE)
myarg <- arg[length(arg)]
myarg <- sub("-","",myarg)
input_id <- as.numeric(myarg)
cat("id = ", input_id)

## simulate data ##
set.seed(input_id) 
side.A = seq(from = 1/54, to = 1-1/54, by = 1/27)
grid.A = expand.grid(side.A, side.A) # grid on the finest resolution
N_A = nrow(grid.A)
D <- as.matrix(dist(grid.A))         # distance matrix of the finest grid
X <- as.matrix(cbind(1, rnorm(nrow(grid.A))))
B <- as.matrix(c(1, 5))
sigma.sq <- 2
tau.sq <- 1
phi <- 3 / 0.6
R <- exp(- phi * D)
w_A <- rmvn(1, rep(0, N_A), sigma.sq * R)
y_A <- rnorm(N_A, X %*% B + w_A, sqrt(tau.sq))

# target units #
plot.centroid =  expand.grid(side.A[seq(2, 27, by = 3)], side.A[seq(2, 27, by = 3)])
D_plot <- rdist(grid.A, plot.centroid) 
A_plot_id <- (D_plot < (sqrt(2)/28 + 0.01)) %*% c(1:81) # plot id on the finest grid
# A_plot_id <- 1:N_A
dt_A <- data.frame(y_A = y_A, w_A = w_A, x_A = X[, 2], plot_id = A_plot_id, 
                   coord.x = grid.A[, 1], coord.y = grid.A[, 2])
dt_A = dt_A %>% arrange(plot_id)
hold_ls = c(20, 21, 22, 29, 31, 38, 40, 47, 49, 56, 57, 58,
            60, 51, 42, 33, 24, 62, 52, 34, 26, 61, 53, 35, 25, 43)
obs_ls = c(1:81)[-hold_ls]
obs_ind <- sapply(dt_A$plot_id, f <- function(x){any(x == obs_ls)})
w_B = dt_A[obs_ind, ] %>% group_by(plot_id) %>% 
  summarize(w_B = mean(w_A)) %>% arrange(plot_id) %>% select(w_B) %>% pull 

# predict region
O_d1 <- rdist(grid.A, t(c(0.2778, 0.69))) 
O_d2 <- rdist(grid.A, t(c(0.2778, 0.31))) 
ind_O <- ((O_d1 + O_d2) < 0.49 & (O_d1 + O_d2) > 0.43)

ind_K1 <- (grid.A[, 1] >0.6 & grid.A[, 1] < 0.68 & grid.A[, 2] < 0.75 & 
             grid.A[, 2] > 0.25)

ind_K2 <- (grid.A[, 1] >0.68 & grid.A[, 1] < 0.9 & grid.A[, 2] < 0.75 & 
             grid.A[, 2] > 0.5) & (grid.A[, 2] < (grid.A[, 1] - 0.1)) &
  (grid.A[, 2] > (grid.A[, 1] - 0.16))

ind_K3 <- (grid.A[, 1] >0.68 & grid.A[, 1] < 0.9 & grid.A[, 2] < 0.5 & 
             grid.A[, 2] > 0.25) & (grid.A[, 2] < (-grid.A[, 1] + 1.15)) &
  (grid.A[, 2] > (-grid.A[, 1] +1.1))

ind_K = ind_K1 | ind_K2 | ind_K3

obs_xy <- as.matrix(dt_A[obs_ind, ] %>% 
                      group_by(plot_id) %>% arrange(plot_id) %>%
                      summarize(x_B = mean(x_A), y_B = mean(y_A)) %>%  
                      mutate(intercept = 1) %>%
                      select(y_B, intercept, x_B))

HX = as.matrix(obs_xy[, -1]); y = obs_xy[, 1]
counts <- dt_A[obs_ind, ] %>% 
  group_by(plot_id) %>% arrange(plot_id) %>%
  summarize(counts = n()) %>% select(counts) %>% pull

dt_A_obs <- dt_A[obs_ind, ]
dt_A_obs <- dt_A_obs %>% arrange(plot_id)
dt_A_obs <- dt_A_obs %>% group_by(plot_id) %>%  
  mutate(n = 1, counts = sum(n)) %>% 
  mutate(weight = n / counts) 
Dh = dt_A_obs %>% group_by(plot_id) %>% # the sum_{i = 1}^{n_a} h^2_{li} 
  summarise(Dh = sum(weight^2)) %>% pull

pred_dta <- data.frame(y_A = y_A[ind_O | ind_K], w_A = w_A[ind_O | ind_K], 
                       x_A = X[ind_O | ind_K, 2], 
                       plot_id = c(ind_O*1 + ind_K*2)[ind_O | ind_K], 
                       coord.x = grid.A[ind_O | ind_K, 1],
                       coord.y = grid.A[ind_O | ind_K, 2])

pred_dta <- pred_dta %>% arrange(plot_id) %>% group_by(plot_id) %>%  
  mutate(n = 1, counts = sum(n)) %>% 
  mutate(weight = n / counts)

HXU <- cbind(1, pred_dta %>% group_by(plot_id) %>%  arrange(plot_id) %>%
               summarize(x_B = mean(x_A)) %>% select(x_B) %>% pull)
counts_u <-  pred_dta %>% group_by(plot_id) %>%  arrange(plot_id) %>%
  summarize(counts = n()) %>% select(counts) %>% pull

DhU = pred_dta %>% group_by(plot_id) %>%  arrange(plot_id) %>%  #the sum_{i = 1}^{n_a} h^2_{li}  for all prediction plots
  summarise(DhU = sum(weight^2)) %>% pull

## the coords of ALS variables 
grid.Aobs <- dt_A[obs_ind, ] %>% arrange(plot_id) %>% 
  select(coord.x, coord.y) # coordinates of predictors over regions with observation
plot_id_A <- dt_A[obs_ind, ] %>% arrange(plot_id) %>%
  select(plot_id) %>% pull
plotid_ind <- c(order(plot_id_A)[!duplicated(plot_id_A)], nrow(grid.Aobs) + 1) # record the first index of each plot
grid.Au <- pred_dta %>% arrange(plot_id) %>% select(coord.x, coord.y) # coordinates of predictors over regions for prediction
pred_id_A <- pred_dta %>% arrange(plot_id) %>% select(plot_id) %>% pull
predid_ind <- c(order(pred_id_A)[!duplicated(pred_id_A)], nrow(grid.Au) + 1) # record the first index of each polygon

na = nrow(grid.Aobs); nb = length(y); p = ncol(HX);

# compute and store the distance matrix #
gamma = 0.6
Dist_M <- rdist(grid.Aobs, grid.Aobs) 


## fit model in stan ##
library(cmdstanr)
library(bayesplot)
file <- file.path(getwd(), "project/stan_code/blowdown_save_RAM_weights_tapering_phi_unif.stan")
# file <- file.path(getwd(), "project/stan_code/blowdown_save_RAM_weights_tapering.stan")
# file <- file.path(getwd(), "project/stan_code/blowdown_save_RAM_weights.stan")
# file <- file.path(getwd(), "project/stan_code/blowdown_flat.stan")
# file <- file.path(getwd(), "project/stan_code/blowdown_save_RAM.stan")
mod <- cmdstan_model(file)

#-------------------------- Set parameters of priors --------------------------#
mu_beta = rep(0, p)     # mean vector in the Gaussian prior of beta
V_beta = diag(p) * 1000    # covariance matrix in the Gaussian prior of beta
## take precison matrix to be zero matrix
ss = 2       # scale parameter in the inverse gamma prior of sigma 
st = 2       # scale parameter in the inverse gamma prior of tau  
ap = 3/500; bp = 3/0.1       # lower and upper bound of uniform prior of phi 

data <- list(na = na, nb = nb, p = p, y = y, HX = HX, Dh = Dh, 
             gridA = grid.Aobs, hA = dt_A_obs$weight,
             plotid_ind = plotid_ind,
             mu_beta = mu_beta, V_beta = V_beta,
             ap = ap, bp = bp, ss = ss, st = st, Dist_M = Dist_M, 
             gamma = gamma)#,
#Dist_X = Dist_X)

fit <- mod$sample(
  data = data,
  seed = 1,
  chains = 4,
  parallel_chains = 4,
  refresh = 100,
  save_warmup = TRUE,
  iter_warmup = 500,
  iter_sampling = 500,
  sig_figs = 18
)

fit_draws <- fit$draws(inc_warmup = FALSE)
print(fit)

# prediction #
pick_sample_id <- seq(5, 2000, by = 10) # pick 200 samples
phi_ls <- c(fit_draws[, , "phi"])[pick_sample_id]
sigmasq_ls <- c(fit_draws[, , "sigmasq"])[pick_sample_id]
tausq_ls <- c(fit_draws[, , "tausq"])[pick_sample_id]
coords_A = grid.Aobs
coords_AU = grid.Au[, -1]
ind_ls_B = plotid_ind
ind_ls_BU = predid_ind
hA = dt_A_obs$weight
hAU = pred_dta$weight

beta_omega_sam <- sample_beta_omega_h_tapering(
  phi_ls, sigmasq_ls, tausq_ls, coords_A, coords_AU, hA, hAU, ind_ls_B, 
  ind_ls_BU, HX, mu_beta, V_beta, gamma, flat_prior = FALSE)

yU_ls <- pred_sample_y(beta_omega_sam$beta_ls, beta_omega_sam$omega_BU_ls, 
                       tausq_ls, HXU, DhU)

# check how many 95% posterior intervals covers the true value
incp_w <- beta_omega_sam$omega_B_ls + rep(1, nb) %*% t(beta_omega_sam$beta_ls[1, ])
qut <- apply(incp_w, 1, f <- function(x){quantile(x, c(0.025, 0.975))})
sum((qut[1,] < (w_B + 1)) & (qut[2,] > (w_B + 1))) / 55 # 98.2%


###### compare to BLOCK approach based on Andy's code #####
library(spBayes)
x <- c(HX[, 2])
n.samples <- 25000
cov.model <- "exponential"

starting.svi <- list("phi"=3/1, "sigma.sq"=sigma.sq, "tau.sq"=tau.sq)
tuning.svi <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)
priors.svi <- list("phi.Unif"=list(ap, bp),
                   "sigma.sq.IG"=list(2, ss),
                   "tau.sq.IG"=c(2, st),
                   "beta.norm"=list(mu_beta, V_beta))

m.1 <- spSVC(y~x, coords=plot.centroid[obs_ls, ], starting=starting.svi, 
             svc.cols=1, tuning=tuning.svi, priors=priors.svi, 
             cov.model=cov.model, n.samples=n.samples, n.omp.threads=1, 
             n.report=5000)

m.1 <- spRecover(m.1, start=floor(0.75*n.samples), thin=2, n.omp.threads=1, 
                 verbose = FALSE)

p.summary <- function(x){
  quantile(x, prob=c(0.5, 0.025, 0.975))
}

m.1.summary <- apply(rbind(t(apply(m.1$p.beta.recover.samples, 2, p.summary)),
                           t(apply(m.1$p.theta.recover.samples, 2, p.summary))),1,format)

HX2 <- as.matrix(dt_A[!obs_ind, ] %>% 
                   group_by(plot_id) %>% arrange(plot_id) %>%
                   summarize(x_B = mean(x_A)) %>%  
                   mutate(intercept = 1) %>%
                   select(intercept, x_B))
## Posterior predictive samples (m^3/ha), joint prediction for all blocks within a given blowdown prediction area.
out <- spPredict(m.1, pred.covars=HX2, pred.coords=plot.centroid[-obs_ls, ], 
                 n.omp.threads=1, joint=TRUE, verbose = FALSE, thin=15)

weight_id <- data.frame(plot_id = A_plot_id, coord.x = grid.A[, 1], coord.y = grid.A[, 2], 
                        pred_id = c(ind_O + ind_K*2), pred_I = c(ind_O + ind_K))
weight_id = weight_id %>% arrange(plot_id)
weight_id %>% glimpse()
weight_sum <- weight_id[!obs_ind, ] %>% arrange(plot_id) %>%
  group_by(plot_id) %>% 
  summarise(pred_n = sum(pred_I), pred_id = max(pred_id)) 
weight_sum %>% glimpse()
n1 = sum(weight_sum$pred_n * as.numeric(weight_sum$pred_id == 1))
n2 = sum(weight_sum$pred_n * as.numeric(weight_sum$pred_id == 2))

yU_O <- colSums(out$p.y.predictive.samples *  
                        (weight_sum$pred_n * as.numeric(weight_sum$pred_id == 1))) / n1

yU_K <- colSums(out$p.y.predictive.samples * 
                        (weight_sum$pred_n * as.numeric(weight_sum$pred_id == 2))) / n2

# save output
output_filename <- paste0("./results/sim/sim_", input_id, ".Rdata")
y_O_true = mean(y_A[ind_O])
y_K_true = mean(y_A[ind_K])

save(yU_ls, yU_O, yU_K, y_O_true, y_K_true, file = output_filename)



# ### compare to benchmark 
# file2 <- file.path(getwd(), "project/stan_code/blowdown_benchmark.stan")
# mod2 <- cmdstan_model(file2)
# 
# data2 <- list(n = nb, p = p, y = y, X = HX, 
#               gridA = plot.centroid[obs_ls, ], 
#               mu_beta = mu_beta, V_beta = V_beta,
#               ap = ap, bp = bp, ss = ss, st = st)
# 
# fit2 <- mod2$sample(
#   data = data2,
#   seed = 1,
#   chains = 4,
#   parallel_chains = 4,
#   refresh = 100,
#   save_warmup = TRUE,
#   iter_warmup = 500,
#   iter_sampling = 500,
#   sig_figs = 18
# )
# 
# fit2_draws <- fit2$draws(inc_warmup = FALSE)
# print(fit2)
# 
# pick_sample_id <- seq(5, 2000, by = 10) # pick 200 samples
# phi_ls2 <- c(fit2_draws[, , "phi"])[pick_sample_id]
# sigmasq_ls2 <- c(fit2_draws[, , "sigmasq"])[pick_sample_id]
# tausq_ls2 <- c(fit2_draws[, , "tausq"])[pick_sample_id]
# coords_A2 = plot.centroid[obs_ls, ]
# coords_AU2 = plot.centroid[-obs_ls, ]
# ind_ls_B2 = 1:(nrow(coords_A2) + 1)
# ind_ls_BU2 = 1:(nrow(coords_AU2) + 1)
# hA2 <- rep(1.0, nrow(coords_A2))
# hAU2 <- rep(1.0, nrow(coords_AU2))
# 
# beta_omega_sam2 <- sample_beta_omega_h(phi_ls2, sigmasq_ls2, tausq_ls2,
#                                        coords_A2, coords_AU2, hA2, hAU2, 
#                                        ind_ls_B2, ind_ls_BU2,
#                                        HX, mu_beta, V_beta, flat_prior = FALSE)
# HX2 <- as.matrix(dt_A[!obs_ind, ] %>% 
#                    group_by(plot_id) %>% arrange(plot_id) %>%
#                    summarize(x_B = mean(x_A)) %>%  
#                    mutate(intercept = 1) %>%
#                    select(intercept, x_B))
# DhU2 <- rep(1, nrow(HX2))
# 
# yU_ls2 <- pred_sample_y(beta_omega_sam2$beta_ls, beta_omega_sam2$omega_BU_ls, 
#                         tausq_ls2, HX2, DhU2)
# 
# # check result
# y_B_U <- dt_A[!obs_ind, ] %>% group_by(plot_id) %>% 
#   summarize(y_B_U = mean(y_A)) %>% arrange(plot_id) %>% select(y_B_U) %>% pull 
# 
# weight_id <- data.frame(plot_id = A_plot_id, coord.x = grid.A[, 1], coord.y = grid.A[, 2], 
#                         pred_id = c(ind_O + ind_K*2), pred_I = c(ind_O + ind_K))
# weight_id = weight_id %>% arrange(plot_id)
# weight_sum <- weight_id[!obs_ind, ] %>% arrange(plot_id) %>%
#   group_by(plot_id) %>% 
#   summarise(pred_n = sum(pred_I), pred_id = max(pred_id)) 
# weight_sum %>% glimpse()
# n1 = sum(weight_sum$pred_n * as.numeric(weight_sum$pred_id == 1))
# n2 = sum(weight_sum$pred_n * as.numeric(weight_sum$pred_id == 2))
# 
# yU_O <- colSums(yU_ls2 * (weight_sum$pred_n * 
#                             as.numeric(weight_sum$pred_id == 1))) / n1
# yU_K <- colSums(yU_ls2 * (weight_sum$pred_n * 
#                             as.numeric(weight_sum$pred_id == 2))) / n2
# 
# # save output
# output_filename <- paste0("./results/sim/sim_", input_id, ".Rdata")
# y_O_true = mean(y_A[ind_O])
# y_K_true = mean(y_A[ind_K])
# 
# save(yU_ls, yU_O, yU_K, y_O_true, y_K_true, file = output_filename)
# 
