rm(list = ls())
library(dplyr)
library(Matrix)
library(fields)
library(ggplot2)
source("./project/utils.R")

## simulate data ##
set.seed(4) #1
side.A = seq(from = 1/72, to = 1-1/72, by = 1/36)
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
# plot.centroid = expand.grid(c(1/16, 3/16, 5/16, 7/16, 9/16, 11/16, 13/16, 15/16), 
#                             c(1/16, 3/16, 5/16, 7/16, 9/16, 11/16, 13/16, 15/16)) # center of 9 plots

plot.centroid =  expand.grid(side.A[seq(3, 36, by = 4)], side.A[seq(3, 36, by = 4)])
D_plot <- rdist(grid.A, plot.centroid) 
A_plot_id <- (D_plot < 0.045) %*% c(1:81) # plot id on the finest grid
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
O_d1 <- rdist(grid.A, t(c(0.3, 0.65))) 
O_d2 <- rdist(grid.A, t(c(0.3, 0.35))) 
ind_O <- ((O_d1 + O_d2) < 0.44 & (O_d1 + O_d2) > 0.39)

ind_K1 <- (grid.A[, 1] >0.6 & grid.A[, 1] < 0.68 & grid.A[, 2] < 0.72 & 
            grid.A[, 2] > 0.28)

ind_K2 <- (grid.A[, 1] >0.68 & grid.A[, 1] < 0.9 & grid.A[, 2] < 0.72 & 
             grid.A[, 2] > 0.5) & (grid.A[, 2] < (grid.A[, 1] - 0.1)) &
  (grid.A[, 2] > (grid.A[, 1] - 0.16))

ind_K3 <- (grid.A[, 1] >0.68 & grid.A[, 1] < 0.9 & grid.A[, 2] < 0.5 & 
             grid.A[, 2] > 0.28) & (grid.A[, 2] < (-grid.A[, 1] + 1.15)) &
  (grid.A[, 2] > (-grid.A[, 1] +1.1))

ind_K = ind_K1 | ind_K2 | ind_K3

## plot the pattern of the source and target units ##
plot(grid.A[, 1], grid.A[, 2])
points(dt_A[obs_ind, "coord.x"], 
       dt_A[obs_ind, "coord.y"], col = "orange", cex = 2, pch = 16)
#points(plot.centroid[obs_ls, 1], plot.centroid[obs_ls, 2], col = "red")
#text(plot.centroid[, 1], plot.centroid[, 2], labels = c(1:81))
#text(grid.A[, 1], grid.A[, 2], labels = c(1:1296))
points(grid.A[c(which(ind_O), which(ind_K)), 1], 
       grid.A[c(which(ind_O), which(ind_K)), 2], col = "blue", cex = 1, pch = 16)
# orange: observed regions. blue: prediction plots

## plot the observed data and the centroid of the observed plots
dt_B_w <- data.frame(w_B = w_B, coords.x = plot.centroid[obs_ls, 1],
                     coords.y = plot.centroid[obs_ls, 2])
p1 <- ggplot(dt_B_w, aes(x = coords.x, y = coords.y, colour = w_B)) +
  geom_point()
p1


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

#C_B <- Block_COV(grid.A, plotid_ind, phi)


## fit model in stan ##
library(cmdstanr)
library(bayesplot)
file <- file.path(getwd(), "project/blowdown_save_RAM_weights.stan")
# file <- file.path(getwd(), "project/blowdown_flat.stan")
# file <- file.path(getwd(), "project/blowdown_save_RAM.stan")
mod <- cmdstan_model(file)

#-------------------------- Set parameters of priors --------------------------#
mu_beta = rep(0, p)     # mean vector in the Gaussian prior of beta
V_beta = diag(p) * 1000    # covariance matrix in the Gaussian prior of beta
## take precison matrix to be zero matrix
ss = 3 * sqrt(2)       # scale parameter in the normal prior of sigma 
st = 3 * sqrt(2)     # scale parameter in the normal prior of tau     
ap = 3; bp = 0.5       # shape and rate parameters in the Gamma prior of phi 

data <- list(na = na, nb = nb, p = p, y = y, HX = HX, Dh = Dh, 
             gridA = grid.Aobs, hA = dt_A_obs$weight,
             plotid_ind = plotid_ind,
             mu_beta = mu_beta, V_beta = V_beta,
             ap = ap, bp = bp, ss = ss, st = st)#,
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
fit$summary()
mcmc_trace(fit$draws("sigmasq"), iter1 = 1) 
mcmc_trace(fit$draws("tausq"), iter1 = 1) 
mcmc_trace(fit$draws("phi"), iter1 = 1) 
n_lf <- fit$sampler_diagnostics(inc_warmup = TRUE)[, , "n_leapfrog__"]
colSums(n_lf)



opt <- mod$optimize(data=data, algorithm='lbfgs', seed = 1,
                    init = function() list(phi = 5.0, sigmasq = 2.0, tausq = 1.0))
#' Check MAP estimates 
opt$draws()
odraws <- opt$draws()

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

beta_omega_sam <- sample_beta_omega_h(phi_ls, sigmasq_ls, tausq_ls,
                                    coords_A, coords_AU, hA, hAU, 
                                    ind_ls_B, ind_ls_BU,
                                    HX, mu_beta, V_beta, flat_prior = FALSE)

yU_ls <- pred_sample_y(beta_omega_sam$beta_ls, beta_omega_sam$omega_BU_ls, 
                       tausq_ls, HXU, DhU)

# check how many 95% posterior intervals covers the true value
incp_w <- beta_omega_sam$omega_B_ls + rep(1, nb) %*% t(beta_omega_sam$beta_ls[1, ])
qut <- apply(incp_w, 1, f <- function(x){quantile(x, c(0.025, 0.975))})
sum((qut[1,] < (w_B + 1)) & (qut[2,] > (w_B + 1))) / 55 # 98.2%

t1 <- rowMeans(beta_omega_sam$omega_B_ls) + rowMeans(beta_omega_sam$beta_ls)[1]
t2 <- dt_A_obs %>% arrange(plot_id) %>% group_by(plot_id) %>%
  summarise(w_B = mean(w_A)) %>% pull
plot(t1, t2 + 1)
abline(a = 0, b = 1)

# check the prediction plots
# O:
hist(yU_ls[1, ])
mean(yU_ls[1, ])
quantile(yU_ls[1, ], c(0.025, 0.975))
mean(y_A[ind_O])


# K:
hist(yU_ls[2, ])
mean(yU_ls[2, ])
quantile(yU_ls[2, ], c(0.025, 0.975))
mean(y_A[ind_K])


