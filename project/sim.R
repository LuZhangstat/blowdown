rm(list = ls())
library(dplyr)
library(Matrix)
library(fields)
library(ggplot2)
source("./project/utils.R")

## simulate data ##
set.seed(1234)
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
plot.centroid = expand.grid(c(1/16, 3/16, 5/16, 7/16, 9/16, 11/16, 13/16, 15/16), 
                            c(1/16, 3/16, 5/16, 7/16, 9/16, 11/16, 13/16, 15/16)) # center of 9 plots
D_plot <- rdist(grid.A, plot.centroid) 
A_plot_id <- (D_plot < 0.05) %*% c(1:64) # plot id on the finest grid
# A_plot_id <- 1:N_A
dt_A <- data.frame(y_A = y_A, w_A = w_A, x_A = X[, 2], plot_id = A_plot_id, 
                   coord.x = grid.A[, 1], coord.y = grid.A[, 2])
dt_A = dt_A %>% arrange(plot_id)
w_B = dt_A %>% filter(plot_id > 0, plot_id < 64) %>% group_by(plot_id) %>% 
  summarize(w_B = mean(w_A)) %>% arrange(plot_id) %>% select(w_B) %>% pull 

## plot the pattern of the source and target units ##
plot(grid.A[, 1], grid.A[, 2])
points(plot.centroid[, 1], plot.centroid[, 2], col = "red")
points(grid.A[which(A_plot_id > 0 & A_plot_id < 64), 1], 
       grid.A[which(A_plot_id > 0 & A_plot_id < 64), 2], col = "orange", cex = 2)
points(grid.A[which(A_plot_id == 64), 1], 
       grid.A[which(A_plot_id == 64), 2], col = "blue", cex = 2)

dt_B_w <- data.frame(w_B = w_B, coords.x = plot.centroid[-64, 1],
                     coords.y = plot.centroid[-64, 2])
p1 <- ggplot(dt_B_w, aes(x = coords.x, y = coords.y, colour = w_B)) +
  geom_point()
p1


obs_xy <- as.matrix(dt_A %>% filter(plot_id > 0, plot_id < 64) %>% 
                      group_by(plot_id) %>% 
                      summarize(x_B = mean(x_A), y_B = mean(y_A)) %>%  
                      mutate(intercept = 1) %>%
                      select(y_B, intercept, x_B))

HX = as.matrix(obs_xy[, -1]); y = obs_xy[, 1]
counts <- dt_A %>% filter(plot_id > 0, plot_id < 64) %>% 
  group_by(plot_id) %>%
  summarize(counts = n()) %>% select(counts) %>% pull
Dh = 1 / counts

HXU <- cbind(1, dt_A %>% filter(plot_id == 64) %>% group_by(plot_id) %>% 
               summarize(x_B = mean(x_A)) %>% select(x_B) %>% pull)
counts_u <- dt_A %>% filter(plot_id == 64) %>% group_by(plot_id) %>% 
  summarize(counts = n()) %>% select(counts) %>% pull
DhU <- 1 / counts_u

## the coords of ALS variables 
grid.A <- dt_A %>% filter(plot_id > 0, plot_id < 64) %>%# select predictors in region Frohn
  select(coord.x, coord.y) # coordinates of predictors over regions with observation
plot_id_A <- dt_A %>% filter(plot_id > 0, plot_id < 64) %>% 
  select(plot_id) %>% pull
plotid_ind <- c(order(plot_id_A)[!duplicated(plot_id_A)], nrow(grid.A) + 1) # record the first index of each plot
grid.Au <- dt_A %>% filter(plot_id == 64) %>% select(coord.x, coord.y) # coordinates of predictors over regions for prediction
pred_id_A <- dt_A %>% filter(plot_id == 64) %>% select(plot_id) %>% pull
predid_ind <- c(order(pred_id_A)[!duplicated(pred_id_A)], nrow(grid.Au) + 1) # record the first index of each polygon

na = nrow(grid.A); nb = length(y); p = ncol(HX);

#C_B <- Block_COV(grid.A, plotid_ind, phi)


## fit model in stan ##
library(cmdstanr)
library(bayesplot)
file <- file.path(getwd(), "project/blowdown_flat.stan")
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
             gridA = grid.A,
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
  iter_warmup = 200,
  iter_sampling = 200,
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
#' Check whether parameters have reasonable values
opt$draws()
odraws <- opt$draws()

# prediction #
phi_ls <- c(fit_draws[, , "phi"])
sigmasq_ls <- c(fit_draws[, , "sigmasq"])
tausq_ls <- c(fit_draws[, , "tausq"])
coords_A = grid.A 
coords_AU = grid.Au
ind_ls_B = plotid_ind
ind_ls_BU = predid_ind
beta_omega_sam <- sample_beta_omega(phi_ls, sigmasq_ls, tausq_ls,
                                    coords_A, coords_AU, ind_ls_B, ind_ls_BU,
                                    HX, mu_beta, V_beta, flat_prior = TRUE)

yU_ls <- pred_sample_y(beta_omega_sam$beta_ls, beta_omega_sam$omega_BU_ls, 
                       tausq_ls, HXU, DhU)

summary(t(beta_omega_sam$beta_ls))
summary(t(beta_omega_sam$omega_B_ls))
plot(rowMeans(beta_omega_sam$omega_B_ls), w_B)
abline(a = 0, b = 1)

plot(rowMeans(beta_omega_sam$omega_B_ls + 
                rep(1, nb) %*% t(beta_omega_sam$beta_ls[1, ])), w_B + 1)
abline(a = 0, b = 1)

# check how many 95% posterior intervals covers the true value
incp_w <- beta_omega_sam$omega_B_ls + rep(1, nb) %*% t(beta_omega_sam$beta_ls[1, ])
qut <- apply(incp_w, 1, f <- function(x){quantile(x, c(0.025, 0.975))})
sum((qut[1,] < (w_B + 1)) & (qut[2,] > (w_B + 1))) / 63
