rm(list = ls())
library(dplyr)
library(Matrix)
library(fields)
library(ggplot2)
source("./project/utils.R")

## read data ##
# we only fit data in Frohn in this example #
plot.y <- read.csv("./data/outcome_y_plot.csv")
plot.y %>% glimpse
plot.y <- plot.y %>% filter(sub.region == "Liesing") #"Frohn") 
plot.y %>% glimpse

plot.grid.x <- read.csv("./data/predictor_x_grid.csv")
plot.grid.x %>% glimpse
plot.grid.x <- plot.grid.x %>% filter(obs.plot.id >=45 & obs.plot.id<=51) #<= 21)
plot.grid.x %>% glimpse


pred.grid.x <- read.csv("./data/prediction_x_grid.csv")
pred.grid.x %>% glimpse
pred.grid.x <- pred.grid.x %>% filter(sub.region == "Liesing") #"Frohn")
# we only predict on three small polygon in this example
# predid_ind <- c(order(pred.grid.x$pred.poly.id)[     # record the first index of each polygon
#   !duplicated(pred.grid.x$pred.poly.id)], nrow(pred.grid.x) + 1)
# pred.grid.x <- pred.grid.x[c(1:(predid_ind[2]-1), 
#                              (predid_ind[3]:(predid_ind[5] - 1))), ]
# pred.grid.x %>% glimpse


## Compute the average height.m over each plot and summary with observation in one dataset
x_obs_bar <- plot.grid.x %>% group_by(obs.plot.id) %>% 
  summarize(height.mean = mean(height.m)) %>% select(obs.plot.id, height.mean) 
obs_xy <- plot.y %>% select(obs.plot.id, sub.region, volume.m3.ha) %>% 
  inner_join(x_obs_bar, by = "obs.plot.id") %>% 
  mutate(intercept = 1) %>%
  # mutate(x.Frohn = height.mean * (sub.region == "Frohn")) %>%
  # mutate(x.Laas = height.mean * (sub.region == "Laas")) %>%
  # mutate(x.Mauthen = height.mean * (sub.region == "Mauthen")) %>%
  mutate(x.Liesing = height.mean * (sub.region == "Liesing")) %>%
  # mutate(x.Ploecken = height.mean * (sub.region == "Ploecken")) %>% 
  # select(volume.m3.ha, intercept, x.Frohn, x.Laas, x.Mauthen, x.Liesing, x.Ploecken)
  select(volume.m3.ha, intercept, x.Liesing)#x.Frohn)

HX = as.matrix(obs_xy[, -1]); y = obs_xy[, 1]
counts <- plot.grid.x %>% group_by(obs.plot.id) %>% summarize(counts = n()) %>%
  select(counts) %>% pull
Dh = 1 / counts

HXU <- cbind(1, pred.grid.x %>% group_by(pred.poly.id) %>% 
               summarize(height.mean = mean(height.m)) %>% select(height.mean) %>% pull)
counts_u <- pred.grid.x %>% group_by(pred.poly.id) %>% summarize(counts = n()) %>%
  select(counts) %>% pull
DhU <- 1 / counts_u

## the coords of ALS variables 
grid.A <- plot.grid.x %>% # select predictors in region Frohn
  select(coord.x, coord.y) # coordinates of predictors over regions with observation
plotid_ind <- c(order(plot.grid.x$obs.plot.id)[      # record the first index of each plot
  !duplicated(plot.grid.x$obs.plot.id)], nrow(grid.A) + 1)       
grid.Au <- pred.grid.x %>% select(coords.x, coords.y) # coordinates of predictors over regions for prediction
predid_ind <- c(order(pred.grid.x$pred.poly.id)[     # record the first index of each polygon
  !duplicated(pred.grid.x$pred.poly.id)], nrow(grid.Au) + 1)

na = nrow(grid.A); nb = length(y); p = ncol(HX);

## visualize the data ##
p1 <- ggplot(plot.y, aes(x = coords.x, y = coords.y, colour = volume.m3.ha)) +
  geom_point()
p1

p2 <- ggplot(plot.grid.x, aes(x = coord.x, y = coord.y, colour = height.m)) + 
  geom_point(size = 0.01) + 
  geom_point(data = plot.y, aes(x=coords.x, y=coords.y), size=3, shape=1,
             color="orange") +
  geom_point(data = pred.grid.x, aes(x = coords.x, y = coords.y), size = 0.01)
p2


## fit model in stan ##
library(cmdstanr)
file <- file.path(getwd(), "project/blowdown_flat.stan")
#file <- file.path(getwd(), "project/blowdown_save_RAM.stan")
mod <- cmdstan_model(file)

#-------------------------- Set parameters of priors --------------------------#
mu_beta = rep(0, p)     # mean vector in the Gaussian prior of beta
V_beta = diag(p) * 1000    # covariance matrix in the Gaussian prior of beta
ss = 3 * sqrt(2)       # scale parameter in the normal prior of sigma 
st = 3 * sqrt(2)     # scale parameter in the normal prior of tau     
ap = 3; bp = 0.5       # shape and rate parameters in the Gamma prior of phi 

data <- list(na = na, nb = nb, p = p, y = y, HX = HX, Dh = Dh, 
             gridA = grid.A / 1000,
             plotid_ind = plotid_ind,
             mu_beta = mu_beta, V_beta = V_beta,
             ap = ap, bp = bp, ss = ss, st = st)

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
n_lf <- fit$sampler_diagnostics(inc_warmup = TRUE)[1:6, 1, "n_leapfrog__"]
sum(n_lf)

opt <- mod$optimize(data=data, algorithm='lbfgs', seed = 321)
#' Check whether parameters have reasonable values
opt$draws()
# around 17.104s per evaluation
# A draws_matrix: 1 draws, and 6 variables
# variable
# draw lp__ sigma    tau phi sigmasq   tausq
# 1 -343    56 0.0056  15    3083 3.1e-05

# compute the target log-density#
phi = 5
t <- proc.time()
C_B <- Block_COV(grid.A, plotid_ind, phi)
proc.time()-t
ll <- log_lik(sigmasq = 1, tausq = 1, y, HX, C_B, mu_beta, V_beta, Dh)
beta_omega_sam <- sample_beta_omega_no_pred(y, HX, C_B, sigmasq = 1, tausq = 1, 
                                            mu_beta, V_beta, Dh)

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



