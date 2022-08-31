## test on generating data on coarser grid
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
plot.y <- plot.y %>% filter(sub.region == "Frohn") 
plot.y %>% glimpse

plot.grid.x <- read.csv("./data/predictor_x_grid.csv")
plot.grid.x %>% glimpse
plot.grid.x <- plot.grid.x %>% filter(obs.plot.id <= 21)
plot.grid.x %>% glimpse


pred.grid.x <- read.csv("./data/prediction_x_grid.csv")
pred.grid.x <- pred.grid.x %>% filter(sub.region == "Frohn")
pred.grid.x %>% glimpse

# generate data on a coarser grid
# create the shared location for every 9 cells
plot.grid.x$coord.x2 <- (plot.grid.x$coord.x + 1 * (plot.grid.x$coord.x %% 3 == 1) - 
                      1 * (plot.grid.x$coord.x %% 3 == 0))
plot.grid.x$coord.y2 <- (plot.grid.x$coord.y + 1 * (plot.grid.x$coord.y %% 3 == 1) - 
                      1 * (plot.grid.x$coord.y %% 3 == 0))

plot.grid.x_lit <- plot.grid.x %>% group_by(obs.plot.id, coord.x2, coord.y2) %>% 
  summarise(height.m = mean(height.m), n = n()) %>%
  select(obs.plot.id, coord.x2, coord.y2, height.m, n) 
plot.grid.x_lit %>% glimpse()


pred.grid.x $coords.x2 <- (pred.grid.x $coords.x + 
                             1 * (pred.grid.x $coords.x %% 3 == 1) - 
                             1 * (pred.grid.x $coords.x %% 3 == 0))
pred.grid.x $coords.y2 <- (pred.grid.x $coords.y + 
                             1 * (pred.grid.x $coords.y %% 3 == 1) - 
                             1 * (pred.grid.x $coords.y %% 3 == 0))
pred.grid.x_lit <- pred.grid.x %>% group_by(pred.poly.id, coords.x2, coords.y2) %>% 
  summarise(height.m = mean(height.m), n = n()) %>%
  select(pred.poly.id, coords.x2, coords.y2, height.m, n) 
pred.grid.x_lit %>% glimpse()


# check grids for observation
test_1 <- plot.grid.x %>% filter(obs.plot.id == 1)
plot(test_1$coord.x, test_1$coord.y, 
     xlim = c(404904, 404943), ylim = c(168376, 168415))

test_1_lit <- plot.grid.x_lit %>% filter(obs.plot.id == 1) 

plot(test_1_lit$coord.x2, test_1_lit$coord.y2, type="n",
     xlim = c(404904, 404943), ylim = c(168376, 168415))
text(test_1_lit$coord.x2, test_1_lit$coord.y2, test_1_lit$n, pos=1, offset = 0.0)


# check grids for prediction 
pre_1 <- pred.grid.x %>% filter(pred.poly.id == 58)
plot(pre_1$coords.x, pre_1$coords.y, 
     xlim = c(min(pre_1$coords.x)-1, max(pre_1$coords.x) + 1), 
     ylim = c(min(pre_1$coords.y)-1, max(pre_1$coords.y) + 1))

pre_1_lit <- pred.grid.x_lit %>% filter(pred.poly.id == 58)

plot(pre_1_lit$coords.x2, pre_1_lit$coords.y2, type="n",
     xlim = c(min(pre_1$coords.x)-1, max(pre_1$coords.x) + 1), 
     ylim = c(min(pre_1$coords.y)-1, max(pre_1$coords.y) + 1))
text(pre_1_lit$coords.x2, pre_1_lit$coords.y2, pre_1_lit$n, pos=1, offset = 0.0)

## Compute the average height.m over each plot and summary with observation in one dataset
x_obs_bar <- plot.grid.x_lit %>% group_by(obs.plot.id) %>% 
  summarize(height.mean = sum(height.m * n) / sum(n)) %>% select(obs.plot.id, height.mean) 
obs_xy <- plot.y %>% select(obs.plot.id, sub.region, volume.m3.ha) %>% 
  inner_join(x_obs_bar, by = "obs.plot.id") %>% 
  mutate(intercept = 1) %>%
  mutate(x.Frohn = height.mean * (sub.region == "Frohn")) %>%
  # mutate(x.Laas = height.mean * (sub.region == "Laas")) %>%
  # mutate(x.Mauthen = height.mean * (sub.region == "Mauthen")) %>%
  # mutate(x.Liesing = height.mean * (sub.region == "Liesing")) %>%
  # mutate(x.Ploecken = height.mean * (sub.region == "Ploecken")) %>% 
  # select(volume.m3.ha, intercept, x.Frohn, x.Laas, x.Mauthen, x.Liesing, x.Ploecken)
  select(volume.m3.ha, intercept, x.Frohn)

HX = as.matrix(obs_xy[, -1]); y = obs_xy[, 1]
plot.grid.x_lit <- plot.grid.x_lit %>% group_by(obs.plot.id) %>% 
  mutate(counts = sum(n)) %>% mutate(weight = n / counts)
plot.grid.x_lit %>% glimpse()


Dh = plot.grid.x_lit %>% group_by(obs.plot.id) %>% 
  summarise(Dh = sum(weight^2)) %>% pull


HXU <- cbind(1, pred.grid.x_lit %>% group_by(pred.poly.id) %>% 
               summarize(height.mean = sum(height.m * n) / sum(n)) %>% 
               select(height.mean) %>% pull)
pred.grid.x_lit <- pred.grid.x_lit %>% group_by(pred.poly.id) %>% 
  mutate(counts = sum(n)) %>% mutate(weight = n / counts)
pred.grid.x_lit %>% glimpse()
DhU = pred.grid.x_lit %>% group_by(pred.poly.id) %>% 
  summarise(DhU = sum(weight^2)) %>% pull

## the coords of ALS variables 
grid.A <- plot.grid.x_lit %>% # select predictors in region Frohn
  select(coord.x2, coord.y2) # coordinates of predictors over regions with observation
plotid_ind <- c(order(plot.grid.x_lit$obs.plot.id)[      # record the first index of each plot
  !duplicated(plot.grid.x_lit$obs.plot.id)], nrow(grid.A) + 1)       
grid.Au <- pred.grid.x_lit %>% select(coords.x2, coords.y2) # coordinates of predictors over regions for prediction
predid_ind <- c(order(pred.grid.x_lit$pred.poly.id)[     # record the first index of each polygon
  !duplicated(pred.grid.x_lit$pred.poly.id)], nrow(grid.Au) + 1)

na = nrow(grid.A); nb = length(y); p = ncol(HX);

## visualize the data ##
p1 <- ggplot(plot.y, aes(x = coords.x, y = coords.y, colour = volume.m3.ha)) +
  geom_point()
p1

p2 <- ggplot(plot.grid.x_lit, aes(x = coord.x2, y = coord.y2)) + 
  geom_point(size = 0.01) + 
  geom_point(data = pred.grid.x_lit, aes(x = coords.x2, y = coords.y2, 
                                         colour = height.m), size = 0.01)
p2


## fit model in stan ##
library(cmdstanr)
library(bayesplot)
file <- file.path(getwd(), "project/blowdown_save_RAM_weights.stan")
mod <- cmdstan_model(file)

#-------------------------- Set parameters of priors --------------------------#
mu_beta = rep(0, p)     # mean vector in the Gaussian prior of beta
V_beta = diag(p) * 1000    # covariance matrix in the Gaussian prior of beta
ss = 3 * sqrt(2)       # scale parameter in the normal prior of sigma 
st = 3 * sqrt(2)     # scale parameter in the normal prior of tau     
ap = 3; bp = 0.5       # shape and rate parameters in the Gamma prior of phi 

data <- list(na = na, nb = nb, p = p, y = y, HX = HX, Dh = Dh, 
             gridA = grid.A[, -1] / 1000, hA = plot.grid.x_lit$weight,
             plotid_ind = plotid_ind,
             mu_beta = mu_beta, V_beta = V_beta,
             ap = ap, bp = bp, ss = ss, st = st)

# fit <- mod$sample(
#   data = data,
#   seed = 1,
#   chains = 4,
#   parallel_chains = 4,
#   refresh = 10,
#   save_warmup = TRUE,
#   iter_warmup = 500,
#   iter_sampling = 500,
#   sig_figs = 18
# )
# 
# fit$summary()
# n_lf <- fit$sampler_diagnostics(inc_warmup = TRUE)[, 1, "n_leapfrog__"]
# sum(n_lf)
# 
# fit$save_object(file = "./results/fit.RDS")

fit <- readRDS("./results/fit.RDS")

fit$summary()
mcmc_trace(fit$draws("sigmasq"), iter1 = 1) 
mcmc_trace(fit$draws("tausq"), iter1 = 1) 
mcmc_trace(fit$draws("phi"), iter1 = 1) 
n_lf <- fit$sampler_diagnostics(inc_warmup = TRUE)[, , "n_leapfrog__"]
colSums(n_lf)
fit_draws <- fit$draws(inc_warmup = FALSE)

# prediction #
pick_sample_id <- seq(5, 2000, by = 10) # pick 200 samples
phi_ls <- c(fit_draws[, , "phi"])[pick_sample_id]
sigmasq_ls <- c(fit_draws[, , "sigmasq"])[pick_sample_id]
tausq_ls <- c(fit_draws[, , "tausq"])[pick_sample_id]
coords_A = grid.A[, -1] /1000
coords_AU = grid.Au[, -1] / 1000
ind_ls_B = plotid_ind
ind_ls_BU = predid_ind
hA = plot.grid.x_lit$weight
hAU = pred.grid.x_lit$weight

# t <- proc.time()
# beta_omega_sam <- sample_beta_omega_h(phi_ls, sigmasq_ls, tausq_ls,
#                                     coords_A, coords_AU, 
#                                     hA, hAU, ind_ls_B, ind_ls_BU,
#                                     HX, mu_beta, V_beta, flat_prior = FALSE)
# proc.time() -t
# save(beta_omega_sam, file = "./results/RDA_recov_sam.RData")
load("./results/RDA_recov_sam.RData")

yU_ls <- pred_sample_y(beta_omega_sam$beta_ls, beta_omega_sam$omega_BU_ls, 
                       tausq_ls, HXU, DhU)

yU_pre <-  data.frame(pred.poly.id = 1:nrow(yU_ls), yU_pm = rowMeans(yU_ls),
                      omega_BU_pm = colMeans(t(beta_omega_sam$omega_BU_ls) + 
                                               beta_omega_sam$beta_ls[1,]))
yB_pre <- data.frame(obs.plot.id = 1:nrow(plot.y), 
                     omega_B_pm = colMeans(t(beta_omega_sam$omega_B_ls) + 
                                             beta_omega_sam$beta_ls[1,]))



## plots ##
obs_plot_dat <- plot.grid.x_lit %>% merge(plot.y[c("obs.plot.id", "volume.m3.ha")], 
                                by = "obs.plot.id", no.dups = FALSE) %>% 
  merge(yB_pre, by = "obs.plot.id")%>%
  select(obs.plot.id, coord.x2, coord.y2, volume.m3.ha, omega_B_pm)
pred_plot_dat <- pred.grid.x_lit %>% merge(yU_pre, by = "pred.poly.id") %>%
  select(pred.poly.id, coords.x2, coords.y2, yU_pm, omega_BU_pm)

library(plyr)
df <- pred_plot_dat
find_hull <- function(df) df[chull(df$coords.x2, df$coords.y2), ]
hulls <- ddply(df, "pred.poly.id", find_hull)


p_o <- ggplot(obs_plot_dat, aes(x = coord.x2, y = coord.y2)) + 
  # geom_polygon(data = hulls, 
  #              aes(x = coords.x2, y = coords.y2))+
  geom_point(size = 0.01, aes(colour = volume.m3.ha)) + xlim(404395, 406938) +
  ylim(168259, 172896) +
  scale_colour_gradient(limits = range(obs_plot_dat$volume.m3.ha,
                                    pred_plot_dat$yU_pm)) 
p_o

p_py <-  ggplot(hulls, aes(x = coords.x2, y = coords.y2)) + 
  #geom_point(size = 0.01) +
  geom_polygon(aes(fill = yU_pm,group = pred.poly.id), col = "white", size = 0.3)+
  geom_point(data = obs_plot_dat,
             size = 0.1, 
             aes(x = coord.x2, y = coord.y2,
                              fill = volume.m3.ha, colour = volume.m3.ha)) +
  geom_point(data = plot.y, aes(x=coords.x, y=coords.y), size=3, shape=1,
             color="orange") +
  scale_fill_gradient(limits = range(obs_plot_dat$volume.m3.ha,
                                       pred_plot_dat$yU_pm))+
  xlim(404395, 406938) + ylim(168259, 172896) 
  #xlim(404395, 404500) + ylim(168259, 168350) 
p_py


p_oW <- ggplot(obs_plot_dat, aes(x = coord.x2, y = coord.y2, colour = omega_B_pm)) + 
  geom_point(size = 0.01) + xlim(404395, 406938) + ylim(168259, 172896) +
  scale_colour_gradient(limits = range(obs_plot_dat$omega_B_pm, 
                                       pred_plot_dat$omega_BU_pm))
p_oW

p_pW <-  ggplot(obs_plot_dat, aes(x = coord.x2, y = coord.y2, colour = omega_B_pm)) + 
  geom_point(size = 0.01) + 
  geom_point(data = plot.y, aes(x=coords.x, y=coords.y), size=3, shape=1, 
             color="orange") +
  geom_point(data = pred_plot_dat, aes(x = coords.x2, y = coords.y2, 
                                       colour = omega_BU_pm), size = 0.01) +
  xlim(404395, 406938) + ylim(168259, 172896) +
  scale_colour_gradient(limits = range(obs_plot_dat$omega_B_pm, 
                                       pred_plot_dat$omega_BU_pm))
p_pW



p_pW <-  ggplot(hulls, aes(x = coords.x2, y = coords.y2)) + 
  #geom_point(size = 0.01) +
  geom_polygon(aes(fill = omega_BU_pm, group = pred.poly.id), 
               col = "white", size = 0.3)+
  geom_point(data = obs_plot_dat,
             size = 0.1, 
             aes(x = coord.x2, y = coord.y2,
                 fill = omega_B_pm, colour = omega_B_pm)) +
  geom_point(data = plot.y, aes(x=coords.x, y=coords.y), size=3, shape=1,
             color="orange") +
  scale_colour_gradient(limits = range(obs_plot_dat$omega_B_pm, 
                                       pred_plot_dat$omega_BU_pm)) +
  xlim(404395, 406938) + ylim(168259, 172896) 
#xlim(404395, 404500) + ylim(168259, 168350) 
p_pW


