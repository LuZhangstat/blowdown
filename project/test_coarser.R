rm(list = ls())
library(Matrix)
library(fields)
library(ggplot2)
source("./project/utils.R")
library(dplyr)

###############
## read data ##
###############
plot.y <- read.csv("./data/outcome_y_plot.csv")
plot.y %>% glimpse

plot.grid.x <- read.csv("./data/predictor_x_grid.csv")
plot.grid.x %>% glimpse


pred.grid.x <- read.csv("./data/prediction_x_grid.csv")
pred.grid.x %>% glimpse

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


pred.grid.x$coords.x2 <- (pred.grid.x $coords.x + 
                             1 * (pred.grid.x $coords.x %% 3 == 1) - 
                             1 * (pred.grid.x $coords.x %% 3 == 0))
pred.grid.x$coords.y2 <- (pred.grid.x $coords.y + 
                             1 * (pred.grid.x $coords.y %% 3 == 1) - 
                             1 * (pred.grid.x $coords.y %% 3 == 0))
pred.grid.x_lit <- pred.grid.x %>% group_by(sub.region, pred.poly.id, coords.x2, 
                                            coords.y2) %>% 
  summarise(height.m = mean(height.m), n = n()) %>% arrange(pred.poly.id) %>%
  select(sub.region, pred.poly.id, coords.x2, coords.y2, height.m, n) 
pred.grid.x_lit %>% glimpse()


## check the coarser grids and the corresponding weights ##
test_1 <- plot.grid.x %>% filter(obs.plot.id == 45)
test_1_lit <- plot.grid.x_lit %>% filter(obs.plot.id == 45) 

plot(test_1$coord.x, test_1$coord.y, 
     xlim = c(min(test_1$coord.x, test_1_lit$coord.x2)-1, 
              max(test_1$coord.x, test_1_lit$coord.x2)+1), 
     ylim = c(min(test_1$coord.y, test_1_lit$coord.y2)-1, 
              max(test_1$coord.y, test_1_lit$coord.y2)+1))

plot(test_1_lit$coord.x2, test_1_lit$coord.y2, type="n",
     xlim = c(min(test_1$coord.x, test_1_lit$coord.x2)-1, 
              max(test_1$coord.x, test_1_lit$coord.x2)+1), 
     ylim = c(min(test_1$coord.y, test_1_lit$coord.y2)-1, 
              max(test_1$coord.y, test_1_lit$coord.y2)+1))
text(test_1_lit$coord.x2, test_1_lit$coord.y2, test_1_lit$n, pos=1, offset = 0.0)

# check for prediction 
pre_1 <- pred.grid.x %>% filter(pred.poly.id == 492)
plot(pre_1$coords.x, pre_1$coords.y, 
     xlim = c(min(pre_1$coords.x)-1, max(pre_1$coords.x) + 1), 
     ylim = c(min(pre_1$coords.y)-1, max(pre_1$coords.y) + 1))

pre_1_lit <- pred.grid.x_lit %>% filter(pred.poly.id == 492)

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
  select(volume.m3.ha, intercept, height.mean)
  # mutate(x.Frohn = height.mean * (sub.region == "Frohn")) %>%
  # mutate(x.Laas = height.mean * (sub.region == "Laas")) %>%
  # mutate(x.Mauthen = height.mean * (sub.region == "Mauthen")) %>%
  # mutate(x.Liesing = height.mean * (sub.region == "Liesing")) %>%
  # mutate(x.Ploecken = height.mean * (sub.region == "Ploecken")) %>% 
  # select(volume.m3.ha, intercept, x.Frohn, x.Laas, x.Mauthen, x.Liesing, x.Ploecken)
  
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


HXU <- cbind(1, pred.grid.x_lit %>% group_by(pred.poly.id) %>% # the sum_{i = 1}^{n_a} h_{li} x(A_i^u) for all predict plots
               summarize(height.mean = sum(height.m * n) / sum(n)) %>% 
               select(height.mean) %>% pull)  # sum_{i = 1}^{n_u} h_i^u 
pred.grid.x_lit <- pred.grid.x_lit %>% group_by(pred.poly.id) %>% 
  mutate(counts = sum(n)) %>% mutate(weight = n / counts)
pred.grid.x_lit %>% glimpse()
DhU = pred.grid.x_lit %>% group_by(pred.poly.id) %>%  #the sum_{i = 1}^{n_a} h^2_{li}  for all prediction plots
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
ss = 10 * sqrt(2)       # scale parameter in the normal prior of sigma 
st = 10 * sqrt(2)     # scale parameter in the normal prior of tau     
ap = 3; bp = 0.5       # shape and rate parameters in the Gamma prior of phi 

data <- list(na = na, nb = nb, p = p, y = y, HX = HX, Dh = Dh, 
             gridA = grid.A[, -1] / 1000, hA = plot.grid.x_lit$weight,
             plotid_ind = plotid_ind,
             mu_beta = mu_beta, V_beta = V_beta,
             ap = ap, bp = bp, ss = ss, st = st)

fit <- mod$sample(
  data = data,
  seed = 1,
  chains = 4,
  parallel_chains = 4,
  refresh = 50,
  save_warmup = TRUE,
  iter_warmup = 500,
  iter_sampling = 500,
  sig_figs = 18
)

fit$summary()
n_lf <- fit$sampler_diagnostics(inc_warmup = TRUE)[, 1, "n_leapfrog__"]
sum(n_lf)

fit$save_object(file = "./results/fit.RDS")

fit <- readRDS("./results/fit.RDS")

## check the fitting results 
fit$summary()
mcmc_trace(fit$draws("sigmasq"), iter1 = 1) 
mcmc_trace(fit$draws("tausq"), iter1 = 1) 
mcmc_trace(fit$draws("phi"), iter1 = 1) 
n_lf <- fit$sampler_diagnostics(inc_warmup = TRUE)[, , "n_leapfrog__"]
colSums(n_lf)
fit_draws <- fit$draws(inc_warmup = FALSE)

################
## prediction ##
################

pick_sample_id <- seq(5, 2000, by = 10) # pick 200 samples
#pick_sample_id <- seq(5, 2000, by = 100) # pick 200 samples
phi_ls <- c(fit_draws[, , "phi"])[pick_sample_id]
sigmasq_ls <- c(fit_draws[, , "sigmasq"])[pick_sample_id]
tausq_ls <- c(fit_draws[, , "tausq"])[pick_sample_id]
coords_A = grid.A[, -1] /1000
coords_AU = grid.Au[, -1] / 1000
ind_ls_B = plotid_ind
ind_ls_BU = predid_ind
hA = plot.grid.x_lit$weight
hAU = pred.grid.x_lit$weight

## recover posterior samples of beta and omega for both observed and unobserved plots ##
t <- proc.time()
beta_omega_sam <- sample_beta_omega_h_quick(phi_ls, sigmasq_ls, tausq_ls,
                                    coords_A, coords_AU,
                                    hA, hAU, ind_ls_B, ind_ls_BU,
                                    HX, mu_beta, V_beta, flat_prior = FALSE)
proc.time() -t
save(beta_omega_sam, file = "./results/RDA_recov_sam.RData")
load("./results/RDA_recov_sam.RData")

## recover posterior samples of predictions ##
yU_ls <- pred_sample_y(beta_omega_sam$beta_ls, beta_omega_sam$omega_BU_ls, 
                       tausq_ls, HXU, DhU)

yU_pre <-  data.frame(pred.poly.id = unique(pred.grid.x_lit$pred.poly.id), 
                      yU_pm = rowMeans(yU_ls),
                      omega_BU_pm = colMeans(t(beta_omega_sam$omega_BU_ls) + 
                                               beta_omega_sam$beta_ls[1,]))

yB_pre <- data.frame(obs.plot.id = unique(plot.grid.x_lit$obs.plot.id), 
                     omega_B_pm = colMeans(t(beta_omega_sam$omega_B_ls) + 
                                             beta_omega_sam$beta_ls[1,]))



#######################
## check the results ##
#######################

## plots ##

obs_plot_dat <- plot.grid.x_lit %>% merge(plot.y[c("obs.plot.id", "volume.m3.ha")], 
                                by = "obs.plot.id") %>% 
  merge(yB_pre, by = "obs.plot.id")%>%
  select(sub.region, obs.plot.id, coord.x2, coord.y2, volume.m3.ha, omega_B_pm)
pred_plot_dat <- pred.grid.x_lit %>% 
  merge(yU_pre, by = "pred.poly.id") %>% mutate(volume.m3.ha = yU_pm) %>%
  select(sub.region, pred.poly.id, coords.x2, coords.y2, volume.m3.ha, omega_BU_pm)
pred_plot_dat %>% glimpse()


yU_pre <- pred.grid.x_lit %>% distinct(sub.region, pred.poly.id) %>% 
  merge(yU_pre, by = "pred.poly.id") %>% arrange() %>%
  select(sub.region, pred.poly.id, yU_pm, omega_BU_pm)
yU_pre %>% glimpse() 


# average response
regions <- unique(pred_plot_dat$sub.region)
region = regions[2]; region  #Liesing

p_o <- ggplot(obs_plot_dat, aes(x = coord.x2, y = coord.y2)) + 
  geom_point(size = 0.01, aes(colour = volume.m3.ha)) + 
  xlim(min(obs_plot_dat$coord.x2[obs_plot_dat$sub.region == region], 
           pred_plot_dat$coords.x2[pred_plot_dat$sub.region == region]) - 1,
       max(obs_plot_dat$coord.x2[obs_plot_dat$sub.region == region], 
           pred_plot_dat$coords.x2[pred_plot_dat$sub.region == region]) + 1) +
  ylim(min(obs_plot_dat$coord.y2[obs_plot_dat$sub.region == region], 
           pred_plot_dat$coords.y2[pred_plot_dat$sub.region == region]) - 1,
        max(obs_plot_dat$coord.y2[obs_plot_dat$sub.region == region], 
            pred_plot_dat$coords.y2[pred_plot_dat$sub.region == region]) + 1) +
  # xlim(410400, 410950) + ylim(171200, 172200) +
  scale_colour_gradient(limits = range(obs_plot_dat$volume.m3.ha,
                                    pred_plot_dat$volume.m3.ha)) 
p_o

p_py <-  ggplot(pred_plot_dat, aes(x = coords.x2, y = coords.y2)) + 
  geom_point(size = 0.01, aes(colour = volume.m3.ha))+
  geom_point(data = obs_plot_dat,
             size = 0.1, 
             aes(x = coord.x2, y = coord.y2,
                              fill = volume.m3.ha, colour = volume.m3.ha)) +
  geom_point(data = plot.y, aes(x=coords.x, y=coords.y), size=3, shape=1,
             color="orange") +
  scale_fill_gradient(limits = range(obs_plot_dat$volume.m3.ha,
                                       pred_plot_dat$volume.m3.ha))+
  xlim(min(obs_plot_dat$coord.x2[obs_plot_dat$sub.region == region], 
           pred_plot_dat$coords.x2[pred_plot_dat$sub.region == region]) - 1,
       max(obs_plot_dat$coord.x2[obs_plot_dat$sub.region == region], 
           pred_plot_dat$coords.x2[pred_plot_dat$sub.region == region]) + 1) +
  ylim(min(obs_plot_dat$coord.y2[obs_plot_dat$sub.region == region], 
           pred_plot_dat$coords.y2[pred_plot_dat$sub.region == region]) - 1,
       max(obs_plot_dat$coord.y2[obs_plot_dat$sub.region == region], 
           pred_plot_dat$coords.y2[pred_plot_dat$sub.region == region]) + 1) 
  #xlim(410400, 410950) + ylim(171200, 172200) 
  
p_py

# average latent w #
p_oW <- ggplot(obs_plot_dat, aes(x = coord.x2, y = coord.y2, colour = omega_B_pm)) + 
  geom_point(size = 0.01) + 
  xlim(min(obs_plot_dat$coord.x2[obs_plot_dat$sub.region == region], 
           pred_plot_dat$coords.x2[pred_plot_dat$sub.region == region]) - 1,
       max(obs_plot_dat$coord.x2[obs_plot_dat$sub.region == region], 
           pred_plot_dat$coords.x2[pred_plot_dat$sub.region == region]) + 1) +
  ylim(min(obs_plot_dat$coord.y2[obs_plot_dat$sub.region == region], 
           pred_plot_dat$coords.y2[pred_plot_dat$sub.region == region]) - 1,
       max(obs_plot_dat$coord.y2[obs_plot_dat$sub.region == region], 
           pred_plot_dat$coords.y2[pred_plot_dat$sub.region == region]) + 1) +
  #xlim(410400, 410950) + ylim(171200, 172200) +
  scale_colour_gradient(limits = range(obs_plot_dat$omega_B_pm, 
                                       pred_plot_dat$omega_BU_pm))
p_oW

p_pW <-  ggplot(obs_plot_dat, aes(x = coord.x2, y = coord.y2, colour = omega_B_pm)) + 
  geom_point(size = 0.01) + 
  geom_point(data = plot.y, aes(x=coords.x, y=coords.y), size=3, shape=1, 
             color="orange") +
  geom_point(data = pred_plot_dat, aes(x = coords.x2, y = coords.y2, 
                                       colour = omega_BU_pm), size = 0.01) +
  xlim(min(obs_plot_dat$coord.x2[obs_plot_dat$sub.region == region], 
           pred_plot_dat$coords.x2[pred_plot_dat$sub.region == region]) - 1,
       max(obs_plot_dat$coord.x2[obs_plot_dat$sub.region == region], 
           pred_plot_dat$coords.x2[pred_plot_dat$sub.region == region]) + 1) +
  ylim(min(obs_plot_dat$coord.y2[obs_plot_dat$sub.region == region], 
           pred_plot_dat$coords.y2[pred_plot_dat$sub.region == region]) - 1,
       max(obs_plot_dat$coord.y2[obs_plot_dat$sub.region == region], 
           pred_plot_dat$coords.y2[pred_plot_dat$sub.region == region]) + 1) +
  #xlim(410400, 410950) + ylim(171200, 172200) +
  scale_colour_gradient(limits = range(obs_plot_dat$omega_B_pm, 
                                       pred_plot_dat$omega_BU_pm))
p_pW


sum(yU_pre$yU_pm[yU_pre$sub.region == region] / 
      DhU[yU_pre$sub.region == region] /1000)

quantile(colSums(yU_ls[yU_pre$sub.region == region, ] * 
                   (1/DhU[yU_pre$sub.region == region])/1000), 
         c(0.025, 0.975))

# "Frohn": 39099.63 (37601.44 41064.57)
# "Laas": 43478.74 (42541.74 44803.45)
# "Mauthen": 34840.27 (33900.24 35563.55)
# "Liesing": 2904.278 (2726.702 3010.066)
# "Ploecken": 37368.63 (36423.37 38572.36)

## zoom in subregions ##
# filter large plots for prediction 
pred_count <- pred_plot_dat %>% group_by(pred.poly.id) %>% 
  summarise(n = n()) %>% select(pred.poly.id, n)
pred_count %>% glimpse()
sum(pred_count$n < 140) / length(pred_count$n) 
# 48% pred plots are smaller than the observed plots
sum(pred_count$n < 140*2) / length(pred_count$n) 
# two-third of pred plots are smaller that 2 times the observed plots

pick_pred_plot_id <- pred_count$pred.poly.id[pred_count$n < 140*2]
pick_small_pred_ind <- sapply(pred_plot_dat$pred.poly.id, 
                              f <- function(x){any(x == pick_pred_plot_id)})

pred_plot_dat_small <- pred_plot_dat[pick_small_pred_ind, ]

# 1. zoom in region Ploecken
region = regions[5]; region  

p_py <-  ggplot(pred_plot_dat_small, aes(x = coords.x2, y = coords.y2)) + 
  geom_point(size = 0.01, aes(colour = volume.m3.ha), 
             show.legend = FALSE)+
  geom_point(data = obs_plot_dat,
             size = 0.01, 
             aes(x = coord.x2, y = coord.y2,
                 #fill = volume.m3.ha, 
                 colour = volume.m3.ha)) +
  geom_point(data = plot.y, aes(x=coords.x, y=coords.y), size=5, shape=1,
             color="orange", stroke = 1) +
  scale_color_gradient(limits = range(obs_plot_dat$volume.m3.ha,
                                     pred_plot_dat$volume.m3.ha),
                      name = expression(m^3/ha))+
  xlim(422300, 423600) + ylim(163223, 164100) + xlab("x") + ylab("y") 
p_py

ggsave("./pics/Ploecken_zoom.png", plot = p_py, 
       width = 6, height = 3.5, units = "in", dpi = 600)

# zoom in Mauthen
region = regions[3]; region  

p_py <-  ggplot(pred_plot_dat_small, aes(x = coords.x2, y = coords.y2)) + 
  geom_point(size = 0.01, aes(colour = volume.m3.ha), 
             show.legend = FALSE)+
  geom_point(data = obs_plot_dat,
             size = 0.01, 
             aes(x = coord.x2, y = coord.y2,
                 #fill = volume.m3.ha, 
                 colour = volume.m3.ha)) +
  geom_point(data = plot.y, aes(x=coords.x, y=coords.y), size=5, shape=1,
             color="orange", stroke = 1) +
  scale_color_gradient(limits = range(obs_plot_dat$volume.m3.ha,
                                      pred_plot_dat$volume.m3.ha),
                       name = expression(m^3/ha))+
  xlim(421000, 422200) + ylim(167500, 168300) + xlab("x") + ylab("y") 
p_py

ggsave("./pics/Mauthen_zoom.png", plot = p_py, 
       width = 6.5, height = 3.5, units = "in", dpi = 600)


# zoom in Lass
region = regions[2]; region  

p_py <-  ggplot(pred_plot_dat_small, aes(x = coords.x2, y = coords.y2)) + 
  geom_point(size = 0.01, aes(colour = volume.m3.ha), 
             show.legend = FALSE)+
  geom_point(data = obs_plot_dat,
             size = 0.01, 
             aes(x = coord.x2, y = coord.y2,
                 #fill = volume.m3.ha, 
                 colour = volume.m3.ha)) +
  geom_point(data = plot.y, aes(x=coords.x, y=coords.y), size=5, shape=1,
             color="orange", stroke = 1) +
  scale_color_gradient(limits = range(obs_plot_dat$volume.m3.ha,
                                      pred_plot_dat$volume.m3.ha),
                       name = expression(m^3/ha))+
  xlim(422800, 424300) + ylim(173500, 174300) + xlab("x") + ylab("y") 
p_py

ggsave("./pics/Laas_zoom_small.png", plot = p_py, 
       width = 6.5, height = 3.0, units = "in", dpi = 600)

p_py2 <-  ggplot(pred_plot_dat, aes(x = coords.x2, y = coords.y2)) + 
  geom_point(size = 0.01, aes(colour = volume.m3.ha), 
             show.legend = FALSE)+
  geom_point(data = obs_plot_dat,
             size = 0.01, 
             aes(x = coord.x2, y = coord.y2,
                 #fill = volume.m3.ha, 
                 colour = volume.m3.ha)) +
  geom_point(data = plot.y, aes(x=coords.x, y=coords.y), size=5, shape=1,
             color="orange", stroke = 1) +
  scale_color_gradient(limits = range(obs_plot_dat$volume.m3.ha,
                                      pred_plot_dat$volume.m3.ha),
                       name = expression(m^3/ha))+
  xlim(422800, 424300) + ylim(173500, 174300) + xlab("x") + ylab("y") 
p_py2

ggsave("./pics/Laas_zoom.png", plot = p_py2, 
       width = 6.5, height = 3.0, units = "in", dpi = 600)


## create a grid for Lass zoom in region ##

p_py_g <-  ggplot(pred_plot_dat, aes(x = coords.x2, y = coords.y2)) + 
  geom_point(size = 0.01, aes(colour = volume.m3.ha), 
             show.legend = FALSE)+
  geom_point(data = obs_plot_dat,
             size = 0.01, 
             aes(x = coord.x2, y = coord.y2,
                 #fill = volume.m3.ha, 
                 colour = volume.m3.ha)) +
  geom_point(data = plot.y, aes(x=coords.x, y=coords.y), size=5, shape=1,
             color="orange", stroke = 1) +
  geom_vline(xintercept = seq(422800, 424300, 36)) + 
  geom_hline(yintercept = seq(173750, 174300, 36)) +
  scale_color_gradient(limits = range(obs_plot_dat$volume.m3.ha,
                                      pred_plot_dat$volume.m3.ha),
                       name = expression(m^3/ha))+
  xlim(422800, 424300) + ylim(173750, 174300) + xlab("x") + ylab("y") 

p_py_g

ggsave("./pics/Laas_zoom_grid.png", plot = p_py_g, 
       width = 7.0, height = 3.0, units = "in", dpi = 600)



