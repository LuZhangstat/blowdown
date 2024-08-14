## real data analysis summary #
rm(list = ls())
library(Matrix)
library(fields)
library(ggplot2)
source("./project/utils.R")
library(dplyr)
library(cmdstanr)
library(bayesplot)

# precomputation #
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

## Compute the average height.m over each plot and summary with observation in one dataset
x_obs_bar <- plot.grid.x_lit %>% group_by(obs.plot.id) %>% 
  summarize(height.mean = sum(height.m * n) / sum(n)) %>% select(obs.plot.id, height.mean) 
obs_xy <- plot.y %>% select(obs.plot.id, sub.region, volume.m3.ha) %>% 
  inner_join(x_obs_bar, by = "obs.plot.id") %>% 
  mutate(intercept = 1) %>%
  # select(volume.m3.ha, intercept, height.mean)
  mutate(x.Frohn = height.mean * (sub.region == "Frohn")) %>%
  mutate(x.Laas = height.mean * (sub.region == "Laas")) %>%
  mutate(x.Mauthen = height.mean * (sub.region == "Mauthen")) %>%
  mutate(x.Liesing = height.mean * (sub.region == "Liesing")) %>%
  mutate(x.Ploecken = height.mean * (sub.region == "Ploecken")) %>%
  select(volume.m3.ha, intercept, x.Frohn, x.Laas, x.Mauthen, x.Liesing, x.Ploecken)

## precomputation for the model fitting ##
HX = as.matrix(obs_xy[, -1]); y = obs_xy[, 1] # HX: the H_{BA} X in the paper
plot.grid.x_lit <- plot.grid.x_lit %>% group_by(obs.plot.id) %>% 
  mutate(counts = sum(n)) %>% mutate(weight = n / counts)
plot.grid.x_lit %>% glimpse()

Dh = plot.grid.x_lit %>% group_by(obs.plot.id) %>% # the sum_{i = 1}^{n_a} h^2_{li} 
  summarise(Dh = sum(weight^2)) %>% pull

HXU <- pred.grid.x_lit %>% group_by(pred.poly.id, sub.region) %>%
  summarize(height.mean = sum(height.m * n) / sum(n)) %>%
  mutate(intercept = 1) %>%
  mutate(x.Frohn = height.mean * (sub.region == "Frohn")) %>%
  mutate(x.Laas = height.mean * (sub.region == "Laas")) %>%
  mutate(x.Mauthen = height.mean * (sub.region == "Mauthen")) %>%
  mutate(x.Liesing = height.mean * (sub.region == "Liesing")) %>%
  mutate(x.Ploecken = height.mean * (sub.region == "Ploecken")) %>%
  select(intercept, x.Frohn, x.Laas, x.Mauthen, x.Liesing, x.Ploecken)

HXU <- as.matrix(HXU[, -1])

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

gamma = 1.8
Dist_M <- rdist(grid.A[, -1] / 1000, grid.A[, -1] / 1000) 

# load the fitted data #
fit <- readRDS("./results/fit_phi_unif.RDS")

## check the fitting results 
fit$summary()
fit_draws <- fit$draws(inc_warmup = FALSE)

################
## prediction ##
################

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

load("./results/RDA_recov_sam_phi_unif.RData")

# recover posterior samples of predictions ##
# set.seed(123)
# yU_ls <- pred_sample_y(beta_omega_sam$beta_ls, beta_omega_sam$omega_BU_ls,
#                        tausq_ls, HXU, DhU)
# save(yU_ls, file = "./results/yU_ls_phi_unif.RData")
load("./results/yU_ls_phi_unif.RData")

yU_pre <-  data.frame(pred.poly.id = unique(pred.grid.x_lit$pred.poly.id), 
                      yU_pm = rowMeans(yU_ls),
                      omega_BU_pm = colMeans(t(beta_omega_sam$omega_BU_ls) + 
                                               beta_omega_sam$beta_ls[1,]))

yB_pre <- data.frame(obs.plot.id = unique(plot.grid.x_lit$obs.plot.id), 
                     omega_B_pm = colMeans(t(beta_omega_sam$omega_B_ls) + 
                                             beta_omega_sam$beta_ls[1,]))


#'Growing stock volume by sub-region and study area wide posterior predictive 
#'distribution median and 95% credible interval

# Area (ha) of prediction #
(pred.grid.x %>% count(sub.region))[, 2]/ 10000

# summary of total Volume #
sub.region_ind <- pred.grid.x_lit %>% group_by(pred.poly.id, sub.region) %>%
  summarize(area = sum(n)/10000) %>%
  mutate(intercept = 1, pred.poly.id = pred.poly.id, sub.region = sub.region) %>%
  select(pred.poly.id, sub.region, area)
summary(1:564-sub.region_ind$pred.poly.id)
volume_yU_ls <- yU_ls * sub.region_ind$area 
# compute the posterior predictive sample of the total volume for each unobserved plot 

region_ind <- (1:564)[!duplicated(sub.region_ind$sub.region)]
combM <- matrix(0, nrow = 5, ncol = nrow(sub.region_ind))
combM[1, 1:(region_ind[2]-1)] = sub.region_ind$area[1:(region_ind[2]-1)]
combM[2, region_ind[2]:(region_ind[3]-1)] = 
  sub.region_ind$area[region_ind[2]:(region_ind[3]-1)]
combM[3, region_ind[3]:(region_ind[4]-1)] = 
  sub.region_ind$area[region_ind[3]:(region_ind[4]-1)]
combM[4, region_ind[4]:(region_ind[5]-1)] = 
  sub.region_ind$area[region_ind[4]:(region_ind[5]-1)]
combM[5, region_ind[5]:nrow(sub.region_ind)] = 
  sub.region_ind$area[region_ind[5]:nrow(sub.region_ind)]
yU_region_sum_ls = combM %*% yU_ls
yU_region_sum <- apply(yU_region_sum_ls, 1, function(x){quantile(x, c(0.5, 0.025, 0.975))})
cbind(unique(sub.region_ind$sub.region), round(t(yU_region_sum), digits = 1))
quantile(t(rep(1, 5)) %*% yU_region_sum_ls, c(0.5, 0.025, 0.975))


