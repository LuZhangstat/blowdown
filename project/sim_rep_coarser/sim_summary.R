rm(list = ls())
library(dplyr)
library(Matrix)
library(fields)
library(ggplot2)
#source("./project/utils.R")

y_O_true_ls <- c()
y_K_true_ls <- c()
y_OK_true_ls <- c()

y_O_m_ls <- c()
y_O_025_ls <- c()
y_O_975_ls <- c()
y_K_m_ls <- c()
y_K_025_ls <- c()
y_K_975_ls <- c()
y_OK_m_ls <- c()
y_OK_025_ls <- c()
y_OK_975_ls <- c()

y_O_2_m_ls <- c()
y_O_2_025_ls <- c()
y_O_2_975_ls <- c()
y_K_2_m_ls <- c()
y_K_2_025_ls <- c()
y_K_2_975_ls <- c()
y_OK_2_m_ls <- c()
y_OK_2_025_ls <- c()
y_OK_2_975_ls <- c()

N_sim = 20
N_ind_O = 126
N_ind_K = 108

## read in results ##
for(i in 1:N_sim){
  result_name <- paste0("./results/sim2/sim_", i, ".Rdata")
  load(result_name)
  y_O_true_ls[i] <- y_O_true
  y_K_true_ls[i] <- y_K_true
  y_OK_true_ls[i] <- y_O_true * N_ind_O + y_K_true * N_ind_K
  
  y_O_m_ls[i] <- mean(yU_ls[1, ])
  y_O_025_ls[i] <- quantile(yU_ls[1, ], 0.025)
  y_O_975_ls[i] <- quantile(yU_ls[1, ], 0.975)
  y_K_m_ls[i] <- mean(yU_ls[2, ])
  y_K_025_ls[i] <- quantile(yU_ls[2, ], 0.025)
  y_K_975_ls[i] <- quantile(yU_ls[2, ], 0.975)
  y_OK_m_ls[i] <- mean(yU_ls[1, ] * N_ind_O + yU_ls[2, ] * N_ind_K) 
  y_OK_025_ls[i] <- quantile(yU_ls[1, ] * N_ind_O + yU_ls[2, ] * N_ind_K, 0.025)
  y_OK_975_ls[i] <- quantile(yU_ls[1, ] * N_ind_O + yU_ls[2, ] * N_ind_K, 0.975)
  
  y_O_2_m_ls[i] <- mean(yU_O)
  y_O_2_025_ls[i] <- quantile(yU_O, 0.025)
  y_O_2_975_ls[i] <- quantile(yU_O, 0.975)
  y_K_2_m_ls[i] <- mean(yU_K)
  y_K_2_025_ls[i] <- quantile(yU_K, 0.025)
  y_K_2_975_ls[i] <- quantile(yU_K, 0.975)
  y_OK_2_m_ls[i] <- mean(yU_O * N_ind_O + yU_K * N_ind_K) 
  y_OK_2_025_ls[i] <- quantile(yU_O * N_ind_O + yU_K * N_ind_K, 0.025)
  y_OK_2_975_ls[i] <- quantile(yU_O * N_ind_O + yU_K * N_ind_K, 0.975)
}

CI_O_dat <- data.frame(
  L_bound = c(y_O_025_ls - y_O_true_ls, y_O_2_025_ls - y_O_true_ls),
  U_bound = c(y_O_975_ls - y_O_true_ls, y_O_2_975_ls - y_O_true_ls),
  sim_id = rep(1:N_sim, 2), true_M = rep(0, 2 * N_sim),
  M_T = c(y_O_m_ls - y_O_true_ls, y_O_2_m_ls - y_O_true_ls),
  approach = rep(c("COS", "Block"), each = N_sim))

CI_O_dat$approach=as.factor(CI_O_dat$approach)

CI_K_dat <- data.frame(
  L_bound = c(y_K_025_ls - y_K_true_ls, y_K_2_025_ls - y_K_true_ls),
  U_bound = c(y_K_975_ls - y_K_true_ls, y_K_2_975_ls - y_K_true_ls),
  sim_id = rep(1:N_sim, 2), true_M = rep(0, 2 * N_sim),
  M_T = c(y_K_m_ls - y_K_true_ls, y_K_2_m_ls - y_K_true_ls),
  approach = rep(c("COS", "Block"), each = N_sim))

CI_K_dat$approach=as.factor(CI_K_dat$approach)

CI_OK_dat <- data.frame(
  L_bound = c(y_OK_025_ls - y_OK_true_ls, y_OK_2_025_ls - y_OK_true_ls),
  U_bound = c(y_OK_975_ls - y_OK_true_ls, y_OK_2_975_ls - y_OK_true_ls),
  sim_id = rep(1:N_sim, 2), true_M = rep(0, 2 * N_sim),
  M_T = c(y_OK_m_ls - y_OK_true_ls, y_OK_2_m_ls - y_OK_true_ls),
  approach = rep(c("COS", "Block"), each = N_sim))

CI_OK_dat$approach=as.factor(CI_OK_dat$approach)


p_O<- ggplot(CI_O_dat, aes(x=sim_id, y=true_M, group=approach, color=approach)) + 
  geom_line(color = "black") +
  geom_point(aes(x=sim_id, y=M_T, group=approach, color=approach),
             position=position_dodge(0.5)) + 
  geom_errorbar(aes(ymin=L_bound, ymax=U_bound), width=.5,
                position=position_dodge(0.5)) + 
  theme_classic() +
  scale_color_manual(values=c('#999999','#E69F00'), name = "") + xlab("Simulation ID") + 
  ylab("CI - true mean")

p_K<- ggplot(CI_K_dat, aes(x=sim_id, y=true_M, group=approach, color=approach)) + 
  geom_line(color = "black") +
  geom_point(aes(x=sim_id, y=M_T, group=approach, color=approach),
             position=position_dodge(0.5)) + 
  geom_errorbar(aes(ymin=L_bound, ymax=U_bound), width=.5,
                position=position_dodge(0.5)) + 
  theme_classic() +
  scale_color_manual(values=c('#999999','#E69F00'), name = "") + xlab("Simulation ID") + 
  ylab("CI - true mean")

p_OK<- ggplot(CI_OK_dat, aes(x=sim_id, y=true_M, group=approach, color=approach)) + 
  geom_line(color = "black") +
  geom_point(aes(x=sim_id, y=M_T, group=approach, color=approach),
             position=position_dodge(0.5)) + 
  geom_errorbar(aes(ymin=L_bound, ymax=U_bound), width=.5,
                position=position_dodge(0.5)) + 
  theme_classic() +
  scale_color_manual(values=c('#999999','#E69F00'), name = "") + xlab("Simulation ID") + 
  ylab("CI - true mean")

p_O
p_K
p_OK

ggsave("./pics/CI_O2.eps", plot = p_O, 
       width = 7.0, height = 3.0, units = "in", dpi = 600)

ggsave("./pics/CI_K2.eps", plot = p_K, 
       width = 7.0, height = 3.0, units = "in", dpi = 600)

# check O+K #
# remove approach and block -> Block
# Simulation

# check all 100 simulations #
for(i in c(1:100)){
  result_name <- paste0("./results/sim2/sim_", i, ".Rdata")
  load(result_name)
  result_name <- paste0("./results/sim2/sim_", i, ".Rdata")
  load(result_name)
  y_O_true_ls[i] <- y_O_true
  y_K_true_ls[i] <- y_K_true
  y_OK_true_ls[i] <- y_O_true * N_ind_O + y_K_true * N_ind_K
  
  y_O_m_ls[i] <- mean(yU_ls[1, ])
  y_O_025_ls[i] <- quantile(yU_ls[1, ], 0.025)
  y_O_975_ls[i] <- quantile(yU_ls[1, ], 0.975)
  y_K_m_ls[i] <- mean(yU_ls[2, ])
  y_K_025_ls[i] <- quantile(yU_ls[2, ], 0.025)
  y_K_975_ls[i] <- quantile(yU_ls[2, ], 0.975)
  y_OK_m_ls[i] <- mean(yU_ls[1, ] * N_ind_O + yU_ls[2, ] * N_ind_K) 
  y_OK_025_ls[i] <- quantile(yU_ls[1, ] * N_ind_O + yU_ls[2, ] * N_ind_K, 0.025)
  y_OK_975_ls[i] <- quantile(yU_ls[1, ] * N_ind_O + yU_ls[2, ] * N_ind_K, 0.975)
  
  y_O_2_m_ls[i] <- mean(yU_O)
  y_O_2_025_ls[i] <- quantile(yU_O, 0.025)
  y_O_2_975_ls[i] <- quantile(yU_O, 0.975)
  y_K_2_m_ls[i] <- mean(yU_K)
  y_K_2_025_ls[i] <- quantile(yU_K, 0.025)
  y_K_2_975_ls[i] <- quantile(yU_K, 0.975)
  y_OK_2_m_ls[i] <- mean(yU_O * N_ind_O + yU_K * N_ind_K) 
  y_OK_2_025_ls[i] <- quantile(yU_O * N_ind_O + yU_K * N_ind_K, 0.025)
  y_OK_2_975_ls[i] <- quantile(yU_O * N_ind_O + yU_K * N_ind_K, 0.975)
}

# fine #
sum((y_O_025_ls <= y_O_true_ls) & 
      (y_O_975_ls >= y_O_true_ls)) 
#0.93
sum((y_K_025_ls <= y_K_true_ls) & 
      (y_K_975_ls >= y_K_true_ls))
#0.95
sum((y_OK_025_ls <= y_OK_true_ls) & 
      (y_OK_975_ls >= y_OK_true_ls))
#0.95 


# Block #
sum((y_O_2_025_ls <= y_O_true_ls) & 
      (y_O_2_975_ls >= y_O_true_ls)) 
#0.97
sum((y_K_2_025_ls <= y_K_true_ls) & 
      (y_K_2_975_ls >= y_K_true_ls)) 
#0.88
sum((y_OK_2_025_ls <= y_OK_true_ls) & 
      (y_OK_2_975_ls >= y_OK_true_ls))
#0.96 

# check Prediction error
MSE_O_COS <- (y_O_m_ls - y_O_true_ls)
MSE_O_block <- (y_O_2_m_ls - y_O_true_ls)
MSE_K_COS <- (y_K_m_ls - y_K_true_ls)
MSE_K_block <- (y_K_2_m_ls - y_K_true_ls)

dta_MSE <- data.frame(MSE = c(MSE_O_block, MSE_O_COS,  MSE_K_block, MSE_K_COS),
                      approach = rep(rep(c(1, 2), each = 100), 2),
                      test = rep(c(1, 2), each = 200))
dta_MSE$approach = factor(dta_MSE$approach, levels = c(1, 2), 
                          labels = c("Block", "COS"))
dta_MSE$test = factor(dta_MSE$test, levels = c(1, 2), 
                      labels = c("O", "K"))

base_plot = dta_MSE %>% ggplot(aes(x=MSE, fill=approach)) + 
  geom_histogram(color="#e9ecef", alpha=0.6, position = 'identity') +
  scale_fill_manual(values=c("#999999", "#E69F00")) +
  theme_bw(base_size = 18) + xlab("Prediction error") + ylab("Count")+
  labs(fill="") + facet_wrap(~ test, nrow = 1) + 
  theme(legend.position="bottom") + 
  geom_vline(aes(xintercept=0),color="black", linetype="dashed", size=1)

base_plot
ggsave("./pics/PE_compar2.png", plot = base_plot, 
       width = 6.5, height = 4.5, units = "in", dpi = 600)

summary(abs(MSE_O_COS)); summary(abs(MSE_O_block))
summary(abs(MSE_K_COS)); summary(abs(MSE_K_block))
median(abs(MSE_O_block))/median(abs(MSE_O_COS)) # 1.090827
median(abs(MSE_K_block))/median(abs(MSE_K_COS)) # 1.12543

#compare CI widths 
summary((y_O_975_ls - y_O_025_ls)/(y_O_2_975_ls - y_O_2_025_ls))
summary((y_K_975_ls - y_K_025_ls)/(y_K_2_975_ls - y_K_2_025_ls))

# rmspe #
RMSPE_O_COS <- sqrt(mean(((y_O_m_ls - y_O_true_ls))^2))
RMSPE_O_block <- sqrt(mean(((y_O_2_m_ls - y_O_true_ls))^2))
RMSPE_K_COS <- sqrt(mean(((y_K_m_ls - y_K_true_ls))^2))
RMSPE_K_block <- sqrt(mean(((y_K_2_m_ls - y_K_true_ls))^2))
round(c(RMSPE_O_COS, RMSPE_O_block, RMSPE_K_COS, RMSPE_K_block), digits = 3)
# 0.282 0.304 0.373 0.434

# empirical coverage probabilities #
block <- c(sqrt(mean(((y_O_2_m_ls - y_O_true_ls))^2)),
           sqrt(mean(((y_K_2_m_ls - y_K_true_ls))^2)),
           median(abs(y_O_2_m_ls - y_O_true_ls)),
           median(abs(y_K_2_m_ls - y_K_true_ls)),
           sum((y_O_2_025_ls <= y_O_true_ls) & 
              (y_O_2_975_ls >= y_O_true_ls))/100,
           sum((y_K_2_025_ls <= y_K_true_ls) & 
                 (y_K_2_975_ls >= y_K_true_ls))/100,
           mean(y_O_2_975_ls - y_O_2_025_ls),
           mean(y_K_2_975_ls - y_K_2_025_ls))
cos <- c(sqrt(mean(((y_O_m_ls - y_O_true_ls))^2)),
         sqrt(mean(((y_K_m_ls - y_K_true_ls))^2)),
         median(abs(y_O_m_ls - y_O_true_ls)),
         median(abs(y_K_m_ls - y_K_true_ls)),
         sum((y_O_025_ls <= y_O_true_ls) & 
               (y_O_975_ls >= y_O_true_ls))/100,
         sum((y_K_025_ls <= y_K_true_ls) & 
               (y_K_975_ls >= y_K_true_ls))/100,
         mean(y_O_975_ls - y_O_025_ls),
         mean(y_K_975_ls - y_K_025_ls))

dt <- t(as.data.frame(round(cbind(block, cos), 2)))
rownames(dt) <- c("Block", "COS")
kable(dt, "latex", col.names = c("RMSPE O", "RMSPE K", "MPE O", "MPE K",
                                 "CI cover O", "CI cover K",  
                                 "CI width O", "CI width K"), 
      booktabs = TRUE)

# \begin{tabular}{lrrrrrrrr}
# \toprule
# & RMSPE O & RMSPE K & MPE O & MPE K & CI cover O & CI cover K & CI width O & CI width K\\
# \midrule
# Block & 0.30 & 0.43 & 0.22 & 0.27 & 0.97 & 0.88 & 1.30 & 1.54\\
# COS & 0.28 & 0.37 & 0.20 & 0.24 & 0.93 & 0.95 & 1.06 & 1.60\\
# \bottomrule
# \end{tabular}

