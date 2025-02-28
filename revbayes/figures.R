####
## BiSSE and HiSSE analyses: Canid diet

###
# load packages
# ggplot2
library(ggplot2)
library(glue)
library(ggridges)
library(forcats)
library(reshape2)

###
# base directory
base_dir <- "/Users/petrucci/Documents/research/canids_chap2/revbayes/"

# read logs
bisse_log <- read.delim(file = paste0(base_dir, "output/bisse.log"), 
                        sep = "\t")
hisse_log <- read.delim(file = paste0(base_dir, "output/hisse.log"), 
                        sep = "\t")
fbdbisse_log <- read.delim(file = paste0(base_dir, "output/fbdbisse.log"), 
                           sep = "\t")
fbdhisse_log <- read.delim(file = paste0(base_dir, "output/fbdhisse.log"), 
                           sep = "\t")

# burnins
bisse_log <- bisse_log[(nrow(bisse_log) / 4):nrow(bisse_log), 5:12]
#hisse_log <- hisse_log[(nrow(hisse_log) / 4):nrow(hisse_log), 5:20]
fbdbisse_log <- fbdbisse_log[(nrow(fbdbisse_log) / 4):nrow(fbdbisse_log), 5:14]
fbdhisse_log <- fbdhisse_log[(nrow(fbdhisse_log) / 4):nrow(fbdhisse_log), 5:24]

# naming
colnames(bisse_log) <- c("lambda0", "lambda1", "mu0", "mu1", "pi0", "pi1", 
                         "q01", "q10")
colnames(hisse_log) <- c("lambda0A", "lambda1A", "lambda0B", "lambda1B",
                         "mu0A", "mu1A", "mu0B", "mu1B",
                         "pi0A", "pi1A", "pi0B", "pi1B",
                         "qAB", "qBA", "q01", "q10")
colnames(fbdbisse_log) <- c("lambda0", "lambda1", "mu0", "mu1", "phi0", "phi1",
                            "pi0", "pi1", "q01", "q10")
colnames(fbdhisse_log) <- c("lambda0A", "lambda1A", "lambda0B", "lambda1B",
                            "mu0A", "mu1A", "mu0B", "mu1B",
                            "phi0A", "phi1A", "phi0B", "phi1B",
                            "pi0A", "pi1A", "pi0B", "pi1B",
                            "qAB", "qBA", "q01", "q10")

###
# BiSSE plots

# get only columns that we care about
bisse_lambda <- bisse_log[, 1:2]
bisse_mu <- bisse_log[, 3:4]

# add iteration column
bisse_lambda$iter <-  bisse_mu$iter <- 1:nrow(bisse_lambda)

# melt
bisse_lambda_final <- melt(bisse_lambda, id.vars = c("iter"))
bisse_mu_final <- melt(bisse_mu, id.vars = c("iter"))

# plot lambda
bisse_ridge_lambda <- ggplot(bisse_lambda_final, 
                              aes(x = value, 
                                  fill = variable,
                                  color = variable)) +
  geom_density(alpha = 0.75) +
  scale_fill_manual(name = "Rate",
                    labels = c(expression(lambda['0']), expression(lambda['1'])),
                    values = c("#E69F00", "#56B4E9")) + 
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  labs(title = "Speciation rate estimates",
       x = "") +
  theme_bw() + 
  guides(color = "none") +
  theme(axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

bisse_ridge_mu <- ggplot(bisse_mu_final, 
                             aes(x = value, 
                                 fill = variable,
                                 color = variable)) +
  geom_density(alpha = 0.75) +
  scale_fill_manual(name = "Rate",
                    labels = c(expression(mu['0']), expression(mu['1'])),
                    values = c("#E69F00", "#56B4E9")) + 
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  labs(title = "Extinction rate estimates",
       x = "") +
  theme_bw() + 
  guides(color = "none") +
  theme(axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

###
# FBD + BiSSE plots

# get only columns that we care about
fbdbisse_lambda <- fbdbisse_log[, 1:2]
fbdbisse_mu <- fbdbisse_log[, 3:4]

# add iteration column
fbdbisse_lambda$iter <-  fbdbisse_mu$iter <- 1:nrow(fbdbisse_lambda)

# melt
fbdbisse_lambda_final <- melt(fbdbisse_lambda, id.vars = c("iter"))
fbdbisse_mu_final <- melt(fbdbisse_mu, id.vars = c("iter"))

# plot lambda
fbdbisse_ridge_lambda <- ggplot(fbdbisse_lambda_final, 
                             aes(x = value, 
                                 fill = variable,
                                 color = variable)) +
  geom_density(alpha = 0.75) +
  scale_fill_manual(name = "Rate",
                    labels = c(expression(lambda['0']), expression(lambda['1'])),
                    values = c("#CB6677", "#89CCEE")) + 
  scale_color_manual(values = c("#CB6677", "#89CCEE")) +
  labs(title = "BISSE speciation rate estimates",
       x = "Rate", y = "Density") +
  theme_bw() + 
  guides(color = "none") +
  theme(axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        #axis.title.y = element_blank(),
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

fbdbisse_ridge_mu <- ggplot(fbdbisse_mu_final, 
                         aes(x = value, 
                             fill = variable,
                             color = variable)) +
  geom_density(alpha = 0.75) +
  scale_fill_manual(name = "Rate",
                    labels = c(expression(mu['0']), expression(mu['1'])),
                    values = c("#CB6677", "#89CCEE")) + 
  scale_color_manual(values = c("#CB6677", "#89CCEE")) +
  labs(title = "BiSSE extinction rate estimates",
       x = "Rate", y = "Density") +
  theme_bw() + 
  guides(color = "none") +
  theme(axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        #axis.title.y = element_blank(),
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

fbdbisse_final <- rbind(fbdbisse_lambda_final, fbdbisse_mu_final)
fbdbisse_final$rate <- c(rep("lambda", nrow(fbdbisse_final) / 2),
                         rep("mu", nrow(fbdbisse_final) / 2))
fbdbisse_final$rate <- factor(fbdbisse_final$rate, 
                              levels = c("lambda", "mu"),
                              labels = c(expression(lambda), expression(mu)))
fbdbisse_final$state <- ifelse(fbdbisse_final$variable %in% c("lambda0", "mu0"),
                               "Hyper", "Non-hyper")

fbdbisse_violin <- ggplot(fbdbisse_final, 
                                aes(x = value, 
                                    fill = variable,
                                    color = variable)) +
  geom_violin(aes(x = state, y = value), alpha = 0.75) +
  facet_wrap(~rate, ncol = 1, labeller = label_parsed) +
  scale_fill_manual(name = "Rate",
                    labels = c(expression(lambda['0']), expression(lambda['1']),
                               expression(mu['0']), expression(mu['1'])),
                    values = c("#CB6677", "#89CCEE", "#CB6677", "#89CCEE")) + 
  scale_color_manual(values = c("#CB6677", "#89CCEE", "#CB6677", "#89CCEE")) +
  scale_x_discrete(label = label_parsed) +
  labs(title = "BISSE speciation rate estimates",
       x = "") +
  theme_bw() + 
  guides(color = "none") +
  theme(#axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

###
# HiSSE plots

# get only columns that we care about
hisse_lambda <- hisse_log[, 1:4]
hisse_mu <- hisse_log[, 5:8]

# add iteration column
hisse_lambda$iter <-  hisse_mu$iter <- 1:nrow(hisse_lambda)

# melt
hisse_lambda_final <- melt(hisse_lambda, id.vars = c("iter"))
hisse_mu_final <- melt(hisse_mu, id.vars = c("iter"))

# plot lambda
hisse_ridge_lambda <- ggplot(hisse_lambda_final, 
                             aes(x = value, 
                                 fill = variable,
                                 color = variable)) +
  geom_density(alpha = 0.75) +
  scale_fill_manual(name = "Rate",
                    labels = c(expression(lambda['0A']), 
                               expression(lambda['1A']),
                               expression(lambda['0B']),
                               expression(lambda['1B'])),
                    values = c("#CB6677", "#89CCEE",  "#44AA99", "#DDCC77")) + 
  scale_color_manual(values = c("#CB6677", "#89CCEE",  "#44AA99", "#DDCC77")) +
  labs(title = "HiSSE speciation rate estimates",
       x = "") +
  theme_bw() + 
  guides(color = "none") +
  theme(axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

hisse_ridge_mu <- ggplot(hisse_mu, 
                             aes(x = value, 
                                 fill = variable,
                                 color = variable)) +
  geom_density(alpha = 0.75) +
  scale_fill_manual(name = "Rate",
                    labels = c(expression(mu['0A']), 
                               expression(mu['1A']),
                               expression(mu['0B']),
                               expression(mu['1B'])),
                    values = c("#CB6677", "#89CCEE",  "#44AA99", "#DDCC77")) + 
  scale_color_manual(values = c("#CB6677", "#89CCEE",  "#44AA99", "#DDCC77")) +
  labs(title = "HiSSE extinctio rate estimates",
       x = "") +
  theme_bw() + 
  guides(color = "none") +
  theme(axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14))

###
# FBD + HiSSE plots

# get only columns that we care about
fbdhisse_lambda <- fbdhisse_log[, 1:4]
fbdhisse_mu <- fbdhisse_log[, 5:8]

# add iteration column
fbdhisse_lambda$iter <-  fbdhisse_mu$iter <- 1:nrow(fbdhisse_lambda)

# melt
fbdhisse_lambda_final <- melt(fbdhisse_lambda, id.vars = c("iter"))
fbdhisse_mu_final <- melt(fbdhisse_mu, id.vars = c("iter"))

# plot lambda
fbdhisse_ridge_lambda <- ggplot(fbdhisse_lambda_final, 
                             aes(x = value, 
                                 fill = variable,
                                 color = variable)) +
  geom_density(alpha = 0.75) +
  scale_fill_manual(name = "Rate",
                    labels = c(expression(lambda['0A']), 
                               expression(lambda['1A']),
                               expression(lambda['0B']),
                               expression(lambda['1B'])),
                    values = c("#CB6677", "#89CCEE",  "#44AA99", "#DDCC77")) + 
  scale_color_manual(values = c("#CB6677", "#89CCEE",  "#44AA99", "#DDCC77")) +
  labs(title = "HiSSE speciation rate estimates",
       x = "") +
  theme_bw() + 
  guides(color = "none") +
  theme(axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

fbdhisse_ridge_mu <- ggplot(fbdhisse_mu_final, 
                         aes(x = value, 
                             fill = variable,
                             color = variable)) +
  geom_density(alpha = 0.75) +
  scale_fill_manual(name = "Rate",
                    labels = c(expression(mu['0A']), 
                               expression(mu['1A']),
                               expression(mu['0B']),
                               expression(mu['1B'])),
                    values = c("#CB6677", "#89CCEE",  "#44AA99", "#DDCC77")) + 
  scale_color_manual(values = c("#CB6677", "#89CCEE",  "#44AA99", "#DDCC77")) +
  labs(title = "HiSSE extinction rate estimates",
       x = "Rate", y = "Density") +
  theme_bw() + 
  guides(color = "none") +
  theme(axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        #axis.title.y = element_blank(),
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

# plot for state A only
fbdhisse_lambda_A <- fbdhisse_lambda_final[fbdhisse_lambda_final$variable %in%
                                             c("lambda0A", "lambda1A"), ]
fbdhisse_mu_A <- fbdhisse_mu_final[fbdhisse_mu_final$variable %in%
                                             c("mu0A", "mu1A"), ]

fbdhisse_ridge_lambda_A <- ggplot(fbdhisse_lambda_A, 
                                aes(x = value, 
                                    fill = variable,
                                    color = variable)) +
  geom_density(alpha = 0.75) +
  scale_fill_manual(name = "Rate",
                    labels = c(expression(lambda['0A']), 
                               expression(lambda['1A'])),
                    values = c("#E69F00", "#56B4E9")) + 
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  labs(title = "Speciation rate estimates",
       x = "") +
  theme_bw() + 
  guides(color = "none") +
  theme(axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

fbdhisse_ridge_mu_A <- ggplot(fbdhisse_mu_A, 
                            aes(x = value, 
                                fill = variable,
                                color = variable)) +
  geom_density(alpha = 0.75) +
  scale_fill_manual(name = "Rate",
                    labels = c(expression(mu['0A']), 
                               expression(mu['1A'])),
                    values = c("#E69F00", "#56B4E9")) + 
  scale_color_manual(values = c("#E69F00", "#56B4E9")) +
  labs(title = "Extinction rate estimates",
       x = "") +
  theme_bw() + 
  guides(color = "none") +
  theme(axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

# plot for state B only
fbdhisse_lambda_B <- fbdhisse_lambda_final[fbdhisse_lambda_final$variable %in%
                                             c("lambda0B", "lambda1B"), ]
fbdhisse_mu_B <- fbdhisse_mu_final[fbdhisse_mu_final$variable %in%
                                     c("mu0B", "mu1B"), ]

fbdhisse_ridge_lambda_B <- ggplot(fbdhisse_lambda_B, 
                                  aes(x = value, 
                                      fill = variable,
                                      color = variable)) +
  geom_density(alpha = 0.75) +
  scale_fill_manual(name = "Rate",
                    labels = c(expression(lambda['0B']), 
                               expression(lambda['1B'])),
                    values = c("#009E73", "#CC79A7")) + 
  scale_color_manual(values = c("#009E73", "#CC79A7")) +
  labs(title = "Speciation rate estimates",
       x = "") +
  theme_bw() + 
  guides(color = "none") +
  theme(axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

fbdhisse_ridge_mu_B <- ggplot(fbdhisse_mu_B, 
                            aes(x = value, 
                                fill = variable,
                                color = variable)) +
  geom_density(alpha = 0.75) +
  scale_fill_manual(name = "Rate",
                    labels = c(expression(mu['0B']), 
                               expression(mu['1B'])),
                    values = c("#009E73", "#CC79A7")) + 
  scale_color_manual(values = c("#009E73", "#CC79A7")) +
  labs(title = "Extinction rate estimates",
       x = "") +
  theme_bw() + 
  guides(color = "none") +
  theme(axis.text.y = element_blank(),  
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        title = element_text(size = 16),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))
