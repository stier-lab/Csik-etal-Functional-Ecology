###**Summary**

# We fit functional responses at the population, treatment, and individual levels using a heirarchical model. In this script, the model is sourced/run, we build Fig4 (functional responses with Bayesian credible interval, and build coefficient plots for the supplement (Fig S1)) 

###**Required Packages:**

    # - tidyverse
    # - here
    # - tidybayes
    # - R2jags
    # - rjags
    # - MCMCvis
    # - cowplot
    # - ggpubr
    # - bbmle
    # - lme4

###**Required Data:**

# (file path: data/foraging/raw)

    # - foraging_assay_data.csv

##############################
# source code to run JAGS model; will automatically import data, packages, and functions
##############################

source(here::here("code", "11b_FR_JAGS_model.R"))

##############################
# this will build a data.frame from the posteriors VERY VERY useful! taken from https://github.com/rpruim/CalvinBayes/blob/master/R/posterior.R 
##############################

df.model <- posterior(model) 

##############################
# pull lobster ids
##############################

ids.a <-names(df.model)[c(1:22, 46, 48:51 )] 
ids.h <-names(df.model)[c(24:45, 47, 52:55 )]

##############################
# basically this produces a list where each slot is a different parameter from the JAGs model as a data.frame. The output in each slot is a data.frame of 4 variables. Basically, the CI function takes a seqence of initial prey densities and then asks the model to predict the posterior distribution of how many prey were killed at that density and extracts the median, and 95% confidence interval. This process in completed for every value in the vector of initial prey densities (N.seq).
##############################

out <- list()
for(i in 1:length(ids.a)){
  out[[i]] <- calculate.CI(a.parameter = ids.a[i], h.parameter = ids.h[i], prob = 0.95)
} 

##############################
# This applies names to each slot in the list.
##############################

names(out) <- c(as.character(sort(unique(foraging_data$lobster_id))), "population_mean", "t11", "t16", "t21", "t26") 

##############################
# this is some funkiness to get a list of data.frames to a single data.frame
##############################

df <- do.call("rbind", out) 

##############################
# this just cleans up the data.frame and gets the id's consistent.
##############################

df <- df %>% mutate(id = as.character(rownames(df))) %>%
  separate(id, into = c("id", "junk"), sep = "[.]") %>% select(-junk) 

##############################
# read in raw data to fiugre otu which lobsters belong in which treatments
##############################

fr_data <- read_csv(here::here("data", "functional_response", "raw", "foraging_assay_data.csv")) %>% 
  mutate(consumption_rate = Killed/24)

formerge <- distinct(fr_data, lobster_id, temp) %>% dplyr::rename(id = lobster_id) %>% mutate(temp = paste("t", temp, sep = "")) 

##############################
# separates individual and treatment level fits for plotting purposes
##############################

df.ind.fits <- df[!df$id %in% c("t11", "t16", "t21", "t26", "population_mean"),] %>% left_join(formerge) 
df.treat.fits <- df[df$id %in% c("t11", "t16", "t21", "t26"), ] %>% dplyr::rename(temp = id) 

##############################
# plot it up and save (Fig 4)
##############################

functional_response_plot <- ggplot() +
  geom_jitter(data = s, aes(x = initial, y = killed/24), alpha = 0.7, shape = 1, width = 0.1) + #shape = lobster_id, , size = 2
  geom_line(data = df.treat.fits, aes(x = N.seq, y = mu, col = factor(temp), group = temp), size = 1.4) +
  geom_ribbon(data = df.treat.fits, aes(x = N.seq, ymin=mu.lower, ymax=mu.upper), linetype=2, alpha=0.1)+
  geom_line(data = df.ind.fits, aes(x = N.seq, y = mu, group = id), alpha = 0.7, color = "darkgray") +
  facet_wrap(~temp, nrow = 1) +
  scale_x_continuous(breaks = seq(0, 60, by = 10)) +
  scale_y_continuous(breaks = seq(0, 1.5, by = 0.25)) +
  scale_color_manual(values = c("lightblue", "lightslategray", "lightcoral", "indianred4")) +
  xlab("Initial Prey Density") +
  ylab(expression(atop("Consumption Rate", paste( "\n(prey consumed" ~ predator^{-1} ~ hr^{-1},")")))) +
  guides(color = guide_legend("Temperature (°C)")) +
  theme_classic(base_size = 12) +
  theme(axis.text = element_text(color = "black", size = 10),
        axis.title = element_text(size = 13),
        axis.title.y = element_text(hjust = 0.5), 
        panel.border = element_rect(colour = "black", fill = NA, size = 0.7),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        strip.text = element_text(size = 12),
        legend.position = "none")

# cowplot::save_plot(here::here("figures", "main_text", "Fig4.pdf"), functional_response_plot, base_width = 10, base_height = 3)

##############################
# plot to look for correlation amongst predictors
##############################

# plot(t.a.1 ~ t.h.1, df.model)
# plot(t.a.2 ~ t.h.2, df.model)
# plot(t.a.3 ~ t.h.3, df.model)
# plot(t.a.4 ~ t.h.4, df.model)
# plot(a.22 ~ h.22, df.model)

##########################################################################################
##########################################################################################
##########################################################################################
# Coefficient Plots
##########################################################################################
##########################################################################################
##########################################################################################

##############################
# get lob ids and associated temperatures
##############################

ids_temps <- read_csv(here::here("data", "metadata", "metadata_lobster_sizes.csv")) %>%
  select(lobster_id, temp) %>%
  rename(id = lobster_id) %>%
  filter(!id %in% c("IV15redo", "GP08", "GP09")) # these lobsters were not usable

##############################
# pull out param estimates from model
##############################

# individual-level (a)
df.ind.a <- as.mcmc(model) %>%
  recover_types(df) %>%
  spread_draws(a[id], h[id]) %>%
  group_by(id) %>%
  median_qi(a)

# individual-level (h)
df.ind.h <- as.mcmc(model) %>%
  recover_types(df) %>%
  spread_draws(a[id], h[id]) %>%
  group_by(id) %>%
  median_qi(h)

# treatment-level (a)
df.treat.a <- as.mcmc(model) %>%
  recover_types(df) %>%
  spread_draws(t.a[tind]) %>%
  median_qi(t.a) %>%
  rename(estimate = t.a) %>%
  mutate(
    temp = case_when(
      tind == "1" ~ "11",
      tind == "2" ~ "16",
      tind == "3" ~ "21",
      tind == "4" ~ "26"
    ),
    parameter = rep("a")
  )

# treatment-level (h)
df.treat.h <- as.mcmc(model) %>%
  recover_types(df) %>%
  spread_draws(t.h[tind]) %>%
  median_qi(t.h) %>%
  rename(estimate = t.h) %>%
  mutate(
    temp = case_when(
      tind == "1" ~ "11",
      tind == "2" ~ "16",
      tind == "3" ~ "21",
      tind == "4" ~ "26"
    ),
    parameter = rep("h")
  )

# pop-level
df.pop <- as.mcmc(model) %>%
  recover_types(df) %>%
  spread_draws(mu.a, mu.h)

##############################
# wrangle and combine INDIVIDUAL param estimates with lob ids and temps
##############################

# individual-level (a)
df.ind.a.temps <- full_join(df.ind.a, ids_temps) %>%
  rename(estimate = a) %>%
  mutate(temp = as.factor(as.character(temp)),
         parameter = rep("a"))

# individual-level (h)
df.ind.h.temps <- full_join(df.ind.h, ids_temps) %>%
  rename(estimate = h) %>%
  mutate(temp = as.factor(as.character(temp)),
         parameter = rep("h"))

# combine dfs for individual-level a and h
df.ind.estimates <- rbind(df.ind.a.temps, df.ind.h.temps) %>%
  mutate(id = fct_relevel(id, c("IV10", "IV11", "IV12", "IV13", "IV14", "MM02",
                                "IV16", "IV17", "IV18", "GP03", "GP05", "GP06",
                                "N14", "N15", "N18", "N19", "L3", "L4",
                                "IV15", "IV19", "GP01", "GP07")))

# create separate df without 11C for inset 
df.ind.estimates.no11 <- df.ind.estimates %>% filter(temp != "11", parameter == "h")

# write.csv(df.ind.estimates, here::here("data", "functional_response", "outputs", "ind_FR_estimates.csv"), row.names = FALSE)

##############################
# wrangle and combine TREATMENT param estimates with lob ids and temps
##############################

# combine dfs for treatment-level a and h
df.treat.estimates <- rbind(df.treat.a, df.treat.h)

# create separate df without 11C for inset 
df.treat.estimates.no11 <- df.treat.estimates %>% filter(temp != "11", parameter == "h")

# write.csv(df.treat.estimates, here::here("data", "functional_response", "outputs", "treat_FR_estimates.csv"), row.names = FALSE)

##############################
# build INDIV coeff plots 
##############################

# primary plot
indiv_plot <- ggplot(df.ind.estimates, aes(y = id, x = estimate)) +
  geom_pointintervalh(aes(color = temp)) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = "Parameter Estimate", y = "Lobster ID", color = "Temperature (°C)") +
  scale_color_manual(values = c("lightslategray", "lightblue", "lightcoral", "indianred4")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 13),
        panel.border = element_rect(colour = "black", fil = NA, size = 0.7),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12))

# inset
indiv_inset_plot <- ggplot(df.ind.estimates.no11, aes(y = id, x = estimate)) +
  geom_pointintervalh(aes(color = temp)) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = "Parameter Estimate", y = "Lobster ID", color = "Temperature (°C)") +
  scale_color_manual(values = c("lightblue", "lightcoral", "indianred4")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 5),
        axis.title = element_text(size = 6),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.7),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")

# add inset to primary plot
indiv_coef_plot <- indiv_plot +
  annotation_custom(ggplotGrob(indiv_inset_plot), xmin = 20, xmax = 150,
                    ymin = 10, ymax = 20)

##############################
# build TREAT coeff plots 
##############################

# primary plot
treat_plot <- ggplot(df.treat.estimates, aes(y = temp, x = estimate)) +
  geom_pointintervalh(aes(color = temp)) +
  facet_wrap(~parameter, scales = "free_x") +
  labs(x = "Parameter Estimate", y = "Lobster ID", color = "Temperature (°C)") +
  scale_color_manual(values = c("lightslategray", "lightblue", "lightcoral", "indianred4")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 13),
        panel.border = element_rect(colour = "black", fil = NA, size = 0.7),
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12))

# inset
treat_inset_plot <- ggplot(df.treat.estimates.no11, aes(y = temp, x = estimate)) +
  geom_pointintervalh(aes(color = temp)) +
  labs(x = "Parameter Estimate", y = "Lobster ID", color = "Temperature (°C)") +
  scale_color_manual(values = c("lightblue", "lightcoral", "indianred4")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 5),
        axis.title = element_text(size = 6),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.7),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = "none")

# add inset to primary plot
treatment_coef_plot <- treat_plot +
  annotation_custom(ggplotGrob(treat_inset_plot), xmin = 6, xmax = 35,
                    ymin = 2, ymax = 4)

# combine plot
coef_plot <- cowplot::plot_grid(treatment_coef_plot, indiv_coef_plot, nrow = 2, labels = "AUTO", align = "vh")

# cowplot::save_plot(here::here("figures", "FigS1.pdf"), coef_plot, base_width = 8, base_height = 8)

##########################################################################################
##########################################################################################
##########################################################################################
# Summary stats
##########################################################################################
##########################################################################################
##########################################################################################

# attack rates
a_summary <- df.ind.estimates %>% 
  filter(parameter == "a") %>% 
  mutate(temp = as.factor(temp)) %>% 
  dplyr::group_by(temp) %>% 
  dplyr::summarize(
    a_avg = round(mean(estimate),3),
    a_sd = round(sd(estimate),3),
    a_se = a_sd/sqrt(length(estimate)),
    a_cv = a_sd/a_avg,
    a_max = max(estimate),
    a_min = min(estimate)
  )

write.csv(a_summary, here::here("data", "functional_response", "outputs", "a_summary_stats.csv"), row.names = FALSE)

# handling time
h_summary <- df.ind.estimates %>% 
  filter(parameter == "h") %>% 
  mutate(temp = as.factor(temp)) %>% 
  dplyr::group_by(temp) %>% 
  dplyr::summarize(
    h_avg = round(mean(estimate),3),
    h_sd = round(sd(estimate),3),
    h_se = h_sd/sqrt(length(estimate)),
    h_cv = h_sd/h_avg,
    h_max = max(estimate),
    h_min = min(estimate)
  )

write.csv(h_summary, here::here("data", "functional_response", "outputs", "h_summary_stats.csv"), row.names = FALSE)
