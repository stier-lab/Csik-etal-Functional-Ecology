

source(here::here("code", "00_libraries.R"))

  
# Treatment level figure


df.treat <- read.csv(here::here("data/foraging/outputs/posteriors_treatments_noGP03.csv")) %>%
  mutate(c_max = 1/t.h*24) %>%
  separate(temp, into = c("junk", "temp_c"), sep = "(?<=[A-Za-z])(?=[0-9])") %>%
  select(-junk) %>%
  mutate(temp_c = as.numeric(temp_c)) %>%
  as_tibble()


boltzmann_constant <- 8.617333262 * 10^-5 #eV K-1

lm_cmax <- lm(log(c_max) ~ I(temp_c/(boltzmann_constant*273.15^2*(1+temp_c/273.15))), data = df.treat)
summary(lm_cmax)

# Bootstrap the CI's... pretty hacky approach.

out <- data.frame(draw = NA, intercept = NA, slope = NA)

for(i in 1:1000){
  temp.df <- df.treat %>%
    group_by(temp_c) %>%
    sample_draws(n = 50)
  
  lm_cmax <- lm(log(c_max) ~ I(temp_c/(boltzmann_constant*273.15^2*(1+temp_c/273.15))), data = temp.df)
  out[i,] <- c(i, lm_cmax$coefficients[1], lm_cmax$coefficients[2]) 
}

slope_ci <- rethinking::PI(out$slope, 0.95)
intercept_ci <- rethinking::PI(out$intercept, 0.95)
median_int <- median(out[,2])
median_slope <- median(out[,3])
mean_slope <- mean(out[, 3])


pred_cmax <- data.frame(x = seq(11,26, by=0.1))
pred_cmax$median <- exp(median_int)*exp(median_slope*pred_cmax$x/(boltzmann_constant*273.15^2*(1+pred_cmax$x/273.15)))
pred_cmax$upper <- exp(intercept_ci[2])*exp(slope_ci[2]*pred_cmax$x/(boltzmann_constant*273.15^2*(1+pred_cmax$x/273.15)))

pred_cmax$lower <- exp(intercept_ci[1])*exp(slope_ci[1]*pred_cmax$x/(boltzmann_constant*273.15^2*(1+pred_cmax$x/273.15)))


background_points <-  df.treat %>%
  group_by(temp_c) %>%
  sample_draws(n = 200)

sum_points <- df.treat %>%
  group_by(temp_c) %>%
  tidybayes::median_qi(c_max)

consumption_plot <- ggplot() + 
  geom_jitter(data = background_points, aes(x = temp_c, y = c_max), color = "black", width = 0.2, alpha = 0.2, size = 3) + 
  geom_line(data = pred_cmax, aes(x = x, y = median), lwd = 1.2)+
  geom_ribbon(data = pred_cmax, aes(x = x, y = median, ymin = lower, ymax = upper), alpha = 0.1)+
  geom_errorbar(data = sum_points, aes(x = temp_c, ymin = .lower, ymax = .upper), width = 0.1, size = 1, color = "black") +
  geom_point(data = sum_points, aes(x =  temp_c,  y = c_max), size = 6, colour = "black") +
  geom_point(data = sum_points, aes(x = temp_c,  y = c_max, color = as.factor(temp_c)), size = 4.5) +
  labs(y = expression(atop("Maximum Consumption Rate", paste("(prey consumed" ~ predator^{-1} ~ "24",~hr^{-1},")"))),
       x = expression(paste("Temperature (", degree, "C)"))) +
  scale_color_manual(values=c("lightslategray", "lightblue", "lightcoral", "indianred4")) +
  theme_classic() + 
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 13),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.7),
        plot.caption = element_text(size = 10, hjust = 0),
        legend.position = "none")+
  scale_x_continuous(limits=c(10, 28),breaks=c(11,16,21,26)) 


cowplot::save_plot(here::here("figures", "main_text", "Fig3.pdf"), consumption_plot, base_width = 5, base_height = 4)



# Individual effects

meta <- read.csv(here::here("data/metadata/", "metadata_lobster_sizes.csv")) %>%
  select(lobster_id, weight_kg, temp) %>%
  rename(id = lobster_id)


df.ind <- read.csv(here::here("data/foraging/outputs/posteriors_individuals_noGP03.csv")) %>%
  mutate(c_max = 1/h*24) %>%
  as_tibble() %>%
  left_join(meta)

# Body mass figure
df.ind %>%
  group_by(id, weight_kg, temp) %>%
  median_qi(c_max) %>%
  ggplot(aes(x = weight_kg, y = c_max))+
  geom_pointinterval(aes(ymin = .lower, ymax = .upper))+
  facet_wrap(~as.factor(temp))
      # Based on this figure I feel that it is safe to assume that body mass isn't having any effect so I don't think its worth correcting for. 

v2 <- df.ind %>%
  group_by(id, temp) %>%
  median_qi(c_max) %>%
  rename(temp_c = temp) %>%
  ggplot(aes(x = temp_c, y = c_max))+
  geom_jitter(size = 3, alpha = 0.5, width = 0.75 )+
  geom_line(data = pred_cmax, aes(x = x, y = median), lwd = 1.2)+
  geom_ribbon(data = pred_cmax, aes(x = x, y = median, ymin = lower, ymax = upper), alpha = 0.1)+
  geom_errorbar(data = sum_points, aes(x = temp_c, ymin = .lower, ymax = .upper), width = 0.4, size = 1, color = "black") +
  geom_point(data = sum_points, aes(x =  temp_c,  y = c_max), size = 6, colour = "black") +
  geom_point(data = sum_points, aes(x = temp_c,  y = c_max, color = as.factor(temp_c)), size = 4.5) +
  labs(y = expression(atop("Maximum consumption rate", paste("(prey consumed" ~ predator^{-1} ~ "24",~hr^{-1},")"))),
       x = expression(paste("Temperature (", degree, "C)"))) +
  scale_color_manual(values=c("lightslategray", "lightblue", "lightcoral", "indianred4")) +
  theme_classic() + 
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 13),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.7),
        plot.caption = element_text(size = 10, hjust = 0),
        legend.position = "none")+
  scale_x_continuous(limits=c(10, 28),breaks=c(11,16,21,26)) 



cowplot::save_plot(here::here("figures", "main_text", "Fig3_v2.pdf"), v2, base_width = 5, base_height = 4)

#-----------------------------------
## Attack rates
#-----------------------------------

sum_points <- df.treat %>%
  group_by(temp_c) %>%
  median_qi(t.a)

plot.a <- df.ind %>%
  group_by(id, temp) %>%
  median_qi(a) %>%
  rename(temp_c = temp) %>%
  ggplot(aes(x = temp_c, y = a))+
  geom_jitter(size = 3, alpha = 0.5, width = 0.75 )+
  geom_errorbar(data = sum_points, aes(x = temp_c, y = t.a, ymin = .lower, ymax = .upper), width = 0.4, size = 1, color = "black") +
  geom_point(data = sum_points, aes(x =  temp_c,  y = t.a), size = 6, colour = "black") +
  geom_point(data = sum_points, aes(x = temp_c,  y = t.a, color = as.factor(temp_c)), size = 4.5) +
  labs(y = "Attack rate",
       x = expression(paste("Temperature (", degree, "C)"))) +
  scale_color_manual(values=c("lightslategray", "lightblue", "lightcoral", "indianred4")) +
  theme_classic() + 
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 13),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.7),
        plot.caption = element_text(size = 10, hjust = 0),
        legend.position = "none")+
  scale_x_continuous(limits=c(10, 28),breaks=c(11,16,21,26)) 

cowplot::save_plot(here::here("figures/supplemental /S_attackrate-by-temp.pdf"), plot.a, base_width = 5, base_height = 4)



#-----------------------------------
## Coefficient plot for individuals
#-----------------------------------

coef_plot <- df.ind %>%
  select(-h) %>%
  mutate(a = a) %>%
  group_by(id, temp) %>%
  pivot_longer(cols = c(a, c_max)) %>%
  group_by(id, temp, name) %>%
  median_qi(value) %>%
  ggplot(aes(x = value, y = forcats::fct_reorder(id, value)))+
  geom_pointinterval(aes(xmin = .lower, xmax = .upper, color = as.factor(temp)))+
  facet_wrap(~name, scales = "free_x")+
  scale_color_manual(values=c("lightslategray", "lightblue", "lightcoral", "indianred4")) +
  labs(x = "", y = "Lobster identity", color = "Temp °C")+
  theme_classic() + 
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 13),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.7),
        plot.caption = element_text(size = 10, hjust = 0), 
        strip.background = element_blank(),
        strip.text.x = element_blank())

cowplot::save_plot(here::here("figures/supplemental /S_coef-individuals.pdf"), coef_plot, base_width = 8, base_height = 6)





