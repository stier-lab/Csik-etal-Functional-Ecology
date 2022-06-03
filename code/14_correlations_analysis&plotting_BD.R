##############################
# load required packages
##############################

source(here::here("code", "00_libraries.R"))

##############################
# load data
##############################

metabolic_rates <- read_csv(here::here("data", "metabolism", "metabolic_traits.csv")) %>%
  rename(id = ID)

df <- read_csv(here::here("data/foraging/outputs/posteriors_individuals_noGP03.csv")) %>%
  mutate(c_max = 1/h*24) %>%
  group_by(id)%>%
  median_qi(c_max, .width = 0.95) %>%
  left_join(metabolic_rates) %>%
  mutate(temp_fac = as.factor(temp), 
         c_max = c_max/BW, 
         .upper = .upper/BW, 
         .lower = .lower/BW) %>%
  filter(id != "GP03")

##############################
# Correlations
##############################

lm_smr <- lm(c_max ~ SMR, df)

lm_mmr <- lm(c_max ~ MMR, df)

lm_aas <- lm(c_max ~ AAS, df)

lm_fas <- lm(c_max ~ FAS, df)



# Some stuff for plotting


points <- filter(df, id %in% c("IV10", "IV19")) %>% select(id, c_max, temp, SMR:FAS)

lm_eqn <- function(mod){
  eq <- substitute(~~italic(r)^2~"="~r2, list(r2 = format(summary(mod)$r.squared, digits = 2)))
  as.character(as.expression(eq));
}

lm_eqn(lm_smr)


##############################
# smr
##############################

smr_corr_plot <- df %>%
ggplot(aes(x = SMR, y = c_max)) + 
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0.002, color = "gray") +
  geom_point(aes(color = temp_fac, shape = temp_fac), size = 5) + 
  geom_smooth(method = "lm", color = "gray36", size = 0.75, level = 0.95) +
  scale_color_manual(values = c("lightslategray", "lightblue", "lightcoral", "indianred4"), name = "Temperature (°C)", labels = c("11", "16", "21", "26")) +
  scale_shape_manual(values = c(15, 16, 17, 18), name = "Temperature (°C)", labels = c("11", "16", "21", "26")) +
  labs(x = expression(atop("Standard metabolic rate", paste("(",mg~O[2]~kg^-1~min^-1,")"))),
       y = expression(atop("Maximum consumption rate", paste("(prey consumed" ~ kg^{-1} ~ "24",~hr^{-1},")")))) +
  scale_x_continuous(breaks = seq(0, 1.5, by = 0.2))+
  scale_y_continuous(breaks = seq(0, 120, by = 40))+
  coord_cartesian(ylim = c(-15, 130))+
  annotate(geom = "point", x = points$SMR, y = points$c_max, pch = 4, size = 5)+
  annotate(geom = "text", x = 1.0, y = 10, label = lm_eqn(lm_smr), parse = T, size = 6)+
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 13),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.7),
        plot.caption = element_text(size = 10, hjust = 0),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 8),
        legend.position = "none")

smr_marg <- ggMarginal(smr_corr_plot, type = "density", size = 4, groupColour = TRUE, groupFill = TRUE)

##############################
# mmr
##############################

mmr_corr_plot <- df %>%
  ggplot(aes(x = MMR, y = c_max)) + 
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0.002, color = "gray") +
  geom_point(aes(color = temp_fac, shape = temp_fac), size = 5) + 
  geom_smooth(method = "lm", color = "gray36", size = 0.75, level = 0.95) +
  scale_color_manual(values = c("lightslategray", "lightblue", "lightcoral", "indianred4"), name = "Temp (°C)") +
  scale_shape_manual(values = c(15, 16, 17, 18), name = "Temp (°C)") +
  labs(x = expression(atop("Maximum metabolic rate", paste("(",mg~O[2]~kg^-1~min^-1,")"))),
       y = expression(atop("Maximum consumption rate", paste("(prey consumed" ~ kg^{-1} ~ "24",~hr^{-1},")")))) +
  scale_x_continuous(breaks = seq(0, 7, by = 1)) +
  scale_y_continuous(breaks = seq(0, 120, by = 40))+
  coord_cartesian(ylim = c(-15, 130))+
  annotate(geom = "point", x = points$MMR, y = points$c_max, pch = 4, size = 5)+
  annotate(geom = "text", x = 5, y = 10, label = lm_eqn(lm_mmr), parse = T, size = 6)+
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 13),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.7),
        plot.caption = element_text(size = 10, hjust = 0),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = c(0.05, 0.95),
        legend.justification = c("left", "top"),
        legend.box.just = "left",
        legend.box.background = element_rect(color = "black", size = 1.1),
        legend.margin = margin(3, 3, 3, 3),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

mmr_marg <- ggMarginal(mmr_corr_plot, type = "density", size = 4, groupColour = TRUE, groupFill = TRUE)

##############################
# aas
##############################

aas_corr_plot <- df %>%
  ggplot(aes(x = AAS, y = c_max)) + 
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0.002, color = "gray") +
  geom_point(aes(color = temp_fac, shape = temp_fac), size = 5) + 
  geom_smooth(method = "lm", color = "gray36", size = 0.75, level = 0.95) +
  scale_color_manual(values = c("lightslategray", "lightblue", "lightcoral", "indianred4"), name = "Temperature (°C)") +
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  scale_size_manual(values = c(4, 4, 4, 5)) +
  labs(x = expression(atop("Absolute aerobic scope", paste("(",mg~O[2]~kg^-1~min^-1,")"))),
       y = expression(atop("Maximum consumption rate", paste("(prey consumed" ~ kg^{-1} ~ "24",~hr^{-1},")")))) +
  scale_x_continuous(breaks = seq(0, 7, by = 1)) +
  scale_y_continuous(breaks = seq(0, 120, by = 40))+
  coord_cartesian(ylim = c(-15, 130))+
  annotate(geom = "point", x = points$AAS, y = points$c_max, pch = 4, size = 5)+
  annotate(geom = "text", x = 4.5, y = 10, label = lm_eqn(lm_aas), parse = T, size = 6)+
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 13),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.7),
        plot.caption = element_text(size = 10, hjust = 0),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "none")

aas_marg <- ggMarginal(aas_corr_plot, type = "density", size = 4, groupColour = TRUE, groupFill = TRUE)

##############################
# fas
##############################

fas_corr_plot <- df %>%
  ggplot(aes(x = FAS, y = c_max)) + 
  geom_errorbar(aes(ymin = .lower, ymax = .upper), width = 0.002, color = "gray") +
  geom_point(aes(color = temp_fac, shape = temp_fac), size = 5) + 
  geom_smooth(method = "lm", color = "gray36", size = 0.75, level = 0.95) +
  scale_color_manual(values = c("lightslategray", "lightblue", "lightcoral", "indianred4"), name = "Temperature (°C)") +
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  scale_size_manual(values = c(4, 4, 4, 5)) +
  scale_y_continuous(breaks = seq(0, 120, by = 40))+
  coord_cartesian(ylim = c(-15, 130))+
  labs(x = expression(atop("Factorial aerobic scope", paste("(MMR/SMR)"))),
       y = expression(atop("Maximum consumption rate", paste("(prey consumed" ~ kg^{-1} ~ "24",~hr^{-1},")")))) +
  annotate(geom = "point", x = points$FAS, y = points$c_max, pch = 4, size = 5)+
  annotate(geom = "text", x = 20, y = 120, label = lm_eqn(lm_fas), parse = T, size = 6)+
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 13),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.7),
        plot.caption = element_text(size = 10, hjust = 0),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = "none",
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

fas_marg <- ggMarginal(fas_corr_plot, type = "density", size = 4, groupColour = TRUE, groupFill = TRUE)

##############################
# combine plots and print
##############################

correlations_plot <- cowplot::plot_grid(smr_marg, mmr_marg, aas_marg, fas_marg, nrow = 2, labels = "AUTO", axis = "tblr", rel_widths = c(0.95, 0.80, 0.95, 0.80))

cowplot::save_plot(here::here("figures", "main_text", "Fig5_updated2_noGP03.pdf"), correlations_plot, base_width = 10, base_height = 10)

cowplot::save_plot(here::here("figures", "main_text", "Fig5_updated2_noGP03.svg"), correlations_plot, base_width = 10, base_height = 10)













