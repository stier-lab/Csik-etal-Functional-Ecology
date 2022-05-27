##############################
# import libraries
##############################

source(here::here("code", "00_libraries.R"))

##############################
# load data
##############################

metabolic_traits <- read_csv(here::here("data", "metabolism", "metabolic_traits.csv"))

smr_sum <- read_csv(here::here("data", "metabolism", "outputs", "smr_summary_stats.csv"))
mmr_sum <- read_csv(here::here("data", "metabolism", "outputs", "mmr_summary_stats.csv"))
aas_sum <- read_csv(here::here("data", "metabolism", "outputs", "aas_summary_stats.csv"))
fas_sum <- read_csv(here::here("data", "metabolism", "outputs", "fas_summary_stats.csv"))
rhr_sum <- read_csv(here::here("data", "heart_rate", "outputs", "rhr_summary_stats.csv"))

# Looks like there is still a mass effect, so need to account for this in the consumption data.

ggplot(metabolic_traits, aes(x = BW, y = SMR))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()+
  geom_smooth(method = "lm")+
  labs(x = "Lobster body size (kg)", y = "Mass-specific SMR")

##############################
# smr
##############################

lm_smr <- lm(log(SMR) ~ temp, data = metabolic_traits)
summary(lm_smr)

pred_smr <- as.data.frame(ggeffects::ggpredict(lm_smr, terms = "temp [11:26]"))


smr_plot <- ggplot() +
  geom_jitter(data = metabolic_traits, aes(x = temp, y = SMR), color = "black", width = 0.4, alpha = 0.2, size = 3) +
  geom_line(data = pred_smr, aes(x = x, y = predicted), lwd = 1.2)+
  geom_ribbon(data = pred_smr, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_errorbar(data = smr_sum, aes(x = temp, ymin = SMR_avg - error, ymax = SMR_avg + error), width = 0.1, size = 1, color = "black") +
  geom_point(data = smr_sum, aes(x = temp,  y = SMR_avg), size = 6, colour = "black") +
  geom_point(data = smr_sum, aes(x = temp,  y = SMR_avg, color = as.factor(temp)), size = 4.5) +
  labs(y = expression(atop("Standard Metabolic Rate", paste("(",mg~O[2]~kg^-1~min^-1,")"))),
       x = "Temperature (°C)") +
  scale_y_continuous(breaks = seq(0, 1.6, by = 0.5)) +
  scale_color_manual(values=c("lightslategray", "lightblue", "lightcoral", "indianred4")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 13),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.7),
        plot.caption = element_text(size = 10, hjust = 0),
        legend.position = "none")+
  scale_x_continuous(limits=c(10, 28),breaks=c(11,16,21,26))


##############################
# mmr panel
##############################

lm_mmr <- lm(MMR ~ temp + I(temp^2), data = metabolic_traits)
summary(lm_mmr)

pred_mmr <- as.data.frame(ggeffects::ggpredict(lm_mmr, terms = "temp [11:26]"))

mmr_plot <- ggplot() +
  geom_jitter(data = metabolic_traits, aes(x = temp, y = MMR), color = "black", width = 0.4, alpha = 0.2, size = 3) +
  geom_line(data = pred_mmr, aes(x = x, y = predicted), lwd = 1.2)+
  geom_ribbon(data = pred_mmr, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_errorbar(data = mmr_sum, aes(x = temp, ymin = MMR_avg - error, ymax = MMR_avg + error), width = 0.1, size = 1, color = "black") +
  geom_point(data = mmr_sum, aes(x = temp,  y = MMR_avg), size = 6, color = "black") +
  geom_point(data = mmr_sum, aes(x = temp,  y = MMR_avg, color = as.factor(temp)), size = 4.5) +
  labs(y = expression(atop("Maximum Metabolic Rate", paste("(",mg~O[2]~kg^-1~min^-1,")"))),
       x = "Temperature (°C)") +
  ylim(1.5, 6) +
  scale_color_manual(values=c("lightslategray", "lightblue", "lightcoral", "indianred4")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 13),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.7),
        plot.caption = element_text(size = 10, hjust = 0),
        legend.position = "none")+
  scale_x_continuous(limits=c(10, 28),breaks=c(11,16,21,26)) 


##############################
# aas panel
##############################

lm_aas <- lm(AAS ~ temp + I(temp^2), data = metabolic_traits)
summary(lm_aas)

pred_aas <- as.data.frame(ggeffects::ggpredict(lm_aas, terms = "temp [11:26]"))

aas_plot <- ggplot() +
  geom_jitter(data = metabolic_traits, aes(x = temp, y = AAS), color = "black", width = 0.4, alpha = 0.2, size = 3) +
  geom_line(data = pred_aas, aes(x = x, y = predicted), lwd = 1.2)+
  geom_ribbon(data = pred_aas, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_errorbar(data = aas_sum, aes(x = temp, ymin = AS_avg - error, ymax = AS_avg + error), width = 0.1, size = 1, color = "black") +
  geom_point(data = aas_sum, aes(x = temp,  y = AS_avg), size = 6, color = "black") +
  geom_point(data = aas_sum, aes(x = temp,  y = AS_avg, color = as.factor(temp)), size = 4.5) +
  labs(y = expression(atop("Absolute Aerobic Scope", paste("(",mg~O[2]~kg^-1~min^-1,")"))),
       x = "Temperature (°C)") +
  ylim(1.5, 5.5) +
  scale_color_manual(values = c("lightslategray", "lightblue", "lightcoral", "indianred4")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 13),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.7),
        plot.caption = element_text(size = 10, hjust = 0),
        legend.position = "none")+
  scale_x_continuous(limits=c(10, 28),breaks=c(11,16,21,26)) 


##############################
# fas panel
##############################

lm_fas <- lm(FAS ~ temp, data = metabolic_traits)
summary(lm_fas)

pred_fas <- as.data.frame(ggeffects::ggpredict(lm_fas, terms = "temp [11:26]"))

fas_plot <- ggplot() +
  geom_jitter(data = metabolic_traits, aes(x = temp, y = FAS), color = "black", width = 0.4, alpha = 0.2, size = 3) +
  geom_line(data = pred_fas, aes(x = x, y = predicted), lwd = 1.2)+
  geom_ribbon(data = pred_fas, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_errorbar(data = fas_sum, aes(x = temp, ymin = FAS_avg - error, ymax = FAS_avg + error), width = 0.1, size = 1, color = "black") +
  geom_point(data = fas_sum, aes(x = temp,  y = FAS_avg), size = 6, color = "black") +
  geom_point(data = fas_sum, aes(x = temp,  y = FAS_avg, color = as.factor(temp)), size = 4.5) +
  labs(y = expression(atop("Factorial Aerobic Scope", paste("(MMR/SMR)"))),
       x = "Temperature (°C)") +
  scale_color_manual(values = c("lightslategray", "lightblue", "lightcoral", "indianred4")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 13),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.7),
        plot.caption = element_text(size = 10, hjust = 0),
        legend.position = "none")+
  scale_x_continuous(limits=c(10, 28),breaks=c(11,16,21,26)) 



#------------------------------------------------------------------------------------------------
## Resting heart rates stuffs
#------------------------------------------------------------------------------------------------

resting_heart_rates <- read_csv(here::here("data", "heart_rate", "outputs", "resting_heart_rates.csv"))

##############################
# add dummy data for 26C
##############################

dummy_data <- c("NA", "26", "NA", "NA", "NA")

# combine data with dummy data
resting_heart_rates2 <- rbind(resting_heart_rates, dummy_data) %>% 
  mutate(rhr = as.numeric(as.character(rhr)))

resting_heart_rates3<-resting_heart_rates2[c(1:9),]
resting_heart_rates3$temp<-as.numeric(resting_heart_rates3$temp)

##############################
# rhr panel
##############################

lm_rhr <- lm(rhr ~ temp, data = resting_heart_rates3)
summary(lm_rhr)

pred_rhr <- as.data.frame(ggeffects::ggpredict(lm_rhr, terms = "temp [11:21]"))

rhr_plot <- ggplot() +
  geom_jitter(data = resting_heart_rates3, mapping = aes(x = temp, y = rhr), alpha = 0.2, width = 0.4, size = 3) +
  geom_line(data = pred_rhr, aes(x = x, y = predicted), lwd = 1.2)+
  geom_ribbon(data = pred_rhr, aes(x = x, y = predicted, ymin = conf.low, ymax = conf.high), alpha = 0.1)+
  geom_errorbar(data = rhr_sum, aes(x = temp, ymin = hr_avg - error, ymax = hr_avg + error), width = 0.1, size = 1, color = "black") +
  geom_point(data = rhr_sum, aes(x = temp,  y = hr_avg), size = 6, color = "black") +
  geom_point(data = rhr_sum, aes(x = temp,  y = hr_avg, color = as.factor(temp)), size = 4.5) +
  labs(y = expression(atop("Resting Heart Rate", paste("(beats",~min^-1,")"))),
       x = "Temperature (°C)") +
  scale_color_manual(values = c("lightslategray", "lightblue", "lightcoral")) +
  theme_classic() +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 13),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.7),
        plot.caption = element_text(size = 10, hjust = 0),
        legend.position = "none")+
  scale_x_continuous(limits=c(10, 28),breaks=c(11,16,21,26)) 


##############################
# combine panels and save plot
##############################

mr_hr_plot <- cowplot::plot_grid(smr_plot, mmr_plot, aas_plot, fas_plot, rhr_plot, nrow = 3, align = "vh", labels = "AUTO")

cowplot::save_plot(here::here("figures", "main_text", "Fig2.pdf"), mr_hr_plot, base_width = 8, base_height = 9)
