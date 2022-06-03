source(here::here("code", "00_libraries.R"))

df.treat <- read_csv("data/foraging/outputs/posteriors_treatments_noGP03.csv") %>%
  mutate(c_max = 1/t.h*24) %>%
  separate(temp, into = c("junk", "temp_c"), sep = "(?<=[A-Za-z])(?=[0-9])") %>%
  select(-junk) %>%
  mutate(temp_c = as.numeric(temp_c)) %>%
  as_tibble()

meta <- read.csv(here::here("data/metadata/", "metadata_lobster_sizes.csv")) %>%
  select(lobster_id, weight_kg, temp) %>%
  rename(id = lobster_id)

df.ind <- read.csv(here::here("data/foraging/outputs/posteriors_individuals_noGP03.csv")) %>%
  mutate(c_max = 1/h*24) %>%
  as_tibble() %>%
  left_join(meta) %>%
  rename(temp_c = temp) %>%
  mutate(c_max = c_max/weight_kg)



foraging_data <- read.csv(here::here("data", "foraging", "raw", "foraging_assay_data.csv")) %>%
  select(temp, lobster_id, Initial, Killed) %>%
  rename(id = lobster_id, initial = Initial, killed = Killed) %>%
  arrange(temp, id, initial) %>%
  filter(id != "GP03") %>%
  mutate(temp = case_when(temp == 11 ~ "11°C", 
                          temp == 16 ~ "16°C", 
                          temp == 21 ~ "21°C", 
                          temp == 26 ~ "26°C"))


treats <- unique(df.treat$temp_c)
ids <- unique(df.ind$id)

library(rethinking)
calculate_CI<- function(a, h, prob = 0.95, name, ...){
  N.vec <- seq(0, 60, length.out = 100)
  p_sim_output <- sapply(N.vec, function(i) a*i / (1 + a*h*i))
  p_mu <- apply(p_sim_output, 2 ,median, na.rm = T)
  p_ci <- t(apply( p_sim_output , 2 , PI, prob = prob))
  
  return(data.frame(N = N.vec, mu = p_mu, mu.lower = p_ci[,1], mu.upper = p_ci[,2], name = name[i]))
}

out <- list()
for(i in 1:length(treats)){
  temp.a <- df.treat$t.a[df.treat$temp_c == treats[i]]
  temp.h <- df.treat$t.h[df.treat$temp_c == treats[i]]
  out[[i]] <- calculate_CI(a = temp.a, h = temp.h, name = treats)
}

pred.treat <- data.frame(do.call(rbind, out)) %>%
  as_tibble() %>%
  rename(temp = name) %>%
  mutate(mu = mu*24, 
         mu.lower = mu.lower*24, 
         mu.upper = mu.upper*24) %>%
  mutate(temp = case_when(temp == 11 ~ "11°C", 
                          temp == 16 ~ "16°C", 
                          temp == 21 ~ "21°C", 
                          temp == 26 ~ "26°C"))

out <- list()
for(i in 1:length(ids)){
  temp.a <- df.ind$a[df.ind$id == ids[i]]
  temp.h <- df.ind$h[df.ind$id == ids[i]]
  out[[i]] <- calculate_CI(a = temp.a, h = temp.h, name = ids)
}

pred.ind <- data.frame(do.call(rbind, out)) %>%
  as_tibble() %>%
  rename(id = name) %>%
  left_join(meta %>% select(id, temp)) %>%
  mutate(mu = mu*24, 
         mu.lower = mu.lower*24, 
         mu.upper = mu.upper*24) %>%
  mutate(temp = case_when(temp == 11 ~ "11°C", 
                          temp == 16 ~ "16°C", 
                          temp == 21 ~ "21°C", 
                          temp == 26 ~ "26°C"))


functional_response_plot <- ggplot() +
  geom_jitter(data = foraging_data, aes(x = initial, y = killed), alpha = 0.7, shape = 1, width = 2) + 
  geom_line(data = pred.ind, aes(x = N, y = mu, group = id), alpha = 0.7, color = "darkgray", size = 1.5) +
  geom_line(data = pred.treat, aes(x = N, y = mu, col = factor(temp), group = factor(temp)), size = 2) +
  geom_ribbon(data = pred.treat, aes(x = N, ymin=mu.lower, ymax=mu.upper), linetype=2, alpha=0.1)+
  facet_wrap(~temp, nrow = 1) +
  scale_x_continuous(breaks = seq(0, 60, by = 10)) +
  scale_color_manual(values = c("lightblue", "lightslategray", "lightcoral", "indianred4")) +
  labs(x = "Initial Prey Density", y = expression(atop("Consumption Rate", paste( "\n(prey consumed" ~ predator^{-1} ~ "24" ~ hr^{-1},")")))) +
  theme_classic(base_size = 12) +
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(size = 13),
        axis.title.y = element_text(hjust = 0.5), 
        panel.border = element_rect(colour = "black", fill = NA, size = 0.7),
        strip.text = element_text(size = 12, hjust = 0, face = "bold"),
        strip.background = element_blank(),
        legend.position = "none")

cowplot::save_plot(here::here("figures", "main_text", "Fig4_updated.pdf"), functional_response_plot, base_width = 10, base_height = 3)
