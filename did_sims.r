

set.seed(12345)
library(tidyverse)
library(parallel)
library(did)
library(bacondecomp)
library(patchwork)


# Construct base sample.
delta1 <- 0.10
delta2 <- 1.00
delta3 <- 0.50

dat <- tibble(pid = c(0, 1, 2, 3, 4),
              group = c(0, 1, 1, 2, 2),
              mint = c(NA, 2, 2, 4, 4)) %>%
  expand_grid(period = 0:5, tfx = 0, d = 0) %>%
  mutate(tfx = ifelse(group == 1 & period %in% c(2,3), delta1, tfx),
         tfx = ifelse(group == 1 & period %in% c(4,5), delta2, tfx),
         tfx = ifelse(group == 2 & period %in% c(4,5), delta3, tfx),
         d = ifelse(group == 1 & period >= 2, 1, d),
         d = ifelse(group == 2 & period >= 4, 1, d)) %>%
  mutate(et = ifelse(!is.na(mint), period - mint, -1)) %>%
  mutate(et = relevel(as_factor(et), ref = "-1"))

dat$y <- dat$tfx + rnorm(n = nrow(dat), mean = 0, sd = 1 / sqrt(10000))


# Plot group time series.
plot_dat <- dat %>%
  group_by(group, period) %>%
  summarize(y = mean(y)) %>%
  mutate(Group = as_factor(group), Time = as_factor(period))
  
p1 <- ggplot(plot_dat, aes(x = Time, y = y, group = Group, color = Group)) +
  geom_point(alpha = 0.5, size = 4) +
  geom_step(alpha = 0.5, size = 4) +
  ylab("Outcome") +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))

p1
ggsave("./dgp.pdf")


# Goodman-Bacon decomposition
bacon_decomp <- bacon(y ~ d, data = dat, 
                      id_var = "pid", 
                      time_var = "period")

att <- mean(filter(dat, d == 1)$tfx)

p2 <- ggplot(bacon_decomp) +
  aes(x = weight, y = estimate, color = factor(type)) +
  labs(x = "Weight", y = "Estimate", color = "Type") +
  geom_point(size = 4) +
  annotate(geom = "curve", x = 0.24, y = -0.2, xend = 0.25, yend = -0.34, 
           curvature = -.3, arrow = arrow(length = unit(2, "mm")),
           size = 1) +
  annotate(geom = "text", x = 0.24, y = -0.2, 
           label = "-ATT(1,4) & -ATT(1,5)", 
           hjust = "right", 
           size = 6) +  
  xlim(c(0.20, 0.25)) +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16)) +
  scale_color_manual(values=c("#7fc97f", "#beaed4", "#fdc086"))
p2

p1 + p2
ggsave("./bacon_decomp.pdf")


# Simulate DID & event study estimates.
sim_func <- function(iter, dat) {
  dat$y <- dat$tfx + rnorm(n = nrow(dat), mean = 0, sd = 1 / sqrt(10000))
  
  d <- lm(y ~ d + as_factor(group) + as_factor(period), dat)$coefficients["d"]
  es <- lm(y ~ et + as_factor(group) + as_factor(period), dat)$coefficients[2:7]
  return(c(d, es))
}

sims <- mclapply(1:1000, sim_func, dat = dat, mc.cores = 8) %>%
  bind_rows()

# DID
did_sim <- sims %>% select(d)
ggplot(did_sim, aes(d)) +
  geom_density(color = "#08589e", fill = "#08589e", size = 2, alpha = 0.5) +
  geom_vline(xintercept = att, size = 2) + 
  geom_vline(xintercept = 0, size = 2, linetype = "dashed") + 
  ylab("Density") +
  xlab("DID Estimate") +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))
ggsave("./did_sims.pdf")


# Event study
es_sim <- sims %>% 
  select(starts_with("et")) %>%
  pivot_longer(cols = everything(),
               names_to = "time",
               names_prefix = "et") %>%
  mutate(time = as.numeric(time)) %>%
  group_by(time) %>%
  summarize(est = mean(value), 
            l95 = quantile(value, probs = 0.025),
            u95 = quantile(value, probs = 0.975)) %>%
  add_row(time = -1, est = 0, u95 = 0, l95 = 0)

ggplot(es_sim) +
  geom_vline(xintercept = -0.5, linetype = "dashed", size = 1) + 
  geom_hline(yintercept = 0, size = 1) + 
  geom_ribbon(aes(x = time, ymin = l95, ymax = u95),
              color = "#08589e", fill = "#08589e", 
              alpha = 0.5, size = 2) +
  geom_line(aes(x = time, y = est), size = 2) +
  geom_point(aes(x = time, y = est), size = 4) +   
  ylab("Estimate") +
  xlab("Event Time") +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 16),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))
ggsave("./es_sims.pdf")





