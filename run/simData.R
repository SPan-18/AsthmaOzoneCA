# rm(list = ls())

library(ggplot2)
library(ggpubr)
library(dplyr)

set.seed(1728)

n_sites <- 50
n_months <- 12
n_avg <- 30

phi_s <- 4
phi_t <- 3/5
nu <- 0.5
sigmasqSp <- 1
nugget <- 1

sites_coords <- cbind(runif(n_sites), runif(n_sites))
# plot(sites_coords)
time_coords <- seq(0, n_months, length.out = n_avg*n_months)

sp_D <- as.matrix(dist(sites_coords))
tm_D <- as.matrix(dist(time_coords))
Vsp <- geoR::matern(sp_D, 1/phi_s, nu)
Vtm <- geoR::matern(tm_D, 1/phi_t, 0.5)
Lsp <- chol(Vsp)
Ltm <- chol(Vtm)
# Lst <- chol(kronecker(Vsp, Vtm))

# tm_trend <- 5 + 2 * sin(0.3 * pi * time_coords - 0.5)
# plot(time_coords, tm_trend)
# trend_df <- data.frame(x = time_coords, y = tm_trend)
# write.csv(trend_df, "trend-df.csv", row.names = FALSE)

# Define periodic kernel
period <- 7
lengthscale <- 10
K <- outer(time_coords, time_coords, function(t1, t2)
  exp(-2 * sin(pi * abs(t1 - t2) / period)^2 / lengthscale^2))
tm_trend2 <- MASS::mvrnorm(1, mu = rep(0, length(time_coords)), Sigma = 10*K)
tm_trend2 <- as.vector(scale(tm_trend2))
tm_trend2 <- 5 + 2*tm_trend2
# plot(time_coords, tm_trend2)
trend_df <- data.frame(x = time_coords, y = tm_trend2)
# write.csv(trend_df, "trend-df.csv", row.names = FALSE)
tm_trend <- tm_trend2

Y_st <- array(dim = c(n_sites, length(time_coords)))
mu_st <- array(dim = c(n_sites, length(time_coords)))
for(i in 1:n_sites){
  mu_st[i, ] = tm_trend
}
w <- matrix(rnorm(n_sites * length(time_coords), 0, sqrt(sigmasqSp)),
            nrow = n_sites, ncol = length(time_coords))

w <- w %*% Ltm
w <- crossprod(Lsp, w)
mu_st <- mu_st + w
Y_st[] <- rnorm(length(mu_st), mean = mu_st, sd = sqrt(nugget))

sim_df <- data.frame(Y = as.vector(Y_st),
                     Time = rep(time_coords, each = n_sites),
                     Month = rep(rep(1:n_months - 0.5, each = n_avg), each = n_sites),
                     Sites = rep(1:n_sites, times = length(time_coords)),
                     s1 = rep(sites_coords[, 1], times = length(time_coords)),
                     s2 = rep(sites_coords[, 2], times = length(time_coords)),
                     w = as.vector(w),
                     mu = as.vector(mu_st))

# write.csv(sim_df, "sim-big.csv", row.names = FALSE)

sim_df_avg <- sim_df %>%
  group_by(Month, Sites, s1, s2) %>%
  summarise(Y_avg = mean(Y, na.rm = TRUE), .groups = "drop")

sim_df_avg <- sim_df_avg %>%
  mutate(t_a = Month - 0.5, t_b = Month + 0.5) %>%
  relocate(t_a, t_b, .after = s2) %>% 
  relocate(Month, .after = Sites)

# write.csv(sim_df_avg, "test-sim-full.csv", row.names = FALSE)

n_delete <- floor(0.1 * nrow(sim_df_avg))
rows_to_delete <- sample(nrow(sim_df_avg), n_delete)
sim_df_avg <- sim_df_avg[-rows_to_delete, ]

# write.csv(sim_df_avg, "test-sim-big.csv", row.names = FALSE)

months_str <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
months_str2 <- paste(months_str, "1", sep = " ")

common_limits <- range(sim_df$Y)
common_breaks <- seq(0, 12, by = 4)
common_colors <- RColorBrewer::brewer.pal(11, "RdYlGn")

snap <- time_coords[121]
snap_text <- "May 1"

qmids <- c(1.5, 4.5, 7.5, 10.5)
qlabels <- c("Q1", "Q2", "Q3", "Q4")
shade_regions <- data.frame(
  xmin = c(0, 6),
  xmax = c(3, 9)
)

p1 <- ggplot() +
  geom_rect(data = shade_regions, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
            fill = "gray90", alpha = 0.5, inherit.aes = FALSE) +
  annotate("text", x = qmids, y = max(sim_df$Y), label = qlabels, 
           hjust = 0.5, vjust = 0.5, size = 3.5, family = "sans") +
  geom_line(data = sim_df, aes(x = Time, y = Y, group = factor(Sites)),
            color = "#abd9e9", alpha = 0.2, linetype = "dotted", linewidth = 0.2) +
  geom_point(data = sim_df, aes(x = Time, y = Y, group = factor(Sites)),
             color = "#abd9e9", alpha = 0.2, size = 0.25) +
  # scale_color_distiller(palette = "RdYlGn", direction = -1, 
  #                       limits = common_limits, breaks = common_breaks,
  #                       oob = scales::squish, guide = "none") +
  labs(color = latex2exp::TeX('$X(s, t)$')) +
  geom_vline(xintercept = snap, color = "gray") +
  annotate("text", x = snap, y = min(sim_df$Y), label = snap_text, 
           hjust = -0.1, vjust = 0.5, size = 3, family = "sans") +
  # geom_line(data = sim_df_avg, aes(x = Month, y = Y_avg, group = factor(Sites)),
  #           color = "#2c7bb6", alpha = 0.3, linetype = "dotted", linewidth = 0.5) +
  geom_point(data = sim_df_avg, aes(x = Month, y = Y_avg, group = factor(Sites)),
             color = "#2c7bb6", alpha = 0.3, size = 0.5) +
  geom_line(data = trend_df, aes(x = x, y = y), 
            color = "black") +
  theme_bw() +
  xlab("Month") + ylab(latex2exp::TeX('$X(s, t)$')) +
  scale_x_continuous(breaks = 1:n_months - 0.5,
                     labels = months_str) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 10))

p2 <- ggplot() +
  geom_line(data = sim_df_avg, aes(x = Month, y = Y_avg, group = factor(Sites)),
            color = "#2c7bb6", alpha = 0.2, linetype = "dotted", linewidth = 0.5) +
  geom_point(data = sim_df_avg, aes(x = Month, y = Y_avg, group = factor(Sites)),
             color = "#2c7bb6", alpha = 0.2, size = 0.5) +
  geom_line(data = trend_df, aes(x = x, y = y), 
            color = "black") +
  theme_bw() +
  xlab("Month") + ylab(latex2exp::TeX('$X(s, I)$')) +
  scale_x_continuous(breaks = 1:n_months - 0.5,
                     labels = months_str) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 10))


sim_df_snap <- sim_df %>%
  filter(Time == snap)

p2.1 <- spStack::surfaceplot(tab = sim_df_snap, 
                             coords_name = c("s1", "s2"), 
                             var_name = "Y", mark_points = TRUE)
p2.1 <- p2.1 +
  scale_fill_distiller(palette = "RdYlGn", direction = -1) +
  # scale_fill_distiller(palette = "RdYlGn", direction = -1, 
  #                      limits = common_limits, breaks = common_breaks,
  #                      oob = scales::squish) +
  xlab(latex2exp::TeX('$s_1$')) + 
  ylab(latex2exp::TeX('$s_2$')) +
  labs(fill = latex2exp::TeX('$X(s, t)$')) +
  annotate("text", x = 1, y = 0, label = snap_text, 
           hjust = 1, vjust = 0, size = 3, family = "sans")

# ggsave("../plots/sim_oz_panel.pdf", height = 3.5, width = 12, units = "in")
# system2(command = "pdfcrop",
#         args    = c("../plots/sim_oz_panel.pdf",
#                     "../plots/sim_oz_panel.pdf"))

