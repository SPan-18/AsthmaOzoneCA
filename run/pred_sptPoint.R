# rm(list = ls())
# 
# library(spStackCOS)
# library(dplyr)
# library(ggplot2)

dat <- sim_df_avg
X <- dat$Y_avg
n <- length(X)
dat$Month2 <- factor(dat$Month + 0.5)
dummy_Xcov <- model.matrix(~ Month2, data = dat)
Xcov <- matrix(dummy_Xcov, nrow = nrow(dummy_Xcov))
sp_coords <- as.matrix(dat[, c("s1", "s2")])
time_coords <- as.matrix(dat[, c("t_a", "t_b")])

phi_s0 <- 4
phi_t0 <- 3/5
nu_s0 <- 0.5
deltasq_0 <- 2

mod1 <- sptLMexactTimeAgg(X = X, Xcov = Xcov, spcoords = sp_coords, 
                          timecoords = time_coords, phi_s = phi_s0, 
                          phi_t = phi_t0, nu_s = nu_s0, 
                          noise_sp_ratio = deltasq_0, n.samples = 1000)

post_beta <- mod1$samples$beta
print(t(apply(post_beta, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))))
post_z <- mod1$samples$z


sp_coords_pred <- sp_coords[which(dat$Month == 0.5), ]
time_coords_pred <- as.matrix(jitter(seq(0.1, 11.9, length.out = 100)))
# time_coords_pred <- as.matrix(sort(runif(50, 0.1, 12)))
sp_pred <- sp_coords_pred[rep(1:nrow(sp_coords_pred), times = nrow(time_coords_pred)), ]
time_pred <- time_coords_pred[rep(1:nrow(time_coords_pred), each = nrow(sp_coords_pred)), ]
Xcov_pred <- matrix(0, length(time_pred), ncol(Xcov))
Xcov_pred[, 1] <- 1
for(i in 1:length(time_pred)){
  id <- floor(time_pred[i])
  if(time_pred[i] - id == 0) Xcov_pred[i, id] <- 1
  else Xcov_pred[i, id + 1] <- 1
}

mod1.pred <- sptPointPredictTimeAgg(mod1, Xcov_pred = Xcov_pred,
                                    spcoords_pred = sp_pred,
                                    timecoords_pred = time_pred)

post_zpred <- mod1.pred$samples$zpred
postpred_y <- mod1.pred$samples$xpred

postpred_yt <- matrix(NA, nrow = length(time_coords_pred), ncol = 3)
for(i in seq_len(length(time_coords_pred))){
  idx <- ((i - 1) * nrow(sp_coords_pred) + 1):(i * nrow(sp_coords_pred))
  postpred_yt[i, ] <- quantile(as.vector(postpred_y[idx, ]), c(0.025, 0.5, 0.975))
}
postpred_yt <- as.data.frame(postpred_yt)
colnames(postpred_yt) <- c("lower", "median", "upper")
postpred_yt$Time <- time_coords_pred

trend_df <- read.csv("trend-df.csv")
trend_df <- trend_df %>%
  filter(x >= 0.5, x <= 11.5)
months_str <- c("Jan", "Feb", "Mar", "Apr", "May", "Jun", 
                "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")

p3 <- ggplot() +
  # geom_line(data = dat, aes(x = Month, y = Y_avg, group = factor(Sites)),
  #           color = "#2c7bb6", alpha = 0.2, linetype = "dotted", linewidth = 0.5) +
  geom_point(data = dat, aes(x = Month, y = Y_avg, group = factor(Sites)),
             color = "#2c7bb6", alpha = 0.2, size = 0.5) +
  geom_line(data = trend_df, aes(x = x, y = y),
            color = "black") +
  geom_ribbon(data = postpred_yt, aes(x = Time, y = median,
                                      ymin = lower, ymax = upper),
              fill = "#d7191c", color = "#d7191c", alpha = 0.1, linewidth = 0.25) +
  theme_bw() +
  xlab("Month") + ylab(latex2exp::TeX('$X$')) +
  scale_x_continuous(breaks = 1:12 - 0.5,
                     labels = months_str) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 10))


pseq <- c(5, 6, 7, 8)
Xcov <- array(dim = c(n, 2*length(pseq) + 1))
for(i in 1:n){
  Xcov[i, 1] <- 1
  for(j in 1:length(pseq)){
    a <- as.numeric(dat[i, "t_a"])
    b <- as.numeric(dat[i, "t_b"])
    p <- pseq[j]
    Xcov[i, 2*j] <- (p/(2*pi*abs(b-a)))*(cos(2*pi*a/p) - cos(2*pi*b/p))
    Xcov[i, 2*j + 1] <- (p/(2*pi*abs(b-a)))*(sin(2*pi*b/p) - sin(2*pi*a/p))
  }
}
sp_coords <- as.matrix(dat[, c("s1", "s2")])
time_coords <- as.matrix(dat[, c("t_a", "t_b")])

mod2 <- sptLMexactTimeAgg(X = X, Xcov = Xcov, spcoords = sp_coords, 
                          timecoords = time_coords, phi_s = phi_s0, 
                          phi_t = phi_t0, nu_s = nu_s0, 
                          noise_sp_ratio = deltasq_0, n.samples = 1000)

post_beta <- mod2$samples$beta
print(t(apply(post_beta, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))))
post_z <- mod2$samples$z

Xcov_pred <- array(dim = c(length(time_pred), ncol(Xcov)))
for(i in 1:length(time_pred)){
  Xcov_pred[i, 1] <- 1
  for(j in 1:length(pseq)){
    Xcov_pred[i, 2*j] <- sin(2*pi*time_pred[i]/pseq[j])
    Xcov_pred[i, 2*j + 1] <- cos(2*pi*time_pred[i]/pseq[j])
  }
}

mod2.pred <- sptPointPredictTimeAgg(mod2, Xcov_pred = Xcov_pred,
                                    spcoords_pred = sp_pred,
                                    timecoords_pred = time_pred)

post_zpred <- mod2.pred$samples$zpred
postpred_y <- mod2.pred$samples$xpred

postpred_yt <- matrix(NA, nrow = length(time_coords_pred), ncol = 3)
for(i in seq_len(length(time_coords_pred))){
  idx <- ((i - 1) * nrow(sp_coords_pred) + 1):(i * nrow(sp_coords_pred))
  postpred_yt[i, ] <- quantile(as.vector(postpred_y[idx, ]), c(0.025, 0.5, 0.975))
}
postpred_yt <- as.data.frame(postpred_yt)
colnames(postpred_yt) <- c("lower", "median", "upper")
postpred_yt$Time <- time_coords_pred

p4 <- ggplot() +
  # geom_line(data = dat, aes(x = Month, y = Y_avg, group = factor(Sites)),
  #           color = "#2c7bb6", alpha = 0.2, linetype = "dotted", linewidth = 0.5) +
  geom_point(data = dat, aes(x = Month, y = Y_avg, group = factor(Sites)),
             color = "#2c7bb6", alpha = 0.2, size = 0.5) +
  geom_line(data = trend_df, aes(x = x, y = y),
            color = "black") +
  geom_ribbon(data = postpred_yt, aes(x = Time, y = median,
                                      ymin = lower, ymax = upper),
              fill = "#d7191c", color = "#d7191c", alpha = 0.1, linewidth = 0.25) +
  theme_bw() +
  xlab("Month") + ylab(latex2exp::TeX('$X$')) +
  scale_x_continuous(breaks = 1:12 - 0.5,
                     labels = months_str) +
  theme(panel.grid = element_blank(),
        strip.text = element_text(size = 10))

ggpubr::ggarrange(p3, p4)

