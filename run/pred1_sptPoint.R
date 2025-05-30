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

trend_df2 <- trend_df %>%
  filter(x >= 0.5, x <= 11.5)

p3 <- ggplot() +
  # geom_line(data = dat, aes(x = Month, y = Y_avg, group = factor(Sites)),
  #           color = "#2c7bb6", alpha = 0.2, linetype = "dotted", linewidth = 0.5) +
  geom_point(data = dat, aes(x = Month, y = Y_avg, group = factor(Sites)),
             color = "#2c7bb6", alpha = 0.2, size = 0.5) +
  geom_line(data = trend_df2, aes(x = x, y = y),
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
