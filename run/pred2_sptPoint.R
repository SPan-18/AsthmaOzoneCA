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
