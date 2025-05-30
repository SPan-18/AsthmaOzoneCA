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
post_z <- mod1$samples$z
Xcov_names <- paste(months_str, "/Jan", sep = "")
Xcov_names[1] <- "Jan"
rownames(post_beta) <- Xcov_names
