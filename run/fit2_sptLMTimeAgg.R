
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
# print(t(apply(post_beta, 1, function(x) quantile(x, c(0.025, 0.5, 0.975)))))
post_z <- mod2$samples$z
