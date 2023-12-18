

# step 1: loading data and subroutines
source("chapter_4_coordinate_transform.R")
source("chapter_4_spatial_estimation_sub.R")

sample_2d <- obs_surface[, c("x", "y", "water_temperature")]

sample_var <- var(sample_2d$water_temperature)
pair_wt <- comb_2(sample_2d$water_temperature)


dist_xy <- as.numeric(dist(sample_2d[, 1:2]))

vd <- data.frame(dist_xy = dist_xy, pair_wt)

head(vd)


# h cloud plot
op <- par(no.readonly = T)
par(mar = c(5,5,3,2))
plot(vd$dist_xy, abs(vd$V1 - vd$V2),
     xlab = "Distance (m)",
     ylab = expression(
       paste(wt[i], " - ", wt[j], " (", i!=j, ")")
       ),
     cex.lab = 2, cex.axis = 1.5,
     pch = 16, las = 1)
abline(v = seq(2000, 20000, 2000), lty = 2, lwd = 1)
par(op)



# conceptual
################################################################################
lag_d_interval <- lag_dist(vd$dist_xy, buffer = 2000, m_factor = 2)


vg_e <- lapply(1:nrow(lag_d_interval),
               function(x){
                 
                 lag_d_id <- which(vd$dist_xy >= lag_d_interval$lag_st[x] &
                                     vd$dist_xy < lag_d_interval$lag_ed[x])
                 mean_lag_d <- mean(vd[lag_d_id, 1], na.rm = T)
                 cov <- sample_var -
                   cov(vd[lag_d_id, 2], vd[lag_d_id, 3], use = "everything")
                 
                 res <- c(lag_d = mean_lag_d, semi_var = cov)
                 return(res)
                 
               }) %>%
  bind_rows %>%
  filter(complete.cases(.)) %>%
  as.data.frame


# experimental variogram
op <- par(no.readonly = T)
par(mar = c(5,5,5,5))
plot(vg_e,
     ylim = c(0, 2),
     xlab = "Distance (m)",
     ylab =  expression(
       paste(wt[i], " - ", wt[j], " (", i!=j, ")")
     ),
     cex = 2, cex.lab = 2, cex.axis = 1.5,
     pch = 16, las = 1)
points(vd$dist_xy, abs(vd$V1 - vd$V2), pch = 16, col = "gray50", cex = 1)
axis(4, cex.axis = 1.5, las = 1)
mtext("semi-variance", side = 4, line = 3, cex = 2)
abline(v = seq(2000, 20000, 2000), lty = 2, lwd = 1)
abline(v = seq(1000, 20000, 2000), lty = 2, col = 2, lwd = 1)

legend("bottomright",
       inset = c(0,1), xpd = TRUE, ncol = 2,
       legend = c(
         expression(paste(wt[i], " - ", wt[j], " (", i!=j, ")")),
         "semi-variance",
         "Separation Distance",
         "Overlapped Separation Distance"),
       col = c("gray50", 1, 1, 2), cex = 1.3,
       pch = c(16, 16, NA, NA), lty = c(NA, NA, 2, 2), bg = "white", bty = "n")

par(op)
################################################################################


# sph model
source("chapter_4_variogram_functions.R")

ex_vgm <- vge(xi = sample_2d, s_var = 1:2, val = 3, valid_n = 3,
              max_range = 1, lag_buffer = 2000, overlap_factor = 2,
              plot = TRUE)

# ignore outliers
ex_vgm <- ex_vgm[-c(1,14),]


# parameter optimization
par0 <- c(0.5, 0.5, 5000)
lower <- c(0.1, 0.1, 5000)
upper <- c(1, 2, 15000)

vg_f_sph <- optim(par0, Sph,
                  d = ex_vgm$lag_s, gamma = ex_vgm$gamma,
                  return = "error",
                  gr=NULL, method = "L-BFGS-B",
                  lower = lower, upper = upper)

# fitted Sph model
new_d <- seq(0, max(ex_vgm$lag_s), length.out = 100)
vg_f_sph_new <- Sph(new_d, pars = vg_f_sph$par)


op <- par(no.readonly = T)
par(mar = c(5,5,3,2))
plot(ex_vgm$lag_s, ex_vgm$gamma,
     xlim = range(new_d), ylim = c(0, max(vg_e$semi_var)),
     xlab = "Separation Distance (m)", ylab = "semi-Variance",
     las = 1,
     type = "p", pch = 16, lwd = 2, col = "gray70",
     cex = 2, cex.main = 2, cex.lab = 2, cex.axis = 1.5)
lines(new_d, vg_f_sph_new, col = 1, lwd = 2)


# exp model
vg_f_exp <- optim(par0, Exp,
                  d = ex_vgm$lag_s, gamma = ex_vgm$gamma,
                  return = "error",
                  gr=NULL, method = "L-BFGS-B",
                  lower = lower, upper = upper)

# fitted Exp model
vg_f_exp_new <- Exp(new_d, pars = vg_f_exp$par)


# gauss model
vg_f_gau <- optim(par0, Gau,
                  d = ex_vgm$lag_s, gamma = ex_vgm$gamma,
                  return = "error",
                  gr=NULL, method = "L-BFGS-B",
                  lower = lower, upper = upper)


vg_f_gau_new <- Gau(new_d, pars = vg_f_gau$par)


# remaining results
lines(new_d, vg_f_exp_new, col = 2, lwd = 2)
lines(new_d, vg_f_gau_new, col = 4, lwd = 2)
legend("bottomright",
       legend = c("Sph", "Exp", "Gau"),
       lty = 1, col = c(1, 2, 4), lwd = 2, cex = 1.5)
par(op)



# wrapped function
vg_fit <- vgt(ex_vgm = ex_vgm,
              vgm_list = list(vgm_s = "Sph"),
              par0 = c(0.5, 0.5, 5000),
              lower = c(0.1, 0.1, 5000),
              upper = c(1, 2, 15000),
              plot = TRUE)


# variogram error 
sph_err <- Sph(d = ex_vgm$lag_s, gamma = ex_vgm$gamma,
               pars = vg_f_sph$par, return = "error")
exp_err <- Exp(d = ex_vgm$lag_s, gamma = ex_vgm$gamma,
               pars = vg_f_sph$par, return = "error")
gau_err <- Gau(d = ex_vgm$lag_s, gamma = ex_vgm$gamma,
               pars = vg_f_sph$par, return = "error")

sph_err; exp_err; gau_err

