

h <- seq(0.1, 1.0, 0.1)
p <- c(0, 1, 2, 3)

inv_0 <- (1/h)^p[1]
lambda_0 <- inv_0/sum(inv_0)
lambda_0


inv_1 <- (1/h)^p[2]
lambda_1 <- inv_1/sum(inv_1)
round(lambda_1, 2)


inv_2 <- (1/h)^p[3]
lambda_2 <- inv_2/sum(inv_2)
round(lambda_2, 2)


inv_3 <- (1/h)^p[4]
lambda_3 <- inv_3/sum(inv_3)
round(lambda_3, 2)



plot(h, lambda_0,
     type = "l", lty = 1, lwd = 2,
     ylim = c(0, 1), ylab = "", xlab = "h (Distance)", las = 1,
     cex.axis = 1.5, cex.lab = 2)

lines(h, lambda_1, lwd = 2, lty = 2)
lines(h, lambda_2, lwd = 2, lty = 3)
lines(h, lambda_3, lwd = 2, lty = 4)

legend("topright", legend = paste0("p = ", c(0:3)),
       lty = 1:4, lwd = 5, cex = 2, seg.len = 5)
