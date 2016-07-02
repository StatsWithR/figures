# posterior ---------------------------------------------------------

posterior <- c(47, 33, 35, 32, 19, 33, 34, 36, 47, 32, 35, 41, 32, 29, 35, 
               25, 32, 36, 20, 47, 37, 32, 35, 25, 37, 40, 36, 38, 40, 35, 49, 
               23, 33, 35, 38, 28, 36, 4, 28, 45, 37, 39, 34, 41, 28, 33, 27, 
               26, 30, 34, 23)

length(posterior) # 51

t(sort(posterior))

mean(posterior) # 33.45098
median(posterior) # 34
table(posterior)[which.max(table(posterior))] # 35

# dotplot of posterior ----------------------------------------------

pdf("posterior.pdf", width = 10, height = 3)
par(mar = c(2, 0, 0, 0), cex.axis = 1.5, cex = 1.5)
BHH2::dotPlot(posterior, pch = 19, xlim = c(0, 50), axes = FALSE)
axis(1, at = seq(from = 0, to = 50, by = 5))
abline(v = mean(posterior), col = "orange", lwd = 4)
abline(v = median(posterior), col = "turquoise4", lwd = 4)
abline(v = 35, col = "pink", lwd = 4)
legend("topleft", col = c("orange", "turquoise4", "pink"), 
       c("mean", "median", "mode"), lty = 1, lwd = 4,
       bty = "n")
dev.off()

# g = 30 ------------------------------------------------------------

g = 30

# L0 for g ----------------------------------------------------------

(abs(sort(posterior)-g) > 1e-6)[c(1,2,3,12,13,14,50,51)]

sum( abs(posterior-g) > 1e-6 )

# L0 ----------------------------------------------------------------

pdf("L0.pdf", width = 10, height = 5)
par(mfrow = c(2, 1), cex.x.axis = 1.5, cex = 1.5, mar = c(2, 4, 0.5, 0), las = 1)
s = seq(0, 50, by = 0.01)
L0 = sapply(s, function(s) sum( abs(posterior - s) > 1e-6 ))
plot(s, L0, ylab = "L0", type = "l", axes = FALSE, xlim = c(0, 50))
axis(2, at = c(44, 46, 48, 50))
#
BHH2::dotPlot(posterior, pch = 19, xlim = c(0, 50), axes = FALSE)
axis(1, at = seq(from = 0, to = 50, by = 5))
abline(v = mean(posterior), col = "orange", lwd = 4)
abline(v = median(posterior), col = "turquoise4", lwd = 4)
abline(v = 35, col = "pink", lwd = 4)
legend("topleft", col = c("orange", "turquoise4", "pink"), 
       c("mean", "median", "mode"), lty = 1, lwd = 4,
       bty = "n")
dev.off()

# L1 for g = 30 -----------------------------------------------------

abs(sort(posterior) - g)[c(1,2,3,12,13,14,50,51)]

sum(abs(sort(posterior) - g))

# L1 ----------------------------------------------------------------

pdf("L1.pdf", width = 10, height = 5)
par(mfrow = c(2, 1), cex.x.axis = 1.5, cex = 1.5, mar = c(2, 4, 0.5, 0), las = 1)
s = seq(0, 50, by = 0.01)
L1 = sapply(s,function(s) sum( abs(posterior - s) ))
plot(s, L1, ylab = "", type = "l", xlim = c(0, 50), axes = FALSE)
axis(2, at = seq(250, 1650, 350))
#
BHH2::dotPlot(posterior, pch = 19, xlim = c(0, 50), axes = FALSE)
axis(1, at = seq(from = 0, to = 50, by = 5))
abline(v = mean(posterior), col = "orange", lwd = 4)
abline(v = median(posterior), col = "turquoise4", lwd = 4)
abline(v = 35, col = "pink", lwd = 4)
legend("topleft", col = c("orange", "turquoise4", "pink"), 
       c("mean", "median", "mode"), lty = 1, lwd = 4,
       bty = "n")
dev.off()

# L2 for g = 30 -----------------------------------------------------

((sort(posterior) - g)^2)[c(1,2,3,12,13,14,50,51)]

sum((sort(posterior) - g)^2)

# L2 ----------------------------------------------------------------

pdf("L2.pdf", width = 10, height = 5)
par(mfrow = c(2, 1), cex.x.axis = 1.5, cex = 1.5, mar = c(2, 4, 0.5, 0), las = 1)
s = seq(0, 50, by = 0.01)
L2 = sapply(s,function(s) sum( (posterior - s)^2 ))
plot(s, L2, ylab = "", type = "l", xlim = c(0, 50), ylim = c(0, 61000), axes = FALSE)
axis(2, at = seq(0, 60000, 20000))
#
BHH2::dotPlot(posterior, pch = 19, xlim = c(0, 50), axes = FALSE)
axis(1, at = seq(from = 0, to = 50, by = 5))
abline(v = mean(posterior), col = "orange", lwd = 4)
abline(v = median(posterior), col = "turquoise4", lwd = 4)
abline(v = 35, col = "pink", lwd = 4)
legend("topleft", col = c("orange", "turquoise4", "pink"), 
       c("mean", "median", "mode"), lty = 1, lwd = 4,
       bty = "n")
dev.off()



