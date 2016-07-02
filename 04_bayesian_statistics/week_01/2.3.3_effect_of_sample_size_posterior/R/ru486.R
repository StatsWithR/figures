# models ------------------------------------------------------------
p <- seq(0.1, 0.9, 0.1)

# priors ------------------------------------------------------------
prior <- c(rep(0.06, 4), 0.52, rep(0.06, 4))

# likelihood --------------------------------------------------------
likelihood <- dbinom(4, size = 20, prob = p)

# posterior ---------------------------------------------------------
numerator <- prior * likelihood
denominator <- sum(numerator)
posterior <- numerator / denominator
sum(posterior)

# plot prior, likelihood, posterior ---------------------------------
pdf("ru486_pri_lik_post.pdf", width = 10, height = 3)
par(mfrow = c(1,3), bg = NA)
barplot(prior, names.arg = p, las = 2, main = "Prior")
barplot(likelihood, names.arg = p, las = 2, main = "Likelihood")
barplot(posterior, names.arg = p, las = 2, main = "Posterior")
dev.off()

# more data: x2 -----------------------------------------------------
likelihood <- dbinom(4*2, size = 20*2, prob = p)
numerator <- prior * likelihood
denominator <- sum(numerator)
posterior <- numerator / denominator

pdf("ru486_pri_lik2n_post.pdf", width = 10, height = 3)
par(mfrow = c(1,3), bg = NA)
barplot(prior, names.arg = p, las = 2, main = "Prior")
barplot(likelihood, names.arg = p, las = 2, main = "Likelihood")
barplot(posterior, names.arg = p, las = 2, main = "Posterior")
dev.off()

# more data: x10 ----------------------------------------------------
likelihood <- dbinom(4*10, size = 20*10, prob = p)
numerator <- prior * likelihood
denominator <- sum(numerator)
posterior <- numerator / denominator

pdf("ru486_pri_lik10n_post.pdf", width = 10, height = 3)
par(mfrow = c(1,3), bg = NA)
barplot(prior, names.arg = p, las = 2, main = "Prior")
barplot(likelihood, names.arg = p, las = 2, main = "Likelihood")
barplot(posterior, names.arg = p, las = 2, main = "Posterior")
dev.off()